#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@File   : scTrace.py
@Author : XinWang
"""

import warnings

warnings.filterwarnings("ignore")
import os
import pandas as pd
from Bio import Phylo
from ChromoPhyloGen.keyCNA import cell_events
from ChromoPhyloGen.events_define import define_CNA_mechnism
from ChromoPhyloGen.estimateDER import estimate_der
from ChromoPhyloGen.tree_class import TreeClass


def extract_parent_child_pairs(tree):
    pairs = []
    for clade in tree.find_clades(order="level"):
        for child in clade.clades:
            pairs.append((clade.name, child.name))
    return pairs


def calc_der(tree, cna_dir, output_dir, prefix):
    parent_child_pairs = extract_parent_child_pairs(tree)
    merge_edges = pd.DataFrame(parent_child_pairs, columns=["p", "c"])
    merge_edges = merge_edges.loc[~merge_edges.p.isin([None]), :]
    cna_profile = pd.read_table(cna_dir)
    cna_profile.index = cna_profile.iloc[:, 0].map(str) + '_' + cna_profile.iloc[:, 1].map(
        str) + '_' + cna_profile.iloc[:, 2].map(str)
    cna_profile = cna_profile.iloc[:, 3::]
    cna_profile = cna_profile.astype('int')
    cna_profile = cna_profile.transpose()
    tc = TreeClass(merge_edges, cna_profile)
    tc.unique_name()
    tc.update_tree()
    tc.modify_length()
    cell_relation = cell_events(tc)
    cell_relation.to_csv(os.path.join(output_dir, prefix + 'cell_info.txt'), sep='\t', index=False)

    tc.node_data.to_csv(os.path.join(output_dir, prefix + 'all_node_data.txt'), sep='\t', index=True)
    with open(os.path.join(output_dir, prefix + 'cell_tree.newick'), 'w') as f:
        f.write('(' + tc.newick() + ');')

    der = estimate_der(tc, cell_relation)
    der.to_csv(os.path.join(output_dir, prefix + 'error_risk_score.txt'), sep='\t', index=True, header=False)


def chromosome_event(node_data,
                     output_dir='./',
                     newick='nonw',
                     prefix='ChromoPhyloGen_',
                     cancer_type='ALL',
                     cores=1,
                     randome_num=1000,
                     verbose=True):
    if verbose:
        print('Scoring the chromosomal rearrangements.')

    if newick == 'none':
        cna_profile = pd.read_table(node_data)
        cna_profile.index = cna_profile.iloc[:, 0].map(str) + '_' + cna_profile.iloc[:, 1].map(
            str) + '_' + cna_profile.iloc[:, 2].map(str)
        cna_profile = cna_profile.iloc[:, 3::]
        cna_profile = cna_profile.astype('int')
        cna_profile = cna_profile.transpose()
        sctc_obj = {'tree': None, 'cnv_data': cna_profile}
    else:
        tree = Phylo.read(newick, 'newick')
        calc_der(tree, node_data, output_dir, prefix)
        all_node_data = pd.read_table(os.path.join(output_dir, prefix + 'all_node_data.txt'), sep='\t', index_col=0)
        tree = Phylo.read(os.path.join(output_dir, prefix + 'cell_tree.newick'), 'newick')
        sctc_obj = {'tree': tree, 'cnv_data': all_node_data}

    res = define_CNA_mechnism(sctc_obj, cancer_type=cancer_type, cores=cores, random_num=randome_num)

    res['mode'].to_csv(f"{output_dir}/{prefix}mode.txt")
    res['orig_res']['rearrange_score']['all_rearrange_score'].to_csv(f"{output_dir}/{prefix}re_score.txt")
    res['orig_res']['rearrange_score']['all_chromothripsis_prop'].to_csv(
        f"{output_dir}/{prefix}chromothripsis_score.txt")
    res['orig_res']['rearrange_score']['all_rearrange_score_pvalue'].to_csv(f"{output_dir}/{prefix}re_score_pvalue.txt")

    if verbose:
        print(f'The resulting file is saved in the {os.path.join(output_dir, prefix)}*:')
        print(f'\t *Detail in : https://github.com/FangWang-SYSU/ChromoPhyloGen')
        print(f'Finished!')
