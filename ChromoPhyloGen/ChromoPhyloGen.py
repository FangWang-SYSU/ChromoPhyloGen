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
from .utils import plot_cn_tree
from .infer_tree import InferTree
from .keyCNA import cell_events
from .events_define import define_CNA_mechnism
from .estimateDER import estimate_der
from Bio import Phylo


def run(cna_dir,
        output_dir,
        prefix='ChromoPhyloGen',
        resolution=1,
        clone_thr=0.9,
        n_neighbors=5,
        min_clone_size=10,
        random_state=1234,
        plot_png=False,
        verbose=True):
    # 1
    if verbose:
        print('Load CNA profile file.')
    if type(cna_dir) is str:
        cna_profile = pd.read_table(cna_dir)
        # cna_profile.index = cna_profile['chr'].map(str) + '_' + cna_profile['start'].map(str) + '_' + cna_profile['end'].map(str)
        cna_profile.index = cna_profile.iloc[:, 0].map(str) + '_' + cna_profile.iloc[:, 1].map(
            str) + '_' + cna_profile.iloc[:, 2].map(str)
        # del cna_profile['chr']
        # del cna_profile['start']
        # del cna_profile['end']
        cna_profile = cna_profile.iloc[:,3::]
        cna_profile = cna_profile.astype('int')
        cna_profile = cna_profile.transpose()
    else:
        cna_profile = cna_dir

    # 2
    if verbose:
        print(f'There are {cna_profile.shape[0]} cells and {cna_profile.shape[1]} position for inference')
        print(f'Inferring...')

    if min_clone_size is None:
        min_clone_size = 0.1 * cna_profile.shape[0]
    ht = InferTree(cna_profile,
                   resolution=resolution,
                   clone_thr=clone_thr,
                   n_neighbors=n_neighbors,
                   min_clone_size=min_clone_size,
                   random_state=random_state)
    tree_obj = ht.run()
    if verbose:
        print('ChromoPhyloGen inferrence finish.')
        print('Saving file.')

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    if plot_png:
        plot_cn_tree(tree_obj.newick(), cna_profile, title_label=prefix,
                     save=os.path.join(output_dir, prefix + 'tree.png'))

    cell_relation = cell_events(tree_obj)
    cell_relation.to_csv(os.path.join(output_dir, prefix + 'cell_info.txt'), sep='\t', index=False)
    tree_obj.node_data.to_csv(os.path.join(output_dir, prefix + 'all_node_data.txt'), sep='\t', index=True)
    with open(os.path.join(output_dir, prefix + 'cell_tree.newick'), 'w') as f:
        f.write('(' + tree_obj.newick() + ');')

    der = estimate_der(tree_obj, cell_relation)
    der.to_csv(os.path.join(output_dir, prefix + 'error_risk_score.txt'), sep='\t', index=True, header=False)

def chromosome_event(output_dir,
                     prefix='ChromoPhyloGen',
                     cancer_type='ALL',
                     cores=1,
                     randome_num=1000,
                     verbose=True):
    if verbose:
        print('Scoring the chromosomal rearrangements.')
    tree = Phylo.read(os.path.join(output_dir, prefix + 'cell_tree.newick'), 'newick')
    all_node_data = pd.read_table(os.path.join(output_dir, prefix + 'all_node_data.txt'), sep='\t', index_col=0)
    sctc_obj = {'tree': tree, 'cnv_data': all_node_data}
    res = define_CNA_mechnism(sctc_obj, cancer_type=cancer_type, cores=cores, random_num=randome_num)

    res['mode'].to_csv(f"{output_dir}/{prefix}mode.txt")
    res['orig_res']['rearrange_score']['all_rearrange_score'].to_csv(f"{output_dir}/{prefix}re_score.txt")
    res['orig_res']['rearrange_score']['all_chromothripsis_prop'].to_csv(
        f"{output_dir}/{prefix}chromothripsis_score.txt")
    res['orig_res']['rearrange_score']['all_rearrange_score_pvalue'].to_csv(f"{output_dir}/{prefix}re_score_pvalue.txt")

    if verbose:
        print(f'The resulting file is saved in the {os.path.join(output_dir, prefix)}*:')
        print(f'\t 1.cell_info.txt: Cell variation information in trace')
        print(f'\t 2.all_node_data.txt: Cell CNA profile, including internal node, name by "virtual_"')
        print(f'\t 3.cell_tree.newick: Single cell trace file，format newick')
        print(f'\t Other: *Detail in : https://github.com/FangWang-SYSU/ChromoPhyloGen')
        print(f'Finished!')
