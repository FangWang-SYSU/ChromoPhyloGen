#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@File   : keyCNA.py
@Author : XinWang
"""

import numpy as np
import pandas as pd


def split_status(pnode, cns0, cns1, tree_obj):
    cns0_data = tree_obj.get_data(cns0)
    cns1_data = tree_obj.get_data(cns1)
    pnod_data = tree_obj.get_data(pnode)

    c0_data, c1_data, p_data = np.array(cns0_data), np.array(cns1_data), np.array(pnod_data)
    _dd = ((c0_data + c1_data) % 2) != 0
    _anu = (((c0_data + c1_data) % 2) == 0) & (c0_data != c1_data)
    sstime = (tree_obj.get_length(cns0) + tree_obj.get_length(cns1)) / 2

    return np.sum(_dd), np.sum(_anu), sstime


def _dd_cn(rdata, cdata):
    _dd = ((cdata + rdata) % 2) != 0
    _dd_gain = _dd & (cdata > rdata)
    _dd_loss = _dd & (cdata < rdata)

    dd_vary_gain_loc, dd_vary_loss_loc = np.sum(_dd_gain), np.sum(_dd_loss)
    dd_cusum_gain_cn = np.sum(np.abs(rdata[_dd_gain] - cdata[_dd_gain]))
    dd_cusum_loss_cn = np.sum(np.abs(rdata[_dd_loss] - cdata[_dd_loss]))

    return dd_vary_gain_loc, dd_vary_loss_loc, dd_cusum_gain_cn, dd_cusum_loss_cn


def cell_status(node, tree_obj):
    node_data = tree_obj.get_data(node)
    root_data = tree_obj.get_data('root')
    ance_data = tree_obj.get_data(tree_obj.get_ancestor(node))

    c_data, r_data, p_data = np.array(node_data), np.array(root_data), np.array(ance_data)
    _sd = c_data == p_data

    # root
    Root_loss_cn, Root_gain_cn, Root_gain_loc, Root_loss_loc = _dd_cn(r_data, c_data)
    Parent_gain_cn, Parent_loss_cn, Parent_gain_loc, Parent_loss_loc = _dd_cn(p_data, c_data)

    return Root_gain_loc, Root_loss_loc, Root_gain_cn, Root_loss_cn, \
           Parent_gain_loc, Parent_loss_loc, Parent_gain_cn, Parent_loss_cn, \
           np.sum(_sd)


def cell_events(tree_obj, copy_rate=None, aneu_rate=None):
    ps = []
    node_length = {}
    q = ['root']
    while q:
        node = q.pop(0)
        node_len = tree_obj.get_length(node)
        if node == 'root':
            cs = [np.nan] * 9
            nl = 0
        else:
            cs = cell_status(node, tree_obj)
            nl = node_length[tree_obj.get_ancestor(node)] + node_len
        node_length[node] = nl
        cns = tree_obj.get_children(node)
        if cns.__len__() != 0:
            cns0, cns1 = cns
            ss = split_status(node, cns0, cns1, tree_obj)
            ps.append([node, *cs, *ss, nl, node_len])
            q.extend([cns0, cns1])
        else:
            ps.append([node, *cs, np.nan, np.nan, np.nan, nl, node_len])

    ps = pd.DataFrame(ps, columns=['name',
                                   'Root_gain_loc', 'Root_loss_loc', 'Root_gain_cn', 'Root_loss_cn',
                                   'Parent_gain_loc', 'Parent_loss_loc', 'Parent_gain_cn', 'Parent_loss_cn',
                                   'Mitosis_copy',
                                   'Mitosis_dd_loc', 'Mitosis_ad_loc', 'Mitosis_time',
                                   'Pseudotime_tree', 'Mitosis_time_next'])

    ps['aneu_rate'] = ps['Mitosis_ad_loc'] / tree_obj.node_data.shape[1]
    ps['copy_rate'] = ps['Mitosis_copy'] / tree_obj.node_data.shape[1]

    if copy_rate is None:
        copy_rate = max(ps.copy_rate.quantile(0.8), 0.5)
    if aneu_rate is None:
        aneu_rate = max(ps.aneu_rate.quantile(0.25), 0.01)
    # define clone
    ps['status'] = 'other'
    ps.loc[ps.aneu_rate >= aneu_rate, 'status'] = 'Aneuploidy'
    ps.loc[ps.copy_rate >= copy_rate, 'status'] = 'Copy'
    return ps


def key_events_graph(tree_obj, copy_rate=None, aneu_rate=None, min_cell_num=None):
    cell_relation = cell_events(tree_obj)

    if copy_rate is None:
        copy_rate = max(cell_relation.copy_rate.quantile(0.8), 0.5)
    if aneu_rate is None:
        aneu_rate = max(cell_relation.aneu_rate.quantile(0.25), 0.01)
    if min_cell_num is None:
        cell_num = len(tree_obj.get_leaves('root'))
        if cell_num <= 10:
            min_cell_num = 1
        else:
            min_cell_num = np.round(cell_relation/5)

    # define clone
    cell_relation['status'] = 'other'
    cell_relation.loc[cell_relation.aneu_rate>=aneu_rate, 'status'] = 'Aneuploidy'
    cell_relation.loc[cell_relation.copy_rate>=copy_rate, 'status'] = 'Copy'

    q = ['root']
    cell_clone = {'root': 'root'}
    clone_num = 0
    while q:
        node = q.pop(0)
        cns = tree_obj.get_children(node)
        if len(cns) != 0:
            for c in cns:
                if len(tree_obj.get_leaves(c)) <= min_cell_num:
                    cell_clone[c] = f'clone_{clone_num}'
                    for k in tree_obj.get_all_children(c):
                        cell_clone[k] = cell_clone[c]
                    clone_num += 1
                else:
                    if cell_relation.loc[cell_relation.name == c, 'status'].values[0] == 'Copy':
                       cell_clone[c] = cell_clone[node]
                    else:
                        cell_clone[c] = f'clone_{clone_num}'
                        clone_num += 1
                    q.append(c)

    clone_tree = tree_obj.edges[['p', 'c']]
    clone_tree.p = clone_tree.p.map(cell_clone)
    clone_tree.c = clone_tree.c.map(cell_clone)
    clone_tree = clone_tree.loc[clone_tree.p!=clone_tree.c, :]
    clone_tree = clone_tree.drop_duplicates()

    return clone_tree, cell_clone, cell_relation
