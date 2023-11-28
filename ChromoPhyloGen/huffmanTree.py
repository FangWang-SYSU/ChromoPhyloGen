#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@File   : huffmanTree.py
@Author : XinWang
"""

import pandas as pd
import numpy as np
import fast_dist
from operator import itemgetter

def add_labels_optim(dist, loci_num, leaves_num, clone_thr=0.9):
    # 0:c1, 1:c2, 2:SD, 3:AD, 4:DD, 5:is_clone 6: is_leaf, 7:AD_SD
    # digitized = np.digitize(dist[:,2], bins)
    dist = np.insert(dist, 5, [2] * dist.shape[0], axis=1)
    dist[dist[:, 2] >= (clone_thr * loci_num), 5] = 0
    dist[(dist[:, 2] < (clone_thr * loci_num)) & ((dist[:, 2] + dist[:, 3]) >= (clone_thr * loci_num)), 5] = 1
    # leaves first
    dist = np.insert(dist, 6, [2] * dist.shape[0], axis=1)
    dist[(dist[:, 0] < leaves_num) | (dist[:, 1] < leaves_num), 6] = 1
    dist[(dist[:, 0] < leaves_num) & (dist[:, 1] < leaves_num), 6] = 0

    dist = np.insert(dist, 7, dist[:, 2] + dist[:, 3], axis=1)
    # dist = np.insert(dist, 7, (dist[:, 2]) // 10, axis=1)
    # dd = loci_num - dist[:, 2] - dist[:, 3]
    # dist = np.insert(dist, 7, dd * dist[:, 3] / (dist[:, 2] + dist[:, 3] + 1e-6) / (dist[:, 3] + 1e-6), axis=1)
    return dist


class TreeObj(object):
    def __init__(self, tree_roots=None, data=None):
        self.tree_roots = tree_roots
        self.edges = pd.DataFrame(columns=['p', 'c'])
        self.node_data = {}
        if data is not None:
            self.node_data = {tree_roots: data}

    def get_leaves(self):
        return list(set(self.edges.c) - set(self.edges.p))

    def parent_cn_optim(self, cn1, cn2):
        c1 = self.node_data.get(cn1)
        c2 = self.node_data.get(cn2)
        return np.array(np.round((np.array(c1) + np.array(c2)) / 2))

    def all_nodes(self):
        if self.edges.shape[0] <= 1:
            return list(self.node_data.keys())
        return self.edges.values.ravel().tolist()



def get_other_set_node(curr, tree_set):
    res = []
    for k, v in tree_set.items():
        a = v.all_nodes()
        if curr not in a:
            res.extend(a)
    return list(set(res))


def get_tree_by_node(name, tree_set):
    for k, v in tree_set.items():
        if name in v.all_nodes():
            return v
    return None


def get_dist_pair(tree_set):
    res = pd.DataFrame(columns=['p', 'c'])
    for k, v in tree_set.items():
        other_set_node = get_other_set_node(k, tree_set)
        res = pd.concat([res, pd.DataFrame(zip([k] * len(other_set_node), other_set_node), columns=['p', 'c'])], axis=0)
    return res


def get_tree_sets_data(tree_set):
    res = {}
    for k, v in tree_set.items():
        res.update(v.node_data)
    return res


def _mergeTree(new_name, obj1, obj2):
    new_obj = TreeObj(new_name)
    new_obj.edges = pd.concat([obj1.edges, obj2.edges])
    new_edge = [[new_name, obj1.tree_roots], [new_name, obj2.tree_roots]]
    new_obj.edges = pd.concat([new_obj.edges, pd.DataFrame(new_edge, columns=['p', 'c'])], axis=0)


    new_obj.node_data[new_name] = np.array(np.round((np.array(obj1.node_data[obj1.tree_roots]) +
                                                     np.array(obj2.node_data[obj2.tree_roots])) / 2))
    new_obj.node_data.update(obj1.node_data)
    new_obj.node_data.update(obj2.node_data)
    return new_obj


def mergeTreeObj(p, k_node1, k_node2, tree_set):
    new_tree = _mergeTree(p, tree_set[k_node1], tree_set[k_node2])
    tree_set.pop(k_node1)
    tree_set.pop(k_node2)
    tree_set[p] = new_tree
    return tree_set


def update_node_data(node, obj):
    if node not in obj.edges.p.values:
        return obj.node_data[node], obj
    c1, c2 = obj.edges.loc[obj.edges.p == node, 'c'].values.ravel().tolist()
    c1_d, obj = update_node_data(c1, obj)
    c2_d, obj = update_node_data(c2, obj)
    obj.node_data[node] = np.array(np.round((np.array(c1_d) + np.array(c2_d)) / 2))
    return obj.node_data[node], obj


def addTreeObjTo(new_name, k_node1, k2_tree, k_node2, tree_sets, update=False):
    k1_tree = tree_sets[k_node1]
    new_obj = TreeObj(k2_tree.tree_roots)
    k2_edges = k2_tree.edges
    k2_edges.loc[k2_edges.c == k_node2, 'c'] = new_name

    new_obj.edges = pd.concat([k1_tree.edges, k2_edges])
    new_edge = [[new_name, k_node2], [new_name, k_node1]]
    new_obj.edges = pd.concat([new_obj.edges, pd.DataFrame(new_edge, columns=['p', 'c'])], axis=0)

    new_obj.node_data[new_name] = np.array(np.round((np.array(k1_tree.node_data[k_node1]) +
                                                     np.array(k2_tree.node_data[k_node2])) / 2))
    new_obj.node_data.update(k1_tree.node_data)
    new_obj.node_data.update(k2_tree.node_data)

    if update:
        _, new_obj = update_node_data(k2_tree.tree_roots, new_obj)

    tree_sets[k2_tree.tree_roots] = new_obj
    tree_sets.pop(k_node1)
    return tree_sets


def drop_parent_nodes(name, edges):
    res = []
    q = [name]
    while q:
        node = q.pop(0)
        if node in edges.c.tolist():
            p, c = edges.loc[edges.c == node, :].iloc[0, :]
            res.append(p)
            q.append(p)
    edges = edges.loc[~edges.p.isin(res), :]
    return edges, res


def subset_TreeObj(name, obj):
    edges = obj.edges
    res = [name]
    q = [name]
    while q:
        node = q.pop(0)
        if node in edges.p.tolist():
            c1, c2 = edges.loc[edges.p == node, 'c'].tolist()
            res.extend([c1, c2])
            q.extend([c1, c2])

    new_obj = TreeObj(name)
    new_obj.edges = edges.loc[edges.p.isin(res), :]
    new_obj.node_data = {k: v for k, v in obj.node_data.items() if k in res}
    return new_obj


def splitTreeObj(k2_tree, k_node2, tree_sets):
    new_edges, drop_res = drop_parent_nodes(k_node2, k2_tree.edges)
    all_node = set(k2_tree.edges.values.ravel())
    new_roots = list(set(new_edges.p) - set(new_edges.c))
    leaves_roots = list(all_node - set(new_edges.values.ravel()) - set(drop_res))
    new_roots += leaves_roots
    for i in new_roots:
        tree_sets[i] = subset_TreeObj(i, k2_tree)
    tree_sets.pop(k2_tree.tree_roots)
    return tree_sets, drop_res


def huffman_tree_optim(profile, init_virtual_num=None, all_dist=None, cell_pos_map=None, clone_thr=0.9):
    # print(profile)
    tree_sets = {cell_pos_map[i]: TreeObj(cell_pos_map[i], profile.loc[i, :].tolist()) for i in profile.index}
    leaves_node = cell_pos_map.tolist()
    if init_virtual_num is None:
        virtual_num = max(cell_pos_map.values) + 1
    else:
        virtual_num = init_virtual_num

    cell_num, seg_num = profile.shape
    all_cell_num = max(cell_pos_map.values) + 1

    if len(tree_sets) == 1:
        node_name = cell_pos_map[profile.index[0]]
        tmp_tree = tree_sets[node_name]
        tmp_tree.edges = pd.DataFrame([[node_name, node_name]], columns=['p', 'c'])
        return tmp_tree

    if all_dist is None:
        all_dist = fast_dist.dist(np.array(profile, dtype=float),
                                  np.array(profile, dtype=float),
                                  np.array(cell_pos_map[profile.index], dtype=float).reshape((1, -1)),
                                  np.array(cell_pos_map[profile.index], dtype=float).reshape((1, -1)))
        all_dist = all_dist[all_dist[:, 0] != all_dist[:, 1], :]
    # all_dist = np.insert(all_dist, 7, all_dist[:, 2] + all_dist[:, 3], axis=1)
    # all_dist = np.insert(all_dist, 8, np.mean(all_dist[:, 5:7], axis=1), axis=1)
    all_dist = add_labels_optim(all_dist,
                                seg_num,
                                all_cell_num,
                                clone_thr=clone_thr)
    dist_num = 2
    # clone_thr = clone_thr * seg_num
    old_dist_dict = {}
    fix_nodes = []
    split_num = {}
    while True:
        all_dist = all_dist[np.lexsort((
            all_dist[:, 4],  # DD
            -all_dist[:, 2],  # SD
            -all_dist[:, 7],  # AD_SD
            all_dist[:, 6],  # is_leaf
        ))]
        # select min
        all_dist = all_dist[all_dist[:, 0] != all_dist[:, 1], :]
        k_node1, k_node2, new_d = all_dist[0, [0, 1, dist_num]]
        all_dist = all_dist[1:, :]
        all_dist = all_dist[~((all_dist[:, 0] == k_node2) & (all_dist[:, 1] == k_node1)), :]
        if k_node2 in tree_sets:
            all_dist = all_dist[all_dist[:, 0] != k_node1, :]
            all_dist = all_dist[all_dist[:, 0] != k_node2, :]
            tree_sets = mergeTreeObj(virtual_num, k_node1, k_node2, tree_sets)
            if k_node2 in leaves_node:
                fix_nodes.extend([k_node1, k_node2])
            old_dist_dict[k_node1] = old_dist_dict[k_node2] = new_d
            # add
            other_set_node = get_other_set_node(virtual_num, tree_sets)
            new_dist_pair = pd.DataFrame(zip([virtual_num] * len(other_set_node), other_set_node), columns=['p', 'c'])
        else:
            k2_tree = get_tree_by_node(k_node2, tree_sets)
            if new_d >= (clone_thr * seg_num):
                old_roots = k2_tree.tree_roots
                k1_tree_nodes = tree_sets[k_node1].all_nodes()
                tree_sets = addTreeObjTo(virtual_num, k_node1, k2_tree, k_node2, tree_sets, update=False)
                fix_nodes.extend([k_node1, k_node2])
                old_dist_dict[virtual_num] = old_dist_dict[k_node2]
                old_dist_dict[k_node1] = old_dist_dict[k_node2] = new_d
                # drop
                all_dist = all_dist[all_dist[:, 0] != k_node1, :]
                all_dist = all_dist[~((all_dist[:, 0] == old_roots) &
                                    np.isin(all_dist[:, 1], k1_tree_nodes)), :]
                # add
                other_root = list(tree_sets.keys())
                other_root.remove(old_roots)
                new_dist_pair = pd.DataFrame(zip(other_root, [virtual_num] * len(other_root)), columns=['p', 'c'])
            else:
                if k_node2 in fix_nodes:
                    continue
                try:
                    k2_parent = k2_tree.edges.loc[k2_tree.edges.c == k_node2, 'p'].to_list()
                    k2_brother = k2_tree.edges.loc[k2_tree.edges.p == k2_parent[0], 'c'].tolist()
                    k2_brother.remove(k_node2)
                    brother_new_d = all_dist[(all_dist[:, 0] == k_node1) & (all_dist[:, 1] == k2_brother[0]), dist_num]
                    new_d = min(brother_new_d[0], new_d)
                except:
                    new_d = new_d
                if new_d <= old_dist_dict[k_node2]:
                    continue

                if split_num.get(k_node2, 0) >= 2:
                    continue
                split_num[k_node2] = split_num.get(k_node2, 0) + 1
                old_roots = list(tree_sets.keys())  # + [k_node2]
                tree_sets, drop_res = splitTreeObj(k2_tree, k_node2, tree_sets)
                for ks in tree_sets.keys():
                    old_dist_dict[ks] = 0
                # drop
                all_dist = all_dist[~np.isin(all_dist[:, 1], drop_res), :]
                all_dist = all_dist[~np.isin(all_dist[:, 0], drop_res), :]
                new_roots = list(set(tree_sets.keys()) - set(old_roots))

                all_dist = all_dist[all_dist[:, 0] != k_node1, :]
                all_dist = all_dist[all_dist[:, 0] != k_node2, :]
                # merge
                tree_sets = mergeTreeObj(virtual_num, k_node1, k_node2, tree_sets)
                old_dist_dict[k_node1] = old_dist_dict[k_node2] = new_d
                # fix_nodes.extend([k_node1, k_node2])
                # add
                new_roots.remove(k_node2)
                other_roots = new_roots + [virtual_num]
                new_dist_pair = pd.DataFrame(columns=['p', 'c'])
                for tmp_i in other_roots:
                    other_set_node = get_other_set_node(tmp_i, tree_sets)
                    new_dist_pair = pd.concat(
                        [new_dist_pair,
                         pd.DataFrame(zip([tmp_i] * len(other_set_node), other_set_node), columns=['p', 'c'])], axis=0)
        #print(len(tree_sets))
        # if len(tree_sets) < 17:
        #     print(111)
        if len(tree_sets) == 1:
            break

        new_dist_pair = new_dist_pair.sort_values('p', ascending=False)
        new_dist_pair = new_dist_pair.loc[~new_dist_pair.apply(lambda x: sorted(-x), 1).duplicated(), :]
        tree_sets_all_data = get_tree_sets_data(tree_sets)
        non_dist = fast_dist.dist2(np.array(itemgetter(*new_dist_pair.p)(tree_sets_all_data), dtype=float).reshape(
            (new_dist_pair.shape[0], -1)),
                                   np.array(itemgetter(*new_dist_pair.c)(tree_sets_all_data), dtype=float).reshape(
                                       (new_dist_pair.shape[0], -1)),
                                   np.array(new_dist_pair.p, dtype=float).reshape((1, -1)),
                                   np.array(new_dist_pair.c, dtype=float).reshape((1, -1))
                                   )
        # non_dist = np.insert(non_dist, 7, non_dist[:, 2] + non_dist[:, 3], axis=1)
        # non_dist = np.insert(non_dist, 8, non_dist[:, 5] * non_dist[:, 6], axis=1)
        non_dist = add_labels_optim(non_dist,
                                    seg_num,
                                    all_cell_num,
                                    clone_thr=clone_thr)

        all_dist = np.vstack((all_dist, non_dist))

        virtual_num += 1
    return tree_sets.popitem()[1]
