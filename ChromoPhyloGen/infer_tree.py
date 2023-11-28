#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@File   : infer_tree.py
@Author : XinWang
"""

import pandas as pd
import snf
from sklearn.cluster import spectral_clustering
from .distance import aneuploidy_dist_optim
from .huffmanTree import huffman_tree_optim
from .tree_class import TreeClass


class InferTree(object):
    def __init__(self,
                 profile,
                 resolution=1,
                 n_neighbors=None,
                 clone_thr=0.9,
                 clone_dist='2',
                 **kwargs):
        self.profile = profile
        self.n_neighbors = n_neighbors
        self.clone_thr = clone_thr
        self.resolution = resolution
        self.min_clone_size = kwargs.get('min_clone_size', 30)  # clone min size
        self.random_state = kwargs.get('random_state', 1994)
        self.cell_num = self.profile.shape[0]
        self.cell_name = self.profile.index

    def run(self):
        # 1. distance
        dist = aneuploidy_dist_optim(self.profile, self.profile)
        SD_dist = pd.DataFrame(dist[:, 2].reshape((self.cell_num, self.cell_num)), columns=self.cell_name,
                               index=self.cell_name)
        AD_dist = pd.DataFrame(dist[:, 3].reshape((self.cell_num, self.cell_num)), columns=self.cell_name,
                               index=self.cell_name)
        DD_dist = pd.DataFrame(dist[:, 4].reshape((self.cell_num, self.cell_num)), columns=self.cell_name,
                               index=self.cell_name)

        # 2. get clone
        self.clone_tree = self.get_clone(SD_dist, AD_dist, DD_dist)

        # print(self.clone_tree)
        # 3. create_tree
        merge_edges = self.create_tree()

        # 4. convert tree
        merge_edges1 = merge_edges.loc[merge_edges.p.isin(merge_edges.loc[merge_edges.p.duplicated(), 'p']), :]
        merge_edges2 = merge_edges.loc[~merge_edges.p.isin(merge_edges.loc[merge_edges.p.duplicated(), 'p']), :]
        for _, r in merge_edges2.iterrows():
            p, c = r['p'], r['c']
            merge_edges1.loc[merge_edges1.c == p, 'c'] = c

        # 5. build tree class
        tc = TreeClass(merge_edges1, self.profile)
        tc.unique_name()
        tc.update_tree()
        tc.modify_length()
        return tc

    def get_clone(self, SD, AD, DD):
        clone_tree = pd.DataFrame(index=self.cell_name)

        def clone(mat, n_neighbors=20, layer=0, min_clone_size=6):
            if mat.shape[0] <= min_clone_size:
                clone_tree.loc[mat.index, f'layer_{layer}'] = mat.index
                return
            affinity_mat = [
                SD.loc[mat.index, mat.index],
                AD.loc[mat.index, mat.index],
                DD.loc[mat.index, mat.index],
            ]
            affinity_mat[2] = snf.make_affinity(affinity_mat[2], metric='euclidean', K=n_neighbors, mu=0.5)
            #affinity_mat = [affinity_mat[0],affinity_mat[1]]
            fused_network = snf.snf(affinity_mat, K=n_neighbors)
            n_clusters = range(2, 6)
            if fused_network.shape[0]<6:
                n_clusters = range(1, fused_network.shape[0])
            best, _ = snf.get_n_clusters(fused_network, n_clusters=n_clusters)
            labels = spectral_clustering(fused_network, n_clusters=best, random_state=self.random_state)
            if f'layer_{layer}' in clone_tree.columns:
                labels += (clone_tree[f'layer_{layer}'].dropna().count() + 1)
            clone_tree.loc[mat.index, f'layer_{layer}'] = labels
            #
            mat['S_clone'] = labels
            for sc, clone_data in mat.groupby(['S_clone']):
                del clone_data['S_clone']
                if clone_data.shape[0] < 30:
                    nn = max(int(clone_data.shape[0] / 4), 2)
                else:
                    nn = 20
                clone(clone_data, n_neighbors=nn, layer=layer + 1, min_clone_size=min_clone_size)
            return
        if self.n_neighbors is not None and self.n_neighbors<self.profile.shape[0]:
            nn = self.n_neighbors
        elif self.profile.shape[0] < 30:
            nn = max(int(self.profile.shape[0] / 4), 2)
        else:
            nn = 20
        clone(self.profile.__deepcopy__(), n_neighbors=nn, layer=1, min_clone_size=self.min_clone_size)
        clone_tree['layer_0'] = 0
        return clone_tree

    def create_tree(self):
        self.ctc = {}
        self.merge_edges = pd.DataFrame(columns=['p', 'c'])

        def build_tree(self, ct, parent_layer=0, parent_clone=0, curr_layer=1):
            clone_roots_data = []
            tmp_clone_tree = ct.loc[ct[f'layer_{parent_layer}'] == parent_clone, :]
            tmp_layer_clones = list(set(tmp_clone_tree[f'layer_{curr_layer}']))
            if len(tmp_layer_clones) == tmp_clone_tree.shape[0]:
                clone_roots_data = self.profile.loc[tmp_layer_clones, :]
                cell_Si = pd.Series(range(len(tmp_layer_clones)), tmp_layer_clones)
            else:
                for tmp_layer in tmp_layer_clones:
                    new_clone_data = build_tree(self, tmp_clone_tree, parent_layer=curr_layer, parent_clone=tmp_layer,
                                                curr_layer=curr_layer + 1)
                    clone_roots_data.append(new_clone_data)
                    # clone_root_names.append(f'layer_{parent_layer}_{tmp_layer}')
                clone_roots_data = pd.DataFrame(clone_roots_data, index=tmp_layer_clones)
                cell_Si = pd.Series(range(len(tmp_layer_clones)), tmp_layer_clones)

            # clone_roots_data.index = range(len(tmp_layer_clones))
            hft = huffman_tree_optim(clone_roots_data.astype('int'),
                                     all_dist=None,
                                     cell_pos_map=cell_Si,
                                     init_virtual_num=None,
                                     clone_thr=self.clone_thr)
            dt = {int(v): k for k, v in dict(cell_Si).items()}
            new_edges = []
            leaves = hft.get_leaves()
            for _, r in hft.edges.iterrows():
                p, c = r['p'], r['c']
                if c == p and hft.edges.shape[0] == 1:
                    c = dt[int(c)]
                    p = f'layer_{parent_layer}_{int(parent_clone)}'
                else:
                    if c in leaves:
                        if type(dt[int(c)]) is str:
                            c = dt[int(c)]
                        else:
                            c = f'layer_{curr_layer}_{int(dt[int(c)])}'
                    else:
                        c = f'layer_{parent_layer}_{int(parent_clone)}_{int(c)}'
                    if p == hft.tree_roots:
                        p = f'layer_{parent_layer}_{int(parent_clone)}'
                    else:
                        p = f'layer_{parent_layer}_{int(parent_clone)}_{int(p)}'
                new_edges.append([p, c])
            # self.ctc[f'layer_{parent_layer}_{int(parent_clone)}'] = hft
            self.merge_edges = pd.concat([self.merge_edges, pd.DataFrame(new_edges, columns=['p', 'c'])], axis=0)

            return hft.node_data.get(hft.tree_roots)

        build_tree(self, self.clone_tree, parent_layer=0, parent_clone=0, curr_layer=1)
        return self.merge_edges
