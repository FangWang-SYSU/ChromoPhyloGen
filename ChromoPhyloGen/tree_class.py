#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@File   : tree_class.py
@Author : XinWang
"""

import numpy as np
import pandas as pd
from .distance import dist

class TreeClass(object):
    def __init__(self, edges=None, node_data=None):
        '''
        :param edges: two columns ['p', 'c']
        :param node_data: pd.DataFrame, rowname is node name
        '''
        self.edges = pd.DataFrame(edges, columns=['p', 'c'])
        self.node_data = node_data
        self.root_length = 0

    def get_leaves(self, name):
        res = []
        q = [name]
        while q:
            node = q.pop(0)
            if node in self.edges.p.values:
                tmp = self.edges.loc[self.edges['p'] == node, :]
                cns = tmp['c'].tolist()
                q.extend(cns)
            else:
                res.append(node)
        return res

    def parent_cn_optim(self, cn1, cn2):
        c1 = self.node_data.loc[cn1, :]
        c2 = self.node_data.loc[cn2, :]
        return np.array(np.round((np.array(c1) + np.array(c2)) / 2))

    def get_all_children(self, name):
        res = []
        q = [name]
        while q:
            node = q.pop(0)
            if node in self.edges.p.values:
                tmp = self.edges.loc[self.edges.p == node, :]
                cns = tmp.c.tolist()
                q.extend(cns)
                res.extend(cns)
        return res

    def get_children(self, name):
        tmp = self.edges.loc[self.edges.p == name, :]
        cns = tmp.c.tolist()
        return cns

    def get_ancestor(self, name):
        tmp = self.edges.loc[self.edges.c == name, :]
        cns = tmp.p.tolist()
        if len(cns) == 0:
            return None
        return cns[0]

    def get_data(self, name):
        # print(name)
        return self.node_data.loc[name, :].tolist()

    def get_length(self, name):
        if name=='root':
            return self.root_length
        return float(self.edges.loc[self.edges.c==name, 'length'])
    # def walk_tree(self):
    #     root = list(set(self.edges.p) - set(self.edges.c))[0]
    #

    def unique_name(self):
        parent_nodes = set(self.edges.p)
        new_name = [f'virtual_{i}' for i in range(len(parent_nodes))]
        name_dict = dict(zip(parent_nodes, new_name))
        self.name_dict = name_dict
        root = list(set(self.edges.p) - set(self.edges.c))[0]
        name_dict[root] = 'root'

        self.edges.p = [name_dict[i] for i in self.edges.p]
        self.edges.c = [name_dict[i] if i in name_dict else i for i in self.edges.c]

    def update_tree(self, node='root'):
        childs = self.get_children(node)
        if len(childs) == 0:
            return self.node_data.loc[node, :].tolist()

        c1 = self.update_tree(childs[0])
        c2 = self.update_tree(childs[1])
        p = np.round((np.array(c1) + np.array(c2)) / 2).tolist()
        self.node_data.loc[node, :] = p

        l1, l2 = dist([p], [c1]), dist([p], [c2])
        self.edges.loc[(self.edges.p==node) & (self.edges.c==childs[0]), 'length'] = l1
        self.edges.loc[(self.edges.p==node) & (self.edges.c==childs[1]), 'length'] = l2
        return p

    def modify_length(self):
        max_l = max(set(self.edges.length) - set((1000000000,)))
        self.edges.length += 1
        self.edges.loc[self.edges.length==1000000000, 'length'] = max_l+1

        rcn = self.node_data.loc['root', :].tolist()
        l = dist([[2] * len(rcn)], [rcn])
        self.root_length = l+1


    def newick(self, node='root'):
        if node == 'root':
            try:
                l = f":{self.root_length}"
            except:
                l = ''
        else:
            try:
                l = f":{self.edges.loc[self.edges.c == node, 'length'].to_list()[0]}"
            except:
                l = ''

        childs = self.get_children(node)
        if len(childs) == 0:
            node_str = f'{node}{l}'
            return node_str

        descendants = ','.join([self.newick(n) for n in childs])
        descendants = f'({descendants}){node}{l}'
        return descendants
