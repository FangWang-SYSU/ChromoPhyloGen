#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@File   : distance.py
@Author : XinWang
"""

import numpy as np
import copy

def matrixbuilder(node):
    matrix = []
    for node1 in node:
        temp = []
        for node2 in node:
            temp.append(dist(node[node1], node[node2]))
        matrix.append(temp)
    return node.keys(), matrix


def dist(node1, node2):
    d = 0
    for i in range(0, len(node1)):
        d = d + disthelper(node1[i], node2[i])
    return d

def disthelper(node1, node2):
    if 0 in node1 or 0 in node2:
        return zerodisthelper(node1, node2)
    return distcalc(node1, node2)


def distcalc(node1, node2):
    assert len(node1) == len(node2)
    if len(node1) == 1:
        return abs(node1[0] - node2[0])
    else:
        d = 0
        newlist = copy.deepcopy(node1)
        for i in range(0, len(node2)):
            newlist[i] -= node2[i]
        while newlist:
            if newlist[0] == 0:
                newlist.pop(0)
            elif newlist[0] > 0:
                k = 0
                for i in range(0, len(newlist)):
                    if newlist[i] > 0:
                        k = i
                    else:
                        break
                for i in range(0, k + 1):
                    newlist[i] -= 1
                d += 1
            elif newlist[0] < 0:
                k = 0
                for i in range(0, len(newlist)):
                    if newlist[i] < 0:
                        k = i
                    else:
                        break
                for i in range(0, k + 1):
                    newlist[i] += 1
                d += 1
        return abs(d)


def zerodisthelper(node1, node2):
    n1 = copy.deepcopy(node1)
    n2 = copy.deepcopy(node2)
    dist = 0
    temp1 = []
    temp2 = []
    while n1:
        x1 = n1.pop()
        x2 = n2.pop()
        if x1 == 0:
            if x2 == 0:
                temp1.append(x1)
                temp2.append(x2)
            else:
                temp1.append(x1)
                temp2.append(0)
                # return 1000000*1000
        else:
            temp1.append(x1)
            temp2.append(x2)
    return distcalc(temp1, temp2)


def split_array(a, b, c, aneu=False):
    _dd1 = np.array(((b + c) % 2) != 0)#~(np.array(((a + b) % 2) == 0) & np.array(a != b))
    _dd2 = ((b+c) / 2) != a
    if aneu:
        _dd = _dd1 | _dd2
    else:
        _dd = _dd1
    d = np.nonzero(_dd != np.roll(_dd, 1))[0]
    if len(d) == 0:
        return [a.tolist()], [b.tolist()], [c.tolist()]
    if d[0] == 0:
        d = d[1:]

    list0 = np.split(a, d)
    list0 = list0[0::2] if _dd[0] else list0[1::2]

    list1 = np.split(b, d)
    list1 = list1[0::2] if _dd[0] else list1[1::2]

    list2 = np.split(c, d)
    list2 = list2[0::2] if _dd[0] else list2[1::2]

    list0 = [i.tolist() for i in list0]
    list1 = [i.tolist() for i in list1]
    list2 = [i.tolist() for i in list2]
    return list0, list1, list2


def updata_ntree_dmed(node, root_data=None, aneu=False):
    # dist([[2, 2, 2, 2, 2, 2, 2, 2]],
    #      [[3, 3, 3, 2, 2, 0, 1, 1]])
    if node.is_leaf:
        return node.data
    cns0, cns1 = node.descendants
    cns0_data = updata_ntree_dmed(cns0)
    cns1_data = updata_ntree_dmed(cns1)
    cns0.comment = '!'.join(list(map(lambda x: str(int(x)), cns0_data)))
    cns1.comment = '!'.join(list(map(lambda x: str(int(x)), cns1_data)))
    if not node.ancestor and root_data is not None:
        node_data = root_data
    else:
        node_data = np.round((np.array(cns0_data) + np.array(cns1_data)) / 2).tolist()
    node.comment = '!'.join(list(map(lambda x: str(int(x)), node_data)))

    a, b, c = np.array(node.data), np.array(cns0.data), np.array(cns1.data)
    list_a0, list_b, list_c = split_array(a, b, c, aneu=aneu)
    # list_a1, list_c = split_array(a, c)
    l0, l1 = dist(list_a0, list_b), dist(list_a0, list_c)
    cns0.length = str(l0)
    cns1.length = str(l1)
    return node_data