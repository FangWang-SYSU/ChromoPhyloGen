#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@File   : utils.py
@Author : XinWang
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from io import StringIO
from Bio import Phylo
import matplotlib.colors as mc
from scipy import signal



def np_fix(x, n=2):
    res = np.asanyarray(np.ceil(x))
    res = np.floor(x, out=res, where=np.greater_equal(x, n))
    #res = np.round(x)
    return res


def load_data(dir):
    if type(dir) == str:
        dataset = pd.read_csv(dir, sep='\t')
    else:
        dataset = dir
    # del dataset['chrompos']
    cnames = dataset.columns.tolist()
    cnames[:2] = ['chr', 'absposstart']
    dataset.columns = cnames
    _map = {'X': 23, 'Y': 24}
    dataset['chr'] = [i if i not in _map else _map[i] for i in dataset['chr']]
    return dataset


def arm_level(data, arm_info=None):
    if arm_info is None:
        arm_info = './data/hg19_chr_arm_info.txt'
    arm_info = pd.read_csv(arm_info, sep='\t')
    breaks = [0] + arm_info.abspos.tolist()
    labels = arm_info['chrom'].astype('str').str.cat(arm_info['arm'], sep='_').tolist()
    new_label = pd.cut(data.start,
           breaks,
           labels=labels)
    data['arm'] = new_label
    del data['chrom']
    del data['start']
    data = np_fix(data.groupby('arm', sort=False).agg('mean')).transpose()
    data = data.dropna(axis=1)
    return data


def step_k(data, k=10 * 1000 * 1000):
    data['group'] = data['absposstart'] // k
    data['group'] = data.apply(lambda x: str(int(x['chr'])) + '_' + str(int(x['group'])), axis=1)
    del data['absposstart']
    del data['chr']

    data = np_fix(data.groupby('group', sort=False).agg('mean')).transpose()
    data = data.dropna(axis=1)
    return data


def plot_cn_tree(newick, profile, title_label=None, save=False):
    if title_label is None:
        title_label = ''
    # print(tmp_cell_order)
    phylo_tree = Phylo.read(StringIO(newick), "newick")
    tmp_cell_order = [i.name for i in phylo_tree.get_terminals()]
    fig = plt.figure(figsize=(16, 8))
    grid = plt.GridSpec(1, 5, wspace=0, hspace=0.5)
    ax1 = plt.subplot(grid[0, 0:2], frame_on=False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    ax1.spines['bottom'].set_visible(False)
    ax1.axis('off')
    ax2 = plt.subplot(grid[0, 2:5])
    plt.title(title_label)
    colors = {'white': 0, '#91FF88': 1, '#C6C6C6': 2, '#FFEB97': 3, '#FCCC00': 4, '#ec9336': 5, '#7d170e': 6, 'black':7}
    # colors2 = ['gray', '#FF9900', 'blue', 'red', 'pink']
    if profile.shape[0] < 30:
        plot_data = profile.loc[tmp_cell_order, :]
    else:
        plot_data = np.array(profile.loc[tmp_cell_order, :])

    sns.heatmap(plot_data,
                cmap=mc.ListedColormap(colors),
                vmin=0, vmax=len(colors), ax=ax2,
                )
    Phylo.draw(phylo_tree,
               do_show=False,
               label_func=lambda x: '',
               branch_labels=lambda x: '',
               axes=ax1)
    if save:
        fig.savefig(save, dpi=300)
    else:
        plt.show()


def plot_cn(data):
    colors = {'white': 0, '#91FF88': 1, '#C6C6C6': 2, '#FFEB97': 3, '#FCCC00': 4, '#ec9336': 5, '#7d170e': 6, 'black':7}
    fig, ax = plt.subplots(1, 1, figsize=(12, 6))
    sns.heatmap(data,
                cmap=mc.ListedColormap(colors), vmin=0, vmax=len(colors),
                ax=ax,
                linewidths=0, linecolor='white', )
    plt.show()


def expand_data(profile, bw=None):
    '''
    :param profile: colnames = [chrom, start, end, cell1, cell2 ... celln], first column is chr, second is start, third is end.
    :param kernel_size: denoise kernel size
    :param bin_width: bin width, default is min segment length
    :return: smooth profile
    '''
    # split bin by chr

    colnames = profile.columns.to_list()
    colnames[:3] = ['chrom', 'start', 'end']
    profile.columns = colnames


    segment_length = profile['end'] - profile['start']
    bin_width = segment_length.min()

    if bw is None or bw > bin_width:
        bw = bin_width

    new_data = []
    new_pos = []
    for idx, sl in enumerate(segment_length):
        cut_num = sl // bw
        row = profile.iloc[idx, :].tolist()
        chrom, start, end = row[:3]
        tmp_cn = row[3:]
        sp_res = np.array_split(range(int(start), int(end+1)), cut_num)
        tmp_pos = [[chrom, sp[0], sp[-1]] for sp in sp_res]
        new_data.extend([tmp_cn]*len(sp_res))
        new_pos.extend(tmp_pos)
    expand_profile = pd.DataFrame(np.hstack([np.array(new_pos), np.array(new_data)]), columns=profile.columns)
    return expand_profile


def denoise(profile, kernel_size=5, bw=None, dup=0.1):
    '''
    :param profile: colnames = [chrom, start, end, cell1, cell2 ... celln], first column is chr, second is start, third is end.
    :param kernel_size: denoise kernel size
    :param bin_width: bin width, default is min segment length
    :return: smooth profile
    '''
    # split bin by chr

    colnames = profile.columns.to_list()
    colnames[:3] = ['chrom', 'start', 'end']
    profile.columns = colnames


    segment_length = profile['end'] - profile['start']

    if bw is None:
        bin_width = segment_length.min()
        bw = bin_width


    new_data = []
    new_pos = []
    for idx, sl in enumerate(segment_length):
        row = profile.iloc[idx, :].tolist()
        chrom, start, end = row[:3]
        tmp_cn = row[3:]

        if sl <= bw:
            new_data.extend([tmp_cn])
            new_pos.append([chrom, start, end])
        else:
            cut_num = sl // bw
            sp_res = np.array_split(range(int(start), int(end+1)), cut_num)
            tmp_pos = [[chrom, sp[0], sp[-1]] for sp in sp_res]
            new_data.extend([tmp_cn]*len(sp_res))
            new_pos.extend(tmp_pos)
    expand_profile = pd.DataFrame(np.hstack([np.array(new_pos), np.array(new_data)]), columns=profile.columns)

    # denoise
    # if cpts is not None:
    #     groupby_col = ['chrom', 'cpts']
    # else:
    #     groupby_col = ['chrom']

    epg = expand_profile.groupby('chrom')
    smooth_profile = pd.DataFrame()
    for c, ep in epg:
        temp_p = ep.iloc[:, 3:].apply(lambda x: signal.medfilt(x, kernel_size), axis=0)
        smooth_profile = pd.concat([smooth_profile, temp_p])

    # drop duplicated sites
    def temp_fun(x):
        vc = x.value_counts() / x.__len__()
        if max(vc) > (1 - dup):
            return True
        else:
            return False
    simply_profile_pos = smooth_profile.apply(temp_fun, axis=1)
    simply_profile = smooth_profile.loc[~simply_profile_pos, :]
    return simply_profile, smooth_profile, expand_profile



