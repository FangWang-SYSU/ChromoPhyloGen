#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@File   : events_define.py
@Author : XinWang
"""
import warnings
warnings.filterwarnings('ignore')

import pandas as pd
import numpy as np
from scipy.stats import norm
from scipy.stats import skew
from collections import Counter
import fast_score

from concurrent.futures import ThreadPoolExecutor, as_completed
import tqdm
import os

def define_WGD(data, TCGA_cancer_name_short='ALL'):
    current_dir = os.path.dirname(__file__)
    file_path = os.path.join(current_dir, './data/cancer_ploity_gmm_param.txt')
    parameter = pd.read_table(file_path, sep='\t')
    if TCGA_cancer_name_short in parameter['cancer_name'].values:
        cell_ploidy = np.mean(data, axis=1)

        tmp_param = parameter.loc[parameter['cancer_name'] == TCGA_cancer_name_short,]

        wgd0_q = tmp_param.iloc[0, 2] * norm.pdf(cell_ploidy, tmp_param['mu_0'], tmp_param['sd_0'])
        wgd1_q = tmp_param.iloc[0, 5] * norm.pdf(cell_ploidy, tmp_param['mu_1'], tmp_param['sd_1'])
        if np.isnan(tmp_param.iloc[0, 8]):
            wgd2_q = np.zeros_like(cell_ploidy)
        else:
            wgd2_q = tmp_param.iloc[0, 8] * norm.pdf(cell_ploidy, tmp_param['mu_2'], tmp_param['sd_2'])

        # select max post
        wgd_res = np.column_stack((wgd0_q, wgd1_q, wgd2_q))
        wgd_res = wgd_res / np.sum(wgd_res, axis=1, keepdims=True)

        # correct
        for i in range(len(cell_ploidy)):
            tmp_pd = cell_ploidy[i]
            if tmp_pd<tmp_param['mu_0'].values[0]:
                wgd_res[i,0] = 1
                wgd_res[i,[1,2]] = 0
            try:
                if tmp_pd > tmp_param['mu_2'].values[0]:
                    wgd_res[i, 2] = 1
                    wgd_res[i, [0,1]] = 0
            except:
                if tmp_pd > tmp_param['mu_1'].values[0]:
                    wgd_res[i, 1] = 1
                    wgd_res[i, [0, 2]] = 0

        label = np.array(["WGD0", "WGD1", "WGD2"])[np.argmax(wgd_res, axis=1)]
        return {'wgd_density': wgd_res, 'wgd_type': label}
    else:
        print('Cannot find cancer name.')
        return None


def cal_lambda(x):
    bkps = np.where(np.diff(x) != 0)[0]
    if len(bkps) == 0:
        return 0
    return np.mean(np.diff(bkps))

def cal_possion_lambda(data, pos):
    def tmp_fun(aa):
        bkps = np.where(np.diff(np.array(aa)) != 0)[0]
        return len(bkps)
    possion_param =[]
    gp = pos.groupby('chr')
    for ix, tmp_pos in gp:
        tmp_x = data.loc[:, tmp_pos['index']]
        bk_num = tmp_x.apply(tmp_fun, 1)
        possion_param.append(np.mean(bk_num))
    return possion_param


def cal_R(cnv_data):
    bkps = np.where(np.diff(cnv_data) != 0)[0]
    if len(bkps) <= 1:
        return 0
    # U
    seg_len = np.abs(np.diff(bkps))
    U = 1/(1+np.exp(-skew(seg_len, nan_policy='omit')))
    if np.isnan(U):
        U = 1e-5
    # A
    A = np.std(cnv_data)
    A = A / np.mean(cnv_data)
    if np.isnan(A):
        A = 0
    if A > 1:
        A = 1
    return U*A


def cal_g(x):
    counter = Counter(x)
    if len(counter) <= 1:
        return -1 #'neutral'

    counter_list = list(counter.items())
    counter_list.sort(key=lambda a: -a[1])
    qmax = counter_list[0][0]
    dci = x - qmax
    g = sum(np.isin(dci, [1,-1])) / sum(dci!=0)
    return g

def permutation_R(seed_num, cna, exp_lam, possion_lam, random_num):
    true_R = cal_R(cna)
    perm_data = fast_score.perm_cnv_pyx(np.array(cna, dtype=float),
                       np.array(exp_lam, dtype=float),
                       np.array(possion_lam, dtype=float),
                       int(random_num),
                       int(seed_num)
                       )
    num_greater = np.sum(np.array(perm_data) >= true_R)
    dc = {'rearrange_score': true_R, 'limit': cal_g(np.array(cna))}
    return {'dc': dc, 'num_greater': num_greater, 'x': seed_num}



def calculate_rearrange_score(all_cell_cnv,
                              random_num=100,
                              cores=1):
    gene_loc = pd.DataFrame([x.split('_') for x in all_cell_cnv.columns], columns=['chr', 'start', 'end'])
    gene_loc['start'] = pd.to_numeric(gene_loc['start'])
    gene_loc['end'] = pd.to_numeric(gene_loc['end'])
    gene_loc['index'] = all_cell_cnv.columns
    gene_loc.index = gene_loc['index']

    exp_lambda = all_cell_cnv.apply(cal_lambda, 1).to_list()
    possion_lambda = cal_possion_lambda(all_cell_cnv, gene_loc)
    chr_name = []
    for i in gene_loc['chr'].unique():
        if i == 'X':
            chr_name.append('23')
        elif i == 'Y':
            chr_name.append('24')
        else:
            chr_name.append(i)

    all_limit_prop = pd.DataFrame()
    all_rearrange_score = pd.DataFrame()
    all_rearrange_score_pvalue = pd.DataFrame()
    for i in chr_name:
        i_gene = gene_loc.loc[gene_loc['chr'] == i, :].copy()
        i_gene.sort_values('start', inplace=True)
        tmp_cnv = all_cell_cnv[i_gene['index']]
        tmp_cnv = np.array(tmp_cnv)
        tmp_res_limit = []
        tmp_res_rearrange = []
        tmp_res_rearrange_pvalue = []
        res = []
        with ThreadPoolExecutor(max_workers=cores) as executor:
            obj_list = []
            for iter_num in range(tmp_cnv.shape[0]):
                obj = executor.submit(permutation_R,
                                      iter_num, tmp_cnv[iter_num, :],
                                      exp_lambda[iter_num],
                                      possion_lambda[int(i)-1],
                                      random_num)
                obj_list.append(obj)

            for future in tqdm.tqdm(as_completed(obj_list), total=tmp_cnv.shape[0], desc=f"Processing chr{i}"):
                tmp_future = future.result()
                res.append(tmp_future)

        cell_order = []
        for icell in range(len(res)):
            dc = res[icell]['dc']
            num_greater = res[icell]['num_greater']
            tmp_res_limit.append(dc['limit'])
            tmp_res_rearrange.append(dc['rearrange_score'])
            tmp_res_rearrange_pvalue.append((num_greater + 1) / (random_num + 1))
            cell_order.append(res[icell]['x'])

        cell_order = np.argsort(cell_order)
        tmp_res_limit = np.array(tmp_res_limit)[cell_order]
        tmp_res_rearrange = np.array(tmp_res_rearrange)[cell_order]
        tmp_res_rearrange_pvalue = np.array(tmp_res_rearrange_pvalue)[cell_order]

        all_limit_prop[i] = tmp_res_limit
        all_rearrange_score[i] = tmp_res_rearrange
        all_rearrange_score_pvalue[i] = tmp_res_rearrange_pvalue

    all_limit_prop.index = all_cell_cnv.index
    all_rearrange_score.index = all_cell_cnv.index
    all_rearrange_score_pvalue.index = all_cell_cnv.index

    all_limit_prop.columns = all_rearrange_score.columns = all_rearrange_score_pvalue.columns = chr_name
    all_limit_prop.index = all_rearrange_score.index = all_rearrange_score_pvalue.index = all_cell_cnv.index

    return {
        'all_chromothripsis_prop': all_limit_prop,
        'all_rearrange_score': all_rearrange_score,
        'all_rearrange_score_pvalue': all_rearrange_score_pvalue
    }


def define_BFB(tree, data):
    bfb_res = {}
    def traverse_tree(clade):
        childs = clade.clades
        if len(childs) == 2:
            c1, c2 = clade.clades
            c1_data = np.array(data.loc[c1.name, :])
            c2_data = np.array(data.loc[c2.name, :])
            p_data = np.array(data.loc[clade.name, :])
            # pos = (np.round((c1_data + c2_data) / 2)) == p_data
            # pos = pos & (c1_data!=c2_data)
            pos = (((c1_data + c2_data) % 2) == 0) & (c1_data!=c2_data)
            bfb_res[c1.name] = {'pos': pos, 'bfb_num': np.sum(pos)}
            bfb_res[c2.name] = {'pos': pos, 'bfb_num': np.sum(pos)}
            traverse_tree(c1)
            traverse_tree(c2)
        else:
            pass
    traverse_tree(tree.root.clades[0])
    return bfb_res

def define_mode(obj, p_thr):

    # WGD
    wgd_label = obj['wgd']['wgd_type']

    # BFB
    BFB_res = {k: v['bfb_num'] for k, v in obj['BFB_res'].items()}

    # chromothripsis
    all_rearrange_score = obj['rearrange_score']['all_rearrange_score']
    all_rearrange_score_pvalue = obj['rearrange_score']['all_rearrange_score_pvalue']
    all_limit_prop = obj['rearrange_score']['all_chromothripsis_prop']

    # is chromothripsis chr
    chr_pos = all_rearrange_score_pvalue<p_thr
    res = pd.DataFrame(chr_pos.sum(axis=1), columns=['chr_num'])

    #
    res['wgd'] = wgd_label
    res['chromothripsis_num'] = ((all_limit_prop > 0.5) & chr_pos).sum(axis=1)
    res['seismic_num'] = ((all_limit_prop <= 0.5) & chr_pos).sum(axis=1)
    res['chromothripsis_score'] = np.mean(np.array(all_rearrange_score), axis=1,
                                 where=np.array((all_limit_prop > 0.5) & chr_pos).tolist())
    res['seismic_score'] = np.mean(np.array(all_rearrange_score), axis=1,
                                   where=np.array((all_limit_prop <= 0.5) & chr_pos).tolist())

    res['BFB'] = [BFB_res[k] if k !='root' else 0 for k in all_rearrange_score.index]

    res.fillna(0, inplace=True)
    return res

def define_CNA_mechnism(sctc_obj, cancer_type='ALL',
                        random_num=100,
                        cores=1):
    all_cell_cnv = sctc_obj['cnv_data']
    all_cell_cnv = all_cell_cnv.loc[~all_cell_cnv.index.str.startswith('virtual'), :]
    all_cell_cnv = all_cell_cnv.loc[~all_cell_cnv.index.str.startswith('root'), :]
    # WGD
    wgd_res = define_WGD(all_cell_cnv, TCGA_cancer_name_short=cancer_type)
    # chromothripsis & chromoplexy
    rearrange_score = calculate_rearrange_score(all_cell_cnv,
                                                random_num=random_num,
                                                cores=cores)
    # BFB
    BFB_res = define_BFB(sctc_obj['tree'], sctc_obj['cnv_data'])

    orig_res = {
        'rearrange_score': rearrange_score,
        'BFB_res': BFB_res,
        'wgd': wgd_res
    }
    mode = define_mode(orig_res, p_thr=1/random_num)
    return {
        'orig_res': orig_res,
        'mode': mode
    }
