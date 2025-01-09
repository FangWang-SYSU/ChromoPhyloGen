#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@File   : cli_main.py
@Author : XinWang
"""

import argparse
from ChromoPhyloGen.main_func import chromosome_event


def main():
    description = "A package for quantifying chromothripsis and seismic events to dissect tumor evolution with single-cell resolution."
    author = 'Author: wangxin, Email: wangx768@mail2.sysu.edu.cn'

    parser = argparse.ArgumentParser(
        description=description,
        epilog=author,
        add_help=True)
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    parser.add_argument("-i", '--cnv_data', help="All node CNA in newick.", required=True)
    parser.add_argument("-o", "--output", default='./', help=u"The output path.", required=True)
    parser.add_argument("-t", '--tree', default='none', help="Newick format file path.")
    parser.add_argument("-p", "--prefix", help=u"Prefix for output file names.")
    parser.add_argument("-r", "--random_num", default=1000, help=u"Random number for creating a null distribution.")
    parser.add_argument("-c", "--cancer_type", default='ALL', help=u"Select a cancer type for estimating WGD. The default is all.")
    parser.add_argument("-n", "--ncores", default=1, help=u"Number of cores required to run copy number variation events.")

    args = parser.parse_args()
    chromosome_event(
        node_data=args.cnv_data,
        newick=args.tree,
        output_dir=args.output,
        prefix=args.prefix,
        cancer_type=args.cancer_type,
        cores=int(args.ncores),
        randome_num=int(args.random_num),
        verbose=True)

