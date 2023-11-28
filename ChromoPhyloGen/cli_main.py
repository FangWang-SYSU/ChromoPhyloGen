#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@File   : cli_main.py
@Author : XinWang
"""


import argparse
from .ChromoPhyloGen import *


def main():
    description = "A package for quantifying chromothripsis and seismic events to dissect tumor evolution with single-cell resolution."
    author = 'Author: wangxin, Email: wangx768@mail2.sysu.edu.cn'

    parser = argparse.ArgumentParser(
        description=description,
        epilog=author,
        add_help=True)
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    parser.add_argument("-I", '--input', help="single-cell copy number profile.", required=True)
    parser.add_argument("-O", "--output", default='scTrace', help=u"The output path.", required=True)
    parser.add_argument("-p", "--prefix", help=u"Prefix for output file names.")
    parser.add_argument("-r", "--resolution", default=1, help=u"Lineage partitioning resolution(default=1).")
    parser.add_argument("-t", "--huffman_split_threshold", default=0.9, help=u"huffman split threshold(default=0.9)")
    parser.add_argument("-n", "--n_neighbors", default=5, help=u"Number of neighbors for creating affinity matrix in SNF(default=5).")
    parser.add_argument("-m", "--min_clone_size", help=u"When min_clone_size is reached, division will no longer continue(default=0.1*cell_number).")
    parser.add_argument("-s", "--scoring", default=1, help=u"Whether to run Scoring the chromosomal rearrangements.")
    parser.add_argument("-R", "--random_num", default=1000, help=u"Random number for creating a null distribution.")
    parser.add_argument("-C", "--cancer_type", default='ALL', help=u"Select a cancer type for estimating WGD. The default is all.")
    parser.add_argument("-c", "--cores", default=1, help=u"Number of cores required to run copy number variation events.")
    parser.add_argument("-d", "--draw", default=0, help=u"Draw tree and CNA heatmap.")

    args = parser.parse_args()

    run(args.input,
        args.output,
        prefix=args.prefix,
        resolution=float(args.resolution),
        clone_thr=float(args.huffman_split_threshold),
        n_neighbors=int(args.n_neighbors),
        min_clone_size=args.min_clone_size,
        plot_png=int(args.draw),
        verbose=True
    )
    if args.scoring != '0':
        chromosome_event(
            args.output,
            prefix=args.prefix,
            cancer_type=args.cancer_type,
            cores=int(args.cores),
            randome_num=int(args.random_num),
            verbose=True)
