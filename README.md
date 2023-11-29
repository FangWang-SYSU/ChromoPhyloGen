Description
===========

This framework infers the phylogeny of tumor cells characterized by complex chromosome rearrangements, distinguishes between chromothripsis and seismic rearrangements and assesses their association with the risk of errors during tumor cell division, enabling the establishment of clonal lineages. 

![avatar](github.png)

System requirements and dependency
==================================
Software package development environment： 
    
    macOS
    Python 3.11.3


This package requires Python version 3.9 or greater. 

In addition, `ChromoPhyloGen` will automatically install the following dependency packages during the installation process. 

`biopython`,`Cython`,`numpy`,`pandas`,`scikit-learn`,`scipy`,`tqdm`,`snfpy`,`matplotlib`,`seaborn`,`statsmodels`

But it is recommended to prepare dependency packages in your environment in advance.


Installation
============
First create a virtual environment for ChromoPhyloGen, but this is not required.
```shell
conda create --name ChromoPhyloGen_env python=3.9
conda activate ChromoPhyloGen_env
```
### 1.From Pypi
You can install the latest release from PyPi, with:
```shell
pip install ChromoPhyloGen
```

### 2.Source code
You can install this package by opening a command terminal and running the following:
```shell
git clone https://github.com/FangWang-SYSU/ChromoPhyloGen.git
cd ChromoPhyloGen
pip install .
```

After the installation is complete, you can use the `ChromoPhyloGen --version`  command to test whether the software was successfully installed.

Usage
=====
```
optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  -I INPUT, --input INPUT
                        single-cell copy number profile.
  -O OUTPUT, --output OUTPUT
                        The output path.
  -p PREFIX, --prefix PREFIX
                        Prefix for output file names.
  -r RESOLUTION, --resolution RESOLUTION
                        Lineage partitioning resolution(default=1).
  -t HUFFMAN_SPLIT_THRESHOLD, --huffman_split_threshold HUFFMAN_SPLIT_THRESHOLD
                        huffman split threshold(default=0.9)
  -n N_NEIGHBORS, --n_neighbors N_NEIGHBORS
                        Number of neighbors for creating affinity matrix in SNF(default=5).
  -m MIN_CLONE_SIZE, --min_clone_size MIN_CLONE_SIZE
                        When min_clone_size is reached, division will no longer continue(default=0.1*cell_number).
  -s SCORING, --scoring SCORING
                        Whether to run Scoring the chromosomal rearrangements.
  -R RANDOM_NUM, --random_num RANDOM_NUM
                        Random number for creating a null distribution.
  -C CANCER_TYPE, --cancer_type CANCER_TYPE
                        Select a cancer type for estimating WGD. The default is all.
  -c CORES, --cores CORES
                        Number of cores required to run copy number variation events.
  -d DRAW, --draw DRAW  Draw tree and CNA heatmap.

Author: wangxin, Email: wangx768@mail2.sysu.edu.cn

```

Input files
===========

The input file of ChromoPhyloGen needs to be an integer copy number spectrum:

    The row is the genome segment,
    the first column is the chromosome, 
    the second column is the genome starting coordinate, 
    the third column is the genome end coordinate,
    and the other columns are the integer copy numbers at the cell level.

  	chr	start	end	cell_1	cell_2	cell_3 ...
    1	100167143	100220943	2	2	2 ...
    1	100504443	100559237	2	1	2 ...
    1	101395562	101451560	2	3	4 ...

>Connection with [inferCNV](https://github.com/broadinstitute/inferCNV):
> 
> To obtain the integer copy number，we propose to identify peaks and infer their intervals, with each interval representing an integer copy number (detail in method).

Examples
============

### Run in command line
The example data `exampleCNA.txt` was included in the `ChromoPhyloGen` package, and you can also download it from [here](https://github.com/FangWang-SYSU/ChromoPhyloGen/blob/main/ChromoPhyloGen/data/exampleCNA.txt)
```shell
#ChromoPhyloGen [-h] [--version] -I INPUT -O OUTPUT [-p PREFIX] [-r RESOLUTION] [-t HUFFMAN_SPLIT_THRESHOLD] [-n N_NEIGHBORS] [-m MIN_CLONE_SIZE] [-s SCORING] [-R RANDOM_NUM] [-C CANCER_TYPE] [-c CORES]
mkdir ChromoPhyloGen_out
ChromoPhyloGen \
    -I ./ChromoPhyloGen/data/exampleCNA1.txt \
    -O ./ChromoPhyloGen_out \
    -p example_ \
    -r 1 \
    -t 1 \
    -n 5 \
    -C ALL \
    -c 8 \
    -d 1
```
> The `-r` parameter is used for lineage partitioning resolution, where a higher value indicates greater precision. 
> The `-t`parameter represents the proportion of subtree splitting during the `Huffman process` and takes values between 0 and 1. 
> A higher value implies a lower probability of splitting two already merged cells.
> Additionally, the estimation of chromosome rearrangement score is influenced by the `-R` parameter, with a larger value leading to longer runtime.

### Run in python
```python
# load package
import ChromoPhyloGen
from ChromoPhyloGen import ChromoPhyloGen as fp
# 1.infer tree
fp.run(cna_dir = ChromoPhyloGen.__path__[0] + '/data/exampleCNA.txt',
       output_dir='ChromoPhyloGen_out',
       prefix='example_',
       resolution=1,
       clone_thr=1,
       n_neighbors=5,
       plot_png=False,
       verbose=True)

# 2.score
fp.chromosome_event(
     'ChromoPhyloGen_out',
     prefix='example_',
     cancer_type='ALL',
     cores=8,
     randome_num=1000,
     verbose=True)
```

Output files
============

| id | FileName | Description |
| :----- | :-----| :---- |
| 1| *cell_info.txt | Cell variation information in trace |
| 2| *all_node_data.txt | Cell CNA profile, including internal node, name by "virtual_" |
| 3| *cell_tree.newick| Single cell trace file，format newick.  |
| 4| *re_score.txt | Chromosomal rearrangements score.|
| 5| *chromothripsis_score.txt | Chromosomal chromothripsis score. |
| 6| *re_socre_pvalue.txt | P-value of chromosomal rearrangements score. |
| 7| *mode.txt | Define mode for each cell. |
| 8| *error_risk_score.txt | Division error risk. |
| 9| *error_risk_score.txt | Division error risk. |

> where * corresponds to the parameter `-p`, indicating the prefix of the output file name。


### 1.*cell_info.txt: Cell variation information in trace.
    name	Root_gain_loc	Root_loss_loc	Root_gain_cn	Root_loss_cn	Parent_gain_loc	Parent_loss_loc	Parent_gain_cn	Parent_loss_cn	Mitosis_copy	Mitosis_dd_loc	Mitosis_ad_loc	Mitosis_time	Pseudotime_tree	Mitosis_time_next	aneu_rate	copy_rate	status
    root										1042.0	163.0	20.0	0.0	40.0	0.013		Aneuploidy
    cell_2	160.0	428.0	428.0	160.0	160.0	428.0	160.0	428.0	11617.0	825.0	478.0	31.0	21.0	21.0	0.039	0.951	Aneuploidy
    cell_3	151.0	629.0	629.0	151.0	151.0	629.0	151.0	629.0	11425.0	1801.0	182.0	32.5	19.0	19.0	0.014	0.936	Aneuploidy
    cell_4	51.0	360.0	360.0	51.0	291.0	410.0	291.0	410.0	11504.0	779.0	189.0	25.0	46.0	25.0	0.015	0.942	Aneuploidy

```
[Root|Parent]_[gain|loss]_[loc|cn]: The number of sites or copies accumulated (gain|loss) relative to the (Root|Parent) node.
Mitosis_copy: The count of genomic segments sharing the same copy number states between the current node and its parent node(D_ss).
Mitosis_dd_loc: The count of genomic segments different copy number states between the current node and its parent node(D_ds).
Mitosis_ad_loc: The count of aneuploidy segregation states between the current node and its parent node(D_as).
Mitosis_time: Branch length of current node.
Pseudotime_tree: Pseudotime of current node in tree.
Mitosis_time_next: The branch length of the next mitosis of the current node.
aneu_rate: Rate of aneuploidy segregation states.
copy_rate: Rate of same copy number states.
status: The current cell mitotic state inferred based on aneu_rate and copy_rate.
```

### 2.*all_node_data.txt: Cell CNA profile, including internal node, name by "virtual_".

    	1_977836_977836	1_1200863_1200863 ...
    cell1	1.0	1.0	...
    cell2	2.0	2.0	...
>Integer copy number profile of all nodes in tree. Rows are cells, columns are genome segments.

### 3.*cell_tree.newick: Single cell trace file，format newick.
>Stores the structural information of the evolutionary tree, including branch length.

### 4.*re_score.txt
    ,1,2,3,...
    cell_1,0.761,0.800,0.793,...
    cell_2,0.88,0.0,0.636,0.659,...
    cell_3,0.88,0.957,0.783,...
>The level of chromosomal rearrangements. Rows are cells, columns are chromosome. If the value is -1, it means that the chromosome has not changed significantly.

### 5.*chromothripsis_score.txt
    ,1,2,3,...
    cell_1,0.001,0.008,0.705,...
    cell_2,0.780,0.005,0.837,...
    cell_3,0.890,0.907,0.463,...
>The level of chromosomal chromothripsis. Rows are cells, columns are chromosome.

### 6.*re_socre_pvalue.txt
>P-value of chromosomal rearrangements (detail in methods). 

### 7.*mode.txt
    ,chr_num,wgd,chromothripsis_num,seismic_num,chromothripsis_score,seismic_score,BFB
    cell_1,19,WGD1,14,5,0.009,0.0149,1368
    cell_2,11,WGD1,7,4,0.001,0.006,234
    cell_3,12,WGD0,6,6,0.002,0.001,209
```
chr_num: The total number of chromosomes in chromothripsis and seismic.
wgd: Cell chromosome WGD type (wgd0 is diploid, wgd1 involves a single whole genome duplication, and wgd2 entails multiple whole genome duplications.).
chromothripsis_num: The total number of chromosomes in chromothripsis.
seismic_num: The total number of chromosomes in seismic.
chromothripsis_score: Average chromothripsis score in chromothripsis chromosome.
seismic_score: Average seismic score in seismic chromosome.
BFB: The number of aneuploidy segregation states.
```

### 8.*error_risk_score.txt
    cell_1	0.24
    cell_2	0.19
    cell_3	0.87
```
first column: Cell id.
second column: error_risk score.
```

### 9.*tree.png
    Optional parameter '-d' or 'plot_png'. 
    If set to 1, it will draw the phylogenetic tree and heatmap of CNA profile. 
    If set to 0, it will not be drawn.

![avatar](example_tree.png)

>Note that drawing requires `matplotlib` and `seaborn` packages

Developer
=========
Xin Wang (wangx768@mail2.sysu.edu.cn)

Draft date
==========
Nov.15, 2023
