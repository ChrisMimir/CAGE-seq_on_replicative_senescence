# CAGE-seq_replicative_senescence
Scripts for publication of CAGE-seq

Please install required packages and change the directory of data to your data location before running the scripts. 

For reproducing most of the figures, please run CAGE-seq_Analysis.R first, then you can reporduce the Extended_fig.R and chip_analysis.R in the same enviroment. 

The processed data of ChIP-seq used in chip_analysis.R can be downloaded in this repository: diff_sites_final_c3.0_cond1.bed and diff_sites_final_c3.0_cond2.bed.And also you can generate by using the command in MACS3_ChIP-seq.txt. 

human_TE_TSS.txt is required for transposable elements analysis(Extended Data Fig 1b)

For reproding the figure 2a, b, Extended Data Figure 2a, b, please run diff_puffin.ipynb. 

For reproducing the figure 4f, running GTEx_aging.R in a separate environment is recommanded.

For reproducing the figure 4a,b, please run CAGE_to_ANANSE.R in a seperate environment, then run run_ananse.sh in the terminal. For the detail please refer to the ANANSE user mannul (https://anansepy.readthedocs.io/en/master/)

For reproducing the Extended Data Figure 2c, run CLIPNET (Extended Data figure 2c).txt first and then transform the data format in diff_puffin.ipynb. 


