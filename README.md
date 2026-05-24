# CAGE-seq_replicative_senescence
Scripts for the publication of CAGE‑seq

Please install all required packages and change the data directory to your own data location before running the scripts.

To reproduce most of the figures, please run CAGE-seq_Analysis.R first. After that, you can reproduce Extended_fig.R and chip_analysis.R in the same environment.

The processed ChIP‑seq data used in chip_analysis.R can be downloaded from this repository: diff_sites_final_c3.0_cond1.bed and diff_sites_final_c3.0_cond2.bed. You may also generate them using the commands provided in MACS3_ChIP-seq.txt.

human_TE_TSS.txt is required for the transposable element analysis (Extended Data Fig. 1b).

To reproduce Figure 2a, 2b and Extended Data Figure 2a, 2b, please run diff_puffin.ipynb.

To reproduce Figure 4f, it is recommended to run GTEx_aging.R in a separate environment.

To reproduce Figure 4a, 4b, please run CAGE_to_ANANSE.R in a separate environment, and then execute run_ananse.sh in the terminal. For details, please refer to the ANANSE user manual:
https://anansepy.readthedocs.io/en/master/

To reproduce Extended Data Figure 2c, run CLIPNET (Extended Data Figure 2c).txt first, and then convert the output format in diff_puffin.ipynb. For the detail usage, please refer to CLIPNET repository (https://github.com/Danko-Lab/clipnet/) and tfmodisco repository (https://github.com/kundajelab/tfmodisco)


