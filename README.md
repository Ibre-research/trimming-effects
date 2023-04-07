# Negligible effects of read trimming on the accuracy of germline variant calling in the human genome

This repository contains data and code pertinent to the analysis of adapter trimming effects on variant calling

![The workflow of the main analysis](./TP_count.png)

The repository contains the following files and folders:

* `preprocessing_scripts` folder contains the scripts used for data analysis (read trimming, alignment, variant calling, and evaluation using the hap.py software;
* `trim_explore.R` is the main script which contains all code used to perform downstream data analysis and figure preparation;
* `trimming_r8.tsv.gx` is the main data file containing all results of the variant calls' evaluation with hap.py;
* `adapter_read_counts.tsv` and `base_counts.tsv` contain the calculated number of reads and bases, respectively, affected by adapter removal.`
