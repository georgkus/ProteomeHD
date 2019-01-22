# About this repository
R scripts and data files related to the manuscript "The human proteome co-regulation map reveals functional relationships between proteins" by Kustatscher <i> et al </i> (2019). The article is available in bioRxiv {ENTER REF HERE}.

To reproduce the results of our paper download the files, set the R working directory to the download folder and execute the appropriate scripts. All scripts were tested in Ubuntu 18.04, but should generally be working on all operating systems. Below is a brief description of each file.


# R scripts
- <b> tune_treeclust.R: </b> This script was used to perform a grid search to optimise treeClust / rpart hyperparameters (serule and cp) against true and false positives pairs annoated in the Reactome_TP_FP_10perc_subset_for_GS.csv file.


# Data files
- <b> ProteomeHD_v1_1.7z: </b> This compressed csv file is ProteomeHD, consisting of 10,323 proteins and 294 SILAC ratios. An exact copy of this file is Supplementary_Table_S1.csv {DOUBLE CHECK TABLE NAME} in the publication.

- <b> Reactome_TP_FP.7z: </b> This is a compressed csv file containing the "gold standard", i.e. a list of true and false positive protein pairs. It was generated as described in the article.

- <b> Reactome_TP_FP_10perc_subset_for_GS.7z: </b> A





