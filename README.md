# About this repository
R scripts and data files related to the manuscript "The human proteome co-regulation map reveals functional relationships between proteins" by Kustatscher <i> et al </i> (2019). The article is available in bioRxiv {ENTER REF HERE}.

To reproduce the results of our paper download the files, set the R working directory to the download folder and execute the appropriate scripts. All scripts were tested in Ubuntu 18.04, but should generally be working on all operating systems. Below is a brief description of each file.


# R scripts
- <b> Coregulation_scores.R: </b> This script creates the co-regulation scores for protein pairs in ProteomeHD. It first calculates treeClust dissimilarities (using hyperparameters optimised with the tune_treeclust.R script) and turns them into similarities (i.e. 1 - dissimilarity). Using the WGCNA package, the script subsequently performs a sigmoid transformation to create an adjacency matrix (using hyperparameters optimised with the tune_sigmoid.R script). Finally, the script applies the topological overlap measure, thereby creating the "co-regulation scores". The script outputs three files:
    (a) <i> treeClust_similarities.csv: </i> Contains the treeClust similarities, the adjacency matrix and the TOM-transformed values.
    (b) <i> coregulation_scores.csv: </i> Contains only the final pairwise coregulation scores (treeClust + sigmoid transformation + TOM). 
    (c) <i> ScaleFreeNess.png: </i> A plot showing that the resulting network is scale free.

- <b> Fig1_coreg_proteins.R </b> This script handles multiple tasks related to Figure 1 of the manuscript. It produces a plot that shows the distribution of the co-regulation scores. It produces a plot showing that co-regulated protein pairs are enriched for proteins from the same complex and subcellular localisation and for enzymes catalyzing consecutive metabolic reactions. Membership of protein complexes and metabolic pathways is based on Reactome annotation (the relavant input file, homo_sapiens.interactions.txt, can be downloaded from the Reactome website). A file containing all reviewed Uniprot proteins together with their subcellular localisation was downloaded from the Uniprot website (<i>uniprot-all.tab</i>) and is available in compressed form in this repository. This script also creates the <i>coregulated_protein_pairs.csv</i> file, which contains all co-regulated protein pairs in the form “Protein X” - “Coregulated with protein X”. This file is part of the manuscript as Supplementary Table 3 {UPDATE NUMBER IF NECESSARY}. Finally, the script produces a barplot showing how many co-regulation partners proteins tend to have (in bins of 0, 1-5, 6-20 etc).

- <b> topGO_enrichment_analysis_BP.R / _CC.R / _MF.R </b> Three almost identical scripts calculating the enrichment of gene ontology terms among groups of co-regulated proteins. The analysis uses the topGO R package. The difference between the three scripts is that they analyse GO Biological Processes (BP), Cellular Component (CC) and Molecular Function (MF). Each script takes the corresponding input file (<i>GO_associations_BP.tsv, GO_associations_CC.tsv or GO_associations_MF.tsv</i>), which are available in this respository in compressed form. These files contain GO associations downloaded from QuickGO using taxon = human, qualifier = NULL and aspect = Biological process (or aspect = Molecular function, or aspect = Cellular component). 

- <b> Coregulation_score_vs_corr_part1.R / _part2.R </b> These two scripts compare the performance of treeClust dissimilarities with three correlation metrics (Spearman's rho, Pearson correlation and biweight midcorrelation). The first script calculates the correlation coefficients for ProteomeHD and saves them in a large, temporary output file called <i> ProHD_correlations.csv </i>. The second script performs a precision - recall analysis using the Reactome gold standard, taking the correlation coefficients from <i> ProHD_correlations.csv </i> and the treeClust and treeClust + TOM similarities from <i> treeClust_similarities.csv: </i> (produced by the Coregulation_scores.R script). 

- <b> tSNE_map.R </b> This script produces the tSNE map of ProteomeHD. It takes the coregulation_scores.csv file as input, which is itself created by the Coregulation_scores.R script (also available as download elsewhere, see see below). Annotation of the map is based on the "tSNE_map_annotation.csv" file, which contains Uniprot-based organelle annotation as well as manual annotation of zoom regions. The script outputs <i> tSNE_coordinates_and_annotation.csv </i>, which contains the tSNE coordinates and the annotations together, and which is also saved as Supplementary Table 4. {DOUBLE CHECK NAME}

- <b> Orthogonal_Methods.R </b> This script compares protein - protein associations detected by co-regulation with those found in STRING, BioGRID, IntAct and BioPlex 2.0. It uses downloaded files from these databases as input, which can be found on the respective websites. These include <i> 9606.protein.links.detailed.v10.5.txt</i> and <i> BIOGRID-ALL-3.4.152.tab2.txt </i>. The IntAct_interactions.csv file is provided here, because it is a pre-processed version of the original data from IntAct (pre-processed only to reduce file size by removing unnecessary information). The file nature22366-s2_Bioplex2.xlsx is Table S1 from the BioPlex 2.0 paper (Huttlin et al, Nature, 2017). Finally, the script uses <i> SwissProt_IDs_HoSa_Nov_18.tab </i>, downloaded from Uniprot on 19.11.2018 and containing their up-to-date IDs for removing old and unused protein IDs. This file is available here.

- <b> Uncharacterised_proteins.R </b> XXXXXXXX SwissProt_HS_AnnotScore_Mass.tab
Cancer_gene_census_Uniprot.csv

- <b> tune_treeclust.R: </b> This script was used to perform a grid search to optimise treeClust / rpart hyperparameters (serule and cp) against true and false positives pairs annoated in the Reactome_TP_FP_10perc_subset_for_GS.csv file. The optimal values turned out to be cp = 0.105 and serule = 1.8, providing a ~10% improvement over default settings.

- <b> tune_sigmoid.R: </b> This script was used to perform a grid search to find the optimal parameters for WGCNA's topological overlap matrix (TOM). These parameters were mu and alpha, which relate to the sigmoid transformation taking place before calculating TOM. The optimal values turned out to be mu = 0.91 and alpha = 37. This step provided a further 10% improvement in performance. Note that this script is designed for execution on a server with 30 free cores.

# Data files
- <b> ProteomeHD_v1_1.7z: </b> This compressed csv file is ProteomeHD, consisting of 10,323 proteins and 294 SILAC ratios. An exact copy of this file is Supplementary_Table_S1.csv {DOUBLE CHECK TABLE NAME} in the publication.

- <b> Reactome_TP_FP.7z: </b> This is a compressed csv file containing the "gold standard", i.e. a list of true and false positive protein pairs. It was generated as described in the article.

- <b> Reactome_TP_FP_10perc_subset_for_GS.7z: </b> This is a compressed csv file containing a subset of the Reactome gold standard (Reactome_TP_FP.7z). It was used to optimise treeClust hyperparameters in a grid search.

- <b> uniprot-all.7z </b> File containing all reviewed Uniprot proteins together with their subcellular location.

- <b> GO_associations_BP.7z / _CC.7z / _MF.7z </b> These files contain GO associations downloaded from QuickGO using taxon = human, qualifier = NULL and aspect = Biological process (or aspect = Molecular function, or aspect = Cellular component).

- <b> tSNE_map_annotation.csv: </b> Used for annotation of the tSNE map by the tSNE_map.R script. Contains Uniprot-based organelle annotation as well as manual annotation of zoom regions.

- <b> IntAct_interactions.7z </b> To reduce file size, the table downloaded from IntAct was pre-filtered to exlude interactions from non-human species, interactors that were not proteins, interactions without Uniprot ID as well as self-interactions and duplicates. Finally, annotation strings were simplified, e.g. by reducing IDs from "uniprotkb:P06748" to "P06748".

- <b> SwissProt_IDs_HoSa_Nov_18.tab </b> Up-to-date Homo Sapiens Uniprot IDs (19.11.2018) for removing old and unused protein IDs from input files.

- <b> SwissProt_HS_AnnotScore_Mass.tab </b> Contains reviewed human Uniprot IDs together with their Uniprot annotation score and molecular weight.

- <b> Cancer_gene_census_Uniprot.csv </b> A file containing annotation from the cancer gene census (COSMIC v81). Originally downloaded from http://cancer.sanger.ac.uk/census. Entrez Gene IDs were converted to reviewed Homo Sapiens Uniprot IDs, ambiguous mappings were removed.

# Data files available elsewhere
- <b> coregulation_scores.csv: </b> This file contains the final pairwise co-regulation scores for all 12,562,578 protein pairs analysed in our study. It exceeds Github's 100 MB limit even in compressed form. It can be generated by simply running the Coregulation_scores.R script, which only takes ProteomeHD_v1_1.csv as input. Alternatively, a copy of the file can be downloaded from either www.proteomeHD.net or the PRIDE repository {ADD IDENTIFIED HERE}.

- <b> 9606.protein.links.detailed.v10.5.txt </b> Pairwise association scores for human proteins in STRING 10.5 (https://string-db.org).  

- <b> BIOGRID-ALL-3.4.152.tab2.txt </b> Pairwise interaction scores from BioGRID (https://thebiogrid.org).

- <b> nature22366-s2_Bioplex2.xlsx </b> This is Table S1 from the BioPlex 2.0 paper (Huttlin et al, Nature, 2017).






