## This script compares pairwise interactions in STRING, BioGRID and IntAct with those based on 
## protein co-regulation

# Load the necessary packages
library(ggplot2); library(data.table); library(VennDiagram); library(readxl); library(egg)

#### Prepare the STRING data ####

# Load the data
STRING <- fread("9606.protein.links.detailed.v10.5.txt")          # Load String databse
STRING[, protein1 := gsub("9606.", "", protein1, fixed = TRUE)]   # Remove species tag
STRING[, protein2 := gsub("9606.", "", protein2, fixed = TRUE)]   # Remove species tag

# It appears that all pairs are annotated as both A <-> B and B <-> A in this file.
# Keep only one of these duplicates (could just to A > B but this is the safer option)
STRING[, Protein_1_sorted := ifelse(protein1 > protein2, protein1, protein2) ]     # Sort the IDs row-wise
STRING[, Protein_2_sorted := ifelse(protein1 < protein2, protein1, protein2) ]     # Sort the IDs row-wise
STRING <- STRING[, .(protein1 = Protein_1_sorted, protein2 = Protein_2_sorted, neighborhood, fusion,
                     cooccurence, coexpression, experimental, database, textmining, combined_score )]
STRING <- unique(STRING)     # Keep only unique pairs

# To convert STRING's Ensembl protein IDs to Uniprot IDs, write them out and convert them using Uniprot
# website's Retrieve tool
write.csv( unique(c(STRING$protein1, STRING$protein2)), "all_string_ids.csv", row.names = FALSE )   # Write out IDs
uniprot_ID_mapping <- fread("uniprot_ID_mapping.tab")                                               # Read converted proteins IDs back in
colnames(uniprot_ID_mapping)[c(1,5)] <- c("Uniprot_ID", "Ensembl_ID")                               # Assign better column names
ambiguous <- names(which(table(uniprot_ID_mapping$Ensembl_ID) > 1))                                           # Ambiguously mapped IDs
uniprot_ID_mapping <- uniprot_ID_mapping[ !Ensembl_ID %in% ambiguous ]                                        # Remove ambiguous mappings
STRING[, uniprot1 := uniprot_ID_mapping$Uniprot_ID[ match(STRING$protein1, uniprot_ID_mapping$Ensembl_ID) ]]  # Add Uniprot ID for protein1
STRING[, uniprot2 := uniprot_ID_mapping$Uniprot_ID[ match(STRING$protein2, uniprot_ID_mapping$Ensembl_ID) ]]  # Add Uniprot ID for protein2
STRING <- STRING[ complete.cases(STRING) ]                                                                    # Delete pairs with missing mapping

# Re-order IDs to avoid pairs annotated as A <-> B and B <-> A for Uniprot IDs
STRING[, Protein_1_sorted := ifelse(uniprot1 > uniprot2, uniprot1, uniprot2) ]           # Sort uniprot1 and uniprot2
STRING[, Protein_2_sorted := ifelse(uniprot1 < uniprot2, uniprot1, uniprot2) ]           # Sort uniprot1 and uniprot2
STRING <- STRING[, .( SimpleID_1 = Protein_1_sorted, SimpleID_2 = Protein_2_sorted,      # Keep relevant columns (rename ID columns for merging later)
                      neighborhood, fusion, cooccurence, coexpression, experimental, database, textmining, combined_score) ]                 
STRING[,.N] == nrow(unique(STRING))                     # Double checking that there are no duplicates. Should be TRUE!
STRING[ SimpleID_1 > SimpleID_2 , .N] == STRING[, .N]   # Double checking that all pairs are in B <-> A order (alphabetically)
setkey(STRING, SimpleID_1, SimpleID_2)                  # Set keys for merging tables later


#### Prepare the BioGRID data ####

# Load the data
BioGRID <- fread("BIOGRID-ALL-3.4.152.tab2.txt")                                           # Load the BioGRID data
BioGRID <- BioGRID[ `Organism Interactor A` == 9606 & `Organism Interactor B` == 9606 ]    # Restrict to human interactions

# To convert Entrez IDs to Uniprot IDs, do the same as for String's Ensembl IDs
Entrez_IDs <- unique(c(BioGRID$`Entrez Gene Interactor A`, BioGRID$`Entrez Gene Interactor B`))   # List the unique IDs 
write.csv(Entrez_IDs, "BioGriD_Entrez_IDs.csv", row.names = FALSE)                                # Write out Entrez IDs and convert them to Uniprot IDs (reviewed) 
# using Uniprot's website Retrieve feature
uniprot_ID_mapping <- fread("entrez_to_uniprot.tab")                                              # Read converted proteins IDs back in
uniprot_ID_mapping <- uniprot_ID_mapping[, .(Entrez = `yourlist:M20170904AAFB7E4D2F1D05654627429E83DA5CCECCF1CD0`, Uniprot_ID = Entry)]  # Rename and select relevant columns (the "yourlist" column is named by Uniprot)
ambiguous <- names(which(table(uniprot_ID_mapping$Entrez) > 1))                                   # Ambiguously mapped IDs
uniprot_ID_mapping <- uniprot_ID_mapping[ !Entrez %in% ambiguous ]                                # Remove ambiguous mappings
BioGRID[, uniprot1 := uniprot_ID_mapping$Uniprot_ID[ match(BioGRID$`Entrez Gene Interactor A`, uniprot_ID_mapping$Entrez) ]]  # Add Uniprot ID for protein1
BioGRID[, uniprot2 := uniprot_ID_mapping$Uniprot_ID[ match(BioGRID$`Entrez Gene Interactor B`, uniprot_ID_mapping$Entrez) ]]  # Add Uniprot ID for protein2
BioGRID <- BioGRID[ complete.cases(BioGRID) ]       # Delete pairs with missing mapping
BioGRID <- BioGRID[ uniprot1 != uniprot2 ]          # Remove self-interactions

# Finalise BioGrid data pre-processing
BioGRID[, Protein_1_sorted := ifelse(uniprot1 > uniprot2, uniprot1, uniprot2) ]   # Sort interactor IDs alphabetically (to keep only B <-> A pairs)
BioGRID[, Protein_2_sorted := ifelse(uniprot1 < uniprot2, uniprot1, uniprot2) ]   # Sort interactor IDs alphabetically (to keep only B <-> A pairs)
BioGRID <- BioGRID[ `Experimental System Type` != "genetic" ]                     # Remove genetic interactions because there aren't enough to do statistics
BioGRID <- BioGRID[, .(SimpleID_1 = Protein_1_sorted, SimpleID_2 = Protein_2_sorted, `Experimental System`)]  # Keep only relevant columns
BioGRID <- unique(BioGRID)                    # Keep only unique rows
setkey(BioGRID, SimpleID_1, SimpleID_2)       # Set keys for merging tables later


#### Prepare the IntAct data ####

# Load the data (to reduce file size, the table downloaded from IntAct was pre-filtered to exlude interactions
# from non-human species, interactors that were not proteins, interactions without Uniprot ID as well as 
# self-interactions and duplicates. Finally, annotation strings were simplified, e.g. by reducing IDs from 
# "uniprotkb:P06748" to "P06748"
IntAct <- fread("IntAct_interactions.csv")
IntAct[, Protein_1_sorted := ifelse(Interactor_A_ID > Interactor_B_ID, Interactor_A_ID, Interactor_B_ID) ]     # Sort the IDs row-wise
IntAct[, Protein_2_sorted := ifelse(Interactor_A_ID < Interactor_B_ID, Interactor_A_ID, Interactor_B_ID) ]     # Sort the IDs row-wise
IntAct <- IntAct[, .(SimpleID_1 = Protein_1_sorted, SimpleID_2 = Protein_2_sorted, Expansion_method,
                     MI_score, Interaction_type, Interaction_detection_method ) ]
IntAct <- unique(IntAct)                  # Keep only unique pairs
setkey(IntAct, SimpleID_1, SimpleID_2)    # Set keys for merging tables later


#### Prepare the BioPlex data ####

# Read and pre-process Bioplex 2.0 (Table S1)
bioPlex <- read_xlsx("nature22366-s2_Bioplex2.xlsx", sheet = 2)
bioPlex <- as.data.table(bioPlex)
bioPlex <- bioPlex[, .(Protein_A = `Uniprot A`, Protein_B = `Uniprot B`)]
bioPlex <- bioPlex[ Protein_A != "UNKNOWN" & Protein_B != "UNKNOWN" ]

# Remove isoform annotation
bioPlex[, Protein_A := gsub("-.+", "", Protein_A)]
bioPlex[, Protein_B := gsub("-.+", "", Protein_B)]

# Sort and remove duplicates
bioPlex[, Protein_A_sorted := ifelse(Protein_A > Protein_B, Protein_A, Protein_B) ] # Sort the IDs row-wise
bioPlex[, Protein_B_sorted := ifelse(Protein_A < Protein_B, Protein_A, Protein_B) ] # Sort the IDs row-wise
bioPlex <- bioPlex[, .(SimpleID_1 = Protein_A_sorted, SimpleID_2 = Protein_B_sorted) ]
bioPlex[,.N] == nrow(unique(bioPlex))                             # There are no duplicates, true or false?
bioPlex$PPI <- "bioplex"
setkey(bioPlex, SimpleID_1, SimpleID_2)                           # Set keys for merging tables later


#### Prepare the treeClust data ####

# Load the coregulation scores
tC <- fread("coregulation_scores.csv")

# We arbitrarily define the highest scoring 0.5% of protein pairs as "co-regulated". Assign those. 
n_coreg <- floor( tC[,.N] * 0.005 )                              # How many pairs are these?
score_cut_off <- tC[ order(-coregulation_score)                  # What score cut-off does that correspond to?           
                     ][ n_coreg , coregulation_score ]
tC[ coregulation_score >= score_cut_off , coregulated := TRUE ]  # Assign the term
tC <- tC[ coregulated == TRUE ]                                  # Remove non-coregulated pairs

# Further prepping as necessary
tC[, SimpleID_1 := gsub(";.+", "", Protein_1)][, SimpleID_1 := gsub("-.+", "", SimpleID_1)]   # Simplify protein_1 IDs
tC[, SimpleID_2 := gsub(";.+", "", Protein_2)][, SimpleID_2 := gsub("-.+", "", SimpleID_2)]   # Simplify protein_1 IDs
tC <- tC[ SimpleID_1 != SimpleID_2 ]                      # Remove self-interactions after ID simplification (isoform removal)
tC <- tC[ !duplicated( tC[,.(SimpleID_1, SimpleID_2)] ) ] # Remove duplicates after ID simplification (isoform removal)
tC[ SimpleID_1 > SimpleID_2 , .N] == tC[, .N]             # Double checking that all pairs are in B <-> A order (alphabetically)
tC <- tC[, .(SimpleID_1, SimpleID_2, coregulated) ]       # Keep only relevant columns
setkey(tC, SimpleID_1, SimpleID_2)                        # Set keys for merging

# Clear workspace
rm( list = ls()[! ls() %in% c("tC", "STRING", "BioGRID", "IntAct", "bioPlex")] )


#### Analysis 1: Compare resource size ####

## How many unique genes / proteins are present in each dataset?

# Some IDs (particularly in IntAct) are not reviewed SwissProt but trEMBL proteins
# To remove those, first download all reviewed human SwissProt IDs and load them here
SwissProtIDs <- fread("SwissProt_IDs_HoSa_Nov_18.tab")  

# Now reduce each dataset on "real" (i.e. reviewed) IDs
 STRING <-  STRING[ SimpleID_1 %in% SwissProtIDs$Entry & SimpleID_2 %in% SwissProtIDs$Entry ]
BioGRID <- BioGRID[ SimpleID_1 %in% SwissProtIDs$Entry & SimpleID_2 %in% SwissProtIDs$Entry ]
 IntAct <-  IntAct[ SimpleID_1 %in% SwissProtIDs$Entry & SimpleID_2 %in% SwissProtIDs$Entry ]
bioPlex <- bioPlex[ SimpleID_1 %in% SwissProtIDs$Entry & SimpleID_2 %in% SwissProtIDs$Entry ]
     tC <-      tC[ SimpleID_1 %in% SwissProtIDs$Entry & SimpleID_2 %in% SwissProtIDs$Entry ]

# Now count the unique genes / proteins
 STRING_prots <-  STRING[ combined_score >= 400, length( unique( c( SimpleID_1, SimpleID_2 )))]
BioGRID_prots <- BioGRID[                      , length( unique( c( SimpleID_1, SimpleID_2 )))]
 IntAct_prots <-  IntAct[                      , length( unique( c( SimpleID_1, SimpleID_2 )))]
bioPlex_prots <- bioPlex[                      , length( unique( c( SimpleID_1, SimpleID_2 )))]
     tC_prots <-      tC[                      , length( unique( c( SimpleID_1, SimpleID_2 )))]

## How many interactions - on average - for these genes?
 avg_int_STRING <- nrow( unique(  STRING[ combined_score >= 400,.(SimpleID_1, SimpleID_2)] )) / STRING_prots
avg_int_BioGRID <- nrow( unique( BioGRID[                      ,.(SimpleID_1, SimpleID_2)] )) / BioGRID_prots
 avg_int_IntAct <- nrow( unique(  IntAct[                      ,.(SimpleID_1, SimpleID_2)] )) / IntAct_prots
avg_int_bioPlex <- nrow( unique( bioPlex[                      ,.(SimpleID_1, SimpleID_2)] )) / bioPlex_prots
     avg_int_tC <- nrow( unique(      tC[                      ,.(SimpleID_1, SimpleID_2)] )) / tC_prots

# Combine these data into one table
resource_size <- data.table(resource = c("STRING", "BioGRID", "IntAct", "BioPlex2", "ProHD"), 
                            n_genes = c(STRING_prots, BioGRID_prots, IntAct_prots, bioPlex_prots, tC_prots ), 
                            avg_n_interactions = c(avg_int_STRING, avg_int_BioGRID, avg_int_IntAct, avg_int_bioPlex, avg_int_tC))     


#### Analysis 2: Compare associations for the overlapping genes - whole resources ####

## (A) STRING vs ProHD coregulation

# Find the proteins that are covered both by treeClust and STRING
tC_proteins <-     tC[, unique(c(SimpleID_1, SimpleID_2)) ]
STRING_proteins <- STRING[, unique(c(SimpleID_1, SimpleID_2)) ]
overlap_tCSTprot <- tC_proteins[ tC_proteins %in% STRING_proteins ]

# Full outer merge for proteins they both have in common 
tC_STRING <- merge(     tC[ SimpleID_1 %in% overlap_tCSTprot & SimpleID_2 %in% overlap_tCSTprot ],
                    STRING[ SimpleID_1 %in% overlap_tCSTprot & SimpleID_2 %in% overlap_tCSTprot ],
                        all = TRUE)

# Get the number of `unique pairs` that are coregulated, STRING-associated, or both
    tC_STRING_coreg <- nrow( unique( tC_STRING[ coregulated == TRUE                        , .(SimpleID_1, SimpleID_2) ]))
tC_STRING_assoc_all <- nrow( unique( tC_STRING[                       combined_score >= 150, .(SimpleID_1, SimpleID_2) ]))
 tC_STRING_both_all <- nrow( unique( tC_STRING[ coregulated == TRUE & combined_score >= 150, .(SimpleID_1, SimpleID_2) ]))

pc_coreg_STR_all <- tC_STRING_both_all / tC_STRING_coreg * 100      # Which % of coregulated protein pairs is also found associated in STRING?
pc_STR_coreg_all <- tC_STRING_both_all / tC_STRING_assoc_all * 100  # Which % of STRING associations also found by co-regulation?


## Same for low-medium confidence interactions in STRING
tC_STRING_assoc_250 <- nrow( unique( tC_STRING[                       combined_score >= 250, .(SimpleID_1, SimpleID_2) ]))
 tC_STRING_both_250 <- nrow( unique( tC_STRING[ coregulated == TRUE & combined_score >= 250, .(SimpleID_1, SimpleID_2) ]))
pc_coreg_STR_250 <- tC_STRING_both_250 / tC_STRING_coreg * 100
pc_STR_coreg_250 <- tC_STRING_both_250 / tC_STRING_assoc_250 * 100


## Same for medium confidence interactions in STRING
tC_STRING_assoc_400 <- nrow( unique( tC_STRING[                       combined_score >= 400, .(SimpleID_1, SimpleID_2) ]))
 tC_STRING_both_400 <- nrow( unique( tC_STRING[ coregulated == TRUE & combined_score >= 400, .(SimpleID_1, SimpleID_2) ]))
pc_coreg_STR_400 <- tC_STRING_both_400 / tC_STRING_coreg * 100
pc_STR_coreg_400 <- tC_STRING_both_400 / tC_STRING_assoc_400 * 100


## Same for medium-high confidence interactions in STRING
tC_STRING_assoc_600 <- nrow( unique( tC_STRING[                       combined_score >= 600, .(SimpleID_1, SimpleID_2) ]))
 tC_STRING_both_600 <- nrow( unique( tC_STRING[ coregulated == TRUE & combined_score >= 600, .(SimpleID_1, SimpleID_2) ]))
pc_coreg_STR_600 <- tC_STRING_both_600 / tC_STRING_coreg * 100
pc_STR_coreg_600 <- tC_STRING_both_600 / tC_STRING_assoc_600 * 100


## Same for highest confidence interactions in STRING
tC_STRING_assoc_900 <- nrow( unique( tC_STRING[                       combined_score >= 900, .(SimpleID_1, SimpleID_2) ]))
 tC_STRING_both_900 <- nrow( unique( tC_STRING[ coregulated == TRUE & combined_score >= 900, .(SimpleID_1, SimpleID_2) ]))
pc_coreg_STR_900 <- tC_STRING_both_900 / tC_STRING_coreg * 100
pc_STR_coreg_900 <- tC_STRING_both_900 / tC_STRING_assoc_900 * 100


## (B) BioGRID vs ProHD coregulation

# Find the proteins that are covered both by treeClust and BioGRID
BioGRID_proteins <- BioGRID[, unique(c(SimpleID_1, SimpleID_2)) ]
overlap_tCBIprot <- tC_proteins[ tC_proteins %in% BioGRID_proteins ]

# Full outer merge for proteins they both have in common 
tC_BioGRID <- merge(     tC[ SimpleID_1 %in% overlap_tCBIprot & SimpleID_2 %in% overlap_tCBIprot ],
                    BioGRID[ SimpleID_1 %in% overlap_tCBIprot & SimpleID_2 %in% overlap_tCBIprot ],
                         all = TRUE)

# Get the number of `unique pairs` that are coregulated, BioGRID-associated, or both
tC_BioGRID_coreg <- nrow( unique( tC_BioGRID[ coregulated == TRUE                                , .(SimpleID_1, SimpleID_2) ]))
tC_BioGRID_assoc <- nrow( unique( tC_BioGRID[                       !is.na(`Experimental System`), .(SimpleID_1, SimpleID_2) ]))
 tC_BioGRID_both <- nrow( unique( tC_BioGRID[ coregulated == TRUE & !is.na(`Experimental System`), .(SimpleID_1, SimpleID_2) ]))

# Which % of coregulated protein pairs is also found associated in BioGRID?
pc_coreg_BIO <- tC_BioGRID_both / tC_BioGRID_coreg * 100

# Which % of BioGRID associations also found by co-regulation?
pc_BIO_coreg <- tC_BioGRID_both / tC_BioGRID_assoc * 100


## (C) IntAct vs ProHD coregulation

# Find the proteins that are covered both by treeClust and IntAct
IntAct_proteins <- IntAct[, unique(c(SimpleID_1, SimpleID_2)) ]
overlap_tCINprot <- tC_proteins[ tC_proteins %in% IntAct_proteins ]

# Full outer merge for proteins they both have in common 
tC_IntAct <- merge(    tC[ SimpleID_1 %in% overlap_tCINprot & SimpleID_2 %in% overlap_tCINprot ],
                       IntAct[ SimpleID_1 %in% overlap_tCINprot & SimpleID_2 %in% overlap_tCINprot ],
                       all = TRUE)

# Get the number of `unique pairs` that are coregulated, IntAct-associated, or both
tC_IntAct_coreg <- nrow( unique( tC_IntAct[ coregulated == TRUE                           , .(SimpleID_1, SimpleID_2) ]))
tC_IntAct_assoc <- nrow( unique( tC_IntAct[                       !is.na(Interaction_type), .(SimpleID_1, SimpleID_2) ]))
 tC_IntAct_both <- nrow( unique( tC_IntAct[ coregulated == TRUE & !is.na(Interaction_type), .(SimpleID_1, SimpleID_2) ]))

# Which % of coregulated protein pairs is also found associated in IntAct?
pc_coreg_INT <- tC_IntAct_both / tC_IntAct_coreg * 100

# Which % of IntAct associations also found by co-regulation?
pc_INT_coreg <- tC_IntAct_both / tC_IntAct_assoc * 100


## Combine these data into one table
overlap_all <- data.table( resource = c("STRING", "STRING_250", "STRING_400", "STRING_600", "STRING_900", "BioGRID", "IntAct"),
                           type = "all",
                           y = c( pc_coreg_STR_all, pc_coreg_STR_250, pc_coreg_STR_400, pc_coreg_STR_600, pc_coreg_STR_900, pc_coreg_BIO, pc_coreg_INT),
                           x = c( pc_STR_coreg_all, pc_STR_coreg_250, pc_STR_coreg_400, pc_STR_coreg_600, pc_STR_coreg_900, pc_BIO_coreg, pc_INT_coreg))


#### Analysis 3: Compare associations for the overlapping genes - breakdown by type ####

## (A) STRING: Get the number of `unique pairs` as before, but broken down by evidence channel
tC_STRING_neighborhood_assoc <- nrow( unique( tC_STRING[                       neighborhood >= 150, .(SimpleID_1, SimpleID_2) ]))
 tC_STRING_neighborhood_both <- nrow( unique( tC_STRING[ coregulated == TRUE & neighborhood >= 150, .(SimpleID_1, SimpleID_2) ]))
pc_coreg_STR_neighborhood <- tC_STRING_neighborhood_both / tC_STRING_coreg * 100               # % coregulated pairs also found associated by neighborhood?
pc_STR_coreg_neighborhood <- tC_STRING_neighborhood_both / tC_STRING_neighborhood_assoc * 100  # % neighborhood associations also found by co-regulation?

tC_STRING_fusion_assoc <- nrow( unique( tC_STRING[                       fusion >= 150, .(SimpleID_1, SimpleID_2) ]))
 tC_STRING_fusion_both <- nrow( unique( tC_STRING[ coregulated == TRUE & fusion >= 150, .(SimpleID_1, SimpleID_2) ]))
pc_coreg_STR_fusion <- tC_STRING_fusion_both / tC_STRING_coreg * 100         # % coregulated pairs also found associated by fusion?
pc_STR_coreg_fusion <- tC_STRING_fusion_both / tC_STRING_fusion_assoc * 100  # % fusion associations also found by co-regulation?

tC_STRING_cooccurence_assoc <- nrow( unique( tC_STRING[                       cooccurence >= 150, .(SimpleID_1, SimpleID_2) ]))
 tC_STRING_cooccurence_both <- nrow( unique( tC_STRING[ coregulated == TRUE & cooccurence >= 150, .(SimpleID_1, SimpleID_2) ]))
pc_coreg_STR_cooccurence <- tC_STRING_cooccurence_both / tC_STRING_coreg * 100              # % coregulated pairs also found associated by cooccurence?
pc_STR_coreg_cooccurence <- tC_STRING_cooccurence_both / tC_STRING_cooccurence_assoc * 100  # % cooccurence associations also found by co-regulation?

tC_STRING_coexpression_assoc <- nrow( unique( tC_STRING[                       coexpression >= 150, .(SimpleID_1, SimpleID_2) ]))
 tC_STRING_coexpression_both <- nrow( unique( tC_STRING[ coregulated == TRUE & coexpression >= 150, .(SimpleID_1, SimpleID_2) ]))
pc_coreg_STR_coexpression <- tC_STRING_coexpression_both / tC_STRING_coreg * 100               # % coregulated pairs also found associated by coexpression?
pc_STR_coreg_coexpression <- tC_STRING_coexpression_both / tC_STRING_coexpression_assoc * 100  # % coexpression associations also found by co-regulation?

tC_STRING_experimental_assoc <- nrow( unique( tC_STRING[                       experimental >= 150, .(SimpleID_1, SimpleID_2) ]))
 tC_STRING_experimental_both <- nrow( unique( tC_STRING[ coregulated == TRUE & experimental >= 150, .(SimpleID_1, SimpleID_2) ]))
pc_coreg_STR_experimental <- tC_STRING_experimental_both / tC_STRING_coreg * 100               # % coregulated pairs also found associated by experimental?
pc_STR_coreg_experimental <- tC_STRING_experimental_both / tC_STRING_experimental_assoc * 100  # % experimental associations also found by co-regulation?

tC_STRING_database_assoc <- nrow( unique( tC_STRING[                       database >= 150, .(SimpleID_1, SimpleID_2) ]))
 tC_STRING_database_both <- nrow( unique( tC_STRING[ coregulated == TRUE & database >= 150, .(SimpleID_1, SimpleID_2) ]))
pc_coreg_STR_database <- tC_STRING_database_both / tC_STRING_coreg * 100           # % coregulated pairs also found associated by database?
pc_STR_coreg_database <- tC_STRING_database_both / tC_STRING_database_assoc * 100  # % database associations also found by co-regulation?

tC_STRING_textmining_assoc <- nrow( unique( tC_STRING[                       textmining >= 150, .(SimpleID_1, SimpleID_2) ]))
 tC_STRING_textmining_both <- nrow( unique( tC_STRING[ coregulated == TRUE & textmining >= 150, .(SimpleID_1, SimpleID_2) ]))
pc_coreg_STR_textmining <- tC_STRING_textmining_both / tC_STRING_coreg * 100             # % coregulated pairs also found associated by textmining?
pc_STR_coreg_textmining <- tC_STRING_textmining_both / tC_STRING_textmining_assoc * 100  # % textmining associations also found by co-regulation?

# Combine these data into one table
STRING_breakdown <- data.table( resource = "STRING",
                                type = c("neighborhood", "fusion", "cooccurence", "coexpression", "experimental", "database", "textmining"),
                           y = c( pc_coreg_STR_neighborhood, pc_coreg_STR_fusion, pc_coreg_STR_cooccurence, pc_coreg_STR_coexpression, pc_coreg_STR_experimental, pc_coreg_STR_database, pc_coreg_STR_textmining),
                           x = c( pc_STR_coreg_neighborhood, pc_STR_coreg_fusion, pc_STR_coreg_cooccurence, pc_STR_coreg_coexpression, pc_STR_coreg_experimental, pc_STR_coreg_database, pc_STR_coreg_textmining))


## (B) BioGRID: Get the number of `unique pairs` as before, but broken down by evidence channel
tC_BioGRID_AffinityCapture_assoc <- nrow( unique( tC_BioGRID[                       `Experimental System` %like% "Affinity Capture", .(SimpleID_1, SimpleID_2) ]))
 tC_BioGRID_AffinityCapture_both <- nrow( unique( tC_BioGRID[ coregulated == TRUE & `Experimental System` %like% "Affinity Capture", .(SimpleID_1, SimpleID_2) ]))
pc_coreg_BIO_AffinityCapture <- tC_BioGRID_AffinityCapture_both / tC_BioGRID_coreg * 100                  # % coregulated pairs also found associated by AffinityCapture?
pc_BIO_coreg_AffinityCapture <- tC_BioGRID_AffinityCapture_both / tC_BioGRID_AffinityCapture_assoc * 100  # % AffinityCapture associations also found by co-regulation?

tC_BioGRID_CoFractionation_assoc <- nrow( unique( tC_BioGRID[                       `Experimental System` %like% "Co-fractionation", .(SimpleID_1, SimpleID_2) ]))
 tC_BioGRID_CoFractionation_both <- nrow( unique( tC_BioGRID[ coregulated == TRUE & `Experimental System` %like% "Co-fractionation", .(SimpleID_1, SimpleID_2) ]))
pc_coreg_BIO_CoFractionation <- tC_BioGRID_CoFractionation_both / tC_BioGRID_coreg * 100                  # % coregulated pairs also found associated by CoFractionation?
pc_BIO_coreg_CoFractionation <- tC_BioGRID_CoFractionation_both / tC_BioGRID_CoFractionation_assoc * 100  # % CoFractionation associations also found by co-regulation?

tC_BioGRID_ReconstitutedComplex_assoc <- nrow( unique( tC_BioGRID[                       `Experimental System` %like% "Reconstituted Complex", .(SimpleID_1, SimpleID_2) ]))
 tC_BioGRID_ReconstitutedComplex_both <- nrow( unique( tC_BioGRID[ coregulated == TRUE & `Experimental System` %like% "Reconstituted Complex", .(SimpleID_1, SimpleID_2) ]))
pc_coreg_BIO_ReconstitutedComplex <- tC_BioGRID_ReconstitutedComplex_both / tC_BioGRID_coreg * 100                  # % coregulated pairs also found associated by ReconstitutedComplex?
pc_BIO_coreg_ReconstitutedComplex <- tC_BioGRID_ReconstitutedComplex_both / tC_BioGRID_ReconstitutedComplex_assoc * 100  # % ReconstitutedComplex associations also found by co-regulation?

tC_BioGRID_TwoHybrid_assoc <- nrow( unique( tC_BioGRID[                       `Experimental System` %like% "Two-hybrid", .(SimpleID_1, SimpleID_2) ]))
 tC_BioGRID_TwoHybrid_both <- nrow( unique( tC_BioGRID[ coregulated == TRUE & `Experimental System` %like% "Two-hybrid", .(SimpleID_1, SimpleID_2) ]))
pc_coreg_BIO_TwoHybrid <- tC_BioGRID_TwoHybrid_both / tC_BioGRID_coreg * 100                  # % coregulated pairs also found associated by TwoHybrid?
pc_BIO_coreg_TwoHybrid <- tC_BioGRID_TwoHybrid_both / tC_BioGRID_TwoHybrid_assoc * 100  # % TwoHybrid associations also found by co-regulation?

# Combine these data into one table
BioGRID_breakdown <- data.table( resource = "BioGRID",
                                type = c("Affinity Capture", "Co-fractionation", "in vitro", "two-hybrid"),
                                y = c( pc_coreg_BIO_AffinityCapture, pc_coreg_BIO_CoFractionation, pc_coreg_BIO_ReconstitutedComplex, pc_coreg_BIO_TwoHybrid),
                                x = c( pc_BIO_coreg_AffinityCapture, pc_BIO_coreg_CoFractionation, pc_BIO_coreg_ReconstitutedComplex, pc_BIO_coreg_TwoHybrid))


#### Plot the results ####

## Resource size

# Set plotting order
resource_size[, resource := factor(resource, levels = c("ProHD", "BioPlex2", "IntAct", "BioGRID", "STRING" )) ]

# Create the two barcharts
p1a <- ggplot(resource_size, aes( x = resource, y = n_genes ))+
         geom_bar( stat = "identity" , fill = "black")+
         xlab("Resource")+
         ylab("Genes covered")+
         coord_flip()+
         theme(panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "black", size=0.25),
               panel.grid = element_blank(), axis.text=element_text(size=5), axis.title=element_text(size=6), 
               axis.ticks.y = element_blank(), axis.ticks.x = element_line(size=0.25))
  
p1b <- ggplot(resource_size, aes( x = resource, y = avg_n_interactions ))+
        geom_bar( stat = "identity" , fill = "black")+
        xlab("Resource")+
        ylab("avg. interactions per gene")+
        coord_flip()+
        theme(panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "black", size=0.25),
              panel.grid = element_blank(), axis.text=element_text(size=5), axis.title=element_text(size=6), 
              axis.ticks.y = element_blank(), axis.ticks.x = element_line(size=0.25))

# Combine & align charts, then save
p1 <- ggarrange(p1a, p1b, nrow = 1)
ggsave("Resource_size.pdf", p1, width = 6, height = 3, units = "cm")


## Whole resources (databases only) 
p2a <- ggplot(overlap_all, aes( x = x, y = y, colour = resource))+
        geom_point()+
        geom_text(aes(label = resource), size=2)+
        scale_x_continuous( limits = c(0,23), breaks=seq(0,50,10))+
        scale_y_continuous( limits = c(0,40), breaks=seq(0,50,10))+
        xlab("[%] PPIs in resource found by co-regulation")+
        ylab("[%] co-regulated protein pairs\nfound by other methods")+
        theme(panel.background = element_rect(fill = "white", colour = NA), panel.border = element_rect(fill = NA, colour = "grey20", size=0.25),
              panel.grid.major = element_line(colour = "grey92", size=0.25), panel.grid.minor = element_line(colour = "grey92", size = 0.25),
              axis.text=element_text(size=6), axis.title=element_text(size=6), axis.ticks = element_line(size=0.25),
              legend.position = "none")

## Breakdown by type of interactions
breakdown <- rbind(STRING_breakdown, BioGRID_breakdown)

p2b <- ggplot(breakdown, aes( x = resource, y = y, colour = resource))+
        geom_point()+
        geom_text(aes(label = type), size=2, hjust = 0)+
        scale_y_continuous( limits = c(0,25), breaks=seq(0,50,5))+
        scale_colour_manual(values = c("firebrick2", "dodgerblue2"))+
        ylab("[%] co-regulated protein pairs\nfound by other methods")+
        theme(panel.background = element_rect(fill = "white", colour = NA), panel.border = element_rect(fill = NA, colour = "grey20", size=0.25),
            panel.grid.major = element_line(colour = "grey92", size=0.25), panel.grid.minor = element_line(colour = "grey92", size = 0.25),
            axis.text=element_text(size=6), axis.title=element_text(size=6), axis.ticks = element_line(size=0.25),
            legend.position = "none", axis.title.x = element_blank())

# Combine & align charts, then save
p2 <- ggarrange(p2a, p2b, nrow = 1)
ggsave("Resource_overlap.pdf", p2, width = 10, height = 4.7, units = "cm")


#### Analysis 4: Compare associations for the overlapping genes - BioPlex 2.0 ####

rm( list = ls()[! ls() %in% c("tC", "STRING", "BioGRID", "IntAct", "bioPlex")] )

## Venn diagram

# Find the proteins that are covered both by treeClust and BioPlex 2.0
tC_proteins <- tC[, unique(c(SimpleID_1, SimpleID_2)) ]               # Genes covered by us
bioPlex_proteins <- bioPlex[, unique(c(SimpleID_1, SimpleID_2)) ]     # Genes covered by Bioplex
   overlap_prots <- tC_proteins[ tC_proteins %in% bioPlex_proteins ]

# Full outer merge for proteins they both have in common 
tC_bioPlex <- merge(     tC[ SimpleID_1 %in% overlap_prots & SimpleID_2 %in% overlap_prots ],
                    bioPlex[ SimpleID_1 %in% overlap_prots & SimpleID_2 %in% overlap_prots ],
                        all = TRUE)

# Get the number of `unique pairs` that are coregulated, bioPlex-associated, or both
    tC_bioPlex_coreg <- nrow( unique( tC_bioPlex[ coregulated == TRUE              , .(SimpleID_1, SimpleID_2) ]))
tC_bioPlex_assoc_all <- nrow( unique( tC_bioPlex[                       !is.na(PPI), .(SimpleID_1, SimpleID_2) ]))
 tC_bioPlex_both_all <- nrow( unique( tC_bioPlex[ coregulated == TRUE & !is.na(PPI), .(SimpleID_1, SimpleID_2) ]))

grid.newpage()
pdf("Venn_Coreg_BioPlex.pdf", width = 3, height = 3)
draw.pairwise.venn( area1 = tC_bioPlex_coreg, area2 = tC_bioPlex_assoc_all, cross.area = tC_bioPlex_both_all,
                    category = c("co-regulated", "BioPlex 2"))
dev.off()


## STRING enrichment; How many of the associations unique to our study are backed up by STRING?

# Which associations are found by us but not BioPlex (for the overlapping genes)?
our_nonBioPlex_PPIs <- unique( tC_bioPlex[ coregulated == TRUE & is.na(PPI) , .(SimpleID_1, SimpleID_2) ])

# How many of our unique PPIs are also captured by STRING?
in_string <- nrow( merge(our_nonBioPlex_PPIs, STRING) )

# How does this compare to a random sample of the possible interactions between the genes we cover?
our_genes <- tC[, unique( c( SimpleID_1, SimpleID_2)) ]                 # For these genes we have detected PPIs
all_possible_combis <- as.data.table(t(combn(our_genes, 2)))            # Which associations *could* there be between these genes
all_possible_combis[, Protein_1_sorted := ifelse(V1 > V2, V1, V2) ]     # Sort the IDs row-wise
all_possible_combis[, Protein_2_sorted := ifelse(V1 < V2, V1, V2) ]     # Sort the IDs row-wise
all_possible_combis <- all_possible_combis[, .(SimpleID_1 = Protein_1_sorted,
                                               SimpleID_2 = Protein_2_sorted)]
random_combis <- all_possible_combis[ sample(.N, our_nonBioPlex_PPIs[,.N] )]   # Randomly sample the same number of pairs as we actually detect
setkey(random_combis, SimpleID_1, SimpleID_2)
ran_in_string <- nrow( merge(random_combis, STRING) )                          # How many random pairs are in STRING

# Plot the results
plot_dt <- data.table( coregulated_pairs = in_string / our_nonBioPlex_PPIs[,.N] * 100 ,
                       random_pairs = ran_in_string / random_combis[,.N] * 100        )
plot_dt <- melt(plot_dt, measure.vars = 1:2 )
plot_dt[, variable := factor(variable, levels = c("random_pairs", "coregulated_pairs")) ]   # Set plotting order

p3 <- ggplot(plot_dt, aes(x = variable, y = value))+
      geom_bar(stat="identity")+
      ylab("protein pairs [%]")+
      coord_flip()+
      theme(panel.background = element_rect(fill = "white", colour = NA), panel.grid = element_blank(),
            axis.text=element_text(size=6), axis.title=element_text(size=6), axis.ticks.x = element_line(size=0.25),
            axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.line.x = element_line(size=0.25))
      
# Save the plot
ggsave("In_string_or_not.pdf", p3, width=4, height=2, units=c("cm"))



