## This script looks at the connectivity of microproteins

# Load the necessary packages
library(ggplot2); library(data.table); library(cowplot)

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


#### Annotate protein size ####

# Load Uniprot annotation data (all reviewed human SwissProt protein IDs where downloaded with 
# Annotation score and molecular weight)
uni <- fread("SwissProt_HS_AnnotScore_Mass.tab")                                  # Load data
uni <- uni[, .(Protein_ID = Entry, Annotation_score = Annotation, Mass)]          # Keep and rename relevant columns
uni[, Mass := gsub(",", "", Mass, fixed = TRUE)][, Mass := as.numeric(Mass)]      # Make Mass a numeric feature
uni[, Annotation_score := as.integer( gsub(" out of 5", "", Annotation_score)) ]  # Simplify scores

     
#### Analysis 1: Connectivity of microproteins vs larger proteins ####

# Subset STRING to very reliable interactions
STRING900 <- STRING[ combined_score > 900 ]

# Get the unique genes / proteins in each dataset
STRING_prots <-  STRING900[ , unique( c( SimpleID_1, SimpleID_2 ))]
BioGRID_prots <- BioGRID[   , unique( c( SimpleID_1, SimpleID_2 ))]
tC_prots <-      tC[        , unique( c( SimpleID_1, SimpleID_2 ))]

# Get STRING900 connectivity
STRING900_connectivity <- integer()
for(i in STRING_prots){ STRING900_connectivity <- c(STRING900_connectivity, STRING900[ SimpleID_1 == i | SimpleID_2 == i , .N ]) }

# Get BioGRID connectivity
BioGRID_connectivity <- integer()
for(i in BioGRID_prots){ BioGRID_connectivity <- c(BioGRID_connectivity, BioGRID[ SimpleID_1 == i | SimpleID_2 == i , .N ]) }

# Get treeClust connectivity
tC_connectivity <- integer()
for(i in tC_prots){ tC_connectivity <- c(tC_connectivity, tC[ SimpleID_1 == i | SimpleID_2 == i , .N ]) }

# Combine data into one table
DT <- rbind( data.table( Protein_ID = tC_prots     , connectivity = tC_connectivity       , type = "tC"),
             data.table( Protein_ID = STRING_prots , connectivity = STRING900_connectivity, type = "STRING"),
             data.table( Protein_ID = BioGRID_prots, connectivity = BioGRID_connectivity  , type = "BioGRID"))

# Annotate with protein size 
DT <- merge(DT, uni, by = "Protein_ID")

# Plot the results
ptC <- ggplot(DT[ type == "tC" ], aes( y = connectivity, fill = Mass > 15000))+
        geom_boxplot( outlier.colour = NA , notch = TRUE , size = 0.25)+
        coord_cartesian( ylim = c(0,280))+
        scale_fill_manual( values = c("cornflowerblue", "firebrick"))+
        ylab("Number of interaction partners")+
        xlab("ProteomeHD")+
        scale_y_continuous( breaks = seq(0,300,50))+
        theme(panel.background = element_blank(), panel.grid = element_blank(),
              panel.border = element_rect(fill = NA, colour = "black", size = 0.25, linetype = "solid"),    
              axis.text.y=element_text(size=5), axis.text.x=element_blank(),
              axis.title=element_text(size=6), axis.ticks.y = element_line(size=0.25), axis.ticks.x = element_blank(),
              axis.line = element_blank(), legend.position = "none")

pSTRING <- ggplot(DT[ type == "STRING" ], aes( y = connectivity, fill = Mass > 15000))+
            geom_boxplot( outlier.colour = NA , notch = TRUE , size = 0.25)+
            coord_cartesian( ylim = c(0,250))+
            scale_fill_manual( values = c("cornflowerblue", "firebrick"))+
            ylab("Number of interaction partners")+
            xlab("STRING")+
            scale_y_continuous( breaks = seq(0,300,50))+
            theme(panel.background = element_blank(), panel.grid = element_blank(),
                  panel.border = element_rect(fill = NA, colour = "black", size = 0.25, linetype = "solid"),    
                  axis.text.y=element_text(size=5), axis.text.x=element_blank(),
                  axis.title=element_text(size=6), axis.ticks.y = element_line(size=0.25), axis.ticks.x = element_blank(),
                  axis.line = element_blank(), legend.position = "none")

pBioGRID <- ggplot(DT[ type == "BioGRID" ], aes( y = connectivity, fill = Mass > 15000))+
              geom_boxplot( outlier.colour = NA , notch = TRUE , size = 0.25)+
              coord_cartesian( ylim = c(0,110))+
              scale_fill_manual( values = c("cornflowerblue", "firebrick"))+
              ylab("Number of interaction partners")+
              xlab("BioGRID")+
              scale_y_continuous( breaks = seq(0,300,25))+
              theme(panel.background = element_blank(), panel.grid = element_blank(),
                    panel.border = element_rect(fill = NA, colour = "black", size = 0.25, linetype = "solid"),    
                    axis.text.y=element_text(size=5), axis.text.x=element_blank(),
                    axis.title=element_text(size=6), axis.ticks.y = element_line(size=0.25), axis.ticks.x = element_blank(),
                    axis.line = element_blank(), legend.position = "none")
            

#### Analysis 2: Experiment types that predict microprotein functions ####

# These are the microproteins in Uniprot
uniprot_MPs <- uni[ Mass <= 15000 , Protein_ID ]

# Subset STRING to keep only interactions that involve at least one microprotein
STRING_mp <- STRING[ SimpleID_1 %in% uniprot_MPs | SimpleID_2 %in% uniprot_MPs ]

# Subset STRING further based on evidence type
STRING_COEX <- STRING_mp[ coexpression > 400 ]
STRING_EXPE <- STRING_mp[ experimental > 400 ]
STRING_DATA <- STRING_mp[     database > 400 ]
STRING_TEXT <- STRING_mp[   textmining > 400 ]

# As above, find out how many interaction partners the different STRING evidences provide per MP
STRING_COEX_MPs <- STRING_COEX[ , unique( c( SimpleID_1, SimpleID_2 ))]   # These are the proteins found in the subsetted STRING version
STRING_EXPE_MPs <- STRING_EXPE[ , unique( c( SimpleID_1, SimpleID_2 ))]
STRING_DATA_MPs <- STRING_DATA[ , unique( c( SimpleID_1, SimpleID_2 ))]
STRING_TEXT_MPs <- STRING_TEXT[ , unique( c( SimpleID_1, SimpleID_2 ))]

STRING_COEX_MPs <- STRING_COEX_MPs[ STRING_COEX_MPs %in% uniprot_MPs ]    # And these are the involved microproteins that I need to query
STRING_EXPE_MPs <- STRING_EXPE_MPs[ STRING_EXPE_MPs %in% uniprot_MPs ] 
STRING_DATA_MPs <- STRING_DATA_MPs[ STRING_DATA_MPs %in% uniprot_MPs ] 
STRING_TEXT_MPs <- STRING_TEXT_MPs[ STRING_TEXT_MPs %in% uniprot_MPs ] 

STRING_COEX_connectivity <- integer()
for(i in STRING_COEX_MPs){ STRING_COEX_connectivity <- c(STRING_COEX_connectivity, STRING_COEX[ SimpleID_1 == i | SimpleID_2 == i , .N ]) }

STRING_EXPE_connectivity <- integer()
for(i in STRING_EXPE_MPs){ STRING_EXPE_connectivity <- c(STRING_EXPE_connectivity, STRING_EXPE[ SimpleID_1 == i | SimpleID_2 == i , .N ]) }

STRING_DATA_connectivity <- integer()
for(i in STRING_DATA_MPs){ STRING_DATA_connectivity <- c(STRING_DATA_connectivity, STRING_DATA[ SimpleID_1 == i | SimpleID_2 == i , .N ]) }

STRING_TEXT_connectivity <- integer()
for(i in STRING_TEXT_MPs){ STRING_TEXT_connectivity <- c(STRING_TEXT_connectivity, STRING_TEXT[ SimpleID_1 == i | SimpleID_2 == i , .N ]) }


## Now the same for BioGrid data ...

# Subset BioGRID to keep only interactions that involve at least one microprotein
BioGRID_mp <- BioGRID[ SimpleID_1 %in% uniprot_MPs | SimpleID_2 %in% uniprot_MPs ]

# Subset BioGRID based on evidence type
BioGRID_AFCA <- BioGRID_mp[ `Experimental System` %like% "Affinity Capture" ]
BioGRID_COFR <- BioGRID_mp[ `Experimental System` %like% "Co-fractionation" ]
BioGRID_RECO <- BioGRID_mp[ `Experimental System` %like% "Reconstituted Complex" ]
BioGRID_Y2HY <- BioGRID_mp[ `Experimental System` %like% "Two-hybrid" ]

# As above, find out how many interaction partners the different BioGRID evidences provide per MP
BioGRID_AFCA_MPs <- BioGRID_AFCA[ , unique( c( SimpleID_1, SimpleID_2 ))]  # These are the proteins found in the subsetted BioGRID version
BioGRID_COFR_MPs <- BioGRID_COFR[ , unique( c( SimpleID_1, SimpleID_2 ))]
BioGRID_RECO_MPs <- BioGRID_RECO[ , unique( c( SimpleID_1, SimpleID_2 ))]
BioGRID_Y2HY_MPs <- BioGRID_Y2HY[ , unique( c( SimpleID_1, SimpleID_2 ))]

BioGRID_AFCA_MPs <- BioGRID_AFCA_MPs[ BioGRID_AFCA_MPs %in% uniprot_MPs ]  # And these are the involved microproteins that I need to query
BioGRID_COFR_MPs <- BioGRID_COFR_MPs[ BioGRID_COFR_MPs %in% uniprot_MPs ] 
BioGRID_RECO_MPs <- BioGRID_RECO_MPs[ BioGRID_RECO_MPs %in% uniprot_MPs ] 
BioGRID_Y2HY_MPs <- BioGRID_Y2HY_MPs[ BioGRID_Y2HY_MPs %in% uniprot_MPs ] 

BioGRID_AFCA_connectivity <- integer()
for(i in BioGRID_AFCA_MPs){ BioGRID_AFCA_connectivity <- c(BioGRID_AFCA_connectivity, BioGRID_AFCA[ SimpleID_1 == i | SimpleID_2 == i , .N ]) }

BioGRID_COFR_connectivity <- integer()
for(i in BioGRID_COFR_MPs){ BioGRID_COFR_connectivity <- c(BioGRID_COFR_connectivity, BioGRID_COFR[ SimpleID_1 == i | SimpleID_2 == i , .N ]) }

BioGRID_RECO_connectivity <- integer()
for(i in BioGRID_RECO_MPs){ BioGRID_RECO_connectivity <- c(BioGRID_RECO_connectivity, BioGRID_RECO[ SimpleID_1 == i | SimpleID_2 == i , .N ]) }

BioGRID_Y2HY_connectivity <- integer()
for(i in BioGRID_Y2HY_MPs){ BioGRID_Y2HY_connectivity <- c(BioGRID_Y2HY_connectivity, BioGRID_Y2HY[ SimpleID_1 == i | SimpleID_2 == i , .N ]) }

## Combine STRING and BioGRID data into one table
DTmp <- rbind( data.table( Protein_ID = STRING_COEX_MPs, connectivity = STRING_COEX_connectivity, source = "STRING", type = "coexpression"),
               data.table( Protein_ID = STRING_EXPE_MPs, connectivity = STRING_EXPE_connectivity, source = "STRING", type = "experimental"),
               data.table( Protein_ID = STRING_DATA_MPs, connectivity = STRING_DATA_connectivity, source = "STRING", type = "database"),
               data.table( Protein_ID = STRING_TEXT_MPs, connectivity = STRING_TEXT_connectivity, source = "STRING", type = "textmining"),
               
               data.table( Protein_ID = BioGRID_AFCA_MPs, connectivity = BioGRID_AFCA_connectivity, source = "BioGRID", type = "Affinity capture"),
               data.table( Protein_ID = BioGRID_COFR_MPs, connectivity = BioGRID_COFR_connectivity, source = "BioGRID", type = "Co-fractionation"),
               data.table( Protein_ID = BioGRID_RECO_MPs, connectivity = BioGRID_RECO_connectivity, source = "BioGRID", type = "In vitro"),
               data.table( Protein_ID = BioGRID_Y2HY_MPs, connectivity = BioGRID_Y2HY_connectivity, source = "BioGRID", type = "Two-hybrid"))

# Plot the results
p1 <- ggplot( DTmp[ source == "STRING"], aes( x = type, y = connectivity))+
        geom_boxplot( outlier.colour = "black" , outlier.shape = 16, outlier.size = 0.25, 
                      notch = TRUE , size = 0.25 , fill = "cornflowerblue")+
        ylab("Number of interaction partners")+
        xlab("Evidence type")+
        scale_y_continuous( breaks = seq(0,500,50))+
        coord_cartesian( ylim = c(0,360))+
        theme(panel.background = element_blank(), panel.grid = element_blank(),
              panel.border = element_rect(fill = NA, colour = "black", size = 0.25, linetype = "solid"),    
              axis.text=element_text(size=5), 
              axis.title=element_text(size=6), axis.ticks.y = element_line(size=0.25), axis.ticks.x = element_blank(),
              axis.line = element_blank(), legend.position = "none")
        

p2 <- ggplot( DTmp[ source == "BioGRID"], aes( x = type, y = connectivity))+
        geom_boxplot( outlier.colour = "black" , outlier.shape = 16, outlier.size = 0.25,
                      notch = TRUE , size = 0.25 , fill = "cornflowerblue")+
        coord_cartesian( ylim = c(0,50) )+
        ylab("Number of interaction partners")+
        xlab("Evidence type")+
        scale_y_continuous( breaks = seq(0,500,10))+
        theme(panel.background = element_blank(), panel.grid = element_blank(),
              panel.border = element_rect(fill = NA, colour = "black", size = 0.25, linetype = "solid"),    
              axis.text=element_text(size=5), 
              axis.title=element_text(size=6), axis.ticks.y = element_line(size=0.25), axis.ticks.x = element_blank(),
              axis.line = element_blank(), legend.position = "none")


#### Combine and print the plots ####

# Use cowplot to combine the plots into a figure
plot3by3 <- plot_grid(ptC, pSTRING, pBioGRID, p1, p2, nrow = 1, align = "hv", rel_widths = c(1,1,1,2,2))

# Save combined plot
save_plot("Microprotein_connectivity.pdf", plot3by3, nrow = 1, base_height = 2.2, base_width = 7.2)


