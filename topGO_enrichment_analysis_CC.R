# Read in the necessary libraries
library(plyr); library(data.table); library(topGO); library(ggplot2); library(grid)

#### Prepare the co-regulation data ####

# Read in ProteomeHD, so we can limit the GO search for proteins for which coregulation has actually been tested
ProHD <- fread("ProteomeHD_v1_1.csv")
ProHD_ratios <- ProHD[, .SD, .SDcols = colnames(ProHD) %like% "Ratio"]      # Limit to data columns
feature_count <- apply( ProHD_ratios, 1, function(x){ sum( !is.na(x)) } )   # Calculate number of features per protein
ProHD <- ProHD[ feature_count >= 95 ,]      # Discard proteins detected in fewer than 95 experiments (as was done to make the co-regulation map)
ProHD_proteins <- unique( ProHD[, SimpleID_1 := gsub(";.+", "", Majority_protein_IDs)][, gsub("-.+", "", SimpleID_1) ] )

# Read in co-regulated proteins (if any)
CP <- fread("coregulated_protein_pairs.csv")
CP[, Coreg_with_protein_X := gsub(";.+", "", Coreg_with_protein_X)][, Coreg_with_protein_X := gsub("-.+", "", Coreg_with_protein_X)]  # Simplify IDs 
CP[, protein_X := gsub(";.+", "", protein_X)][, protein_X := gsub("-.+", "", protein_X)]                                              # Simplify IDs 
CP <- CP[ Coreg_with_protein_X != protein_X ]                             # Remove coregulation pairs coming from different isoforms of same protein
CP <- unique(CP)                                                          # Remove duplicate coregulation pairs 
CP <- dlply( CP, "protein_X", function(x){ unique( c( x$Coreg_with_protein_X, x$protein_X )) })


#### GO enrichment analysis: CC ####

# Read in gene-to-GO associations. To download these from www.ebi.ac.uk/QuickGO/GAnnotation, click on "Filter"
# and select taxon = human, qualifier = NULL, aspect = Biological Process
GO <- fread("GO_associations_CC.tsv")
GO <- GO[ ID %in% ProHD_proteins , .(ID = ID, GO_ID = `GO ID`)]    # Limit to proteins in ProHD and the important columns
GO <- unique(GO)                                                   # Discard duplicate annotations
geneID2GO <- dlply( GO, "ID", function(x){as.character(x$GO_ID)})  # GO annotation for all proteins that are in coregulation analysis and for which GO annotation is available

# Turn co-regulated protein groups into allGenes input for topGO
# allGenes_list is a list of named factors, each of which can be used for topGO's allGenes parameter. The names are all proteins in
# ProHD for which coregulation has been tested using treeClust AND for which we have GO annotations (the "universe"). 
# The factor levels indicate whether proteins are in a co-regulation group or not
# This will be a list of such factors, each element corresponding to one coregulation group
geneNames <- names(geneID2GO)
allGenes_list <- lapply(CP, function(x){  geneList <- factor(as.integer(geneNames %in% x))
names(geneList) <- geneNames
return(geneList) })

# Discard coregulation groups that are too small for GO term enrichment testing
large_enough <- sapply( allGenes_list, function(x){ sum( x == 1) }) >= 11   # Group size 11 means at least 10 proteins coregulated with protein X
allGenes_list <- allGenes_list[ large_enough ]     

# Create pilot topGOdata object (which will be updated with new allGenes in each iteration below)
GOdata <- new("topGOdata", ontology = "CC", allGenes = allGenes_list[[1]],
              annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)

# Screen every coregulation group for GO-BP enrichment
pvalues_CC <- data.frame(GO_term = usedGO(GOdata))
for(i in 1:length(allGenes_list)){
  GOdata <- updateGenes(GOdata, allGenes_list[[i]])                           # Update the pilot GOdata object with the real coregulation group to be tested
  result <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")     # Calculate enrichment considering graph structure
  current_pvalue <- as.data.frame( score(result))                             # Get the pvalues for enrichment of all GO terms for this gene subset
  colnames(current_pvalue) <- names(allGenes_list)[i]                         # Attach the name of protein X
  pvalues_CC <- merge(pvalues_CC, current_pvalue, by.x="GO_term", by.y="row.names")
  print(i)
}


#### Plot the results ####

below_0.01 <- function(x){ sum( x < 0.01 )     > 0 }   # Is at least 1 GO term significantly enriched?
below_0.001 <- function(x){ sum( x < 0.001 )   > 0 }
below_0.0001 <- function(x){ sum( x < 0.0001 ) > 0 }

theme_pie <- theme(panel.background=element_blank(), panel.grid.major=element_blank(), axis.ticks=element_blank(),
                   axis.text=element_blank(), axis.title=element_blank())

# Turn to data.table
CC <- as.data.table( pvalues_CC )

# Calculations and plotting for CC
CC[, GO_term := NULL]                      # Remove GO_term column to leave only p-values per neighbourhood
CC_0.01 <- apply(CC, 2, below_0.01)        # Get GO terms enriched with p < 0.01
CC_0.001 <- apply(CC, 2, below_0.001)      # Get GO terms enriched with p < 0.001
CC_0.0001 <- apply(CC, 2, below_0.0001)    # Get GO terms enriched with p < 0.0001

CC <- melt( data.frame( p0.01 = CC_0.01, p0.001 = CC_0.001, p0.0001 = CC_0.0001), measure.vars = c("p0.01", "p0.001", "p0.0001"))

CC <- ggplot(CC, aes(x=variable, fill=value))+
  geom_bar() + coord_polar(theta="y", direction = -1)+
  theme_pie + scale_fill_manual(guide = guide_legend(title = "Co-regulation groups"), values = c("seagreen2", "royalblue4" ),
                                labels = c("No GO term enriched", "GO terms enriched"))

ggsave("GO_CC_pie.pdf", CC, width = 10, height = 5, units = c("cm"))








