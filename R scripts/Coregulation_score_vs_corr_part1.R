## This script calculates pair-wise correlation scores using various metrics

library(data.table); library(reshape2); library(ggplot2); library(treeClust); library(WGCNA)

# Load ProteomeHD
prohd <- fread("ProteomeHD_v1_1.csv")
prohd <- data.frame(prohd, row.names = 1)
prohd <- prohd[,-c(1:3)]

# Select only proteins with at least 95 ratios
feature_count <- apply(prohd, 1, function(x){ sum(!is.na(x))} )
prohd_ratios_min95 <- prohd[feature_count >= 95,]

# Median normalise
prohd_ratios_min95 <- t(prohd_ratios_min95)
my_medians <- apply(prohd_ratios_min95, 1, median, na.rm=TRUE)
prohd_ratios_min95 <- sweep(prohd_ratios_min95, 1, my_medians, FUN="-")

# Run Pearson's correlation
PCCs <- stats::cor( prohd_ratios_min95, method="pearson", use="pairwise.complete.obs")
PCCs <- as.data.table( melt( as.matrix( PCCs  )))       # Convert to a long data table
PCCs <- PCCs[, .( Protein_1 = as.character(Var1),       # Change colnames and types
                  Protein_2 = as.character(Var2), 
                        PCC = value ) ]
PCCs <- PCCs[ Protein_1 > Protein_2 ]                   # Remove self-comparisons and duplicate pairs

# Run Spearman's correlation
RHOs <- stats::cor( prohd_ratios_min95, method="spearman", use="pairwise.complete.obs" )
RHOs <- as.data.table( melt( as.matrix( RHOs  )))       # Convert to a long data table
RHOs <- RHOs[, .( Protein_1 = as.character(Var1),       # Change colnames and types
                  Protein_2 = as.character(Var2), 
                  RHO = value ) ]
RHOs <- RHOs[ Protein_1 > Protein_2 ]                   # Remove self-comparisons and duplicate pairs

# Run bicor correlation
BICs <- bicor( prohd_ratios_min95, use="pairwise.complete.obs")
BICs <- as.data.table( melt( as.matrix( BICs  )))       # Convert to a long data table
BICs <- BICs[, .( Protein_1 = as.character(Var1),       # Change colnames and types
                  Protein_2 = as.character(Var2), 
                  BIC = value ) ]
BICs <- BICs[ Protein_1 > Protein_2 ]                   # Remove self-comparisons and duplicate pairs

# Merge the three datasets
setkey(PCCs, Protein_1, Protein_2)
setkey(RHOs, Protein_1, Protein_2)
setkey(BICs, Protein_1, Protein_2)

DT <- merge(PCCs, RHOs)
DT <- merge(DT, BICs)

# Write out result
fwrite(DT, "ProHD_correlations.csv")

