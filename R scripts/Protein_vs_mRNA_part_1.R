## This script compares the performance of protein and mRNA expression data for gene function prediction

# Load required libraries
library(ggplot2); library(treeClust); library(data.table); library(ROCR); library(scales); library(VennDiagram)

#### mRNA vs proteins from the same samples: PR curve ####

# Load the RNA - protein comparison data (no normalisation necessary)
df <- read.csv("BattleSILAC_PickrellRPKM.csv", stringsAsFactors = FALSE)
rownames(df) <- df$Majority.protein.IDs
RPKM <- df[, grep("RPKM_", colnames(df))]
SILAC <- df[, grep("SILAC_", colnames(df))]

# Use treeClust to learn a dissimilarity matrix
set.seed(42)
tc_dist_RPKM  <- treeClust.dist( RPKM, d.num = 2, verbose = TRUE)
tc_dist_SILAC <- treeClust.dist(SILAC, d.num = 2, verbose = TRUE)

# Turn into one table with pairwise treeClust dissimilarities
tc_dt_RPKM  <- as.data.table( melt( as.matrix( tc_dist_RPKM  )))     # Convert to a long data table
tc_dt_SILAC <- as.data.table( melt( as.matrix( tc_dist_SILAC )))     

tc_dt_RPKM  <-  tc_dt_RPKM[, .( Protein_1 = as.character(Var1), Protein_2 = as.character(Var2), tC_dissimilarity_RPKM = value )]     # Change colnames and types
tc_dt_SILAC <- tc_dt_SILAC[, .( Protein_1 = as.character(Var1), Protein_2 = as.character(Var2), tC_dissimilarity_SILAC = value )]

 tc_dt_RPKM <-  tc_dt_RPKM[ Protein_1 > Protein_2 ]                 # Remove self-comparisons and duplicate pairs
tc_dt_SILAC <- tc_dt_SILAC[ Protein_1 > Protein_2 ]   

setkey( tc_dt_RPKM, Protein_1, Protein_2)                           # Set key for merging
setkey(tc_dt_SILAC, Protein_1, Protein_2)

Battle_dataset <- merge( tc_dt_RPKM, tc_dt_SILAC )                  # Merge into one table

# For the barchart later, save the number of experiments
Battle_N_exp <- melt(data.table( mRNA = ncol(RPKM), protein = ncol(SILAC)), measure.vars = c("mRNA", "protein"))

# Clear workspace
rm( list = ls()[! ls() %in% c("Battle_dataset", "Battle_N_exp")] )

# Load gold standard of true and false positives (based on Reactome)
TP_FP_pairs <- fread("Reactome_TP_FP.csv")   # Note that these pairs are already sorted such that Protein_1 > Protein_2
setkey(TP_FP_pairs, Protein_1, Protein_2)    # Set keys for merging

# Merge data with gold standard
Battle_dataset <- merge(Battle_dataset, TP_FP_pairs)

# Calculate the number of genes used for the analysis
N_genes_battle <- Battle_dataset[, length( unique( c( Protein_1, Protein_2 ))) ]

# Add a randomised (scrambled) classifier
Battle_dataset[, Random := sample( tC_dissimilarity_SILAC ) ]   

# Perform precision recall analysis using ROCR package
  RPKM <- 1 - Battle_dataset$tC_dissimilarity_RPKM    # Get the (inverse) scores for ranking
 SILAC <- 1 - Battle_dataset$tC_dissimilarity_SILAC
Random <- 1 - Battle_dataset$Random
labels <- Battle_dataset$Class                        # Get the class labels

pred <- prediction( predictions = list( RPKM,   SILAC,  Random),
                         labels = list( labels, labels, labels),
                    label.ordering = c("FP", "TP"))
perf <- performance(pred, measure = "prec", x.measure = "rec")

# Extract the data for a PR curve
 RPKM <- data.table(       Recall = perf@x.values[[1]],
                        Precision = perf@y.values[[1]],
                          Measure = "mRNA")
SILAC <- data.table(       Recall = perf@x.values[[2]],
                        Precision = perf@y.values[[2]],
                          Measure = "protein")
Random <- data.table(      Recall = perf@x.values[[3]],
                        Precision = perf@y.values[[3]],
                          Measure = "random")
pre_tp_dt <- rbindlist( list(RPKM, SILAC, Random ))

pre_tp_dt <- pre_tp_dt[ Recall > 0.005 ]                     # Drop low recall points because their pretty much random
pre_tp_dt <- pre_tp_dt[ sample( pre_tp_dt[,.N], 2000000 ) ]  # Randomly downsample to speed up loading times in Inkscape

p1 <- ggplot(pre_tp_dt, aes(x = Recall, y = Precision, colour = Measure))+
      geom_line( size = 0.25 )+
      annotate("text", label = paste("N genes =", N_genes_battle), size = 2, x = 0.7, y = 0.7)+
      scale_colour_manual(values = c("#00A0BB", "#EC008C", "grey50"))+
      scale_x_continuous(limits=c(0,1), breaks=seq(0,1,0.2))+
      scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.2))+
      theme(panel.background = element_blank(), axis.text=element_text(size=5), axis.title=element_text(size=6),
            axis.ticks = element_line(size=0.25), axis.line = element_line(colour="black", size=0.25),
            legend.position = "none", plot.margin = rep(unit(0,"null"),4))

ggsave("mRNA_protein_Battle_PR.png", p1, width=6.2, height=5.2, units=c("cm"), dpi = 600)


#### mRNA vs proteins from the same samples: coverage ####

# I previously used Battle et al's and Pickrell et al's original tables to calculate the number of genes (ENSG IDs)
# covered by both, and simply hard-code these here
N_mRNA <- 18302
N_prot <- 4381
N_overlap <- 4338

# Print the Venn Diagram
grid.newpage()
pdf("mRNA_protein_Battle_coverage.pdf")
draw.pairwise.venn(area1 = N_mRNA, area2 = N_prot, cross.area = N_overlap, category = c("mRNA", "protein"))
dev.off()

# Barplot to show number of experiments
p2 <- ggplot(Battle_N_exp, aes( x = variable, y = value, fill = variable))+
      geom_bar(stat="identity")+
      ylab("# experiments")+
      scale_fill_manual(values = c("#00A0BB", "#EC008C"))+
      scale_y_continuous(breaks = c(0,30,60), limits = c(0,60), expand = c(0,0))+
      theme(panel.background = element_blank(), axis.text=element_text(size=5), axis.title=element_text(size=6),
            axis.ticks = element_line(size=0.25), axis.line = element_line(colour="black", size=0.25), 
            legend.position = "none")

# Print barplot
ggsave("N_exp_bars_Battle.pdf", p2, width = 2, height = 2, units = "cm")


