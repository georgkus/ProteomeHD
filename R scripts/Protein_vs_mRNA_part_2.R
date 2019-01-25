## This script compares the performance of protein and mRNA expression data for gene function prediction

# Load required libraries
library(ggplot2); library(treeClust); library(data.table); library(ROCR); library(scales); library(VennDiagram)

#### Load and prep the data ####

# Load the ProteomeHD correlations (keep only Pearson's correlations for this)
ProHD_cor <- fread("ProHD_correlations.csv")
ProHD_cor <- ProHD_cor[, .( Protein_1, Protein_2, ProHD_PCC = PCC )]
setkey(ProHD_cor, Protein_1, Protein_2)                 # Set key for merging

# Load the coregulation scores
tC <- fread("coregulation_scores.csv")
setkey(tC, Protein_1, Protein_2)                        # Set keys for merging

# Merge and simplify table
tC <- merge(tC, ProHD_cor)                                                                    # Merge with PCCs from ProteomeHD
tC[, SimpleID_1 := gsub(";.+", "", Protein_1)][, SimpleID_1 := gsub("-.+", "", SimpleID_1)]   # Simplify protein_1 IDs
tC[, SimpleID_2 := gsub(";.+", "", Protein_2)][, SimpleID_2 := gsub("-.+", "", SimpleID_2)]   # Simplify protein_2 IDs
tC <- tC[ SimpleID_1 != SimpleID_2 ]                      # Remove self-interactions after ID simplification (isoform removal)
tC <- tC[ !duplicated( tC[,.(SimpleID_1, SimpleID_2)] ) ] # Remove duplicates after ID simplification (isoform removal)
tC[ SimpleID_1 > SimpleID_2 , .N] == tC[, .N]   # Double checking that all pairs are in B <-> A order (alphabetically)
tC <- tC[, .(SimpleID_1, SimpleID_2, coregulation_score, ProHD_PCC)]                          # Keep only relevant columns
tC <- tC[ complete.cases(tC) ]                  # Remove rows with missing values (some PCCs weren't calculated)
setkey(tC, SimpleID_1, SimpleID_2)                                                            # Set keys for merging

# Clear workspace and memory
rm( list = ls()[! ls() %in% c("tC")] )                                  
gc()    

# Load the pre-processed STRING data. These are PCCs obtained by STRING for microarrays and RNA-seq data in GEO. The original
# file belongs to STRING so I won't make it available for download.
STRING <- fread("RNA_STRING_Uniprot.csv")
STRING[, STR_coexpr_score := NULL ]                      # Remove unnecessary column
STRING <- STRING[ complete.cases(STRING) ]               # Consider only protein pairs for which there are correlations in all three channels
STRING[, RNA_PCC_mean := rowMeans(.SD, na.rm = TRUE ),   # Get the mean PCC for 1- and 2-channel microarrays and RNA-seq data
         .SDcols = c("ch1_PCC", "ch2_PCC", "seq_PCC") ]  # (this gives the highest accuracy for mRNA data)

STRING[, ch1_PCC := NULL ]                               # Remove unnecessary columns
STRING[, ch2_PCC := NULL ]
STRING[, seq_PCC := NULL ]

STRING[ SimpleID_1 > SimpleID_2 , .N] == STRING[, .N]   # Double checking that all pairs are in B <-> A order (alphabetically)
setkey(STRING, SimpleID_1, SimpleID_2)                  # Set key for merging

rm( list = ls()[! ls() %in% c("tC", "STRING")] )        # Clear workspace and memory                         
gc() 

# Get unique genes covered by the relevant datasets and find intersection
 STRING_genes <- STRING[, unique( c( SimpleID_1, SimpleID_2 )) ]
     tC_genes <- tC[, unique(c(SimpleID_1, SimpleID_2))]     

# Load gold standard of true and false positives (based on Reactome)
TP_FP_pairs <- fread("Reactome_TP_FP.csv")              # Note that these pairs are already sorted such that Protein_1 > Protein_2
names(TP_FP_pairs) <- c("SimpleID_1", "SimpleID_2", "Class")
setkey(TP_FP_pairs, SimpleID_1, SimpleID_2)             # Set keys for merging

# Merge data
DT <- merge( STRING, tC )        
DT <- merge( DT, TP_FP_pairs )

rm( list = ls()[! ls() %in% c("DT", "STRING_genes", "tC_genes")] )   # Clear workspace
gc()  

# How many genes are left in the comparison?
N_genes <- DT[, length( unique( c( SimpleID_1, SimpleID_2 ))) ]


#### Precision - Recall analysis ####

# Append a random classifier
DT[, Random := sample(coregulation_score) ]

# Perform precision recall analysis using ROCR package
coregulation_score <- DT$coregulation_score
ProHD_PCC <- DT$ProHD_PCC
RNA_PCC_mean <- DT$RNA_PCC_mean
Random <- DT$Random
labels <- DT$Class

pred <- prediction( predictions = list( coregulation_score,   ProHD_PCC,  RNA_PCC_mean, Random),
                         labels = list( labels,               labels,     labels,       labels),
                    label.ordering = c("FP", "TP"))
perf <- performance(pred, measure = "prec", x.measure = "rec")

# Extract the data for a PR curve
coregulation_score <- data.table( Recall = perf@x.values[[1]],
                               Precision = perf@y.values[[1]],
                                 Measure = "ProHD + treeClust + TOM")
ProHD_PCC <- data.table(          Recall = perf@x.values[[2]],
                               Precision = perf@y.values[[2]],
                                 Measure = "ProHD + PCC")
RNA_PCC_mean <- data.table(       Recall = perf@x.values[[3]],
                               Precision = perf@y.values[[3]],
                                 Measure = "mRNA PCC")
Random <- data.table(             Recall = perf@x.values[[4]],
                               Precision = perf@y.values[[4]],
                                 Measure = "random")

pre_tp_dt <- rbindlist( list( coregulation_score, ProHD_PCC, RNA_PCC_mean, Random ))

pre_tp_dt <- pre_tp_dt[ Recall > 0.005 ]                     # Drop low recall points because their pretty much random

p1 <- ggplot(pre_tp_dt[ Measure != "ProHD + PCC"], aes(x = Recall, y = Precision, colour = Measure))+
      geom_line( size = 0.25 )+
      geom_line( data = pre_tp_dt[ Measure == "ProHD + PCC"], size = 0.25 , linetype = "dashed")+
      annotate("text", label = paste("N genes =", N_genes), size = 2, x = 0.7, y = 0.7)+
      scale_colour_manual(values = c("#00A0BB", "#EC008C", "#EC008C", "grey50"))+
      scale_x_continuous(limits=c(0,1), breaks=seq(0,1,0.2))+
      scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.2))+
      theme(panel.background = element_blank(), axis.text=element_text(size=5), axis.title=element_text(size=6),
            axis.ticks = element_line(size=0.25), axis.line = element_line(colour="black", size=0.25),
            legend.position = "none", plot.margin = rep(unit(0,"null"),4))

ggsave("mRNA_protein_STRING_PR.png", p1, width=6.2, height=5.2, units=c("cm"), dpi = 600)


#### Plot coverage, i.e. protein numbers ####
grid.newpage()
pdf("mRNA_protein_STRING_coverage.pdf")
draw.pairwise.venn(area1 = length(STRING_genes), area2 = length(tC_genes),
                   cross.area = length(intersect(STRING_genes, tC_genes)),
                   category = c("mRNA", "protein"))
dev.off()

# Barplot to show number of experiments. I counted the number of GEO experiments going into STRING previously,
# and ProteomeHD has 294 variables
cont_table <- melt(data.table( mRNA = 1112784, protein = 294 ), measure.vars = c("mRNA", "protein"))

p2 <- ggplot(cont_table, aes( x = variable, y = value, fill = variable))+
        geom_bar(stat="identity")+
        ylab("# experiments")+
        scale_fill_manual(values = c("#00A0BB", "#EC008C"))+
        scale_y_log10(breaks = c(100,10000,1000000), expand = c(0,0))+
        theme(panel.background = element_blank(), axis.text=element_text(size=5), axis.title=element_text(size=6),
              axis.ticks = element_line(size=0.25), axis.line = element_line(colour="black", size=0.25), 
              legend.position = "none")

# Print barplot
ggsave("N_exp_bars_STRING.pdf", p2, width = 2, height = 2, units = "cm")

