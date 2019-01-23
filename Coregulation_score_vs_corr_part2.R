## A script to compare the distance/correlation metrics to treeClust using REACTOME test set and precision-recall curves

# Load the required packages
library(data.table); library(reshape2); library(ggplot2); library(ROCR); library(scales)

#### Prepare the data ####

# Read in the data
   treeclust <- fread("treeClust_similarities.csv")
correlations <- fread("ProHD_correlations.csv")

# Prepare the tables for merging
treeclust <- treeclust[, .(Protein_1, Protein_2, treeClust_sim = tc_sim, treeClust_tom = tc_tom) ]
setkey(treeclust, Protein_1, Protein_2)
setkey(correlations, Protein_1, Protein_2)

# Merge into one table
DT <- merge(treeclust, correlations)

# Remove rows with missing values (only 19 rows)
DT <- DT[ complete.cases(DT) ]

# Remove isoform information (which can't be mapped to gold standard)
DT[, SimpleID_1 := gsub(";.+", "", Protein_1)][, SimpleID_1 := gsub("-.+", "", SimpleID_1)]  # Simplify protein 1 IDs of DT
DT[, SimpleID_2 := gsub(";.+", "", Protein_2)][, SimpleID_2 := gsub("-.+", "", SimpleID_2)]  # Simplify protein 2 IDs of DT
DT <- DT[ SimpleID_1 != SimpleID_2 ]                      # Remove self-interactions after ID simplification (isoform removal)
DT <- DT[ !duplicated( DT[,.(SimpleID_1, SimpleID_2)] ) ] # Remove duplicates after ID simplification (isoform removal)
DT[ SimpleID_1 > SimpleID_2 , .N] == DT[, .N]             # Double checking that all pairs are in B <-> A order (alphabetically)
DT <- DT[, .(SimpleID_1, SimpleID_2, treeClust_sim,       # Keep only relevant columns
             treeClust_tom, PCC, RHO, BIC)]                          
setkey(DT, SimpleID_1, SimpleID_2)                        # Set keys for merging


#### Compare metrics using precision recall curves ####

# Load gold standard of true and false positives (based on Reactome)
TP_FP_pairs <- fread("Reactome_TP_FP.csv")                     # Note that these pairs are already sorted such that Protein_1 > Protein_2
names(TP_FP_pairs) <- c("SimpleID_1", "SimpleID_2", "Class")   # But rename them to fit DT
setkey(TP_FP_pairs, SimpleID_1, SimpleID_2)                    # Set keys for merging

# Merge gold standard with data
DT <- merge(DT, TP_FP_pairs)

# Sample FPs down to get 5% TP
TP <- DT[ Class == "TP" ]            # All true positives
FP <- DT[ Class == "FP" ]            # All false positives
n_FP <- TP[, .N] * 19                # Number of false positives we need
FP <- FP[ sample( FP[,.N] , n_FP) ]  # Restrict to random subset of false positives
DT <- rbindlist( list( TP, FP))      # Downsample DT

# Add a randomised classifier
DT[, Random := sample( treeClust_tom ) ]   

# Get precision recall data using ROCR package
treeClust_sim <- DT$treeClust_sim
treeClust_tom <- DT$treeClust_tom
PCC <- DT$PCC
RHO <- DT$RHO
BIC <- DT$BIC
Random <- DT$Random
labels <- DT$Class

pred <- prediction( predictions = list(treeClust_sim, treeClust_tom, PCC,    RHO,    BIC,    Random),
                         labels = list(labels,        labels,        labels, labels, labels, labels),
                    label.ordering = c("FP", "TP"))
perf <- performance(pred, measure = "prec", x.measure = "rec")

# Make the precision recall plot
treeClust_sim = data.table(    Recall = perf@x.values[[1]],
                            Precision = perf@y.values[[1]],
                              Measure = "treeClust")
treeClust_tom = data.table(    Recall = perf@x.values[[2]],
                            Precision = perf@y.values[[2]],
                              Measure = "treeClust + TOM")
PCC = data.table(              Recall = perf@x.values[[3]],
                            Precision = perf@y.values[[3]],
                              Measure = "PCC")
RHO = data.table(              Recall = perf@x.values[[4]],
                            Precision = perf@y.values[[4]],
                              Measure = "Rho")
BIC = data.table(              Recall = perf@x.values[[5]],
                            Precision = perf@y.values[[5]],
                              Measure = "bicor")
Random = data.table(           Recall = perf@x.values[[6]],
                            Precision = perf@y.values[[6]],
                              Measure = "random classifier")

pre_rec_dt <- rbindlist( list( treeClust_sim, treeClust_tom, PCC, RHO, BIC, Random) )
pre_rec_dt <- pre_rec_dt[ Recall > 0.005 ]                      # Drop low recall points because their pretty much random
pre_rec_dt <- pre_rec_dt[ sample( pre_rec_dt[,.N], 2000000 ) ]  # Randomly downsample to speed up loading times in Inkscape

# Set plotting order
pre_rec_dt[, Measure := factor( Measure, levels = c("treeClust", "treeClust + TOM", "PCC", "Rho", "bicor", "random classifier")) ]

# Plot the result
p1 <- ggplot(pre_rec_dt, aes(x = Recall, y = Precision, colour = Measure))+
        geom_line( size = 0.25 )+
        scale_colour_manual(values = c("royalblue1", "navy", "lightseagreen", "violetred2", "mediumorchid", "grey50"))+
        scale_x_continuous(limits=c(0,1), breaks=seq(0,1,0.2))+
        scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.2))+
        theme(panel.background = element_blank(), axis.text=element_text(size=5), axis.title=element_text(size=6),
              axis.ticks = element_line(size=0.25), axis.line = element_line(colour="black", size=0.25),
              legend.position = "none", plot.margin = rep(unit(0,"null"),4))

ggsave("TreeClust_vs_Cor.png", p1, width=4.6, height=4.6, units=c("cm"), dpi = 600)



