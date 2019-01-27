## This script drops data randomly from ProteomeHD and looks at the impact this has on performance, using a precision - recall analysis.
## It is designed to be executed on a larger server, otherwise it will take a long time...

# Load necessary libraries
library(treeClust); library(ggplot2); library(data.table); library(ROCR);  library(scales); library(parallel)

# Force data.table to use a single core (in order to avoid a bug in reading / writing parallelising)
setDTthreads(1)

#### Drop random data points across the ProteomeHD data matrix ####

# Define necessary functions
count_features <- function(x){ sum( !is.na(x) ) }                  # To calculate number of SILAC ratios per protein

# Load and prepare ProteomeHD
ProHD <- read.csv("ProteomeHD_v1_1.csv", stringsAsFactors=FALSE)
rownames(ProHD) <- ProHD$Majority_protein_IDs       # Set protein IDs as rownames
ProHD <- ProHD[, grep("Ratio", colnames(ProHD))]    # Keep only columns with SILAC ratios

# Restrict to the working subset of proteins
feature_count <- apply(ProHD, 1, count_features)    # Count number of SILAC ratios per protein
ProHD <- ProHD[ feature_count >= 95 ,]              # Discard proteins detected in fewer than 95 experiments

# Create datasets where 5%, 10% or 15% of data points have been dropped at random (3 replicates each)
number_of_datapoints_in_ProHD <- sum( !is.na(ProHD) )
nr <- nrow(ProHD)
nc <- ncol(ProHD)

pc00 <- floor( number_of_datapoints_in_ProHD/100*0  )  # Include this as "full ProHD" - just to make sure the whole procedure doesn't mess up the data
pc05 <- floor( number_of_datapoints_in_ProHD/100*5  )
pc10 <- floor( number_of_datapoints_in_ProHD/100*10 )
pc15 <- floor( number_of_datapoints_in_ProHD/100*15 )


# One replicate 0% data
 ProHD_logical <- as.logical(!is.na(ProHD))  # Turn ProHD into a logical vector, where TRUE means there is a value at this position
TRUE_positions <- which(ProHD_logical)       # Get the indices (positions) of the TRUEs
pc00_s1 <- sample(TRUE_positions, pc00)      # Get 0% of these TRUE positions at random
ProHD_logical[ pc00_s1 ] <- FALSE            # Set those randomly chosen TRUE positions to FALSE
pc00_s1 <- matrix(ProHD_logical, nrow = nr, ncol = nc)   # Create a ProHD-sized logical matrix, TRUE means there should be a value
ProHD_pc00_s1 <- ProHD                       # Create a copy of ProHD to be blanked out
ProHD_pc00_s1[ !pc00_s1 ] <- NA              # Turn the random subset plus the original NAs into NAs


# Three replicates of 5% data
 ProHD_logical <- as.logical(!is.na(ProHD))  # Turn ProHD into a logical vector, where TRUE means there is a value at this position
TRUE_positions <- which(ProHD_logical)       # Get the indices (positions) of the TRUEs
pc05_s1 <- sample(TRUE_positions, pc05)      # Get 5% of these TRUE positions at random
ProHD_logical[ pc05_s1 ] <- FALSE            # Set those randomly chosen TRUE positions to FALSE
pc05_s1 <- matrix(ProHD_logical, nrow = nr, ncol = nc)   # Create a ProHD-sized logical matrix, TRUE means there should be a value
ProHD_pc05_s1 <- ProHD                       # Create a copy of ProHD to be blanked out
ProHD_pc05_s1[ !pc05_s1 ] <- NA              # Turn the random subset plus the original NAs into NAs

 ProHD_logical <- as.logical(!is.na(ProHD))  # Turn ProHD into a logical vector, where TRUE means there is a value at this position
TRUE_positions <- which(ProHD_logical)       # Get the indices (positions) of the TRUEs
pc05_s2 <- sample(TRUE_positions, pc05)      # Get 5% of these TRUE positions at random
ProHD_logical[ pc05_s2 ] <- FALSE            # Set those randomly chosen TRUE positions to FALSE
pc05_s2 <- matrix(ProHD_logical, nrow = nr, ncol = nc)   # Create a ProHD-sized logical matrix, TRUE means there should be a value
ProHD_pc05_s2 <- ProHD                       # Create a copy of ProHD to be blanked out
ProHD_pc05_s2[ !pc05_s2 ] <- NA              # Turn the random subset plus the original NAs into NAs

 ProHD_logical <- as.logical(!is.na(ProHD))  # Turn ProHD into a logical vector, where TRUE means there is a value at this position
TRUE_positions <- which(ProHD_logical)       # Get the indices (positions) of the TRUEs
pc05_s3 <- sample(TRUE_positions, pc05)      # Get 5% of these TRUE positions at random
ProHD_logical[ pc05_s3 ] <- FALSE            # Set those randomly chosen TRUE positions to FALSE
pc05_s3 <- matrix(ProHD_logical, nrow = nr, ncol = nc)   # Create a ProHD-sized logical matrix, TRUE means there should be a value
ProHD_pc05_s3 <- ProHD                       # Create a copy of ProHD to be blanked out
ProHD_pc05_s3[ !pc05_s3 ] <- NA              # Turn the random subset plus the original NAs into NAs


# Three replicates of 10% data
 ProHD_logical <- as.logical(!is.na(ProHD))  # Turn ProHD into a logical vector, where TRUE means there is a value at this position
TRUE_positions <- which(ProHD_logical)       # Get the indices (positions) of the TRUEs
pc10_s1 <- sample(TRUE_positions, pc10)      # Get 10% of these TRUE positions at random
ProHD_logical[ pc10_s1 ] <- FALSE            # Set those randomly chosen TRUE positions to FALSE
pc10_s1 <- matrix(ProHD_logical, nrow = nr, ncol = nc)   # Create a ProHD-sized logical matrix, TRUE means there should be a value
ProHD_pc10_s1 <- ProHD                       # Create a copy of ProHD to be blanked out
ProHD_pc10_s1[ !pc10_s1 ] <- NA              # Turn the random subset plus the original NAs into NAs

 ProHD_logical <- as.logical(!is.na(ProHD))  # Turn ProHD into a logical vector, where TRUE means there is a value at this position
TRUE_positions <- which(ProHD_logical)       # Get the indices (positions) of the TRUEs
pc10_s2 <- sample(TRUE_positions, pc10)      # Get 10% of these TRUE positions at random
ProHD_logical[ pc10_s2 ] <- FALSE            # Set those randomly chosen TRUE positions to FALSE
pc10_s2 <- matrix(ProHD_logical, nrow = nr, ncol = nc)   # Create a ProHD-sized logical matrix, TRUE means there should be a value
ProHD_pc10_s2 <- ProHD                       # Create a copy of ProHD to be blanked out
ProHD_pc10_s2[ !pc10_s2 ] <- NA              # Turn the random subset plus the original NAs into NAs

 ProHD_logical <- as.logical(!is.na(ProHD))  # Turn ProHD into a logical vector, where TRUE means there is a value at this position
TRUE_positions <- which(ProHD_logical)       # Get the indices (positions) of the TRUEs
pc10_s3 <- sample(TRUE_positions, pc10)      # Get 10% of these TRUE positions at random
ProHD_logical[ pc10_s3 ] <- FALSE            # Set those randomly chosen TRUE positions to FALSE
pc10_s3 <- matrix(ProHD_logical, nrow = nr, ncol = nc)   # Create a ProHD-sized logical matrix, TRUE means there should be a value
ProHD_pc10_s3 <- ProHD                       # Create a copy of ProHD to be blanked out
ProHD_pc10_s3[ !pc10_s3 ] <- NA              # Turn the random subset plus the original NAs into NAs


# Three replicates of 15% data
 ProHD_logical <- as.logical(!is.na(ProHD))  # Turn ProHD into a logical vector, where TRUE means there is a value at this position
TRUE_positions <- which(ProHD_logical)       # Get the indices (positions) of the TRUEs
pc15_s1 <- sample(TRUE_positions, pc15)      # Get 15% of these TRUE positions at random
ProHD_logical[ pc15_s1 ] <- FALSE            # Set those randomly chosen TRUE positions to FALSE
pc15_s1 <- matrix(ProHD_logical, nrow = nr, ncol = nc)   # Create a ProHD-sized logical matrix, TRUE means there should be a value
ProHD_pc15_s1 <- ProHD                       # Create a copy of ProHD to be blanked out
ProHD_pc15_s1[ !pc15_s1 ] <- NA              # Turn the random subset plus the original NAs into NAs

 ProHD_logical <- as.logical(!is.na(ProHD))  # Turn ProHD into a logical vector, where TRUE means there is a value at this position
TRUE_positions <- which(ProHD_logical)       # Get the indices (positions) of the TRUEs
pc15_s2 <- sample(TRUE_positions, pc15)      # Get 15% of these TRUE positions at random
ProHD_logical[ pc15_s2 ] <- FALSE            # Set those randomly chosen TRUE positions to FALSE
pc15_s2 <- matrix(ProHD_logical, nrow = nr, ncol = nc)   # Create a ProHD-sized logical matrix, TRUE means there should be a value
ProHD_pc15_s2 <- ProHD                       # Create a copy of ProHD to be blanked out
ProHD_pc15_s2[ !pc15_s2 ] <- NA              # Turn the random subset plus the original NAs into NAs

 ProHD_logical <- as.logical(!is.na(ProHD))  # Turn ProHD into a logical vector, where TRUE means there is a value at this position
TRUE_positions <- which(ProHD_logical)       # Get the indices (positions) of the TRUEs
pc15_s3 <- sample(TRUE_positions, pc15)      # Get 15% of these TRUE positions at random
ProHD_logical[ pc15_s3 ] <- FALSE            # Set those randomly chosen TRUE positions to FALSE
pc15_s3 <- matrix(ProHD_logical, nrow = nr, ncol = nc)   # Create a ProHD-sized logical matrix, TRUE means there should be a value
ProHD_pc15_s3 <- ProHD                       # Create a copy of ProHD to be blanked out
ProHD_pc15_s3[ !pc15_s3 ] <- NA              # Turn the random subset plus the original NAs into NAs


# Control: Did I drop the correct number of data points?
sum( !is.na(ProHD) )          == number_of_datapoints_in_ProHD - pc00
sum( !is.na(ProHD_pc00_s1) )  == number_of_datapoints_in_ProHD - pc00
sum( !is.na(ProHD_pc05_s1) )  == number_of_datapoints_in_ProHD - pc05
sum( !is.na(ProHD_pc05_s2) )  == number_of_datapoints_in_ProHD - pc05
sum( !is.na(ProHD_pc05_s3) )  == number_of_datapoints_in_ProHD - pc05
sum( !is.na(ProHD_pc10_s1) )  == number_of_datapoints_in_ProHD - pc10
sum( !is.na(ProHD_pc10_s2) )  == number_of_datapoints_in_ProHD - pc10
sum( !is.na(ProHD_pc10_s3) )  == number_of_datapoints_in_ProHD - pc10
sum( !is.na(ProHD_pc15_s1) )  == number_of_datapoints_in_ProHD - pc15
sum( !is.na(ProHD_pc15_s2) )  == number_of_datapoints_in_ProHD - pc15
sum( !is.na(ProHD_pc15_s3) )  == number_of_datapoints_in_ProHD - pc15

# Control: is ProHD and ProHD_pc00_s1 the same (does restoring work?)
sum(ProHD == ProHD_pc00_s1, na.rm = TRUE) == sum( !is.na(ProHD) )


#### Random data drops: use treeClust to learn dissimilarities ####

# Do the learning for all 10 datasets
tC_pc00_s1 <- treeClust.dist(ProHD_pc00_s1, d.num = 2, control = treeClust.control(parallelnodes = 40))

tC_pc05_s1 <- treeClust.dist(ProHD_pc05_s1, d.num = 2, control = treeClust.control(parallelnodes = 40))
tC_pc05_s2 <- treeClust.dist(ProHD_pc05_s2, d.num = 2, control = treeClust.control(parallelnodes = 40))
tC_pc05_s3 <- treeClust.dist(ProHD_pc05_s3, d.num = 2, control = treeClust.control(parallelnodes = 40))

tC_pc10_s1 <- treeClust.dist(ProHD_pc10_s1, d.num = 2, control = treeClust.control(parallelnodes = 40))
tC_pc10_s2 <- treeClust.dist(ProHD_pc10_s2, d.num = 2, control = treeClust.control(parallelnodes = 40))
tC_pc10_s3 <- treeClust.dist(ProHD_pc10_s3, d.num = 2, control = treeClust.control(parallelnodes = 40))

tC_pc15_s1 <- treeClust.dist(ProHD_pc15_s1, d.num = 2, control = treeClust.control(parallelnodes = 40))
tC_pc15_s2 <- treeClust.dist(ProHD_pc15_s2, d.num = 2, control = treeClust.control(parallelnodes = 40))
tC_pc15_s3 <- treeClust.dist(ProHD_pc15_s3, d.num = 2, control = treeClust.control(parallelnodes = 40))

# Get molten, pairwise treeClust dissimilarities
tC_pc00_s1 <- as.data.table( melt( as.matrix(tC_pc00_s1)))

tC_pc05_s1 <- as.data.table( melt( as.matrix(tC_pc05_s1)))
tC_pc05_s2 <- as.data.table( melt( as.matrix(tC_pc05_s2)))
tC_pc05_s3 <- as.data.table( melt( as.matrix(tC_pc05_s3)))

tC_pc10_s1 <- as.data.table( melt( as.matrix(tC_pc10_s1)))
tC_pc10_s2 <- as.data.table( melt( as.matrix(tC_pc10_s2)))
tC_pc10_s3 <- as.data.table( melt( as.matrix(tC_pc10_s3)))

tC_pc15_s1 <- as.data.table( melt( as.matrix(tC_pc15_s1)))
tC_pc15_s2 <- as.data.table( melt( as.matrix(tC_pc15_s2)))
tC_pc15_s3 <- as.data.table( melt( as.matrix(tC_pc15_s3)))

tC_pc00_s1 <- tC_pc00_s1[, .( Protein_1 = as.character(Var1), Protein_2 = as.character(Var2), tC_pc00_s1 = value ) ]

tC_pc05_s1 <- tC_pc05_s1[, .( Protein_1 = as.character(Var1), Protein_2 = as.character(Var2), tC_pc05_s1 = value ) ]
tC_pc05_s2 <- tC_pc05_s2[, .( Protein_1 = as.character(Var1), Protein_2 = as.character(Var2), tC_pc05_s2 = value ) ]
tC_pc05_s3 <- tC_pc05_s3[, .( Protein_1 = as.character(Var1), Protein_2 = as.character(Var2), tC_pc05_s3 = value ) ]

tC_pc10_s1 <- tC_pc10_s1[, .( Protein_1 = as.character(Var1), Protein_2 = as.character(Var2), tC_pc10_s1 = value ) ]
tC_pc10_s2 <- tC_pc10_s2[, .( Protein_1 = as.character(Var1), Protein_2 = as.character(Var2), tC_pc10_s2 = value ) ]
tC_pc10_s3 <- tC_pc10_s3[, .( Protein_1 = as.character(Var1), Protein_2 = as.character(Var2), tC_pc10_s3 = value ) ]

tC_pc15_s1 <- tC_pc15_s1[, .( Protein_1 = as.character(Var1), Protein_2 = as.character(Var2), tC_pc15_s1 = value ) ]
tC_pc15_s2 <- tC_pc15_s2[, .( Protein_1 = as.character(Var1), Protein_2 = as.character(Var2), tC_pc15_s2 = value ) ]
tC_pc15_s3 <- tC_pc15_s3[, .( Protein_1 = as.character(Var1), Protein_2 = as.character(Var2), tC_pc15_s3 = value ) ]

tC_pc00_s1 <- tC_pc00_s1[ Protein_1 > Protein_2 ]    # Removes self-comparisons and duplicate pairs

tC_pc05_s1 <- tC_pc05_s1[ Protein_1 > Protein_2 ]    # Removes self-comparisons and duplicate pairs
tC_pc05_s2 <- tC_pc05_s2[ Protein_1 > Protein_2 ]    # Removes self-comparisons and duplicate pairs
tC_pc05_s3 <- tC_pc05_s3[ Protein_1 > Protein_2 ]    # Removes self-comparisons and duplicate pairs

tC_pc10_s1 <- tC_pc10_s1[ Protein_1 > Protein_2 ]    # Removes self-comparisons and duplicate pairs
tC_pc10_s2 <- tC_pc10_s2[ Protein_1 > Protein_2 ]    # Removes self-comparisons and duplicate pairs
tC_pc10_s3 <- tC_pc10_s3[ Protein_1 > Protein_2 ]    # Removes self-comparisons and duplicate pairs

tC_pc15_s1 <- tC_pc15_s1[ Protein_1 > Protein_2 ]    # Removes self-comparisons and duplicate pairs
tC_pc15_s2 <- tC_pc15_s2[ Protein_1 > Protein_2 ]    # Removes self-comparisons and duplicate pairs
tC_pc15_s3 <- tC_pc15_s3[ Protein_1 > Protein_2 ]    # Removes self-comparisons and duplicate pairs

setkey(tC_pc00_s1, Protein_1, Protein_2)             # Set keys for merging

setkey(tC_pc05_s1, Protein_1, Protein_2)             # Set keys for merging
setkey(tC_pc05_s2, Protein_1, Protein_2)             # Set keys for merging
setkey(tC_pc05_s3, Protein_1, Protein_2)             # Set keys for merging

setkey(tC_pc10_s1, Protein_1, Protein_2)             # Set keys for merging
setkey(tC_pc10_s2, Protein_1, Protein_2)             # Set keys for merging
setkey(tC_pc10_s3, Protein_1, Protein_2)             # Set keys for merging

setkey(tC_pc15_s1, Protein_1, Protein_2)             # Set keys for merging
setkey(tC_pc15_s2, Protein_1, Protein_2)             # Set keys for merging
setkey(tC_pc15_s3, Protein_1, Protein_2)             # Set keys for merging

DT1 <- merge(tC_pc00_s1, tC_pc05_s1)                 # Combine the treeClust distances   
DT1 <- merge(       DT1, tC_pc05_s2)                 # Combine the treeClust distances
DT1 <- merge(       DT1, tC_pc05_s3)                 # Combine the treeClust distances
DT1 <- merge(       DT1, tC_pc10_s1)                 # Combine the treeClust distances
DT1 <- merge(       DT1, tC_pc10_s2)                 # Combine the treeClust distances
DT1 <- merge(       DT1, tC_pc10_s3)                 # Combine the treeClust distances
DT1 <- merge(       DT1, tC_pc15_s1)                 # Combine the treeClust distances
DT1 <- merge(       DT1, tC_pc15_s2)                 # Combine the treeClust distances
DT1 <- merge(       DT1, tC_pc15_s3)                 # Combine the treeClust distances

DT1[, SimpleID_1 := gsub(";.+", "", Protein_1)][, SimpleID_1 := gsub("-.+", "", SimpleID_1)]  # Simplify protein IDs
DT1[, SimpleID_2 := gsub(";.+", "", Protein_2)][, SimpleID_2 := gsub("-.+", "", SimpleID_2)]  # Simplify protein IDs

setkey(DT1, SimpleID_1, SimpleID_2)                   # Re-set keys for merging


#### Random data drops: Precision recall analysis ####

# Get list of true and false positives (in same Reactome pathway or not) 
TP_FP_pairs <- fread("Reactome_TP_FP.csv")
names(TP_FP_pairs) <- c("SimpleID_1", "SimpleID_2", "Class")
setkey(TP_FP_pairs, SimpleID_1, SimpleID_2)    

# Merge labels with data
DT1 <- merge(DT1, TP_FP_pairs)
sum(is.na(DT1)) == 0     # There are no NAs in DT1, correct?

# Sample down to 10% TP
TP <- DT1[ Class == "TP" ]           # All true positives
FP <- DT1[ Class == "FP" ]           # All false positives
n_FP <- TP[, .N] * 9                 # Number of false positives we need
FP <- FP[ sample( FP[,.N] , n_FP) ]  # Restrict to random subset of false positives
DT1 <- rbindlist( list( TP, FP))     # Downsample DT1

# Create a random (scrambled) classifier
DT1[, Random := sample(tC_pc00_s1) ]   

# Get precision recall data using ROCR package
tC_pc00_s1 <- 1-DT1$tC_pc00_s1       # Get "inverse" treeClust distances (lower distance ~ higher score)

tC_pc05_s1 <- 1-DT1$tC_pc05_s1       # Get "inverse" treeClust distances (lower distance ~ higher score)
tC_pc05_s2 <- 1-DT1$tC_pc05_s2       # Get "inverse" treeClust distances (lower distance ~ higher score)
tC_pc05_s3 <- 1-DT1$tC_pc05_s3       # Get "inverse" treeClust distances (lower distance ~ higher score)

tC_pc10_s1 <- 1-DT1$tC_pc10_s1       # Get "inverse" treeClust distances (lower distance ~ higher score)
tC_pc10_s2 <- 1-DT1$tC_pc10_s2       # Get "inverse" treeClust distances (lower distance ~ higher score)
tC_pc10_s3 <- 1-DT1$tC_pc10_s3       # Get "inverse" treeClust distances (lower distance ~ higher score)

tC_pc15_s1 <- 1-DT1$tC_pc15_s1       # Get "inverse" treeClust distances (lower distance ~ higher score)
tC_pc15_s2 <- 1-DT1$tC_pc15_s2       # Get "inverse" treeClust distances (lower distance ~ higher score)
tC_pc15_s3 <- 1-DT1$tC_pc15_s3       # Get "inverse" treeClust distances (lower distance ~ higher score)

Random <- DT1$Random                 # Get the random classifier
labels <- DT1$Class                  # Get the class labels

pred <- prediction( predictions = list(tC_pc00_s1, tC_pc05_s1, tC_pc05_s2, tC_pc05_s3, tC_pc10_s1, tC_pc10_s2, 
                                       tC_pc10_s3, tC_pc15_s1, tC_pc15_s2, tC_pc15_s3, Random),
                    labels = list(labels, labels, labels, labels, labels,
                                  labels, labels, labels, labels, labels, labels),
                    label.ordering = c("FP", "TP"))
perf <- performance(pred, measure = "prec", x.measure = "rec")


#### Random data drops: Plotting ####

# (1) modified precision recall plot, with number of recall pairs (not rate) at log10
tC_pc00_s1 <- data.table(Recall = pred@tp[[1]],
                         Precision = perf@y.values[[1]],
                         Measure = "none")

tC_pc05_s1 <- data.table(Recall = pred@tp[[2]],
                         Precision = perf@y.values[[2]],
                         Measure = "5% (sample 1)")
tC_pc05_s2 <- data.table(Recall = pred@tp[[3]],
                         Precision = perf@y.values[[3]],
                         Measure = "5% (sample 2)")
tC_pc05_s3 <- data.table(Recall = pred@tp[[4]],
                         Precision = perf@y.values[[4]],
                         Measure = "5% (sample 3)")

tC_pc10_s1 <- data.table(Recall = pred@tp[[5]],
                         Precision = perf@y.values[[5]],
                         Measure = "10% (sample 1)")
tC_pc10_s2 <- data.table(Recall = pred@tp[[6]],
                         Precision = perf@y.values[[6]],
                         Measure = "10% (sample 2)")
tC_pc10_s3 <- data.table(Recall = pred@tp[[7]],
                         Precision = perf@y.values[[7]],
                         Measure = "10% (sample 3)")

tC_pc15_s1 <- data.table(Recall = pred@tp[[8]],
                         Precision = perf@y.values[[8]],
                         Measure = "15% (sample 1)")
tC_pc15_s2 <- data.table(Recall = pred@tp[[9]],
                         Precision = perf@y.values[[9]],
                         Measure = "15% (sample 2)")
tC_pc15_s3 <- data.table(Recall = pred@tp[[10]],
                         Precision = perf@y.values[[10]],
                         Measure = "15% (sample 3)")

    Random <- data.table(Recall = pred@tp[[11]],
                         Precision = perf@y.values[[11]],
                         Measure = "random")

pre_tp_dt <- rbindlist(  list(tC_pc00_s1, tC_pc05_s1, tC_pc05_s2, tC_pc05_s3, tC_pc10_s1, tC_pc10_s2, 
                              tC_pc10_s3, tC_pc15_s1, tC_pc15_s2, tC_pc15_s3, Random))

pre_tp_dt[ Recall == 1000 , mean(Precision) , by = Measure]  # Show the precision at a recall of 1,000 pairs

pre_tp_dt <- pre_tp_dt[ Recall > 500 ]                       # Drop low recall points because their pretty much random
pre_tp_dt <- pre_tp_dt[ sample( pre_tp_dt[,.N], 5000000 ) ]  # Randomly downsample to speed up loading times in Inkscape

p1 <- ggplot(pre_tp_dt, aes(x = Recall, y = Precision, colour = Measure))+
      geom_step(size=0.25)+
      scale_colour_manual(values = c("blueviolet", "blueviolet", "blueviolet", "darkcyan", "darkcyan", "darkcyan", "deeppink", "deeppink", "deeppink", "black", "grey50"))+
      scale_x_log10(breaks=c(1000,10000,100000), labels = comma)+
      annotation_logticks(sides="b", size=0.25, mid=unit(0.1,"cm"), long=unit(0.1,"cm"))+
      scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.2))+
      #geom_vline(xintercept = 1000, size=0.25, linetype = "dashed")+
      theme(panel.background = element_blank(), axis.text=element_text(size=5), axis.title=element_text(size=6),
            axis.ticks = element_line(size=0.25), axis.line = element_line(colour="black", size=0.25),
            legend.position = "none")

p1
ggsave("PR_random_drop.png", p1, dpi = 600,
       width=7, height=7, units=c("cm"))


# (2) Classic precision recall plot
tC_pc00_s1 <- data.table(Recall = perf@x.values[[1]],
                         Precision = perf@y.values[[1]],
                         Measure = "none")

tC_pc05_s1 <- data.table(Recall = perf@x.values[[2]],
                         Precision = perf@y.values[[2]],
                         Measure = "5% (sample 1)")
tC_pc05_s2 <- data.table(Recall = perf@x.values[[3]],
                         Precision = perf@y.values[[3]],
                         Measure = "5% (sample 2)")
tC_pc05_s3 <- data.table(Recall = perf@x.values[[4]],
                         Precision = perf@y.values[[4]],
                         Measure = "5% (sample 3)")

tC_pc10_s1 <- data.table(Recall = perf@x.values[[5]],
                         Precision = perf@y.values[[5]],
                         Measure = "10% (sample 1)")
tC_pc10_s2 <- data.table(Recall = perf@x.values[[6]],
                         Precision = perf@y.values[[6]],
                         Measure = "10% (sample 2)")
tC_pc10_s3 <- data.table(Recall = perf@x.values[[7]],
                         Precision = perf@y.values[[7]],
                         Measure = "10% (sample 3)")

tC_pc15_s1 <- data.table(Recall = perf@x.values[[8]],
                         Precision = perf@y.values[[8]],
                         Measure = "15% (sample 1)")
tC_pc15_s2 <- data.table(Recall = perf@x.values[[9]],
                         Precision = perf@y.values[[9]],
                         Measure = "15% (sample 2)")
tC_pc15_s3 <- data.table(Recall = perf@x.values[[10]],
                         Precision = perf@y.values[[10]],
                         Measure = "15% (sample 3)")

Random <- data.table(Recall = perf@x.values[[11]],
                     Precision = perf@y.values[[11]],
                     Measure = "random")

pre_rec_dt <-  rbindlist(  list(tC_pc00_s1, tC_pc05_s1, tC_pc05_s2, tC_pc05_s3, tC_pc10_s1, tC_pc10_s2, 
                                tC_pc10_s3, tC_pc15_s1, tC_pc15_s2, tC_pc15_s3, Random))

pre_rec_dt <- pre_rec_dt[ Recall > 0.005 ]                      # Drop low recall points because their pretty much random
pre_rec_dt <- pre_rec_dt[ sample( pre_rec_dt[,.N], 5000000 ) ]  # Randomly downsample to speed up loading times in Inkscape

p2 <- ggplot(pre_rec_dt, aes(x = Recall, y = Precision, colour = Measure))+
      geom_step(size=0.25)+
      scale_colour_manual(values = c("blueviolet", "blueviolet", "blueviolet", "darkcyan", "darkcyan", "darkcyan", "deeppink", "deeppink", "deeppink", "black", "grey50"))+
      scale_x_continuous(limits=c(0,1), breaks=seq(0,1,0.2))+
      scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.2))+
      theme(panel.background = element_blank(), axis.text=element_text(size=5), axis.title=element_text(size=6),
            axis.ticks = element_line(size=0.25), axis.line = element_line(colour="black", size=0.25),
            legend.position = "none")

p2
ggsave("PR_random_drop_2.png", p2, dpi = 600,
       width=7, height=7, units=c("cm"))

