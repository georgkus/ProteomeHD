# Load the required libraries
library(data.table); library(ggplot2); library(egg); library(dbscan)


#### Co-regulation score distribution (incl. histogram plot) ####

tC <- fread("coregulation_scores.csv")                                # Load the coregulation score (treeClust + TOM)
tC[, SimpleID_1 := gsub(";.+", "", Protein_1)][, SimpleID_1 := gsub("-.+", "", SimpleID_1)]  # Simplify protein 1 IDs 
tC[, SimpleID_2 := gsub(";.+", "", Protein_2)][, SimpleID_2 := gsub("-.+", "", SimpleID_2)]  # Simplify protein 2 IDs 
tC[ SimpleID_1 >= SimpleID_2 , .N] == tC[, .N]           # Test: are all protein pairs already sorted alphabetically?

# We arbitrarily define the highest scoring 0.5% of protein pairs as "co-regulated". 
n_coreg <- floor( tC[,.N] * 0.005 )                         # How many are these?
score_cut_off <- tC[ order(-coregulation_score)             # What score cut-off does that correspond to?           
                     ][ n_coreg , coregulation_score ]
  
# Plot a histogram
p1 <- ggplot(tC, aes( coregulation_score) )+
        geom_histogram( binwidth = 0.001 , boundary = score_cut_off, fill = "grey80", colour = NA )+
        geom_histogram( data = tC[ coregulation_score >= score_cut_off ],
                        binwidth = 0.001 , boundary = score_cut_off, fill = "magenta", colour = NA )+
        xlab("Co-regulation score")+
        ylab("Protein pairs")+
        geom_vline( xintercept = score_cut_off, colour = "magenta" , size = 0.25)+
        scale_x_continuous( breaks = seq(0,  0.4, 0.01))+
        theme(panel.background=element_blank(), panel.grid.major=element_blank(), panel.border=element_rect(fill=NA, colour="black", size=0.25),
              axis.text=element_text(size=5), axis.title=element_text(size=6), axis.ticks = element_line(size=0.25))

p1a <- p1 + coord_cartesian( xlim = c(0, 0.05) )+
            scale_y_continuous( breaks = seq(0, 1e+8, 1e+6))
  
p1b <- p1 + coord_cartesian( xlim = c(0, 0.05), ylim = c( 0, 20000))+
            scale_y_continuous( breaks = seq(0, 20000, 10000))

p1c <- ggarrange(p1a, p1b, nrow = 2 , heights = c(2,1))

ggsave("Score_distribution.pdf", p1c, width = 4.5, height = 6, units=c("cm"))


#### Co-regulation ~ function: prep protein complex & consecutive reaction data from Reactome ####

# Read human interactions dataset from Reactome and get necessary columns and rows
# The file homo_sapiens.interactions.txt is available directly from Reactome
DT <- fread("homo_sapiens.interactions.txt", header = FALSE)
DT <- DT[, .( Protein_1 = gsub("UniProt:", "", V1), Protein_2 = gsub("UniProt:", "", V4), Interaction_type = V7 ) ]
DT <- DT[ Protein_1 != Protein_2 ]   # Remove self-interactions

# Order protein 1 and 2 row-wise, so that A <-> B and B <-> A pairs can be detected as duplicates
DT[, Protein_1_sorted := ifelse(Protein_1 > Protein_2, Protein_1, Protein_2) ]
DT[, Protein_2_sorted := ifelse(Protein_1 < Protein_2, Protein_1, Protein_2) ]
DT <- DT[, .( Protein_1 = Protein_1_sorted, Protein_2 = Protein_2_sorted, Interaction_type) ]

# Get interaction pairs in which we are interested
# NR = neighbouring reaction pairs; CO = subunits of same protein complex
NR <- DT[ Interaction_type == "neighbouring_reaction" , paste(Protein_1, Protein_2, sep="_") ]   # Get neighbouring_reactions pairs
CO <- DT[ Interaction_type == "direct_complex" |
          Interaction_type == "indirect_complex"      , paste(Protein_1, Protein_2, sep="_") ]   # Get complex pairs
NR <- unique(NR)  # Remove duplicates
CO <- unique(CO)  # Remove duplicates


#### Co-regulation ~ function: prep subcellular location data from Uniprot ####

# Read in a file downloaded from Uniprot, containing all human SwissProt proteins and their subcellular annotations
uniprot <- fread("uniprot-all.tab")

# Use informative column names
colnames(uniprot) <- c("Uniprot_ID", "Protein_name", "Gene_name", "Localisation")

# Remove irrelavant rows and terms
uniprot <- uniprot[ Localisation != "" ]                                               # Remove proteins without localisation annotation
uniprot <- uniprot[ !grepl("Isoform", Localisation, ignore.case = TRUE) ]              # Remove proteins which contain isoform-specific localisation data
uniprot[, Localisation := gsub("SUBCELLULAR LOCATION: ", "", Localisation) ]           # Remove irrelevant term
uniprot[, Localisation := gsub(" Note=.+", "", Localisation) ]                         # Remove irrelevant term (the notes are always at the end of the string)
uniprot[, Localisation := gsub(" \\{.+\\}", "", Localisation) ]                        # Remove evidence codes
uniprot[, Localisation := gsub("; Single-pass.+membrane protein", "", Localisation) ]  # Remove irrelevant descriptions
uniprot[, Localisation := gsub("; Multi-pass membrane protein", "", Localisation) ]    # Remove irrelevant descriptions
uniprot[, Localisation := gsub("; Peripheral membrane protein", "", Localisation) ]    # Remove irrelevant descriptions
uniprot[, Localisation := gsub("; Lipid-anchor", "", Localisation) ]                   # Remove irrelevant descriptions
uniprot[, Localisation := gsub("; Lipid-anchor, GPI-anchor", "", Localisation) ]       # Remove irrelevant descriptions
uniprot[, Localisation := gsub(", GPI-anchor", "", Localisation) ]                     # Remove irrelevant descriptions
uniprot[, Localisation := gsub("Transmembrane protein: ", "", Localisation) ]          # Remove irrelevant descriptions
uniprot[, Localisation := gsub("Ubiquitin: ", "", Localisation) ]                      # Remove irrelevant descriptions
uniprot <- uniprot[ Localisation != "Membrane." ]                                      # Remove proteins which are only annotated as "Membrane" because we don't know which membrane that is

# Remove proteins which have more than one subcellular localisation
uniprot <- uniprot[ !grepl(".+\\..+\\.", Localisation) ]                               # Distinct localisations are separated by full stop

# Remove localisations with only 1 (cannot make pairs) annotation
exclude_loc <- names( which( table(uniprot$Localisation) == 1 ))
uniprot <- uniprot[ !Localisation %in% exclude_loc ]

# Combine proteins with same subcellular localisation in a pair-wise manner
starting_pairs <- aggregate(uniprot$Uniprot_ID, list(uniprot$Localisation), paste, collapse=",")
pairs <- data.frame(Protein_1=character(), Protein_2=character(), Subcellular_location=character())  # Set up empty df for the for loop
for (i in 1:nrow(starting_pairs) ) {
  x <- unlist( strsplit( starting_pairs$x[i], ",", fixed=TRUE) )
  y <- as.data.frame( t( combn(x, 2)))
  y$group <- starting_pairs$Group.1[i]
  colnames(y) <- c("Protein_1", "Protein_2", "Subcellular_location")
  pairs <- rbind(pairs,y)
}

# Sort the protein pairs alphabetically so that A <-> B and B <-> A duplicates could be removed
pairs <- as.data.table(pairs)
pairs <- pairs[, .( Protein_1 = as.character(Protein_1), Protein_2 = as.character(Protein_2) ) ]
pairs[, Protein_1_sorted := ifelse(Protein_1 > Protein_2, Protein_1, Protein_2) ]
pairs[, Protein_2_sorted := ifelse(Protein_1 < Protein_2, Protein_1, Protein_2) ]
pairs <- pairs[, .( Protein_1 = Protein_1_sorted, Protein_2 = Protein_2_sorted) ]
SU <- unique( paste( pairs$Protein_1, pairs$Protein_2, sep="_") )


#### Co-regulation ~ function: Plot percentage of functionally related protein pairs among coregulated / noncoregulated proteins  ####

# Find the co-regulated protein pairs and a same-sized random sample of the non-co-regulated ones
coregulated_pairs <-         tC[ coregulation_score >= score_cut_off , paste(SimpleID_1, SimpleID_2, sep="_") ]
noncoregula_pairs <- sample( tC[ coregulation_score <  score_cut_off , paste(SimpleID_1, SimpleID_2, sep="_") ] , length(coregulated_pairs) )

# Find the percentages per category
df1 <- data.frame( value = c( sum( coregulated_pairs %in% CO ) / (length(coregulated_pairs)/100),
                              sum( coregulated_pairs %in% NR ) / (length(coregulated_pairs)/100),
                              sum( coregulated_pairs %in% SU ) / (length(coregulated_pairs)/100) ),
                   type = c("complex subunits", "consecutive reaction", "subcellular location"),
                   variable = "co-regulated")

df2 <- data.frame( value = c( sum( noncoregula_pairs %in% CO ) / (length(noncoregula_pairs)/100),
                              sum( noncoregula_pairs %in% NR ) / (length(noncoregula_pairs)/100),
                              sum( noncoregula_pairs %in% SU ) / (length(noncoregula_pairs)/100) ),
                   type = c("complex subunits", "consecutive reaction", "subcellular location"),
                   variable = "non co-regulated")

df <- rbind(df1, df2)

# Create a barchart with the coregulated and non-coregulated protein pairs
p2 <- ggplot(df, aes(x = type, y = value, fill = variable))+
      geom_bar(stat="identity", position="dodge", colour="black")+
      scale_fill_manual(values = c("blue", NA))+
      xlab("delete this placeholder")+
      ylab("Protein pairs [%]")+
      theme(panel.background=element_blank(), panel.grid.major=element_blank(), panel.border=element_rect(fill=NA, colour="black", size=0.25),
            axis.text=element_text(size=5), axis.title=element_text(size=6), axis.ticks=element_line(size=0.25),
            legend.position="right", legend.title=element_blank(), legend.text = element_text(size=6))

ggsave("Coreg_funct_enrich_bars.pdf", p2, width = 6, height = 3.4, units=c("cm"))


#### Supplementary table: co-regulated protein pairs #### 

# Clear memory
rm( list = ls()[! ls() %in% c("tC", "score_cut_off")] )
gc()

# Turn co-regulation score back into a "distance" metric and define distance cut-off
tC[, coreg_distance := ( 1 - coregulation_score) ]
dist_cut_off <- ( 1 - score_cut_off )

# Turn the melted pairwise table back into a dist object
tC_m <- dcast( data = rbind( tC[, .(Protein_1, Protein_2, coreg_distance)],                             # This steps creates a "redundant" table...
                             tC[, .(Protein_1 = Protein_2, Protein_2 = Protein_1, coreg_distance)]),    # ... containing both A <-> B and B <-> pairs
                             Protein_1 ~ Protein_2 , value.var = "coreg_distance")                      # And this casts them into a matrix-shaped data.table
tC_m <- as.data.frame(tC_m)                 # Turn into data.frame
rownames(tC_m) <- tC_m$Protein_1            # Create rownames
tC_m$Protein_1 <- NULL                      # Drop original name column
tC_m <- as.matrix(tC_m)                     # Turn into numeric matrix
tc_dist <- as.dist(tC_m)                    # Turn into dist object
protein_IDs <- attr(tc_dist, "Labels")      # Extract protein IDs from dist object

# For each protein, calculate its co-regulation partners
NN <- frNN(tc_dist, dist_cut_off)                      # Fixed radius nearest neighbour search from dbscan package
NN <- lapply(NN$id, function(x){ protein_IDs[x] })     # Turn indexes into protein IDs
number_coreg <- sapply(NN, length)                     # Save the data (incl. not co-regulated proteins) for the plot later
NN <- as.data.table( melt(NN) )                        # Turn list into data.frame (drops proteins which are not co-regulated with any other protein)

# Rearrange and rename table
names(NN) <- c("Coreg_with_protein_X", "protein_X")
NN <- NN[, .(protein_X, Coreg_with_protein_X)]

# Write out table
fwrite(NN, "coregulated_protein_pairs.csv")


#### Co-regulation ~ function: Number of co-regulated proteins ####

plot.dt <- data.table(bins = cut(number_coreg, c(0,1,6,21,51,101,501), right = FALSE,
                                 labels = c("0", "1-5", "6-20", "21-50", "51-100", ">100")),
                      number_coreg = number_coreg)
bin_plot <- ggplot(plot.dt, aes(x = bins))+
              geom_bar(stat="count")+
              xlab("Number of co-regulated proteins")+
              ylab("Number of proteins")+
              theme(panel.background=element_blank(), panel.grid.major=element_blank(), panel.border=element_rect(fill=NA, colour="black", size=0.25),
                    axis.text=element_text(size=6), axis.title=element_text(size=7))

ggsave("Number_coreg_partners.pdf", bin_plot, width = 5, height = 3, units=c("cm"))
