## In this script we analyse the distribution of uncharacterised proteins and how ProteomeHD can be used
## to annotated them

# Load the required libraries
library(data.table); library(ggplot2); library(readxl)

#### Prepare the data ####

# Load coregulation scores and simplify IDs
tC <- fread("coregulation_scores.csv")
tC[, SimpleID_1 := gsub(";.+", "", Protein_1)][, SimpleID_1 := gsub("-.+", "", SimpleID_1)]  # Simplify protein 1 IDs 
tC[, SimpleID_2 := gsub(";.+", "", Protein_2)][, SimpleID_2 := gsub("-.+", "", SimpleID_2)]  # Simplify protein 2 IDs 
tC[  SimpleID_1 >= SimpleID_2 , .N] == tC[, .N]          # Test: are all protein pairs already sorted alphabetically?

n_coreg <- floor( tC[,.N] * 0.005 )                              # We arbitrarily define the highest scoring 0.5% of protein pairs as "co-regulated". Assign those. 
score_cut_off <- tC[ order(-coregulation_score)                  # What score cut-off does that correspond to?           
                     ][ n_coreg , coregulation_score ]
tC[ coregulation_score >= score_cut_off , coregulated := TRUE ]  # Assign the term

# Load Uniprot annotation data (all reviewed human SwissProt protein IDs where downloaded with Annotation score and molecular weight)
uni <- fread("SwissProt_HS_AnnotScore_Mass.tab")                                # Load data
uni <- uni[, .(Protein_ID = Entry, Annotation_score = Annotation, Mass)]        # Keep and rename relevant columns
uni[, Mass := gsub(",", "", Mass, fixed = TRUE)][, Mass := as.numeric(Mass)]    # Make Mass a numeric feature
uni[, Annotation_score := gsub(" out of 5", "/5", Annotation_score) ]           # Simplify scores


#### Annotation score distribution of the proteins in the co-regulation map ####

protein_isoforms <- unique( c( tC$Protein_1,  tC$Protein_2))   # Get all proteins (incl. isoforms) in co-regulation map
protein_genes <- gsub(";.+", "", protein_isoforms)             # Remove isoform information
protein_genes <- gsub("-.+", "", protein_genes)
counting_dt <- data.table(protein_isoforms, protein_genes)     # Combine into data table
counting_dt <- merge(counting_dt, uni, by.x = "protein_genes", # Add annotation score to our proteins
                                       by.y = "Protein_ID")
uncharacterised <- counting_dt[ Annotation_score != "4/5" &    # Define uncharacterised genes / proteins
                                Annotation_score != "5/5" ]    # (as those having an annotation score of 3 or lower)

length(unique(uncharacterised$protein_genes))          # This is the number of uncharacterised genes
length(unique(uncharacterised$protein_isoforms))       # This is the number of uncharacterised proteins, incl. isoforms


# Pie chart for uncharacterised content of co-regulation map
p1 <- ggplot(counting_dt, aes(x = 1, fill = Annotation_score))+
             geom_bar()+
             ylab("# proteins")+
             scale_fill_manual(values = c("orange4", "orange3", "orange", "dodgerblue1", "dodgerblue3"))+
             theme(panel.background = element_rect(fill = "white", colour = NA), panel.grid = element_blank(),
                   axis.text.y=element_text(size=6), axis.text.x = element_blank(), axis.title=element_text(size=6), axis.ticks.y = element_line(size=0.25),
                   axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.line.y = element_line(size=0.25),
                   legend.position="right", legend.text = element_text(size=6))

  
p1
ggsave("Annotation_scores.pdf", p1, width = 8, height=5, units = "cm")


#### Connectivity of uncharacterised proteins to well-characterised proteins ####

# Shortlist co-regulation events of uncharacterised proteins with characterised proteins
unc_DT <- tC[ (Protein_1 %in% uncharacterised$protein_isoforms | Protein_2 %in% uncharacterised$protein_isoforms ) &  # Pairs where one of the two proteins is uncharacterised
             !(Protein_1 %in% uncharacterised$protein_isoforms & Protein_2 %in% uncharacterised$protein_isoforms ) &  # But not both are uncharacterised (i.e. one partner is characterised)
              coregulated == TRUE ]                                                                                   # The two proteins are co-regulated

# Count number of fully characterised partners for each uncharacterised protein
unc_connectivity_counts <- integer()
for(i in uncharacterised$protein_isoforms){
    temp <- unc_DT[ Protein_1 == i | Protein_2 == i , .N]
    unc_connectivity_counts <- c(unc_connectivity_counts, temp)
}

length(unc_connectivity_counts)  # Just making sure the length corresponds to the number of uncharacterised proteins
mean(unc_connectivity_counts)    # What's the mean number of characterised co-regulation partners?


#### Connectivity of "cancer gene census" proteins to well-characterised proteins ####

# Load the cancer gene cencus converted to Uniprot IDs
CGC <- fread("Cancer_gene_census_Uniprot.csv")
CGC <- CGC[ Uniprot_ID %in% tC$SimpleID_1 | Uniprot_ID %in% tC$SimpleID_2 ,   # Get the subset that is in the co-regulation map
            .(Uniprot_ID, Role = `Role in Cancer`) ]

# Shortlist co-regulation events of CGC genes with other proteins
CGC_DT <- tC[ (SimpleID_1 %in% CGC$Uniprot_ID | SimpleID_2 %in% CGC$Uniprot_ID) &  # Pairs where one of the two proteins is in CGC
              coregulated == TRUE  ]                                               # The two proteins are co-regulated

# Count number of fully characterised partners for each CGC protein
CGC_connectivity_counts <- integer()
for(i in CGC$Uniprot_ID){
  temp_1 <- CGC_DT[ SimpleID_1 == i | SimpleID_2 == i , c(SimpleID_1, SimpleID_2)]  # Shortlist all co-regulation partners
  temp_2 <- temp_1[ temp_1 != i ]                                                   # Remove the current gene i itself 
  temp_3 <- sum(!temp_2 %in% uncharacterised$protein_genes)                         # How many of the partners are NOT uncharacterised
  CGC_connectivity_counts <- c(CGC_connectivity_counts, temp_3)
}

length(CGC_connectivity_counts)  # Just making sure the length corresponds to the number of CGC proteins
mean(CGC_connectivity_counts)    # What's the mean number of characterised co-regulation partners?


#### Connectivity of DisGeNET proteins to well-characterised proteins ####

# Load and prep DisGeNET data (these files can be downloaded from their website)
dgn <- fread("curated_gene_disease_associations.tsv")                      # Gene - disease associations
mapping <- fread("mapa_geneid_4_uniprot_crossref.tsv")                     # ID mapping file to Uniprot
dgn$Uniprot_ID <- mapping$UniProtKB[ match(dgn$geneId, mapping$GENEID) ]   # Append Uniprot IDs

# Get the DisGeNET genes that are in the co-regulation map
dgn_genes <- dgn[ Uniprot_ID %in% tC$SimpleID_1 | Uniprot_ID %in% tC$SimpleID_2 , unique(Uniprot_ID) ]

# Shortlist co-regulation events of DisGeNET genes with other proteins
dgn_DT <- tC[ (SimpleID_1 %in% dgn_genes | SimpleID_2 %in% dgn_genes) &     # Pairs where one of the two proteins is in DisGeNET
               coregulated == TRUE ]                                        # The two proteins are co-regulated

# Count number of fully characterised partners for each DisGeNET protein
dgn_connectivity_counts <- integer()
for(i in dgn_genes){
  temp_1 <- dgn_DT[ SimpleID_1 == i | SimpleID_2 == i , c(SimpleID_1, SimpleID_2)]  # Shortlist all co-regulation partners
  temp_2 <- temp_1[ temp_1 != i ]                                                   # Remove the current gene i itself 
  temp_3 <- sum(!temp_2 %in% uncharacterised$protein_genes)                         # How many of the partners are NOT uncharacterised
  dgn_connectivity_counts <- c(dgn_connectivity_counts, temp_3)
}

length(dgn_connectivity_counts)  # Just making sure the length corresponds to the number of dgn proteins
mean(dgn_connectivity_counts)    # What's the mean number of characterised co-regulation partners?


#### Make the connectivity barchart (uncharacterised, CGC, DisGeNET) ####

# Combine the data
combined_dist <- rbind( data.table( connectivity_count = unc_connectivity_counts, type = "Uncharacterized genes"),
                        data.table( connectivity_count = CGC_connectivity_counts, type = "Cancer gene census"),
                        data.table( connectivity_count = dgn_connectivity_counts, type = "DisGeNET"))

combined_dist[, bin := ifelse(connectivity_count == 0, "0",
                         ifelse(connectivity_count >= 1 & connectivity_count <= 5,  "1-5", 
                           ifelse(connectivity_count >= 6 & connectivity_count <= 20, "6-20",
                            ifelse(connectivity_count >= 21 & connectivity_count <= 50, "21-50", ">50"))))]

combined_dist[, bin := factor(bin, levels = c("0", "1-5", "6-20", "21-50", ">50"))]

# Make the plot
p2 <- ggplot(combined_dist, aes(x = type, fill = bin))+
        geom_bar(position = "fill")+
        scale_fill_brewer(palette = "Blues")+
        ylab("Fraction of genes")+
        annotate("text", x = "Uncharacterized genes", y = 0.5, label = sum(unc_connectivity_counts >= 1))+
        annotate("text", x = "Cancer gene census", y = 0.5, label = sum(CGC_connectivity_counts >= 1))+
        annotate("text", x = "DisGeNET", y = 0.5, label = sum(dgn_connectivity_counts >= 1))+
        scale_y_continuous(breaks=seq(0,1,0.1), expand = c(0,0), limits = c(0, 0.55))+
        theme(panel.background = element_rect(fill = "white", colour = NA), panel.grid = element_blank(),
              axis.text=element_text(size=6), axis.title=element_text(size=6), axis.ticks.y = element_line(size=0.25),
              axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.line.y = element_line(size=0.25),
              legend.position="right", legend.text = element_text(size=6))

p2
ggsave("Connectivity.pdf", p2, width = 5.5, height=4.3, units = "cm")


#### Microproteins (proteins < 15 kDa) per annotation score ####

# Get the percentage of microproteins per annotation category
pc_microproteins <- uni[, sum(Mass < 15000) / (.N/100) , by = Annotation_score ]

# Prep for plotting
pc_microproteins$Annotation_score <- gsub("/5", "", pc_microproteins$Annotation_score, fixed = TRUE)
pc_microproteins$V2 <- 100-pc_microproteins$V1
pc_microproteins <- pc_microproteins[, .(Annotation_score, Micro = V1, Macro = V2) ]
pc_microproteins <- melt(pc_microproteins, measure.vars = c("Micro", "Macro"))
pc_microproteins[, variable := factor(variable, levels = c("Macro", "Micro"))]

p3 <- ggplot(pc_microproteins, aes(x = Annotation_score, y = value, fill = variable))+
        geom_bar( stat = "identity" )+
        ylab("Proteins [%]")+
        xlab("Uniprot annotation score")+
        scale_y_continuous(breaks = seq(0,50,10), expand = c(0,0))+
        coord_cartesian(ylim=c(0,50))+
        theme(panel.background = element_rect(fill = "white", colour = NA), panel.border = element_rect(fill = NA, colour = "grey20", size=0.25),
              panel.grid.major = element_line(colour = "grey92", size=0.25), panel.grid.minor = element_line(colour = "grey92", size = 0.25),
              axis.text=element_text(size=6), axis.title=element_text(size=6), axis.ticks = element_line(size=0.25),
              axis.ticks.x = element_blank(), legend.position = "none")

p3
ggsave("Microprotein_annot_scores.pdf", p3, width = 3, height=4.75, units = "cm")


#### Percentage of microproteins in SwissProt, ProHD and BioPlex ####

# Load BioPlex 2.0 data and prep them for analysis
bio <- read_xlsx("nature22366-s2_Bioplex2.xlsx", sheet = 2)
bio <- as.data.table(bio)
bio <- bio[, .(Protein_A = `Uniprot A`, Protein_B = `Uniprot B`)]
bio <- bio[ Protein_A != "UNKNOWN" & Protein_B != "UNKNOWN" ]

bio[, Protein_A := gsub("-.+", "", Protein_A)]    # Remove isoform annotation
bio[, Protein_B := gsub("-.+", "", Protein_B)]

bio[, Protein_A_sorted := ifelse(Protein_A > Protein_B, Protein_A, Protein_B) ] # Sort the IDs row-wise
bio[, Protein_B_sorted := ifelse(Protein_A < Protein_B, Protein_A, Protein_B) ] # Sort the IDs row-wise
bio <- bio[, .(SimpleID_1 = Protein_A_sorted, SimpleID_2 = Protein_B_sorted) ]
bio[,.N] == nrow(unique(bio))                                 # There are no duplicates, true or false?
bio$PPI <- "Bioplex"


# Get numbers of microproteins among uncharacterised proteins
  bio_genes <- bio[                   , unique( c( SimpleID_1, SimpleID_2 )) ]     # Genes that have PPIs in BioPlex 2.0
coreg_genes <- tC[ coregulated == TRUE, unique( c( SimpleID_1, SimpleID_2 )) ]     # Genes that have co-regulation partners

uncharacterised_uniprot <- uni[ Annotation_score != "4/5" & Annotation_score != "5/5" ]   # All uncharacterised proteins

all_uni <- uncharacterised_uniprot[, .(uni = .N) , by = .(is_micro = Mass < 15000)]  # No. uncharact. genes > or < 15 kDa, in Uniprot overall
crg_uni <- uncharacterised_uniprot[ Protein_ID %in% coreg_genes,
                                     .(crg = .N) , by = .(is_micro = Mass < 15000)]  # No. uncharact. genes > or < 15 kDa, among those in co-regulation map
bio_uni <- uncharacterised_uniprot[ Protein_ID %in% bio_genes,
                                     .(bio = .N) , by = .(is_micro = Mass < 15000)]  # No. uncharact. genes > or < 15 kDa, among those that have PPIs in BioPlex2

# Create contingency tables
all_crg <- merge(all_uni, crg_uni)                   # Uniprot overall vs co-regulation
all_crg <- all_crg[ order(-is_micro) ]               # Ensure that microproteins (TRUE) are in top row
all_crg <- as.matrix( all_crg[, is_micro := NULL ])  # Drop the category column and turn into matrix

all_bio <- merge(all_uni, bio_uni)                   # Uniprot overall vs BioPlex
all_bio <- all_bio[ order(-is_micro) ]               # Ensure that microproteins (TRUE) are in top row
all_bio <- as.matrix( all_bio[, is_micro := NULL ])  # Drop the category column and turn into matrix

crg_bio <- merge(crg_uni, bio_uni)                   # Co-regulation vs BioPlex
crg_bio <- crg_bio[ order(-is_micro) ]               # Ensure that microproteins (TRUE) are in top row
crg_bio <- as.matrix( crg_bio[, is_micro := NULL ])  # Drop the category column and turn into matrix

# Perform Fisher's exact test to see which differences are significant
P_all_crg <- fisher.test(all_crg, alternative="greater")$p.value
P_all_bio <- fisher.test(all_bio, alternative="greater")$p.value
P_crg_bio <- fisher.test(crg_bio, alternative="greater")$p.value

# Get the percentages
pc <- data.table( `SwissProt`         = all_uni[ is_micro == TRUE , uni ] / all_uni[, sum(uni) ] * 100,
                  `Co-regulation map` = crg_uni[ is_micro == TRUE , crg ] / crg_uni[, sum(crg) ] * 100,
                  `BioPlex 2.0`       = bio_uni[ is_micro == TRUE , bio ] / bio_uni[, sum(bio) ] * 100)
pc <- melt(pc, measure.vars = 1:3)

# Plot the results
p4 <- ggplot(pc, aes(x = variable, y = value))+
        geom_bar(stat="identity")+
        ylab("[%] microproteins")+
        annotate("text", x = 1.5, y = 15, label = signif(P_all_crg, 3), size = 1)+
        annotate("text", x = 2.5, y = 15, label = signif( P_crg_bio, 3), size = 1)+      
        scale_y_continuous(limits=c(0,20), breaks = seq(0,20,5))+
        theme(panel.background = element_rect(fill = "white", colour = NA), panel.border = element_rect(fill = NA, colour = "grey20", size=0.25),
              panel.grid.major = element_line(colour = "grey92", size=0.25), panel.grid.minor = element_line(colour = "grey92", size = 0.25),
              axis.text=element_text(size=6), axis.title=element_text(size=6), axis.ticks = element_line(size=0.25),
              axis.title.x = element_blank(), panel.grid.major.x = element_blank())
p4
ggsave("Uncharacterised_microproteins_per_db.pdf", p4, width = 3, height=4, units = "cm")


