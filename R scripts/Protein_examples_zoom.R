## In this script we analyse the distribution of four example proteins in the coregulation map: An uncharacterised microprotein,
## an uncharacterised protein without coregulation partners above threshold and two multifunctional proteins

# Load the required libraries
library(data.table); library(ggplot2); library(gridExtra); library(grid)

#### Prepare the data ####

# Load coregulation scores and annotate arbitrary threshold
tC <- fread("coregulation_scores.csv")
n_coreg <- floor( tC[,.N] * 0.005 )                              # We arbitrarily define the highest scoring 0.5% of protein pairs as "co-regulated". Assign those. 
score_cut_off <- tC[ order(-coregulation_score)                  # What score cut-off does that correspond to?           
                     ][ n_coreg , coregulation_score ]
tC[ coregulation_score >= score_cut_off , coregulated := TRUE ]  # Assign the term

# Load the tSNE map coordinates
SNE <- fread("tSNE_coordinates_and_annotation.csv")

# Add Gene names to the SNE map
ProHD <- fread("ProteomeHD_v1_1.csv")
ProHD <- ProHD[, .(Majority_protein_IDs, Gene_names)]
SNE <- merge(SNE, ProHD, by.x = "ProteinID", by.y = "Majority_protein_IDs")

# A function to "weight" the tSNE map by the treeClust score for a target protein
make_custom_SNE <- function(x){
  d1 <- tC[ Protein_1 == x , .(Target = Protein_1, Partner = Protein_2, coregulated) ]
  d2 <- tC[ Protein_2 == x , .(Target = Protein_2, Partner = Protein_1, coregulated) ]
  d3 <- rbind(d1, d2)
  d4 <- merge(SNE, d3, by.x = "ProteinID", by.y = "Partner", all.x = TRUE)
  return( d4[, .(ProteinID, coregulated) ] )
}


#### Get coregulation partners for the selected proteins of interest #### 

# Use the custom function to assign coregulation partners
TMEM256 <- make_custom_SNE( "Q8N2U0" )
HEATR5B <- make_custom_SNE( "Q9P2D3;Q9P2D3-3" )
  DDX3X <- make_custom_SNE( "O00571;O00571-2" )
    PHB <- make_custom_SNE( "P35232" )

# Modify colnames
names(TMEM256)[ names(TMEM256) == "coregulated" ] <- "coreg_with_TMEM256"
names(HEATR5B)[ names(HEATR5B) == "coregulated" ] <- "coreg_with_HEATR5B"
names(DDX3X)[ names(DDX3X) == "coregulated" ] <- "coreg_with_DDX3X"
names(PHB)[ names(PHB) == "coregulated" ] <- "coreg_with_PHB"

# Merge into one table, together with tSNE coordinates
SNE <- merge( SNE, TMEM256 )
SNE <- merge( SNE, HEATR5B )
SNE <- merge( SNE, DDX3X )
SNE <- merge( SNE, PHB )


#### Uncharacterised proteins: TMEM256 and HEATR5B #### 

# Set TMEM256 zoom region
TMEM256_zoom_x <- c(   5,  11 )
TMEM256_zoom_y <- c( -48, -42 )

# Set HEATR5B zoom region
HEATR5B_zoom_x <- c( 35, 37)
HEATR5B_zoom_y <- c( 3, 6)

# Global map with the zoomed regions annotated
p1 <- ggplot(SNE, aes( x = tSNE_x_dim, y = tSNE_y_dim ))+
       geom_point(shape = 16, size = 0.1, colour = "grey50")+
       geom_point(data = SNE[ coreg_with_TMEM256 == TRUE ], shape = 16, size = 0.2, alpha = 0.5, colour = "royalblue4")+
       geom_point(data = SNE[ coreg_with_HEATR5B == TRUE ], shape = 16, size = 0.2, alpha = 0.5, colour = "red")+
       annotate("rect", xmin = TMEM256_zoom_x[1], xmax = TMEM256_zoom_x[2], ymin = TMEM256_zoom_y[1], ymax = TMEM256_zoom_y[2], colour="black", fill=NA, size=0.25)+
       annotate("rect", xmin = HEATR5B_zoom_x[1], xmax = HEATR5B_zoom_x[2], ymin = HEATR5B_zoom_y[1], ymax = HEATR5B_zoom_y[2], colour="black", fill=NA, size=0.25)+
       theme(panel.background=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(),
             panel.border=element_rect(fill=NA, colour="black", size=0.25), panel.grid.major=element_blank(),
             axis.title=element_blank(), legend.position = "none", plot.margin = unit( c(0.5, 0.5, 0.5, 0.5), "mm"))
p1


## Map zoom for TMEM256

# In this zoom I want to show the enrichment for GO term "mitochondrion inner membrane" (GO:0005743). Download list of UniProt IDs
# associated with this term from QuickGO and load here
mito_IMM <- fread("QuickGO_mito_IMM.tsv")

# Simplify protein IDs for assigning subcellular location
SNE[, SimpleID := gsub(";.+", "", ProteinID) ][, SimpleID := gsub("-.+", "", SimpleID) ]

pTMEM256 <- ggplot(SNE, aes( x = tSNE_x_dim, y = tSNE_y_dim))+
              geom_point(size = 0.6, shape = 16, colour = "grey50" )+
              geom_point(data = SNE[ coreg_with_TMEM256 == TRUE], shape = 16, size = 0.6, colour = "royalblue4")+
              geom_point(data = SNE[ SimpleID %in% mito_IMM$`GENE PRODUCT ID` ], shape = 21, size = 3, colour = "darkorange")+
              geom_point(data = SNE[ Gene_names == "TMEM256" ], shape=16, size = 2, colour = "royalblue4")+
              geom_text( data = SNE[ Gene_names == "TMEM256" ], aes(label = Gene_names), size = 3, colour="royalblue4")+              
              xlim( TMEM256_zoom_x )+
              ylim( TMEM256_zoom_y )+
              theme(panel.background=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(),
                    panel.border=element_rect(fill=NA, colour="black", size=0.25), panel.grid.major=element_blank(),
                    axis.title=element_blank(), legend.position = "top", plot.margin = unit( c(0.5, 0.5, 0.5, 0.5), "mm"))
pTMEM256


# Map zoom for HEATR5B
pHEATR5B <- ggplot(SNE, aes( x = tSNE_x_dim, y = tSNE_y_dim))+
              geom_point(size = 0.6, shape = 16, colour = "grey50" )+
              geom_point(data = SNE[ coreg_with_HEATR5B == TRUE], shape = 16, size = 0.6, colour = "red")+
              geom_point(data = SNE[ Gene_names == "HEATR5B" ], shape=16, size = 2, colour = "red")+
              geom_text( aes(label = Gene_names), size = 3 , hjust = -0.1)+              
              xlim( HEATR5B_zoom_x )+
              ylim( HEATR5B_zoom_y )+
              theme(panel.background=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(),
                    panel.border=element_rect(fill=NA, colour="black", size=0.25), panel.grid.major=element_blank(),
                    axis.title=element_blank(), legend.position = "top", plot.margin = unit( c(0.5, 0.5, 0.5, 0.5), "mm"))
pHEATR5B


# Output combined plot
p1b <- arrangeGrob(pHEATR5B, pTMEM256)
pUnc <- arrangeGrob(p1, p1b, nrow = 1)
grid.draw(pUnc)

ggsave("TMEM256_HEATR5B_plot.pdf", pUnc, width=9, height=4.5, units=c("cm"))


#### Multifunctional proteins: DDX3X AND prohibitin ####

# Global map with the co-regulated proteins shown
p2 <- ggplot(SNE, aes( x = tSNE_x_dim, y = tSNE_y_dim ))+
      geom_point(shape = 16, size = 0.1, colour = "grey50")+
      geom_point(data = SNE[ coreg_with_DDX3X == TRUE ], shape = 16, size = 0.2, alpha = 0.5, colour = "orangered")+
      geom_point(data = SNE[ coreg_with_PHB == TRUE ], shape = 16, size = 0.2, alpha = 0.5, colour = "deepskyblue1")+
      theme(panel.background=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(),
            panel.border=element_rect(fill=NA, colour="black", size=0.25), panel.grid.major=element_blank(),
            axis.title=element_blank(), legend.position = "none", plot.margin = unit( c(0.5, 0.5, 0.5, 0.5), "mm"))
p2

# Save plot
ggsave("DDX3X_PHB.pdf", p2, width=4.5, height=4.5, units=c("cm"))
