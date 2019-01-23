# Load necessary libraries
library(data.table); library(Rtsne); library(ggplot2); library(gridExtra)

#### Prepare co-regulation scores ####

# Load the coregulation scores (treeClust + TOM)
DT <- fread("coregulation_scores.csv")                                

# Turn co-regulation score back into a "distance" metric and log2-transform for better tSNE performance
DT[, coreg_distance := (1 - log2(coregulation_score)) ]

# Turn the melted pairwise table back into a dist object
DTm <- dcast( data = rbind( DT[, .(Protein_1, Protein_2, coreg_distance)],                             # These steps create a "redundant" table...
                            DT[, .(Protein_1 = Protein_2, Protein_2 = Protein_1, coreg_distance)]),    # ... containing both A <-> B and B <-> pairs
               Protein_1 ~ Protein_2 , value.var = "coreg_distance")                      # And this casts them into a matrix-shaped data.table
DTm <- as.data.frame(DTm)                 # Turn into data.frame
rownames(DTm) <- DTm$Protein_1            # Create rownames
DTm$Protein_1 <- NULL                     # Drop original name column
DTm <- as.dist( as.matrix( DTm ))         # Turn into numeric matrix then dist object
protein_IDs <- attr(DTm, "Labels")        # Extract protein IDs from dist object


#### Create tSNE map ####
set.seed(123)
SNE <- Rtsne(DTm, is_distance = TRUE, theta = 0.0, perplexity = 50, max_iter = 1500, verbose = TRUE)
SNE <- as.data.table( SNE$Y )
SNE[, ID := protein_IDs ]

pSNE <- ggplot(SNE, aes(x = V1, y = V2))+
          geom_point(shape = 16, size=0.2, alpha=0.8)+
          theme(panel.background=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(),
                panel.border=element_rect(fill=NA, colour="black", size=0.25), panel.grid.major=element_blank(),
                axis.title=element_blank(), plot.margin = unit( c(0.5, 0.5, 0.5, 0.5), "mm"))


#### Organelle annotation on tSNE map ####

# Read in annotation file
my_annotations <- fread("tSNE_map_annotation.csv")
my_annotations[ Organelle == "" , Organelle := NA ]                    # Replace empty strings with NAs
my_annotations[ Manual_annotation == "" , Manual_annotation := NA ]

# Merge with the map coordinates
SNE <- merge(SNE, my_annotations, by = "ID", all.x = TRUE )

# Create organelle plot
pORG <- ggplot(SNE, aes(x=V1, y=V2))+
         geom_point(data = SNE[ Organelle == "Cytoplasm"],     shape = 16, size=0.2, alpha=0.5, colour = "springgreen3")+
         geom_point(data = SNE[ Organelle == "Mitochondrion"], shape = 16, size=0.2, alpha=0.5, colour = "magenta")+
         geom_point(data = SNE[ Organelle == "ER"],            shape = 16, size=0.2, alpha=0.5, colour = "yellow3")+
         geom_point(data = SNE[ Organelle == "Nucleus"],       shape = 16, size=0.2, alpha=0.5, colour = "mediumblue")+
         geom_point(data = SNE[ Organelle == "Nucleolus"],     shape = 16, size=0.2, alpha=0.5, colour = "deepskyblue3")+
         geom_point(data = SNE[ Organelle == "Secreted"],      shape = 16, size=0.2, alpha=0.5, colour = "firebrick")+
         geom_point(data = SNE[ Organelle == "Ribosome"],      shape = 16, size=0.2, alpha=0.5, colour = "orange2")+
         annotate("rect", xmin=-60.5, xmax=-34.5, ymin=5, ymax=31, colour="black", fill=NA, linetype="dashed", size=0.25)+
         annotate("rect", xmin=5, xmax=27, ymin=-62, ymax=-42, colour="black", fill=NA, linetype="dashed", size=0.25)+
         annotate("rect", xmin=45, xmax=67, ymin=-23, ymax=3, colour="black", fill=NA, linetype="dashed", size=0.25)+
         theme(panel.background=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(),
               panel.border=element_rect(fill=NA, colour="black", size=0.25), panel.grid.major=element_blank(),
               axis.title=element_blank(), legend.position = "none", plot.margin = unit( c(0.5, 0.5, 0.5, 0.5), "mm"))


#### Sub-panel annotation on tSNE map ####

# Zoom area 1
pZ1 <- ggplot(SNE, aes(x=V1, y=V2))+
        geom_point(size = 0.1, shape = 16, colour = "grey50" )+
        geom_point(data = SNE[ Manual_annotation == "nucleolar rRNA processing"], shape = 16, size = 0.6, colour = "firebrick")+
        geom_point(data = SNE[ Manual_annotation == "mRNA splicing"],             shape = 16, size = 0.6, colour = "gold")+
        geom_point(data = SNE[ Manual_annotation == "EJC"],                       shape = 16, size = 0.6, colour = "blue")+
        geom_point(data = SNE[ Manual_annotation == "hnRNP"],                     shape = 16, size = 0.6, colour = "orange3")+
        geom_point(data = SNE[ Manual_annotation == "Nuclear pore"],              shape = 16, size = 0.6, colour = "green")+
        geom_point(data = SNE[ Manual_annotation == "Sm protein"],                shape = 16, size = 0.6, colour = "magenta")+
        geom_point(data = SNE[ Manual_annotation == "Exosome"],                   shape = 16, size = 0.6, colour = "purple")+
        scale_x_continuous( limits = c(-60.5, -34.5) , expand = c(0,0))+
        scale_y_continuous( limits = c(5, 31) , expand = c(0,0))+
        theme(panel.background=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(),
              panel.border=element_rect(fill=NA, colour="black", size=0.25), panel.grid.major=element_blank(),
              axis.title=element_blank(), legend.position = "none", plot.margin = unit( c(0.5, 0.5, 0.5, 0.5), "mm"))

# Zoom area 2
pZ2 <- ggplot(SNE, aes(x=V1, y=V2))+
        geom_point(size = 0.1, shape = 16, colour = "grey50" )+
        geom_point(data = SNE[ Manual_annotation == "ATP synthase"], shape = 16, size = 0.6, colour = "firebrick")+
        geom_point(data = SNE[ Manual_annotation == "complex I"],    shape = 16, size = 0.6, colour = "gold")+
        geom_point(data = SNE[ Manual_annotation == "complex II"],   shape = 16, size = 0.6, colour = "blue")+
        geom_point(data = SNE[ Manual_annotation == "complex III"],  shape = 16, size = 0.6, colour = "orange3")+
        geom_point(data = SNE[ Manual_annotation == "complex IV"],   shape = 16, size = 0.6, colour = "green")+
        geom_point(data = SNE[ Manual_annotation == "VDAC"],         shape = 16, size = 0.6, colour = "magenta")+
        geom_text(data = SNE[ Manual_annotation == "PEX11B"], aes(label = Manual_annotation), size = 0.2)+
        geom_text(data = SNE[ Manual_annotation == "ADP translocase"], aes(label = Manual_annotation), size = 0.2)+
        geom_text(data = SNE[ Manual_annotation == "Pi carrier"], aes(label = Manual_annotation), size = 0.2)+
        geom_text(data = SNE[ Manual_annotation == "SirT3"], aes(label = Manual_annotation), size = 0.2)+
        geom_text(data = SNE[ Manual_annotation == "ATPIF1"], aes(label = Manual_annotation), size = 0.2)+
        scale_x_continuous( limits = c(  5,  27) , expand = c(0,0))+
        scale_y_continuous( limits = c(-62, -42) , expand = c(0,0))+
        theme(panel.background=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(),
              panel.border=element_rect(fill=NA, colour="black", size=0.25), panel.grid.major=element_blank(),
              axis.title=element_blank(), legend.position = "none", plot.margin = unit( c(0.5, 0.5, 0.5, 0.5), "mm"))

# Zoom area 3
pZ3 <- ggplot(SNE, aes(x=V1, y=V2))+
        geom_point(size = 0.1, shape = 16, colour = "grey50" )+
        geom_point(data = SNE[ Manual_annotation == "Actins"],                shape = 16, size = 0.6, colour = "purple")+      
        geom_point(data = SNE[ Manual_annotation == "Actin regulation"],      shape = 16, size = 0.6, colour = "firebrick")+
        geom_point(data = SNE[ Manual_annotation == "Arp2/3 complex"],        shape = 16, size = 0.6, colour = "orangered")+
        geom_point(data = SNE[ Manual_annotation == "Rho GTPase"],            shape = 16, size = 0.6, colour = "dodgerblue")+
        geom_point(data = SNE[ Manual_annotation == "Rho GTPase regulation"], shape = 16, size = 0.6, colour = "firebrick")+
        geom_point(data = SNE[ Manual_annotation == "Myosins"],               shape = 16, size = 0.6, colour = "gold")+
        geom_point(data = SNE[ Manual_annotation == "Tropomyosin"],           shape = 16, size = 0.6, colour = "gold3")+      
        geom_point(data = SNE[ Manual_annotation == "G proteins"],            shape = 16, size = 0.6, colour = "green")+
        geom_point(data = SNE[ Manual_annotation == "Integrins"],             shape = 16, size = 0.6, colour = "forestgreen")+
        scale_x_continuous( limits = c( 45,  67) , expand = c(0,0))+
        scale_y_continuous( limits = c(-23,   3) , expand = c(0,0))+
        theme(panel.background=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(),
              panel.border=element_rect(fill=NA, colour="black", size=0.25), panel.grid.major=element_blank(),
              axis.title=element_blank(), legend.position = "none", plot.margin = unit( c(0.5, 0.5, 0.5, 0.5), "mm"))


#### Combine the plots and write out files ####

# Combine the subplots first
subplot <- arrangeGrob(pORG, pZ2, pZ1, pZ3)

# Combine with the SNE plot
final_plot <- arrangeGrob( pSNE, subplot, nrow = 1 )

# Save the plot
ggsave("tSNE_map.pdf", final_plot, width=18.3, height=9.15, units=c("cm"))

# Output the coordinates and the annotation as supplementary file
SNE <- SNE[, .(ProteinID = ID, tSNE_x_dim = V1, tSNE_y_dim = V2, Organelle_annotation = Organelle, Detailed_manual_annotation = Manual_annotation)]
fwrite(SNE, "tSNE_coordinates_and_annotation.csv")

