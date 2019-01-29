## This script outputs some simple stats about ProteomeHD v1.1

library(data.table); library(ggplot2); library(gridExtra)

# Load ProteomeHD
prohd <- fread("ProteomeHD_v1_1.csv")                             # Read in fast using data.table
prohd <- data.frame(prohd, row.names = "Majority_protein_IDs")    # Convert to data.frame; use protein IDs as rownames
prohd <- prohd[, grep("Ratio", colnames(prohd)) ]                 # Keep only columns with SILAC ratios

#### How many proteins per SILAC ratio? ####
pro_per_ratio <- apply(prohd, 2, function(x){ sum( !is.na(x) ) })                # Count proteins
pro_per_ratio <- data.table( Experiment = factor( names( pro_per_ratio )),       # Format as data.table
                             proteins_per_ratio = pro_per_ratio )    
pro_per_ratio[, Experiment := reorder(Experiment, proteins_per_ratio) ]          # Sort experiments by N detected proteins

p1 <- ggplot( pro_per_ratio, aes( x = Experiment, y = proteins_per_ratio ) )+
        geom_bar( stat = "identity" )+
        geom_hline( yintercept = mean(pro_per_ratio$proteins_per_ratio) , size = 0.25, linetype = "dashed", colour = "black")+
        annotate("text", label = paste( "mean =", round( mean(pro_per_ratio$proteins_per_ratio), 0)), x = 20, y = 4100, size = 2.5)+
        xlab("294 experiments (SILAC ratios) in ProteomeHD")+
        ylab("Number of quantified proteins")+
        scale_y_continuous( breaks = c(0,2000,4000,6000), limits = c(0,6100), expand = c(0,0))+
        theme(panel.background = element_blank(), panel.grid = element_blank(), 
              axis.text.x = element_blank(), axis.text.y=element_text(size=5), axis.title=element_text(size=6),
              axis.ticks.y = element_line(size=0.25), axis.ticks.x = element_blank(),
              axis.line = element_line(colour="black", size=0.25))

#### How many SILAC ratios per protein? ####
ratios_per_pro <- apply(prohd, 1, function(x){ sum( !is.na(x) ) })               # Count SILAC ratios per protein
ratios_per_pro <- data.table( ProteinID = names( ratios_per_pro ),               # Format as data.table
                              ratios_per_protein = ratios_per_pro )

p2 <- ggplot( ratios_per_pro, aes( ratios_per_protein , fill = ratios_per_protein >= 95 ))+
        geom_histogram( binwidth = 5, boundary = 4)+
        geom_vline( xintercept = ratios_per_pro[, mean(ratios_per_protein)], size = 0.25, linetype = "dashed", colour = "black" )+
        annotate("text", label = paste( "mean =", round( ratios_per_pro[, mean(ratios_per_protein)] , 0)), x = 120, y = 400, size = 2.5)+
        geom_vline( xintercept = ratios_per_pro[ ratios_per_protein >= 95 , mean(ratios_per_protein)], size = 0.25, linetype = "dashed", colour = "blue" )+
        annotate("text", label = paste( "mean =", round( ratios_per_pro[ ratios_per_protein >= 95 , mean(ratios_per_protein)] , 0)), x = 200, y = 400, size = 2.5)+
        xlab("SILAC ratios per protein")+
        ylab("Number of quantified proteins")+
        scale_fill_manual(values = c("red", "blue"))+
        scale_y_continuous( breaks = c(0,200,400,600), expand = c(0,0))+
        scale_x_continuous( breaks = seq(0,300,50), expand = c(0,0))+
        theme(panel.background = element_blank(), panel.grid = element_blank(), 
              axis.text=element_text(size=5), axis.title=element_text(size=6),
              axis.ticks = element_line(size=0.25),
              axis.line = element_line(colour="black", size=0.25), legend.position = "none")


#### Cumulative distribution: Which fraction of proteins have been detected in which fraction of experiments? ####

# Assess fraction of proteins at these ratio cut-offs
exp_pc_cutoffs <- seq( 10, 90, 10)                                             

all_counts <- ratios_per_pro[                         , .N, ratios_per_protein ]   # Get all counts
map_counts <- ratios_per_pro[ ratios_per_protein >= 95, .N, ratios_per_protein ]   # Get counts for proteins that are in coregulation map

all_prots <- all_counts[, sum(N) ]                                                 # N proteins in ProHD
map_prots <- map_counts[, sum(N) ]                                                 # N proteins in map

cumulative_fractions <- data.frame( exp_pc_cutoffs ,                            # Initialise result table
                                    pc_all_proteins = NA ,
                                    pc_map_proteins = NA)         

# For each percentage of ratios, count the percentage of proteins
for(i in 1:length(exp_pc_cutoffs)){                                             
  cumulative_fractions[i, "pc_all_proteins"] <- all_counts[ ratios_per_protein >= (294/100 * exp_pc_cutoffs[i]) , sum(N) ] / (all_prots / 100)
  cumulative_fractions[i, "pc_map_proteins"] <- map_counts[ ratios_per_protein >= (294/100 * exp_pc_cutoffs[i]) , sum(N) ] / (map_prots / 100)
}

# Melt for plotting
cumulative_fractions <- melt( cumulative_fractions , id.vars = "exp_pc_cutoffs" )

p3 <- ggplot( cumulative_fractions, aes( x = exp_pc_cutoffs, y = value, fill = variable))+
        facet_grid(~variable)+
        geom_bar(stat = "identity")+
        xlab("[%] experiments")+
        ylab("[%] proteins")+
        scale_fill_manual(values = c("grey50", "blue"))+
        scale_y_continuous( breaks = seq(0,100,10), limits = c(0,100), expand = c(0,0))+
        scale_x_continuous( breaks = exp_pc_cutoffs, limits = c(5,95), expand = c(0,0))+
        theme(panel.background = element_blank(), panel.grid.major.x = element_blank(),
              panel.grid.major.y = element_line(size = 0.25, linetype = "dashed", colour = "grey50"),
              axis.text=element_text(size=5), axis.title=element_text(size=6),
              axis.ticks = element_line(size=0.25), legend.position = "none",
              axis.line = element_line(colour="black", size=0.25))


#### Combine & print plot ####

p <- arrangeGrob(p1, p2, p3)

ggsave("ProHD_stats.pdf", p, width = 12, height = 12, units = "cm")

