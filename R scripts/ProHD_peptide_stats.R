## This script outputs ProteomeHD stats relating to protein coverage

# Load the required libraries
library(ggplot2); library(data.table); library(cowplot)

# Limit data.table threads to avoid known bug when reading / writing large tables
setDTthreads(1)

#### Load and pre-process the data ####

# Load the unprocessed proteinGroups file from MaxQuant (888 MB, available in compressed form, extract before running this script)
DT <- fread("proteinGroups.txt")

# Remove proteins that don't match quality criteria
DT <- DT[ `Only identified by site` != "+" ]
DT <- DT[ Reverse != "+" ]
DT <- DT[ `Potential contaminant` != "+" ]

# Remove proteins that have been detected in less than 4 experiments
ratio_count_cols <- names(DT)[ grepl("Ratio H/L count g", names(DT)) | grepl("Ratio M/L count g", names(DT)) ]    # These are the ratio counts for the SILAC ratios we used for ProteomeHD
ratio_counts <- apply( DT[, ratio_count_cols, with = FALSE ], 1, function(x){ sum(x > 0) } )                      # How often was each protein detected with at least 1 ratio count?
DT[ ratio_counts >= 95 , in_coreg_map := TRUE ]                                                                   # Assign which proteins have also been used for coregulation analysis
DT[ is.na(in_coreg_map), in_coreg_map := FALSE ]
DT <- DT[ ratio_counts >= 4 ]                                                                                     # These are the proteins we kept for ProteomeHD


#### Number of peptides per protein ####

# There are several different ways to look at this. I want to show both the average number of peptides per protein in each
# of the 294 experiments, as well as how many peptides accumulate per protein over these experiments
peptide_cols <- unique( gsub("Ratio ./L count ", "Peptides ", ratio_count_cols))                                              # Number of peptides (distinct peptide sequences) in each experiment
avg_peps_per_exp_all <- apply( DT[, peptide_cols, with = FALSE ], 2, function(x){ mean( x[ x > 0 ] )})                        # Average number of peptides across the detected proteins in each of the 264 experiments (note: 30 of these had M/L ratios, which are not counted as separate exps by MQ)
avg_peps_per_exp_cor <- apply( DT[ in_coreg_map == TRUE, peptide_cols, with = FALSE ], 2, function(x){ mean( x[ x > 0 ] )})   # The same for proteins in the coreg map

# Create a table for plotting
pep_plot_dt <- merge( data.table( exp = names(avg_peps_per_exp_all), avg_peps_per_exp_all = avg_peps_per_exp_all), 
                      data.table( exp = names(avg_peps_per_exp_cor), avg_peps_per_exp_cor = avg_peps_per_exp_cor) )

# Create plot labels
mean_in_map_pep_label <- paste("mean =", round( DT[ in_coreg_map == TRUE, mean( Peptides )], 2) )
   mean_all_pep_label <- paste("mean =", round( DT[, mean( Peptides )]                     , 2) )

# Plot the results: distribution per experiment and cumulative line
p1 <- ggplot( pep_plot_dt, aes( x = exp))+
        geom_bar(stat = "identity", aes( y = avg_peps_per_exp_cor ), fill = "midnightblue")+
        geom_bar(stat = "identity", aes( y = avg_peps_per_exp_all ), fill = "lightsteelblue")+
        # geom_hline( yintercept = DT[ in_coreg_map == TRUE, mean( Peptides )], size = 0.25, linetype = "dashed", colour = "midnightblue")+      # Since the cumulative mean is shown in the next panel (p2) already, I dropped it from here
        # geom_hline( yintercept = DT[, mean( Peptides )], size = 0.25, linetype = "dashed", colour = "lightsteelblue")+
        # annotate("text", x = 50, y = 39, label = mean_in_map_pep_label, size = 3, colour = "midnightblue")+      
        # annotate("text", x = 50, y = 27, label = mean_all_pep_label, size = 3, colour = "lightsteelblue")+      
        xlab("Experiments in ProteomeHD")+
        ylab("Avg. peptides per protein")+
        scale_y_continuous( breaks = seq(0,50,5), limits = c(0,15), expand = c(0,0))+
        theme(panel.background = element_blank(), panel.grid = element_blank(), 
              axis.text.x = element_blank(), axis.text.y=element_text(size=5), axis.title=element_text(size=6),
              axis.ticks.y = element_line(size=0.25), axis.ticks.x = element_blank(),
              axis.line = element_line(colour="black", size=0.25))

# Plot the distribution of peptides per protein (cumulative)
p2 <- ggplot(DT, aes(x = Peptides))+
        geom_histogram( binwidth = 1, fill = "lightsteelblue")+
        geom_histogram( data = DT[ in_coreg_map == TRUE ], binwidth = 1, fill = "midnightblue")+
        geom_vline(xintercept = DT[                     , mean(Peptides) ], size = 0.25, linetype = "dashed", colour = "lightsteelblue")+
        geom_vline(xintercept = DT[ in_coreg_map == TRUE, mean(Peptides) ], size = 0.25, linetype = "dashed", colour = "midnightblue")+
        annotate("text", x = 50, y = 350, label = mean_in_map_pep_label, size = 3, colour = "midnightblue")+      
        annotate("text", x = 20, y = 350, label = mean_all_pep_label, size = 3, colour = "lightsteelblue")+      
        scale_x_continuous( limits = c(0,150), breaks = c(1, seq(25,1000,25)), expand = c(0,0))+
        scale_y_continuous( limits = c(0,400), breaks = seq(0,1000,100), expand = c(0,0))+
        xlab("Number of identified peptides")+
        ylab("Protein count")+
        theme(panel.background = element_blank(), panel.grid = element_blank(), axis.text=element_text(size=5),
              axis.title=element_text(size=6), axis.ticks = element_line(size=0.25), 
              axis.line = element_line(colour="black", size=0.25), legend.position = "none")

# Plot the same cumulative distribution of peptides per MICROprotein
mean_in_map_pep_label_mp <- paste("mean =", round( DT[ `Mol. weight [kDa]` <= 15 & in_coreg_map == TRUE, mean(Peptides) ], 2) )
   mean_all_pep_label_mp <- paste("mean =", round( DT[ `Mol. weight [kDa]` <= 15                       , mean(Peptides) ], 2) )
 
p3 <- ggplot( DT[ `Mol. weight [kDa]` <= 15 ], aes(x = Peptides))+
        geom_histogram( binwidth = 1, fill = "lightsteelblue")+
        geom_histogram( data = DT[ `Mol. weight [kDa]` <= 15 & in_coreg_map == TRUE ], binwidth = 1, fill = "midnightblue")+
        geom_vline(xintercept = DT[ `Mol. weight [kDa]` <= 15                       , mean(Peptides) ], size = 0.25, linetype = "dashed", colour = "lightsteelblue")+
        geom_vline(xintercept = DT[ `Mol. weight [kDa]` <= 15 & in_coreg_map == TRUE, mean(Peptides) ], size = 0.25, linetype = "dashed", colour = "midnightblue")+
        annotate("text", x = 20, y = 50, label = mean_in_map_pep_label_mp, size = 3, colour = "midnightblue")+      
        annotate("text", x = 10, y = 50, label = mean_all_pep_label_mp, size = 3, colour = "lightsteelblue")+      
        scale_x_continuous( limits = c(0,50), breaks = c(1, seq(10,100,10)), expand = c(0,0))+
        scale_y_continuous( limits = c(0,65), breaks = seq(0,1000,20), expand = c(0,0))+
        xlab("Number of identified peptides")+
        ylab("Microprotein count")+
        theme(panel.background = element_blank(), panel.grid = element_blank(), axis.text=element_text(size=5),
              axis.title=element_text(size=6), axis.ticks = element_line(size=0.25), 
              axis.line = element_line(colour="black", size=0.25), legend.position = "none")


#### Number of ratio counts per protein ####

# Get the relevant data as above
avg_RC_per_exp_all <- apply( DT[                     , ratio_count_cols, with = FALSE ], 2, function(x){ mean( x[ x > 0 ] )})                        
avg_RC_per_exp_cor <- apply( DT[ in_coreg_map == TRUE, ratio_count_cols, with = FALSE ], 2, function(x){ mean( x[ x > 0 ] )})

# Create a table for plotting
RC_plot_dt <- merge( data.table( exp = names(avg_RC_per_exp_all), avg_RC_per_exp_all = avg_RC_per_exp_all), 
                     data.table( exp = names(avg_RC_per_exp_cor), avg_RC_per_exp_cor = avg_RC_per_exp_cor) )

# Plot the results: distribution per experiment and cumulative line
# Don't show means / medians here because they are too large
p4 <- ggplot( RC_plot_dt, aes( x = exp))+
        geom_bar(stat = "identity", aes( y = avg_RC_per_exp_cor ), fill = "midnightblue")+
        geom_bar(stat = "identity", aes( y = avg_RC_per_exp_all ), fill = "lightsteelblue")+
        xlab("Experiments in ProteomeHD")+
        ylab("Avg. ratio counts per protein")+
        scale_y_continuous( breaks = seq(0,50,10), limits = c(0,51), expand = c(0,0))+
        theme(panel.background = element_blank(), panel.grid = element_blank(), 
              axis.text.x = element_blank(), axis.text.y=element_text(size=5), axis.title=element_text(size=6),
              axis.ticks.y = element_line(size=0.25), axis.ticks.x = element_blank(),
              axis.line = element_line(colour="black", size=0.25))

# Create plot labels (include medians because there is a large difference here)
  mean_in_map_rat_label <- paste("mean =", round( DT[ in_coreg_map == TRUE,   mean(`Ratio H/L count`) ], 2) )
     mean_all_rat_label <- paste("mean =", round( DT[                     ,   mean(`Ratio H/L count`) ], 2) )
median_in_map_rat_label <- paste("median =", round( DT[ in_coreg_map == TRUE, median(`Ratio H/L count`) ], 2) )
   median_all_rat_label <- paste("median =", round( DT[                     , median(`Ratio H/L count`) ], 2) )

# Plot the distribution of ratio counts per protein (cumulative)
p5 <- ggplot(DT, aes(x = `Ratio H/L count`))+
        geom_histogram(                                    binwidth = 50, boundary = 0, fill = "lightsteelblue")+
        geom_histogram( data = DT[ in_coreg_map == TRUE ], binwidth = 50, boundary = 0, fill = "midnightblue")+
        geom_vline(xintercept = DT[                     , mean(`Ratio H/L count`) ], size = 0.25, linetype = "dashed", colour = "lightsteelblue")+
        geom_vline(xintercept = DT[ in_coreg_map == TRUE, mean(`Ratio H/L count`) ], size = 0.25, linetype = "dashed", colour = "midnightblue")+
        geom_vline(xintercept = DT[                     , median(`Ratio H/L count`) ], size = 0.25, linetype = "dashed", colour = "lightsteelblue")+
        geom_vline(xintercept = DT[ in_coreg_map == TRUE, median(`Ratio H/L count`) ], size = 0.25, linetype = "dashed", colour = "midnightblue")+
        annotate("text", x = 2500, y = 2000, label = mean_in_map_rat_label, size = 3, colour = "midnightblue")+      
        annotate("text", x = 1500, y = 2000, label = mean_all_rat_label   , size = 3, colour = "lightsteelblue")+      
        annotate("text", x = 1000, y = 2000, label = median_in_map_rat_label, size = 3, colour = "midnightblue")+      
        annotate("text", x = 300 , y = 2000, label = median_all_rat_label   , size = 3, colour = "lightsteelblue")+      
        scale_x_continuous( limits = c(0,4000), breaks = seq(0,10000,1000), expand = c(0,0))+
        scale_y_continuous( breaks = seq(0,3000,1000), expand = c(0,0))+
        xlab("Redundant peptide observations used\nfor quantitation (ratio counts)")+
        ylab("Protein count")+
        theme(panel.background = element_blank(), panel.grid = element_blank(), axis.text=element_text(size=5),
              axis.title=element_text(size=6), axis.ticks = element_line(size=0.25), 
              axis.line = element_line(colour="black", size=0.25), legend.position = "none")

# Create plot labels for MICROproteins (include medians because there is a large difference here)
  mean_in_map_rat_label_mp <- paste("mean =", round( DT[ `Mol. weight [kDa]` <= 15 & in_coreg_map == TRUE,   mean(`Ratio H/L count`) ], 2) )
     mean_all_rat_label_mp <- paste("mean =", round( DT[ `Mol. weight [kDa]` <= 15                       ,   mean(`Ratio H/L count`) ], 2) )
median_in_map_rat_label_mp <- paste("median =", round( DT[ `Mol. weight [kDa]` <= 15 & in_coreg_map == TRUE, median(`Ratio H/L count`) ], 2) )
   median_all_rat_label_mp <- paste("median =", round( DT[ `Mol. weight [kDa]` <= 15                       , median(`Ratio H/L count`) ], 2) )

# Plot the same cumulative distribution of peptides per MICROprotein
p6 <- ggplot(DT[ `Mol. weight [kDa]` <= 15 ], aes(x = `Ratio H/L count`))+
        geom_histogram(                                                                 binwidth = 50, boundary = 0, fill = "lightsteelblue")+
        geom_histogram( data = DT[  `Mol. weight [kDa]` <= 15 & in_coreg_map == TRUE ], binwidth = 50, boundary = 0, fill = "midnightblue")+
        geom_vline(xintercept = DT[ `Mol. weight [kDa]` <= 15                       , mean(`Ratio H/L count`) ], size = 0.25, linetype = "dashed", colour = "lightsteelblue")+
        geom_vline(xintercept = DT[ `Mol. weight [kDa]` <= 15 & in_coreg_map == TRUE, mean(`Ratio H/L count`) ], size = 0.25, linetype = "dashed", colour = "midnightblue")+
        geom_vline(xintercept = DT[ `Mol. weight [kDa]` <= 15                       , median(`Ratio H/L count`) ], size = 0.25, linetype = "dashed", colour = "lightsteelblue")+
        geom_vline(xintercept = DT[ `Mol. weight [kDa]` <= 15 & in_coreg_map == TRUE, median(`Ratio H/L count`) ], size = 0.25, linetype = "dashed", colour = "midnightblue")+
        annotate("text", x = 1500, y = 125, label = mean_in_map_rat_label_mp  , size = 3, colour = "midnightblue")+      
        annotate("text", x = 800 , y = 125, label = mean_all_rat_label_mp     , size = 3, colour = "lightsteelblue")+      
        annotate("text", x = 700 , y = 125, label = median_in_map_rat_label_mp, size = 3, colour = "midnightblue")+      
        annotate("text", x = 200 , y = 125, label = median_all_rat_label_mp   , size = 3, colour = "lightsteelblue")+      
        scale_x_continuous( limits = c(0,3000), breaks = seq(0,10000,1000), expand = c(0,0))+
        scale_y_continuous( breaks = seq(0,300,50), expand = c(0,0))+
        xlab("Redundant peptide observations used\nfor quantitation (ratio counts)")+
        ylab("Microprotein count")+
        theme(panel.background = element_blank(), panel.grid = element_blank(), axis.text=element_text(size=5),
              axis.title=element_text(size=6), axis.ticks = element_line(size=0.25), 
              axis.line = element_line(colour="black", size=0.25), legend.position = "none")


#### Sequence coverage per protein ####

# Get relevant data (as above)
seq_cov_cols <- unique( gsub("Ratio ./L count ", "Sequence coverage ", ratio_count_cols))                                    
seq_cov_cols <- paste( seq_cov_cols, " [%]", sep = "")
avg_sqc_per_exp_all <- apply( DT[                     , seq_cov_cols, with = FALSE ], 2, function(x){ mean( x[ x > 0 ] )})   
avg_sqc_per_exp_cor <- apply( DT[ in_coreg_map == TRUE, seq_cov_cols, with = FALSE ], 2, function(x){ mean( x[ x > 0 ] )})   

# Create a table for plotting
sqc_plot_dt <- merge( data.table( exp = names(avg_sqc_per_exp_all), avg_sqc_per_exp_all = avg_sqc_per_exp_all), 
                      data.table( exp = names(avg_sqc_per_exp_cor), avg_sqc_per_exp_cor = avg_sqc_per_exp_cor) )

# Create plot labels
mean_in_map_seq_label <- paste("mean =", round( DT[ in_coreg_map == TRUE, mean( `Sequence coverage [%]` )], 2) )
   mean_all_seq_label <- paste("mean =", round( DT[                     , mean( `Sequence coverage [%]` )], 2) )

# Plot the results: distribution per experiment and cumulative line
p7 <- ggplot( sqc_plot_dt, aes( x = exp))+
        geom_bar(stat = "identity", aes( y = avg_sqc_per_exp_cor ), fill = "midnightblue")+
        geom_bar(stat = "identity", aes( y = avg_sqc_per_exp_all ), fill = "lightsteelblue")+
        # geom_hline( yintercept = DT[ in_coreg_map == TRUE, mean( `Sequence coverage [%]` )], size = 0.25, linetype = "dashed", colour = "midnightblue")+       # Since the cumulative means are shown in p8 already, drop them from here
        # geom_hline( yintercept = DT[                     , mean( `Sequence coverage [%]` )], size = 0.25, linetype = "dashed", colour = "lightsteelblue")+
        # annotate("text", x = 50, y = 60, label = mean_in_map_seq_label, size = 3, colour = "midnightblue")+      
        # annotate("text", x = 50, y = 45, label = mean_all_seq_label   , size = 3, colour = "lightsteelblue")+      
        xlab("Experiments in ProteomeHD")+
        ylab("Avg. protein sequence coverage in %")+
        scale_y_continuous( breaks = seq(0,65,5), limits = c(0,35), expand = c(0,0))+
        theme(panel.background = element_blank(), panel.grid = element_blank(), 
              axis.text.x = element_blank(), axis.text.y=element_text(size=5), axis.title=element_text(size=6),
              axis.ticks.y = element_line(size=0.25), axis.ticks.x = element_blank(),
              axis.line = element_line(colour="black", size=0.25))

# Plot the distribution of peptides per protein (cumulative)
p8 <- ggplot(DT, aes(x = `Sequence coverage [%]`))+
        geom_histogram(                                    binwidth = 2.5, fill = "lightsteelblue")+
        geom_histogram( data = DT[ in_coreg_map == TRUE ], binwidth = 2.5, fill = "midnightblue")+
        geom_vline(xintercept = DT[                     , mean(`Sequence coverage [%]`) ], size = 0.25, linetype = "dashed", colour = "lightsteelblue")+
        geom_vline(xintercept = DT[ in_coreg_map == TRUE, mean(`Sequence coverage [%]`) ], size = 0.25, linetype = "dashed", colour = "midnightblue")+
        annotate("text", x = 60, y = 450, label = mean_in_map_seq_label, size = 3, colour = "midnightblue")+      
        annotate("text", x = 45, y = 450, label = mean_all_seq_label   , size = 3, colour = "lightsteelblue")+      
        scale_x_continuous( limits = c(0,100), breaks = seq(0,100,25), expand = c(0,0))+
        scale_y_continuous( limits = c(0,500), breaks = seq(0,500,100), expand = c(0,0))+
        xlab("% sequence coverage in ProteomeHD")+
        ylab("Protein count")+
        theme(panel.background = element_blank(), panel.grid = element_blank(), axis.text=element_text(size=5),
              axis.title=element_text(size=6), axis.ticks = element_line(size=0.25), 
              axis.line = element_line(colour="black", size=0.25), legend.position = "none")

# Plot the same cumulative distribution of peptides per MICROprotein
mean_in_map_seq_label_mp <- paste("mean =", round( DT[ `Mol. weight [kDa]` <= 15 & in_coreg_map == TRUE, mean( `Sequence coverage [%]` )], 2) )
   mean_all_seq_label_mp <- paste("mean =", round( DT[ `Mol. weight [kDa]` <= 15                       , mean( `Sequence coverage [%]` )], 2) )

p9 <- ggplot(DT[ `Mol. weight [kDa]` <= 15 ], aes(x = `Sequence coverage [%]`))+
        geom_histogram(                                                                binwidth = 2.5, fill = "lightsteelblue")+
        geom_histogram( data = DT[ `Mol. weight [kDa]` <= 15 & in_coreg_map == TRUE ], binwidth = 2.5, fill = "midnightblue")+
        geom_vline(xintercept = DT[ `Mol. weight [kDa]` <= 15                       , mean(`Sequence coverage [%]`) ], size = 0.25, linetype = "dashed", colour = "lightsteelblue")+
        geom_vline(xintercept = DT[ `Mol. weight [kDa]` <= 15 & in_coreg_map == TRUE, mean(`Sequence coverage [%]`) ], size = 0.25, linetype = "dashed", colour = "midnightblue")+
        annotate("text", x = 75, y = 35, label = mean_in_map_seq_label_mp, size = 3, colour = "midnightblue")+      
        annotate("text", x = 60, y = 35, label = mean_all_seq_label_mp   , size = 3, colour = "lightsteelblue")+      
        scale_x_continuous( limits = c(0,100), breaks = seq(0,100,25), expand = c(0,0))+
        scale_y_continuous( limits = c(0,40), breaks = seq(0,40,10), expand = c(0,0))+
        xlab("% sequence coverage in ProteomeHD")+
        ylab("Microprotein count")+
        theme(panel.background = element_blank(), panel.grid = element_blank(), axis.text=element_text(size=5),
              axis.title=element_text(size=6), axis.ticks = element_line(size=0.25), 
              axis.line = element_line(colour="black", size=0.25), legend.position = "none")


#### Combine and output stats plot ####

# Use cowplot to combine the plots into a figure
plot3by3 <- plot_grid( p1, p2, p3, p4, p5, p6, p7, p8, p9, align = "hv", nrow = 3, rel_widths = c(2.2,2,1.5) )

# Save combined plot
save_plot("ProHD_peptide_stats.pdf", plot3by3, nrow = 3, ncol = 3, base_height = 1.43, base_width = 2.31)


#### Annotate the microprotein MP68 with its peptides ####

## The following script was run on our server to narrow down the 10 GB evidence file to the unique entries
## that are relevant for this analysis (a compressed copy of this file is available)
    # library(data.table)
    # setDTthreads(1)
    # evi <- fread("evidence.txt")
    # mysubset <- evi[, .(Sequence, `Protein group IDs`, `Peptide ID`, Experiment)]
    # mysubset <- unique(mysubset)
    # fwrite(mysubset, "evidence_subset.csv", showProgress = TRUE)

# Load the "evidence"
evidence_subset <- fread("evidence_subset.csv")
names(evidence_subset) <- gsub(" ", "_", names(evidence_subset))

# Print peptide info for the microproteins MP68
DT[ `Gene names` == "MP68" , id ]
DT[ `Gene names` == "MP68" , Peptides ]
DT[ `Gene names` == "MP68" , `Majority protein IDs` ]
DT[ `Gene names` == "MP68" , `Mol. weight [kDa]` ]
DT[ `Gene names` == "MP68" , `Ratio H/L count` ]
DT[ `Gene names` == "MP68" , `Sequence coverage [%]` ]

# Find the peptides associated with this protein
MP68_peptides <- evidence_subset[ Protein_group_IDs ==  DT[ `Gene names` == "MP68" , id ] ]

# Print peptide sequence and occurrence
MP68_peptides[, .N, Sequence]



