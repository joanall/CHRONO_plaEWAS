#==============================================================================#
#                                                                              #
#               MATERNAL DIET MEALTIME AND PLACENTA DNA METHYLATION            #
#                             Analysis DMPs                                    #
#                                                                              #                                  
#                   Camille Lassale & Joana Llauradó Pont                      #         
#                               January 2025                                   #
#                                                                              #
#==============================================================================#

# DESCRIPTION: Read the the EWAS results after proceesing them and 
# make volcano and manthattna plot, as well as sumamry tbales. 



# [INPUT]: RDS File (list of files) storing a dataframe where rows are CpGs and columns: 
# - CpGs: chr/string CpG Id
# - Effect coefficient of the exposure of interest
# - SE  standard error of the mode
# - pvalue 


# OUTPUT: Volcano and Manhattan plot + CSV top (cutoff) CpGs



#==============================================================================#
#                                                                              #
#                              USER INPUT                                      #
#                                                                              #
#==============================================================================#

suggestive_cutoff <- 1E-05
array <- 'EPIC'


# Options: nighttime_fasting, first_meal, eating_occasions, eating_jetlag, last_meal (older versions)---- 

# Subfolder names "main": exposure_model
#  'model_1' 'model_2' 'model_3' (adj.first meal)
# for Nightfast, first meal and last meal ->  'model_1_sep_exp' 'model_2_sep_exp' 'model_3_sep_exp'}

# Subfolder names "sens": exposure_model
# if exposure jetlag: 
# 'model_noGDM' 'model_ethnicity' 'model_noconf'  'early_eaters', 'short_fasting'

# if exposure nighttime: 
# 'model_noGDM_adj_first_meal' 'model_ethnicity' 'model_noPC' 'model_noconf'  'early_eaters'

# if exposure first or last meal or eating occasions: 
# 'model_noGDM_adj_firstlast' 'model_ethnicity' 'model_noPC' 'model_noconf'  'early_eaters'

# Input folder 
# main_ewas_folder <- "results/2_quality_control/Bisc_14012025" #05022025
# input_ewas_file <- file.path(main_ewas_folder,  
# #"eating_occasions_model_noGDM_sep_exp/eating_occasions_Bisc_model_noGDM_sep_exp_05022025_adjusted.rds")
# "first_meal_model_1_sep_exp/first_meal_Bisc_sep_exp_14012025_adjusted.rds")

# --------------------------------------------------------------------------------------

# Input  
input_ewas_file <- "results/2_quality_control/Main_13062025/eating_jetlag_model_3_sep_exp/eating_jetlag_model_3_sep_exp_adjusted.rds"

# Output folder 
results_folder <- "results/3_analysis_dmps/Main_13062025/eating_jetlag_model_3_sep_exp/"


#==============================================================================#
#                       INIT  WORKING ENVIRONMENT                              #----
#==============================================================================#
library(ggplot2)
library(dplyr)

if (!dir.exists(results_folder)) {dir.create(results_folder, recursive = T)}
results_folder

#==============================================================================#
#                                                                              #
#                               ANALYSIS EWAS                                  #----
#                                                                              #
#==============================================================================#

cat("\n----------------------  ANALYSIS EWAS  ---------------------------------\n")

# Load the EWAS result and get df 
file <- input_ewas_file
cat('Reading file', basename(file),'\n')
results <- readRDS(file)

#------------------------------------------------------------------------------#
#                             GET TOPS CPGs                                    #----
#------------------------------------------------------------------------------#

cat("\n------------------------- Top CpGs  ----------------------------------\n")


top <- subset(results, results$pvalue <= suggestive_cutoff)

top <- top [, c('CpG', 'Effect', 'SE', 'pvalue', 'bf')]

save_to <- file.path(results_folder,  
                  sub(".rds", paste0("_", suggestive_cutoff , ".csv"), 
                  basename(file)))

if (nrow(top) != 0) {write.table(top, save_to, sep=";", row.names = FALSE)}


if ('pval_bacon' %in% colnames(results)) {

  top <- subset(results, results$pval_bacon <= suggestive_cutoff)
  
  save_to <- file.path(results_folder, 
                    sub(".rds", paste0(suggestive_cutoff, "_bacon.csv"),
                    basename(file)))
  
  if (nrow(top) != 0) {write.table(top, save_to, sep=";")}
  
}



#------------------------------------------------------------------------------#
#                             GET VOLCANO PLOTS                                #----
#------------------------------------------------------------------------------#

# Based on/Modified from  PACEplacenta analysis

# Plot coefficients of EWAS and pvalues 
volcano_plot <- ggplot(results, aes(Effect, -log10(pvalue))) +
  geom_point(alpha=0.6) +
  
  geom_hline(yintercept=-log10(suggestive_cutoff), col="#78c679") + # suggestive, green
  geom_hline(yintercept=-log10(0.05/nrow(results)), col="#e7298a") + # Bonferroni, magenta
  
  # Annotate CpG below the threshold
  geom_text(data=subset(results, pvalue < suggestive_cutoff), 
            aes(label=CpG), 
            vjust=0.5,  # Increase vertical adjustment
            hjust=-0.1,  # Increase horizontal adjustment
            size=5, 
            color="red") +  # Customize color and size as needed
  
  # Styling 
  theme_bw() +
  
  # Labeling
  xlab("Coefficient") + ylab(expression(-log[10](pval))) +
  
  # Optional title
  # ggtitle(title) +
  theme(axis.title=element_text(face="bold", size=12),
        plot.title=element_text(face="bold", size=18, hjust=0),
        plot.margin=margin(20, 20, 20, 20)) # Increase margin more to avoid cutoff

# Save with increased resolution and specific dimensions
save_to <- file.path(results_folder, sub(".rds", "_volcano.png", basename(file)))
save_to
ggsave(save_to, plot = volcano_plot, width = 12, height = 10, dpi = 300)

if ('pval_bacon' %in% colnames(results)) {
  
  # Plot coefficients of EWAS and pvalues 
  volcano_plot <- ggplot(results, aes(Effect, -log10(pval_bacon))) +
    geom_point(alpha=0.6) +
    
    geom_hline(yintercept=-log10(suggestive_cutoff), col="#78c679") + # suggestive, green
    geom_hline(yintercept=-log10(0.05/nrow(results)), col="#e7298a") + # Bonferroni, magenta
    
    # Annotate CpG below the threshold
    geom_text(data=subset(results, pval_bacon < suggestive_cutoff), 
              aes(label=CpG), 
              vjust=0.5,  # Increase vertical adjustment
              hjust=-0.1,  # Increase horizontal adjustment
              size=5, 
              color="red") +  # Customize color and size as needed
    
    # Styling 
    theme_bw() +
    
    # Labeling
    xlab("Coefficient") + ylab(expression(-log[10](pval))) +
    
    # Optional title
    # ggtitle(title) +
    theme(axis.title=element_text(face="bold", size=12),
          plot.title=element_text(face="bold", size=18, hjust=0),
          plot.margin=margin(20, 20, 20, 20)) # Increase margin more to avoid cutoff
  
  # Save with increased resolution and specific dimensions
  save_to <- file.path(results_folder, sub(".rds", "_volcano_bacon.png", basename(file)))
  save_to
  ggsave(save_to, plot = volcano_plot, width = 12, height = 10, dpi = 300)
  
}




#------------------------------------------------------------------------------#
#                             GET MANHATTAN PLOTS                              #----
#------------------------------------------------------------------------------#

cat("\n------------------------ Manhattan Plot  --------------------------------\n")

# Get annotation
if (array == "450K") {
  library('IlluminaHumanMethylation450kanno.ilmn12.hg19')
  annotation <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
} else if (array == "EPIC") {
  library('IlluminaHumanMethylationEPICanno.ilm10b2.hg19')
  annotation <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
}

annotation <- as.data.frame(annotation)

# Add annotation to results by merging on row names = CpGs
results <- merge(results, annotation, by.x = 'row.names', by.y = 'row.names')

# Get numeric for chromosome (e.g. from "chr6" to 6)
results$chromosome <- as.numeric(gsub('chr', "", results$chr))

# Add a column for -log10 p-values
results$logp <- -log10(results$pvalue)

# Prepare data to make plot
don <- results %>%
  # Compute chromosome size
  group_by(chromosome) %>%
  summarise(chr_len = max(as.numeric(pos))) %>%
  # Calculate cumulative position of each chromosome
  mutate(tot = cumsum(as.numeric(chr_len)) - as.numeric(chr_len)) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(results, by = c("chromosome" = "chromosome")) %>%
  # Add a cumulative position of each SNP
  arrange(chromosome, pos) %>%
  mutate(BPcum = as.numeric(pos) + tot)

axisdf <- don %>%
  group_by(chromosome) %>%
  summarize(center = (max(BPcum, na.rm = TRUE) + min(BPcum, na.rm = TRUE)) / 2)

# Define thresholds
bonferroni_threshold <- -log10(0.05 / nrow(results))
suggestive_threshold <- -log10(suggestive_cutoff)

# Create the Manhattan plot
manhattan_plot <- ggplot(don, aes(x = BPcum, y = logp)) +
  geom_point(aes(color = as.factor(chromosome)), alpha = 0.8, size = 1.3) +
  scale_color_manual(values = rep(c("grey", "black"), 22)) +
  scale_x_continuous(label = axisdf$chromosome, breaks = axisdf$center) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(don$logp) + 2)) +
  geom_hline(yintercept = bonferroni_threshold, color = "#e7298a", linetype = "dashed") +
  geom_hline(yintercept = suggestive_threshold, color = "#78c679", linetype = "dashed") +
  geom_text(data = subset(don, logp > suggestive_threshold), aes(label = CpG), vjust = -1, size = 2) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  labs(x = "Chromosome", y = "-log10(p)")

# Save the plot
save_to <- paste0(results_folder, '/', sub(".rds", "_manhattan.png", basename(file)))
ggsave(save_to, plot = manhattan_plot, width = 14, height = 7, dpi = 300)


if ('pval_bacon' %in% colnames(results)) {
  # results already annotated,just redo manhattan with pal_bacon 
  
  # Add a column for -log10 p-values
  results$logp <- -log10(results$pval_bacon)
  
  # Prepare data to make plot
  don <- results %>%
    # Compute chromosome size
    group_by(chromosome) %>%
    summarise(chr_len = max(as.numeric(pos))) %>%
    # Calculate cumulative position of each chromosome
    mutate(tot = cumsum(as.numeric(chr_len)) - as.numeric(chr_len)) %>%
    select(-chr_len) %>%
    # Add this info to the initial dataset
    left_join(results, by = c("chromosome" = "chromosome")) %>%
    # Add a cumulative position of each SNP
    arrange(chromosome, pos) %>%
    mutate(BPcum = as.numeric(pos) + tot)
  
  axisdf <- don %>%
    group_by(chromosome) %>%
    summarize(center = (max(BPcum, na.rm = TRUE) + min(BPcum, na.rm = TRUE)) / 2)
  
  # Define thresholds
  bonferroni_threshold <- -log10(0.05 / nrow(results))
  suggestive_threshold <- -log10(suggestive_cutoff)
  
  # Create the Manhattan plot
  manhattan_plot <- ggplot(don, aes(x = BPcum, y = logp)) +
    geom_point(aes(color = as.factor(chromosome)), alpha = 0.8, size = 1.3) +
    scale_color_manual(values = rep(c("grey", "black"), 22)) +
    scale_x_continuous(label = axisdf$chromosome, breaks = axisdf$center) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, max(don$logp) + 2)) +
    geom_hline(yintercept = bonferroni_threshold, color = "#e7298a", linetype = "dashed") +
    geom_hline(yintercept = suggestive_threshold, color = "#78c679", linetype = "dashed") +
    geom_text(data = subset(don, logp > suggestive_threshold), aes(label = CpG), vjust = -1, size = 2) +
    theme_bw() +
    theme(
      legend.position = "none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    ) +
    labs(x = "Chromosome", y = "-log10(p)")
  
  # Save the plot
  save_to <- paste0(results_folder, '/', sub(".rds", "_manhattan_bacon.png", basename(file)))
  ggsave(save_to, plot = manhattan_plot, width = 14, height = 7, dpi = 300)
  
  
  
}



