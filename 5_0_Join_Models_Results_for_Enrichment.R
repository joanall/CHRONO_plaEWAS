#==============================================================================#
#                                                                              #
#               MATERNAL DIET MEALTIME AND PLACENTA DNA METHYLATION            #
#                             Functional Enrichment                            #
#                                                                              #                                  
#                   Camille Lassale & Joana Llauradó Pont                      #         
#                               January 2025                                   #
#                                                                              #
#==============================================================================#

# Description: Join all mealtime related CpG-genes associated with placental DNAm 
# to later on perform functiona enrichment analysis. 

# [INPUT]:

# [OUTPUT]:


#==============================================================================#
#                       INIT  WORKING ENVIRONMENT                              #
#==============================================================================#
results_path <- "results/5_functional_analysis_enrichment/Bisc_13062025"
if (!dir.exists(results_path)) {dir.create(results_path, recursive = TRUE)}

#==============================================================================#
#                          DEFINE PATHS                                        #
#==============================================================================#

# ewas_files <- c('results/2_quality_control/Bisc_05022025/nighttime_fasting_model_3_sep_exp/nighttime_fasting_Bisc_model_3_sep_exp_05022025_adjusted.rds',
#                 'results/2_quality_control/Bisc_05022025/nighttime_fasting_model_3/nighttime_fasting_Bisc_model_3_05022025_adjusted.rds',  # adjsuted first meal
#                 'results/2_quality_control/Bisc_05022025/first_meal_model_3_sep_exp/first_meal_Bisc_model_3_sep_exp_05022025_adjusted.rds', 
#                 'results/2_quality_control/Bisc_05022025/first_meal_model_3_firstlast/first_meal_Bisc_model_3_firstlast_05022025_adjusted.rds',
#                 'results/2_quality_control/Bisc_05022025/last_meal_model_3_sep_exp/last_meal_Bisc_model_3_sep_exp_05022025_adjusted.rds',
#                 'results/2_quality_control/Bisc_05022025/last_meal_model_3_firstlast/last_meal_Bisc_model_3_firstlast_05022025_adjusted.rds',
#                 'results/2_quality_control/Bisc_05022025/eating_occasions_model_3_sep_exp/eating_occasions_Bisc_model_3_sep_exp_05022025_adjusted.rds',
#                 'results/2_quality_control/Bisc_05022025/eating_occasions_model_3_firstlast/eating_occasions_Bisc_model_3_firstlast_05022025_adjusted.rds',
#                 'results/2_quality_control/Bisc_05022025/eating_jetlag_model_3/eating_jetlag_Bisc_model_3_05022025_adjusted.rds')

# Models 1 & 2
# ewas_files <- c('results/2_quality_control/Bisc_14012025/nighttime_fasting_model_3_sep_exp/nighttime_fasting_Bisc_model_3_sep_exp_14012025_adjusted.rds',
#                 'results/2_quality_control/Bisc_14012025/nighttime_fasting_model_3/nighttime_fasting_Bisc_model_3_14012025_adjusted.rds',  # adjsuted first meal
#                 'results/2_quality_control/Bisc_14012025/first_meal_model_3_sep_exp/first_meal_Bisc_model_3_sep_exp_14012025_adjusted.rds', 
#                 'results/2_quality_control/Bisc_14012025/first_meal_model_3_firstlast/first_meal_Bisc_model_3_firstlast_14012025_adjusted.rds',
#                 'results/2_quality_control/Bisc_14012025/last_meal_model_3_sep_exp/last_meal_Bisc_model_3_sep_exp_14012025_adjusted.rds',
#                 'results/2_quality_control/Bisc_14012025/last_meal_model_3_firstlast/last_meal_Bisc_model_3_firstlast_14012025_adjusted.rds',
#                 'results/2_quality_control/Bisc_14012025/eating_occasions_model_3_sep_exp/eating_occasions_Bisc_model_3_sep_exp_14012025_adjusted.rds',
#                 'results/2_quality_control/Bisc_14012025/eating_occasions_model_3_firstlast/eating_occasions_Bisc_model_3_firstlast_14012025_adjusted.rds',
#                 'results/2_quality_control/Bisc_14012025/eating_jetlag_model_3/eating_jetlag_Bisc_model_3_14012025_adjusted.rds')

ewas_files <- c('results/2_quality_control/Main_13062025/last_meal_model_1_adjusted/last_meal_model_1_adjusted_adjusted.rds', 
                'results/2_quality_control/Main_13062025/last_meal_model_1_sep_exp/last_meal_model_1_sep_exp_adjusted.rds') 
                
                # 'results/2_quality_control/Main_13062025/first_meal_model_1_adjusted/first_meal_model_1_adjusted_adjusted.rds', 
                # 'results/2_quality_control/Main_13062025/first_meal_model_1_sep_exp/first_meal_model_1_sep_exp_adjusted.rds', 
                # 
                # 'results/2_quality_control/Main_13062025/nighttime_fasting_model_1_adjusted/nighttime_fasting_model_1_adjusted_adjusted.rds', 
                # 'results/2_quality_control/Main_13062025/nighttime_fasting_model_1_sep_exp/nighttime_fasting_model_1_sep_exp_adjusted.rds', 
                # 
                # 'results/2_quality_control/Main_13062025/eating_occasions_model_1_adjusted/eating_occasions_model_1_adjusted_adjusted.rds', 
                # 'results/2_quality_control/Main_13062025/eating_occasions_model_1_sep_exp/eating_occasions_model_1_sep_exp_adjusted.rds', 
                # 
                # 'results/2_quality_control/Main_13062025/eating_jetlag_model_1_sep_exp/eating_jetlag_model_1_sep_exp_adjusted.rds')
                
                
pval_cutoff <- 1E-04
array <- 'EPIC'


#==============================================================================#
#                       JOIN MEALTIME DATAFRAMES                               #
#==============================================================================#

# Load and combine the dataframes
dfs <- lapply(ewas_files, readRDS)


# Option with csv (dmrs)
#dfs <- lapply(ewas_files, function(file) read.table(file, sep = ";",  header = T))
if (array == "450K") {
  library('IlluminaHumanMethylation450kanno.ilmn12.hg19')
  annotation <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
} else if (array == "EPIC") {
  library('IlluminaHumanMethylationEPICanno.ilm10b2.hg19')
  annotation <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
}
annotation <- as.data.frame(annotation)

# Take only "CpG, pval" columns and subset by p-value cutoff using base R
filtered_dfs <- lapply(dfs, function(df) {
  
  # Decide which p-value column to use
  pval_col <- if ("pval_bacon" %in% colnames(df)) "pval_bacon" else "pvalue"
  
  # Select and filter
  to_select <- c("CpG", pval_col)
  filtered_df <- df[df[[pval_col]] < pval_cutoff, to_select, drop = FALSE]
  
  # Merge gene annotation using CpG as key
  if (nrow(filtered_df) > 0) {
    filtered_df$UCSC_RefGene_Name <- annotation[filtered_df$CpG, "UCSC_RefGene_Name"]
  } else {
    filtered_df$UCSC_RefGene_Name <- character(0)  # Avoid error on empty data frame
  }
  
  return(filtered_df)
})


# Combine and keep unique CpGs (here we don't care about pvals)
combined_df <- do.call(rbind, filtered_dfs)
combined_df <- combined_df[!duplicated(combined_df$CpG), ]

# Add fake columns "coef", "bf", "pval" (pval cutoff)
set.seed(123)  # For reproducibility of random numbers
combined_df$Effect <- runif(nrow(combined_df))
combined_df$SE <- runif(nrow(combined_df))
combined_df$bf <- runif(nrow(combined_df))
combined_df$pvalue <- pval_cutoff


# Save the combined dataframe to an RDS file
#output_file <- paste0(results_folder, "/", cohort_name, date, "_cpgs_mealtime_model1_adjusted_", pval_cutoff, ".rds")
#save_to <- file.path(results_path, "mealtime_cpgs_toenrich_model_3_1E06.rds")
save_to <- file.path(results_path, "last_meal_cpgs_toenrich_model_1_1E04.rds")
saveRDS(combined_df, save_to)

# Save from one of the files allcpgs to set the backround
all_cpgs <- data.frame(rownames(dfs[[1]]))
colnames(all_cpgs) <- "CpG" 
#save_to <- file.path(results_path, "mealtime_cpgs_all_model_3_1E06.rds")
save_to <- file.path(results_path, "last_meal_cpgs_all_model_1_1E05.rds")
saveRDS(all_cpgs, save_to)
