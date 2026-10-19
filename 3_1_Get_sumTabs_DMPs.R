#==============================================================================#
#                                                                              #
#               MATERNAL DIET MEALTIME AND PLACENTA DNA METHYLATION            #
#                             Analysis DMPs                                    #
#                                                                              #                                  
#                   Camille Lassale & Joana Llauradó Pont                      #         
#                               January 2025                                   #
#                                                                              #
#==============================================================================#


# INPUT: RDS or CSV with Processed EWAS results. If not dataframe with ewas results 
# where rownames are CpGs ans the columnames are (and should be for the code tow ork):
# Effect, SE, tvalue, pvalue, CpG, fdr, bf

# OUTPUT: Summary tables grouped by exposure and/or model 


# Note: Not prepared to run all together, go to each chunk, speicfi the files
# and run the summ tables 

#==============================================================================#
#                                                                              #
#                       INIT  WORKING ENVIRONMENT                              #----
#                                                                              #
#==============================================================================#

# rm(list=ls())
#.rs.restartR()
# gc()

print(paste0('YOUR WORKING DIRECTORY IS:', getwd()))

# Root folder -> Z:/analyses/PACE/mat_veg_diet/BISC"
# Root folder -> /PROJECTES/INMA_OMICS/analyses/PACE/mat_veg_diet/BISC

# Load libraries 
library(R.utils)
library(httr)
library(rvest)
library(data.table)
source("script/Define_functions.R")

#------------------------------------------------------------------------------#
#                           Input from User                                    #
#------------------------------------------------------------------------------#

suggestive_cutoff <- 1E-05



# Output folder 
results_newfolder  <- "results/3_analysis_dmps/SumTables_13062025/"

cat("OUTPUT FOLDER:\n",results_newfolder, "\n")

# Create folder (if it doesn't exist) to save outputs (main ones)
if (!dir.exists(results_newfolder)) {dir.create(results_newfolder)}



#==============================================================================#
#                                                                              #
#                       Summary Table Model Main Analysis                      #----
#                                                                              #
#==============================================================================#  
# [NOTE]: Change filenames each time
# [CHOOSE EXPOSURE]: nighttime_fasting , eating_occasions first_meal, eating_jetlag
exposure_name <- "nighttime_fasting"

#file_name <- paste0(exposure_name, cohort_name, )  
input_files_main <- c(
'results/2_quality_control/Main_13062025/nighttime_fasting_model_1_sep_exp/nighttime_fasting_model_1_sep_exp_adjusted.rds', 
'results/2_quality_control/Main_13062025/nighttime_fasting_model_2_sep_exp/nighttime_fasting_model_2_sep_exp_adjusted.rds', 
'results/2_quality_control/Main_13062025/nighttime_fasting_model_3_sep_exp/nighttime_fasting_model_3_sep_exp_adjusted.rds') 



# Read all files into a list
list_ewas_results_main <- lapply(input_files_main, function(file) readRDS(file))

# Extract top CpGs for each file base on suggestive cutoff
top_cpgs_main <- lapply(list_ewas_results_main, 
                   function(ewas) {subset(ewas, pvalue <= suggestive_cutoff)$CpG})

# Get a unique list of all top CpGs
unique_top_cpgs_main <- unique(unlist(top_cpgs_main))

cat('Number of suggestive CpGs:', length(unique_top_cpgs_main))


# Create a summary table for each EWAS result with the specified columns
summary_tables <- lapply(seq_along(list_files), function(i) {
  ewas <- list_files[[i]]
  summtab <- ewas[ewas$CpG %in% list_cpgs, c('CpG', 'Effect', 'SE', 'pvalue', 'bf')]
  colnames(summtab) <- c('CpG', paste0('coef_model', i), 
                         paste0('se_model', i), 
                         paste0('pval_model', i), 
                         paste0('bf_model', i))
  return(summtab)
})

# Merge all summary tables by CpG
combined_summary_table <- Reduce(function(x, y) merge(x, y, by = "CpG", all = TRUE), summary_tables)

# Order by CpG
combined_summary_table <- combined_summary_table[order(combined_summary_table$CpG), ]

write.table(combined_summary_table, file = new_filename, sep = ";", row.names = FALSE)
# Save to CSV file 
save_to <- paste0(results_newfolder, "nigghtime_fasting_sep_exp.csv")
  
  
get_summary_table(unique_top_cpgs_main,input_files_main,save_to)



#==============================================================================#
#                                                                              #
#                       Summary Table M1 to Sensitivity Models                 #----
#                                                                              #
#==============================================================================# 

# FROM M1 TO SENS 
# [CHOOSE EXPOSURE]: nighttime_fasting , eating_occasions first_meal, eating_jetlag
exposure_name <- "eating_occasions"

input_files_sens <- c(paste0("results/2_quality_control/Bisc_05022025/", 
                             exposure_name, '_model_ethnicity_sep_exp/', 
                             exposure_name, "_Bisc_model_ethnicity_sep_exp_05022025_adjusted.rds"),
                      
                      paste0("results/2_quality_control/Bisc_05022025/", 
                             exposure_name, '_model_noGDM_sep_exp/', 
                             exposure_name, "_Bisc_model_noGDM_sep_exp_05022025_adjusted.rds"),
                      
                      paste0("results/2_quality_control/Bisc_05022025/", 
                             exposure_name, '_model_noconf_sep_exp/', 
                             exposure_name, "_Bisc_model_noconf_sep_exp_05022025_adjusted.rds"),
                      
                      paste0("results/2_quality_control/Bisc_05022025/",
                             exposure_name, '_early_eaters_sep_exp/',
                             exposure_name, "_Bisc_early_eaters_sep_exp_05022025_adjusted.rds"), 
                      
                      paste0("results/2_quality_control/Bisc_05022025/",
                             exposure_name, '_short_fasting_sep_exp/',
                             exposure_name, "_Bisc_short_fasting_sep_exp_05022025_adjusted.rds"))

# 1) GET Top CpGs M1 
# statistic top unique cpg min + sens from model 1, (we use it as reference)
ewas <- readRDS(paste0("results/2_quality_control/Bisc_14012025/", 
                        exposure_name, '_model_1_sep_exp/', 
                        exposure_name, "_Bisc_model_1_sep_exp_14012025_adjusted.rds"))


top_cpgs_m1 <- subset(ewas, pvalue <= suggestive_cutoff)$CpG

summtab_m1 <- ewas[ewas$CpG %in% top_cpgs_m1, c('CpG', 'Effect', 'SE', 'pvalue', 'bf')]
colnames(summtab_m1) <- c('CpG', paste0('coef_m1'), 
                          paste0('se_m1'), 
                          paste0('pval_m1'), 
                          paste0('bf_m1'))


# 2) Get STAITSITS for TOP CPGs in MOdel 1 from sensitivities
list_ewas_results_sens <- lapply(input_files_sens, function(file) readRDS(file))

if (length(top_cpgs_m1) != 0) {
  
  summary_tables_sens <- lapply(seq_along(list_ewas_results_sens), function(i) {
    ewas <- list_ewas_results_sens[[i]]
    summtab <- ewas[ewas$CpG %in% top_cpgs_m1, c('CpG', 'Effect', 'SE', 'pvalue', 'bf')]
    colnames(summtab) <- c('CpG', paste0('effect_sens', i), 
                           paste0('se_sens', i), 
                           paste0('pval_sens', i), 
                           paste0('bf_sens', i))
    return(summtab)
  })
  
  # Merge all summary tables by CpG
  combined_summary_table_sens <- Reduce(function(x, y) merge(x, y, by = "CpG", all = TRUE), summary_tables_sens)
  
  # attach m1 
  final_summary_table_sens <- merge(summtab_m1, combined_summary_table_sens, by='CpG')
  
  # Order by CpG
  final_summary_table_sens <- final_summary_table_sens[order(final_summary_table_sens$CpG), ]
  
  # Save to CSV file 
  save_to <- paste0(results_newfolder, "sumtab_", exposure_name, "_sep_exp_cpg_m1_sensmodels.csv")
  
  write.table(final_summary_table_sens, file = save_to, sep = ";", row.names = FALSE)
  
}





#==============================================================================#
#                                                                              #
#                       Summary Table Sensitivity Models to M1                 #----
#                                                                              #
#==============================================================================# 

# # FROM SENS TO M1 
# # [CHOOSE EXPOSURE]: nighttime_fasting , eating_occasions first_meal, eating_jetlag
# exposure_name <- "first_meal"
# 
# input_files_sens <- c(paste0("results/2_quality_control/Bisc_14012025/", 
#                                 exposure_name, '_model_ethnicity_firstlast/', 
#                                 exposure_name, "_Bisc_model_ethnicity_firstlast_14012025_adjusted.rds"),
#                       
#                       paste0("results/2_quality_control/Bisc_14012025/", 
#                                 exposure_name, '_model_noGDM_firstlast/', 
#                                 exposure_name, "_Bisc_model_noGDM_firstlast_14012025_adjusted.rds"),
#                     
#                       paste0("results/2_quality_control/Bisc_14012025/", 
#                                 exposure_name, '_model_noconf_firstlast/', 
#                                 exposure_name, "_Bisc_model_noconf_firstlast_14012025_adjusted.rds"),
#                       
#                       paste0("results/2_quality_control/Bisc_14012025/",
#                                 exposure_name, '_early_eaters_adj_first_meal/',
#                                 exposure_name, "_Bisc_early_eaters_adj_first_meal_14012025_adjusted.rds"), 
#                       
#                       paste0("results/2_quality_control/Bisc_05022025/",
#                                 exposure_name, '_short_fasting_adj_first_meal/',
#                                 exposure_name, "_Bisc_early_short_fasting_adj_first_meal_05022025_adjusted.rds"))
# 
#                
#                       
# # Read all files into a list
# list_ewas_results_sens <- lapply(input_files_sens, function(file) readRDS(file))
# 
# # Extract top CpGs for each file base on suggestive cutoff
# top_cpgs_sens <- lapply(list_ewas_results_sens, 
#                         function(ewas) {subset(ewas, pvalue <= suggestive_cutoff)$CpG})
# 
# 
# # Get a unique list of all top CpGs
# unique_top_cpgs_sens <- unique(unlist(top_cpgs_sens))
# 
# cat('Number of sens suggestive CpGs:', length(unique_top_cpgs_sens), exposure_name)
# 
# #unique_top_cpgs_main_sens <- c(unique_top_cpgs_main, unique_top_cpgs_sens)
# 
# if (length(unique_top_cpgs_sens) != 0) {
#   # Sum table fortop Cpgs in the main model but with statistics from sensitivitesd models
#   # + top cpGs in main model 
#   summary_tables_sens <- lapply(seq_along(list_ewas_results_sens), function(i) {
#     ewas <- list_ewas_results_sens[[i]]
#     summtab <- ewas[ewas$CpG %in% unique_top_cpgs_sens, c('CpG', 'Effect', 'SE', 'pvalue', 'bf')]
#     colnames(summtab) <- c('CpG', paste0('effect_sens', i), 
#                                   paste0('se_sens', i), 
#                                   paste0('pval_sens', i), 
#                                   paste0('bf_sens', i))
#     return(summtab)
#   })
#   
#   # Merge all summary tables by CpG
#   combined_summary_table_sens <- Reduce(function(x, y) merge(x, y, by = "CpG", all = TRUE), summary_tables_sens)
#   
#   # Order by CpG
#   combined_summary_table_sens <- combined_summary_table_sens[order(combined_summary_table_sens$CpG), ]
#   
#   
#   # statistic top unique cpg mian + sens from model 1, (we use it as reference)
#   ewas <- readRDS(file.path(input_folder, 
#                             paste0(exposure_name, '_model_1_firstlast'), 
#                             paste0(exposure_name, "_Bisc_model_1_firstlast_14012025_adjusted.rds")))
# 
#   
#   summtab_m1_sens <- ewas[ewas$CpG %in% unique_top_cpgs_sens, c('CpG', 'Effect', 'SE', 'pvalue', 'bf')]
#   colnames(summtab_m1_sens) <- c('CpG', paste0('coef_m1'), 
#                            paste0('se_m1'), 
#                            paste0('pval_m1'), 
#                            paste0('bf_m1'))
#   
#   # Order by CpG
#   summtab_m1_sens <- summtab_m1_sens[order(summtab_m1_sens$CpG), ]
#   
#   # attach m1 
#   final_summary_table_sens <- merge(summtab_m1_sens, combined_summary_table_sens, by='CpG')
#   
#   # Order by CpG
#   final_summary_table_sens <- final_summary_table_sens[order(final_summary_table_sens$CpG), ]
#   
#   
#   # Save to CSV file 
#   save_to <- file.path(results_newfolder, paste0( "sumtab_", 
#                     exposure_name, "adj_firstlast_meal_cpg_sensmodels_m1.csv"))
#   write.table(final_summary_table_sens, file = save_to, sep = ";", row.names = FALSE)
# }
# 








