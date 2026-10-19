#==============================================================================#
#                                                                              #
#               MATERNAL DIET MEALTIME AND PLACENTA DNA METHYLATION            #
#                             Quality Control                                  #
#                                                                              #                                  
#                   Camille Lassale & Joana Llauradó Pont                      #         
#                               Janaury 2025                                   #
#                                                                              #
#==============================================================================#

# DESCRIPTION: Read the the dafraframes and get plots of quality control and data analisis
#   - Compute lambdas, and adjusted  
#   - QQ plot with Lambda and CI 
#     - If lambda is inflated, apply bacon to pvalues  
#   - Distribution SE and p vals


# [INPUT]: RDS File (list of files) storing a dataframe where rows are CpGs and columns: 

# [OUTPUT]: QQ plot, SE and PVALUE distribution plots and newd df.RDS with adj values


#==============================================================================#
#                                                                              #
#                                  USER INPUT                                  #----
#                                                                              #
#==============================================================================#

# Options: nighttime_fasting , first_meal, eating_occasions, eating_jetlag, last_meal
# exposure_name <- "nighttime_fasting"


# Input & output folder - path 
input_ewas_file  <- "results/1_raw_results/Sensitivity_13062025/eating_jet"
results_folder <- "results/2_quality_control/Sensitivity_13062025/eating_occasions_model_noGDM_adjusted"
if (!dir.exists(results_folder)) {dir.create(results_folder, recursive = T)}



#==============================================================================#
#                                                                              #
#                       INIT  WORKING ENVIRONMENT                              #----
#                                                                              #
#==============================================================================#

# rm(list=ls())
# gc()

print(paste0('YOUR WORKING DIRECTORY IS:', getwd()))


# Load libraries 
library(readr)
library(qqman)
library(bacon)
#library(plotflow) 
library(ggplot2)
library(stats)

source("script/Define_functions.R")

# Create folder (if it doesn't exist) to save outputs (main ones)
if (!dir.exists(results_folder)) {dir.create(results_folder, recursive = TRUE)}




#==============================================================================#
#                                                                              #
#                             QUALITY CONTROL                                  #----
#                                                                              #
#==============================================================================#

cat("\n----------------------  QUALITY CONTROL   ----------------------------\n")

# Load the EWAS result and get df 
cat('Reading file', basename(input_ewas_file),'\n')
results <- readRDS(input_ewas_file)
input_ewas_file


#------------------------------------------------------------------------------#
#                        LAMBDA AND  QQ PLOTS                                  #----
#------------------------------------------------------------------------------#

cat("\n------- LAMBDA and  QQ  Plots  -----------------\n")
# Get lambda and CI 
lambda_ci <- as.data.frame(t(inflation(results$pvalue)))
colnames(lambda_ci) <- c('lambda','lower_ci','upper_ci')
save_to <- file.path(results_folder, sub(".rds", "_lambda.csv", 
                                         basename(input_ewas_file)))
save_to
write.table(lambda_ci, save_to, sep =";", row.names = FALSE)

# Get QQ PLOT with lambda
qq_plot <- qq_plot2(results, 'pvalue')
qq_lambda <- add_lambda_to_plot2(lambda_ci, qq_plot)


# Save
save_to <- file.path(results_folder, "QQplot_lambda.png")
ggsave(save_to, plot = qq_lambda, width = 6, height = 6, dpi = 600)

save_to <- file.path(results_folder, "QQplot_lambda.rds")
saveRDS(qq_lambda, file = save_to)


if (lambda_ci[1] > 1.1) {
  print('Lambda above 1.1, applying bacon to correct pvalues')
  
  # Apply bacon adjustment to EWAS results
  bc <- bacon(effectsizes = results$Effect,
              standarderrors = results$SE) 

  # Adjust pblaues based on bacon
  adjusted_pvalues <- bacon::pval(bc)
  
  # Add to results  adjusted p values 
  results$pval_bacon <- adjusted_pvalues
  
  # Calculate again lambda 
  corrected_lambda_ci <- as.data.frame(t(inflation(results$pval_bacon)))
  colnames(corrected_lambda_ci) <- c('lambda','lower_ci', 'upper_ci' )
  print(corrected_lambda_ci)
  save_to <- file.path(results_folder, sub(".rds", "_corrected_lambda.csv", basename(input_ewas_file)))
  write.table(corrected_lambda_ci, save_to, sep=";", row.names = FALSE)
  
  #Get QQ PLOT with lambda
  qq_plot <- qq_plot2(results, 'pval_bacon')
  qq_lambda <- add_lambda_to_plot2(corrected_lambda_ci, qq_plot)
  # Save
  save_to <- file.path(results_folder, "QQplot_corrected_lambda.png")
  ggsave(save_to, plot = qq_lambda, width = 6, height = 6, dpi = 600)
  
  save_to <- file.path(results_folder, "QQplot_corrected_lambda.rds")
  saveRDS(qq_lambda, file = save_to)
  
}







#------------------------------------------------------------------------------#
#                        MULTIPLE TESTING CORRECTION                           #----
#------------------------------------------------------------------------------#

cat("\n--------- Multiple Testing Correction  ------------------\n")


# P value adjustment
results$fdr <- p.adjust(results$pvalue, method = 'fdr')
results$bf <- p.adjust(results$pvalue, method = 'bonferroni')

# Order descending based on pval 
results <- results[order(-results$pvalue), ]

if ('pval_bacon' %in% colnames(results)) {
  results$fdr_bacon <- p.adjust(results$pval_bacon, method= 'fdr')
  results$bf_bacon <- p.adjust(results$pval_bacon, method= 'bonferroni')
  
  results <- results[order(-results$pval_bacon), ]
}

# Save df 
saveRDS(results, file.path(results_folder, sub(".rds", "_adjusted.rds", 
                                               basename(input_ewas_file))))
write.table(results,  
            file.path(results_folder, sub(".rds", "_adjusted.csv", 
                                          basename(input_ewas_file))), 
            sep=";", row.names = FALSE)






#------------------------------------------------------------------------------#
#                         DISTRIBUTION PLOTS                                   #----
#------------------------------------------------------------------------------#
# 
# cat("\n--------------------- Distribution Plots   ----------------------------\n")
# 
# 
# # Get P values distribution plot 
# save_to <- file.path(results_folder, "pvals_distplot.png")
# 
# png(save_to, type="cairo")
# 
# plot_distribution(results$pvalue, 
#                   main = paste('p-values of', exposure_name), 
#                   xlab = 'p-value')
# 
# dev.off()
# 
# if ('pval_bacon' %in% colnames(results)) {
#   save_to <- file.path(results_folder, "corrected_pvals_distribution_plot.png")
#   png(save_to, type="cairo")
#   plot_distribution(results$pval_bacon, main = paste('corrected p-values of', exposure_name), xlab = 'p-value')
#   dev.off()
# 
# }
# 
# # Get SE distribution plot 
# 
# save_to <- file.path(results_folder, "SE_distplot.png")                                         
# png(save_to, type="cairo")
# plot_distribution(results$pvalue, main = paste('Standard Errors of', exposure_name), xlab = 'SE')
# dev.off()



#------------------------------------------------------------------------------#
#                           TO CHECK IF OUTLIERS                               #----
#------------------------------------------------------------------------------#


# qq_plot <- qq_norm_plot(results, 'pvalue')
# qq_lambda <- add_lambda_to_plot(results$pvalue, qq_plot)
# lambda_ci
# # lambda  lower_ci  upper_ci
# # 1 0.9251366 0.9016693 0.9486038
# 
# # check outlier 
# outlier_cpg <- results %>%
#   filter(pvalue == min(pvalue, na.rm = TRUE)) %>%
#   select(CpG, pvalue)
# 
# print(outlier_cpg) 
# # 
# # CpG      pvalue
# # 1 cg02972356 9.41852e-07
# 
# 
# 
# # Remove the outlier CpG from the dataset
# results_without_outlier <- results %>%
#   filter(CpG != outlier_cpg$CpG)
# 
# 
# 
# # Get lambda and CI 
# lambda_ci <- as.data.frame(t(inflation(results_without_outlier$pvalue)))
# colnames(lambda_ci) <- c('lambda','lower_ci','upper_ci')
# print(lambda_ci)
# 
# qq_plot <- qq_norm_plot(results_without_outlier, 'pvalue')
# qq_lambda <- add_lambda_to_plot(results_without_outlier$pvalue, qq_plot)
# # lambda  lower_ci  upper_ci
# # 1 0.9251354 0.9016741 0.9485967
# 
# 
# # P value adjustment
# results_without_outlier$fdr <- p.adjust(results_without_outlier$pvalue, method = 'fdr')
# results_without_outlier$bf <- p.adjust(results_without_outlier$pvalue, method = 'bonferroni')
# 
# # Order descending based on pval 
# results_without_outlier <- results_without_outlier[order(-results_without_outlier$pvalue), ]



