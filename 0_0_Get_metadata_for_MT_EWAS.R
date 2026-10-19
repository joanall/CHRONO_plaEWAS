#==============================================================================#
#                                                                              #
#                MATERNAL DIET MEALTIME AND PLACENTA DNA METHYLATION           #
#                              EWAS Analysis                                   #
#                                                                              #                                  
#                   Camille Lassale & Joana Llauradó Pont                      #         
#                              April 2024                                      #
#                                                                              #
#==============================================================================#
# [NOTE]: This script won't be published in the github.
# Description: 

# [INPUT]:
# - RDS dataframe with diet data, colnames food groups and rownames samples
# - RDS dataframe with cell type data, colnames -- and --- 
# - RDS dataframe with colnames as covariables and Nrows samples

# [OUTPUT]:



#==============================================================================#
#                       INIT  WORKING ENVIRONMENT                              #----
#==============================================================================#


current_date <- "test"
print(paste0('YOUR WORKING DIRECTORY IS:', getwd()))
cat('CURRENT DATE IS : ', current_date)

# rm(list=ls())
#.rs.restartR()
# gc()

# Load libraries 
library(dplyr)
library(tidyverse) #quicker with large data 
library(lubridate)
library(mice)
library(gtsummary)



#==============================================================================#
#                                                                              #
#                        GET DATA IN ONE DATAFRAME                             #----
#                                                                              #
#==============================================================================#

# Load data covariables data
covariables_data <- readRDS("db/general_covars_data/covariables_data_02012025.rds") # N = 1080 
# Add weight gain mothers 
load("db/general_covars_data/weight_gained_jllaurado.RData") # named df 

# add weighted bars to covariables data by id 
covariables_data <- merge(covariables_data, df, by.x= 'id_mother', by.y = "id_bisc")

# Load diet data
diet_data <- readRDS("db/diet_data/diet_data_02012025.rds") # N = 809 

# Load sleep data (File Ctrl-C + Ctrl-V from DATA CLEANING FOLDER)
sleep_data <- readRDS("db/sleep_data/final_sleep_data_bisc_12w_duration_31012025.rds")  # N = 1080 (with NAs)

# Load mealtime data (File Ctrl-C + Ctrl-V from DATA CLEANING FOLDER)
mealtime_data <- readRDS("db/mealtime_data/final_mealtime_data_format_duration_11012025.rds") # N = 795
mealtime_data$id_mother <- as.double(mealtime_data$id_mother) 


# Load other variables from metadata in dnam (study center, dna contmaiantione etc)
covars_meth_data <- as.data.frame(readRDS("db/methylation_data/covars_meth_BISC.RDS"))  #N = 565

# Load CT data 
ct_data <- readRDS("db/methylation_data/celltype_BISC.RDS")  # N = 565 
ct_data$Basename <- rownames(ct_data)

# Dictionary of ids  -> Basename 
ids <- read.csv('db/methylation_data/ids_dictionary.csv', sep=";")
ids <- ids[,c('id_mother', 'Basename')]
ids_dict <- tibble::deframe(ids)
ct_data$id_mother <- as.double(names(ids_dict)[ids_dict == rownames(ct_data)])

# Merging 
df_list <- list(covariables_data, diet_data , mealtime_data, sleep_data, ct_data, covars_meth_data)
metadata <- df_list %>% purrr::reduce(full_join, by='id_mother') # 1080 x 56
saveRDS(metadata, file = paste0("db/joined_metadata_incomplete_cases_",current_date,".rds"))

rm(covariables_data, df, diet_data , mealtime_data, sleep_data, ct_data, 
   covars_meth_data, ids, df_list)


#==============================================================================#
#                                                                              #
#                           SLEEP VARIABLES                                    #-----
#                                                                              #
#==============================================================================#

#------------------------------------------------------------------------------#
#                       Impute sleep variables                                 #----
#------------------------------------------------------------------------------#

# Variables to be used for imputation
impute_vars <- c('sleep_duration', 'psqi_score', 'maternal_age', 'maternal_bmi',
                 'maternal_educational_level', 'maternal_smoking', 'parity',
                 'maternal_physical_activity', 'gestational_age', 'Sex',
                 'child_ethnicity', 'dna_contamination', 'gestational_diabetes',
                 'energy_intake', 'pregnancy_complications', 'study_center',
                 'covid_confinement', 'first_meal_wday', 'first_meal_wend',
                 'last_meal_wday', 'last_meal_wend', 'eating_occasions_wday',
                 'eating_occasions_wend')

# Subset the data to include only the variables needed for imputation
impute_data <- metadata[, impute_vars]

# Create a predictor matrix where only sleep_duration and sleep_score are imputed
predictorMatrix <- make.predictorMatrix(impute_data)
predictorMatrix[, c('sleep_duration', 'psqi_score')] <- 1
predictorMatrix[c('sleep_duration', 'psqi_score'), c('sleep_duration', 'psqi_score')] <- 0

# Perform the imputation using the mice package
imputed_data <- mice(impute_data, method = 'pmm', predictorMatrix = predictorMatrix, m = 5)
saveRDS(imputed_data, "imputed_data_05022025.rds")

#------------------------------------------------------------------------------#
#                            Plots & Statistics                                #----
#------------------------------------------------------------------------------#

# Density plot for original and imputed sleep_duration
pdf(paste0("db/sleep_data/density_plot_sleep_duration",current_date,".pdf"))
densityplot(imputed_data, ~ sleep_duration)
dev.off()

# Density plot for original and imputed sleep_score
pdf(paste0("db/sleep_data/density_plot_sleep_score",current_date,".pdf"))
densityplot(imputed_data, ~ psqi_score)
dev.off()

# Trace plot for sleep_duration
pdf(paste0("db/sleep_data/trace_plot_sleep_duration",current_date,".pdf"))
stripplot(imputed_data, pch = 20, cex = 1.2, sleep_duration ~ .imp)
dev.off()

# Trace plot for sleep_score
pdf(paste0("db/sleep_data/trace_plot_sleep_score",current_date,".pdf"))
stripplot(imputed_data, pch = 20, cex = 1.2, psqi_score ~ .imp)
dev.off()

# Summary statistics for original and imputed data
summary_original <- summary(impute_data)
completed_data <- complete(imputed_data)
summary_imputed <- summary(completed_data)

write.csv(summary_original, 
          paste0("db/sleep_data/summary_original_sleep_impu",current_date,".csv"))
write.csv(summary_imputed, 
          paste0("db/sleep_data/summary_imputed_sleep_impu",current_date,".csv"))





# Replace the original sleep_duration and sleep_score in the metadata with the imputed values
metadata$sleep_duration <- completed_data$sleep_duration
metadata$psqi_score <- completed_data$psqi_score

save_to <-  paste0("db/sleep_data/metadata_plaEWAS_incomplete_cases_sleepimputed",current_date,".rds")
saveRDS(metadata, save_to)

rm(impute_data, imputed_data, predictorMatrix)

# Note it chanegs everutime because it's a mulitple imputation - dim 1080
mean(as.numeric(metadata$sleep_duration)/3600) # 8.954473
mean(as.numeric(metadata$sleep_duration)/3600) # 8.942537
# 8.930042
# 8.965276


#==============================================================================#
#                                                                              #
#                              FINAL DATASET                                   #----
#                                                                              #
#==============================================================================#

# Note: We do complete. cases with the averages as is what we use in the rlm models
# and later on we add the individual mealtime and sleep variables so we can do the
# distribution tables. 

#------------------------------------------------------------------------------#
#                          Complete Cases Data                                 #----
#------------------------------------------------------------------------------#

# Define variables to use in EWAS MT
general_covars <- c( "Basename","maternal_age", "maternal_bmi", 
                     "maternal_educational_level", "maternal_smoking", "parity",
                     "maternal_physical_activity", "gestational_age", "Sex",
                     "dna_contamination", "study_center",
                     'gestational_diabetes',"pregnancy_complications", "covid_confinement", 
                     "child_ethnicity") 


mealtime_vars <- c('first_meal_wday','first_meal_wend',
                   'last_meal_wday','last_meal_wend',
                   'eating_occasions_wday', 'eating_occasions_wend',
                   'first_meal', 'last_meal','eating_occasions')

sleep_vars <- c("sleep_duration", "psqi_score")

diet_vars <- c('cereals', 'fruit_nuts', 'meat', 'legumes','fish',
               'legumes','olive_oil', 'vegetables', 'dairy',"energy_intake")

cell_type_vars <- c('Trophoblasts','Stromal', 'Hofbauer', 'Endothelial','nRBC', 'Syncytiotrophoblast')

all_covars <-  c(general_covars,  sleep_vars, mealtime_vars, diet_vars, cell_type_vars)

final_mealtime <-  metadata[ ,all_covars] # 1080 
cc <- final_mealtime

final_mealtime <- final_mealtime[complete.cases(final_mealtime),] #389 or 366
cat('FINAL SAMPLE SIZE: ', nrow(final_mealtime))

# Note: From 1080 participants we have: 
# - 809 with diet data 
# - 795 with mealtime data 
# - 565 with plaDNAm

# - Complete cases considering intermediate and averages mealtime vars -> N = 389
# - Complete cases with considering only averages mealtime vars -> N = 395
# - Complete cases considering sleep vars -> N = 272


# Checks 
# Original Sample = 1080
# Samples with DNAm + diet = 
# Samples with DNAm + mealtime = 
# Samples with diet + mealtime = 
# Samples with DNAm + diet + mealtime = 565

# Step 1: Create logical vectors for non-missing values
has_dnam <- !is.na(metadata$Trophoblasts)
has_diet <- !is.na(metadata$cereals)
has_mealtime <- !is.na(metadata$eating_occasions)

# Step 2: Count original sample size
original_sample_count <- nrow(metadata)

# Step 3: Compute intersections
dnam_diet <- sum(has_dnam & has_diet)
dnam_mealtime <- sum(has_dnam & has_mealtime)
diet_mealtime <- sum(has_diet & has_mealtime)
dnam_diet_mealtime <- sum(has_dnam & has_diet & has_mealtime)

# Step 4: Print results
cat("Original Sample =", original_sample_count, "\n")
cat("Samples with DNAm + diet =", dnam_diet, "\n")
cat("Samples with DNAm + mealtime =", dnam_mealtime, "\n")
cat("Samples with diet + mealtime =", diet_mealtime, "\n")
cat("Samples with DNAm + diet + mealtime =", dnam_diet_mealtime, "\n")


mean(as.numeric(final_mealtime$sleep_duration)/3600) # 8.804462

# From 
# complete_metadata_with_dna <- metadata[complete.cases(metadata[, cell_type_vars]), ]

# Add wake_up and go_bed -> They have missings but not the sleep_duration, which is the
# one we use for the the models. 


sleep_vars <- c("wakeup_wday", "wakeup_wend",  
                "gobed_wday", "gobed_wend",  
                "sleep_duration_wday", "sleep_duration_wend", 
                "gest_weight_gain_m_preg")

final_mealtime <- merge(final_mealtime, metadata[,c('Basename', sleep_vars)], by = "Basename")
print(colnames(final_mealtime))

save_to <- paste0("db/metadata_plaEWAS_MT_incomplete_cases2_wg",current_date,".rds")
saveRDS(final_mealtime, save_to) 


#==============================================================================#
#                                                                              #
#                           CHECK SUBSAMPLE - SAMPLE                           #----
#                                                                              #
#==============================================================================#


#all_metadata <- readRDS("db/final_data/metadata_plaEWAS_incomplete_cases_sleepimputed.rds") #1080 x 55


# Select the required variables
sel_vars <- c("Basename", "maternal_age", "maternal_educational_level", 
              "maternal_physical_activity", "maternal_smoking", 
              "maternal_ethnicity")

all_metadata <- metadata[,sel_vars]


# Get ids from the subsample 
#sub_metadata <- readRDS("db/final_data/metadata_plaEWAS_MT_complete_cases_02012025.rds")
sub_metadata <- final_mealtime
ids_analytical <- sub_metadata$Basename


# Create new column where 0 = not included and 1 = is included  (ids_analytical)
all_metadata$sample <- ifelse(all_metadata$Basename %in% ids_analytical, 1, 0)

# Remove the 'Basename' column
all_metadata <- all_metadata[, !(names(all_metadata) %in% 'Basename')]



# Perform the Chi-square test if there is a singificant difference between the sub_metadata
# and the exlcuded data. do the Chi square tes twith onlu the sel_vars
# Put all result sin a table 

# Perform Chi-square tests and add p-values to the summary table
p_values <- sapply(sel_vars[-1], function(var) {
  if (is.factor(all_metadata[[var]])) {
    test <- chisq.test(table(all_metadata[[var]], all_metadata$sample))
    return(test$p.value)
  } else { return(NA)}
})

# Create summary table with p-values
summary_table <- all_metadata %>%
  tbl_summary(
    by = sample, # Variable to split the data (columns)
    missing_text = "DK/M", # How you want missings to be named
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})", # Statistics for continuous variables to show [in this case: median (quartile 1, quartile 2)]
      all_categorical() ~ "{n} ({p}%)" # Statistics for categorical variables to show [in this case: N (percentage %)]
    )
  ) %>%
  add_p(
    test = list(
      maternal_educational_level ~ "chisq.test",
      maternal_physical_activity ~ "chisq.test",
      maternal_smoking ~ "chisq.test",
      maternal_ethnicity ~ "chisq.test"
    )
  )

# Print summary table with p-values
print(summary_table)

# Convert summary table to a tibble (data frame)
summary_table_df <- as_tibble(summary_table)

# Save the summary table to a CSV file
save_to <- paste0("db/chi_square_results_",current_date,".csv")
write.csv(summary_table_df, save_to, row.names = FALSE)







