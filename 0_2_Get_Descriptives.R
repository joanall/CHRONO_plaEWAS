#==============================================================================#
#                                                                              #
#                MATERNAL DIET MEALTIME AND PLACENTA DNA METHYLATION           #
#                              EWAS Analysis                                   #
#                                                                              #                                  
#                   Camille Lassale & Joana Llauradó Pont                      #         
#                              April 2024                                      #
#                                                                              #
#==============================================================================#

# Description: Get the descriptive tables of the cohort used in this analysis.

# [INPUT]: RDS dataframe with Nrow samples and colnames covariables included in the EWAS
# [OUTPUT]:



#==============================================================================#
#                             DEFINE PARAMETERS                                #
#==============================================================================#

metadata_path <- "db/metadata_plaEWAS_MT_complete_cases_with_scores_05022025.rds"
cohort_specific_variables <- c("study_center", "covid_confinement")

# Directory where outputs/results will be stored. 
results_folder <- "results/0_descriptives"
current_date <- "test"
cohort_name <- "Bisc"



#==============================================================================#
#                       INIT  WORKING ENVIRONMENT                              #
#==============================================================================#


print(paste0('YOUR WORKING DIRECTORY IS:', getwd()))

# Load libraries 
library(dplyr)
library(lubridate) # mealtime data 
library(tableone) # basic descriptives
library(psych) # basic descriptives

source("script/Define_functions.R")

# Create results folder
out_folder <- file.path(results_folder, paste0(cohort_name, "_", current_date))
if (!dir.exists(out_folder)) {dir.create(out_folder, recursive=TRUE)}


#==============================================================================#
#                                                                              #
#                          GET DESCRIPTIVES                                    #----
#                                                                              #
#==============================================================================#

# Load data
final_metadata <- readRDS(metadata_path)  #389

#------------------------------------------------------------------------------#
#                       BASELINE CHARACTERISTICS (COVARIABLES)                 #----
#------------------------------------------------------------------------------#

# Summary statistics (TableOne): 
# - continuous variables:  means and standard deviations for continuous variables   
# - categorical variables: proportions (%) 

final_metadata$first_meal_num <- as.numeric(final_metadata$first_meal/3600) # maybe as.numeric
final_metadata$nighttime_fasting_num <- as.numeric(final_metadata$nighttime_fasting/3600)
final_metadata$eating_jetlag_num <- as.numeric(final_metadata$eating_jetlag/3600)
final_metadata$last_meal_num <- as.numeric(final_metadata$last_meal/3600) # maybe as.numeric
final_metadata$sleep_duration_num <- as.numeric(final_metadata$sleep_duration/3600) # maybe as.numeric


general_covars <- c( "maternal_age", "maternal_bmi", "maternal_educational_level", 
                     "maternal_smoking", "parity","maternal_physical_activity", 
                     "gestational_age", "Sex","child_ethnicity", "dna_contamination", 
                     "gestational_diabetes", "energy_intake", 
                     "sleep_duration_num", "psqi_score", 
                     
                     'eating_occasions', 'first_meal_num', 
                     'last_meal_num',"eating_jetlag_num", "nighttime_fasting_num","rMEDp", 
                     "gest_weight_gain_m_preg", 
                     
                     'sleep_midpoint_imp', 'fasting_before_sleep', 'fasting_after_sleep') 

# Add cohort specific covariables 
if (!is.null(cohort_specific_variables)){general_covars <- c(general_covars,cohort_specific_variables)}

# Get descriptive table one
descriptive_object <- CreateTableOne(vars = general_covars, data = final_metadata, includeNA=T) 
descriptive_table <- print(descriptive_object,  noSpaces = TRUE, varLabels=T)

# Save descriptive table
file_path <- file.path(out_folder, paste0("descriptives_table1_v2",current_date,".csv"))
write.table(descriptive_table , file = file_path, sep=";")

stop()

#------------------------------------------------------------------------------#
#                         MEALTIME AND SLEEP CHARACTERISTICS                   #----
#------------------------------------------------------------------------------#
# Subset mealtime data 
mealtime_vars <- c('first_meal_wday','last_meal_wday','first_meal_wend','last_meal_wend', 
                   'eating_occasions_wday',  'eating_occasions_wend',
                   'eating_occasions', 'first_meal', 'last_meal',
                   "nighttime_fasting_wday","nighttime_fasting_wend","nighttime_fasting",
                   "eating_window_wday", "eating_window_wend", "eating_midpoint_wday",
                   "eating_midpoint_wend","eating_jetlag", 
                   
                   "wakeup_wday", "wakeup_wend", "gobed_wday", "gobed_wend",
                   "sleep_duration_wday", "sleep_duration_wend","sleep_duration", 
                   "sleep_midpoint_imp",
                   "fasting_before_sleep_wday", "fasting_before_sleep_wend", "fasting_before_sleep",
                   "fasting_after_sleep_wday", "fasting_after_sleep_wend", "fasting_after_sleep")

mealtime_complete <- final_metadata[,mealtime_vars]

# Convert to character
mealtime_complete_chr <- as.data.frame(lapply(mealtime_complete, as.character))
descriptive_object <- CreateTableOne(vars = mealtime_vars, data = mealtime_complete_chr , includeNA=T) 
descriptive_table <- print(descriptive_object,  noSpaces = TRUE, varLabels=T)

file_path <- file.path(out_folder, paste0("descriptive_mealtime_sleep_",current_date,".csv"))
write.table(descriptive_table , file = file_path, sep=";")


# Check distribution of mealtime 
# Init df to store quantiles 
mealtime_distribute <- data.frame(matrix(ncol = 18, nrow = 0))

for (var_name in mealtime_vars) {
  x <- show_values(mealtime_complete,var_name)
  mealtime_distribute <- rbind(mealtime_distribute, x)
}

file_path <- file.path(out_folder, paste0("distribution_mealtime_sleep_",current_date,".csv"))
write.table(mealtime_distribute,  file_path, sep = ";", row.names = FALSE)

# Convert to from: 
# Food_group	    N	    Min	              P10	                P20	
# first_meal_wday	389	18000s (~5 hours)	25200s (~7 hours)	27000s (~7.5 hours)	

# To this: 
# Variable	              Min	P10	   P20	P25	   P30	P33.3	P40	   P50	P60	
# First Meal on Weekday	5:00	7:00	7:30	8:00	8:00	8:00	8:00	8:00	9:00	

convert_to_time <- function(seconds) {
  if (is.na(seconds)) {return(NA)}
  
  hours <- floor(seconds / 3600)
  minutes <- floor((seconds %% 3600) / 60)
  sprintf("%02d:%02d", hours, minutes)
}

mealtime_distribute_h <- mealtime_distribute

# Convert time-based variables to hours:minutes format
for (col in colnames(mealtime_distribute_h)) {
    mealtime_distribute_h[[col]] <- sapply(mealtime_distribute_h[[col]], function(x) {
      time_value <- as.numeric(sub("s.*", "", x))  # Extract the numeric part
      if (!is.na(time_value)) {
        convert_to_time(time_value)  # Convert to hh:mm format
      } else {
        x  # Return the original value if conversion isn't applicable
      }
    })
}

# Add eating occasions avriables and remove "s" + float 2units
mealtime_distribute_h[mealtime_distribute_h$Food_group == "eating_occasions_wday",] <- mealtime_distribute[mealtime_distribute$Food_group == "eating_occasions_wday",]
mealtime_distribute_h[mealtime_distribute_h$Food_group == "eating_occasions_wend",] <- mealtime_distribute[mealtime_distribute$Food_group == "eating_occasions_wend",]
mealtime_distribute_h[mealtime_distribute_h$Food_group == "eating_occasions",] <- mealtime_distribute[mealtime_distribute$Food_group == "eating_occasions",]

# Save
file_path <- file.path(out_folder, paste0("distribution_mealtime_sleep_hours_",current_date,".csv"))
write.table(mealtime_distribute_h,  file_path, sep = ";", row.names = FALSE)




#------------------------------------------------------------------------------#
#                   TABLE ONE FOR EARLY VS. NON-EARLY EATERS                   #----
#------------------------------------------------------------------------------#

# Generate Table One for early vs. non-early eaters
descriptive_object <- CreateTableOne(vars = general_covars, 
                                     strata = "first_meal_2c", 
                                     data = final_metadata, 
                                     includeNA = TRUE,
                                     test = TRUE) 
descriptive_table <- print(descriptive_object, noSpaces = TRUE, varLabels = TRUE)

# Save descriptive table
file_path <- file.path(out_folder, paste0("descriptives_table1_EarlyEat_",current_date,".csv"))
write.table(descriptive_table, file = file_path, sep = ";")


# Generate distribution based on early / late eaters 

# Initialize a data frame to store the distributions for both groups
mealtime_distribute_early_late <- data.frame(matrix(ncol = 18, nrow = 0))

# Define the groups
groups <- c("early", "non_early")

# Loop through each group and variable to calculate quantiles
for (group in groups) {
  for (var_name in mealtime_vars) {
    # Subset the data based on the group
    group_data <- mealtime_complete[final_metadata$first_meal_2c == group, ]
    
    # Calculate quantiles using the show_quantile_values function
    x <- show_values(group_data, var_name)
    
    # Add the group name to the output
    x <- cbind(Group = group, x)
    
    # Combine the results
    mealtime_distribute_early_late <- rbind(mealtime_distribute_early_late, x)
  }
}

file_path <- file.path(out_folder, paste0("distribution_mealtime_EarlyEat_",current_date,".csv"))
write.table(mealtime_distribute_early_late,  file_path, sep = ";", row.names = FALSE)


# CONVERT TO HOURS 

# Create a copy for conversion
mealtime_distribute_h <- mealtime_distribute_early_late

# Convert time-based variables to hours:minutes format
for (col in colnames(mealtime_distribute_h)) {
  mealtime_distribute_h[[col]] <- sapply(mealtime_distribute_h[[col]], function(x) {
    time_value <- as.numeric(sub("s.*", "", x))  # Extract the numeric part
    if (!is.na(time_value)) {
      convert_to_time(time_value)  # Convert to hh:mm format
    } else {
      x  # Return the original value if conversion isn't applicable
    }
  })
}

# Reassign eating occasions rows with formatted values
mealtime_distribute_h[mealtime_distribute_h$Group == "early" & mealtime_distribute_h$Food_group == "eating_occasions_wday", ] <- 
  mealtime_distribute_early_late[mealtime_distribute_early_late$Group == "early" & mealtime_distribute_early_late$Food_group == "eating_occasions_wday", ]

mealtime_distribute_h[mealtime_distribute_h$Group == "early" & mealtime_distribute_h$Food_group == "eating_occasions_wend", ] <- 
  mealtime_distribute_early_late[mealtime_distribute_early_late$Group == "early" & mealtime_distribute_early_late$Food_group == "eating_occasions_wend", ]

mealtime_distribute_h[mealtime_distribute_h$Group == "early" & mealtime_distribute_h$Food_group == "eating_occasions", ] <- 
  mealtime_distribute_early_late[mealtime_distribute_early_late$Group == "early" & mealtime_distribute_early_late$Food_group == "eating_occasions", ]

mealtime_distribute_h[mealtime_distribute_h$Group == "non_early" & mealtime_distribute_h$Food_group == "eating_occasions_wday", ] <- 
  mealtime_distribute_early_late[mealtime_distribute_early_late$Group == "non_early" & mealtime_distribute_early_late$Food_group == "eating_occasions_wday", ]

mealtime_distribute_h[mealtime_distribute_h$Group == "non_early" & mealtime_distribute_h$Food_group == "eating_occasions_wend", ] <- 
  mealtime_distribute_early_late[mealtime_distribute_early_late$Group == "non_early" & mealtime_distribute_early_late$Food_group == "eating_occasions_wend", ]

mealtime_distribute_h[mealtime_distribute_h$Group == "non_early" & mealtime_distribute_h$Food_group == "eating_occasions", ] <- 
  mealtime_distribute_early_late[mealtime_distribute_early_late$Group == "non_early" & mealtime_distribute_early_late$Food_group == "eating_occasions", ]


# Save
file_path <- file.path(out_folder, paste0("distribution_mealtime_hours_EarlyEat_",current_date,".csv"))
write.table(mealtime_distribute_h,  file_path, sep = ";", row.names = FALSE)






#------------------------------------------------------------------------------#
#                         TABLE ONE FOR SHORT VS. LONG FASTING                 #----
#------------------------------------------------------------------------------#

# Generate Table One for early vs. non-early eaters
descriptive_object <- CreateTableOne(vars = general_covars, 
                                     strata = "nighttime_fasting_2c", 
                                     data = final_metadata, 
                                     includeNA = TRUE,
                                     test = TRUE) 
descriptive_table <- print(descriptive_object, noSpaces = TRUE, varLabels = TRUE)

# Save descriptive table
file_path <- file.path(out_folder, paste0("descriptives_table1_ShortFasting_",current_date,".csv"))
write.table(descriptive_table, file = file_path, sep = ";")


# Generate distribution based on early / late eaters 

# Initialize a data frame to store the distributions for both groups
mealtime_distribute_short_long <- data.frame(matrix(ncol = 18, nrow = 0))

# Define the groups
groups <- c("short", "long")

# Loop through each group and variable to calculate quantiles
for (group in groups) {
  for (var_name in mealtime_vars) {
    # Subset the data based on the group
    group_data <- mealtime_complete[final_metadata$nighttime_fasting_2c == group, ]
    
    # Calculate quantiles using the show_quantile_values function
    x <- show_values(group_data, var_name)
    
    # Add the group name to the output
    x <- cbind(Group = group, x)
    
    # Combine the results
    mealtime_distribute_short_long <- rbind(mealtime_distribute_short_long, x)
  }
}

file_path <- file.path(out_folder, paste0("distribution_mealtime_ShortFasting_",current_date,".csv"))
write.table(mealtime_distribute_short_long,  file_path, sep = ";", row.names = FALSE)


# CONVERT TO HOURS 

# Create a copy for conversion
mealtime_distribute_h <- mealtime_distribute_short_long

# Convert time-based variables to hours:minutes format
for (col in colnames(mealtime_distribute_h)) {
  mealtime_distribute_h[[col]] <- sapply(mealtime_distribute_h[[col]], function(x) {
    time_value <- as.numeric(sub("s.*", "", x))  # Extract the numeric part
    if (!is.na(time_value)) {
      convert_to_time(time_value)  # Convert to hh:mm format
    } else {
      x  # Return the original value if conversion isn't applicable
    }
  })
}

# Reassign eating occasions rows with formatted values
mealtime_distribute_h[mealtime_distribute_h$Group == "short" & mealtime_distribute_h$Food_group == "eating_occasions_wday", ] <- 
  mealtime_distribute_short_long[mealtime_distribute_short_long$Group == "short" & mealtime_distribute_short_long$Food_group == "eating_occasions_wday", ]

mealtime_distribute_h[mealtime_distribute_h$Group == "short" & mealtime_distribute_h$Food_group == "eating_occasions_wend", ] <- 
  mealtime_distribute_short_long[mealtime_distribute_short_long$Group == "short" & mealtime_distribute_short_long$Food_group == "eating_occasions_wend", ]

mealtime_distribute_h[mealtime_distribute_h$Group == "short" & mealtime_distribute_h$Food_group == "eating_occasions", ] <- 
  mealtime_distribute_short_long[mealtime_distribute_short_long$Group == "short" & mealtime_distribute_short_long$Food_group == "eating_occasions", ]

mealtime_distribute_h[mealtime_distribute_h$Group == "long" & mealtime_distribute_h$Food_group == "eating_occasions_wday", ] <- 
  mealtime_distribute_short_long[mealtime_distribute_short_long$Group == "long" & mealtime_distribute_short_long$Food_group == "eating_occasions_wday", ]

mealtime_distribute_h[mealtime_distribute_h$Group == "long" & mealtime_distribute_h$Food_group == "eating_occasions_wend", ] <- 
  mealtime_distribute_short_long[mealtime_distribute_short_long$Group == "long" & mealtime_distribute_short_long$Food_group == "eating_occasions_wend", ]

mealtime_distribute_h[mealtime_distribute_h$Group == "long" & mealtime_distribute_h$Food_group == "eating_occasions", ] <- 
  mealtime_distribute_short_long[mealtime_distribute_short_long$Group == "long" & mealtime_distribute_short_long$Food_group == "eating_occasions", ]


# Save
file_path <- file.path(out_folder, paste0("distribution_mealtime_hours_SortFasting_",current_date,".csv"))
write.table(mealtime_distribute_h,  file_path, sep = ";", row.names = FALSE)



#------------------------------------------------------------------------------#
#                         CORRELATION MEALTIMe - WGW                           #----
#------------------------------------------------------------------------------#

library(corrplot)

mealtime_vars <- c('first_meal_wday','last_meal_wday','first_meal_wend','last_meal_wend', 
                   'eating_occasions_wday',  'eating_occasions_wend',
                   'eating_occasions', 'first_meal', 'last_meal',
                   "nighttime_fasting_wday","nighttime_fasting_wend","nighttime_fasting",
                   "eating_window_wday", "eating_window_wend", "eating_midpoint_wday",
                   "eating_midpoint_wend","eating_jetlag", 
                   
                   "wakeup_wday", "wakeup_wend", "gobed_wday", "gobed_wend",
                   "sleep_duration_wday", "sleep_duration_wend","sleep_duration")


# Subset relevant data
vars_to_corr <- c("gest_weight_gain_m_preg", mealtime_vars)
corr_data <- final_metadata[ , vars_to_corr]

# Remove rows with NA
corr_data <- na.omit(corr_data)

# Compute correlation matrix
corr_matrix <- cor(corr_data, use = "complete.obs")

# Extract only correlations of mealtime vars with gestational weight gain
cor_with_gwg <- corr_matrix["gest_weight_gain_m_preg", -which(colnames(corr_matrix) == "gest_weight_gain_m_preg")]

# Create matrix for plotting
cor_matrix_plot <- matrix(cor_with_gwg, nrow = 1)
colnames(cor_matrix_plot) <- names(cor_with_gwg)
rownames(cor_matrix_plot) <- "GWG"

# Plot correlation
corrplot(cor_matrix_plot, method = "color", tl.cex = 0.7, tl.col = "black", is.corr = FALSE, cl.pos = "n")


##### more cool plot 

# Select relevant variables
vars_to_corr <- c("gest_weight_gain_m_preg", mealtime_vars)
corr_data <- final_metadata[, vars_to_corr]

# Remove rows with missing values
corr_data <- na.omit(corr_data)

# Compute correlation matrix
corr_matrix <- cor(corr_data, use = "complete.obs")

# Function to compute p-values for correlation matrix
cor.mtest <- function(mat, conf.level = 0.95) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], conf.level = conf.level)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

# Compute p-values
p.mat <- cor.mtest(corr_data)

# Extract correlation and p-values for GWG
cor_with_gwg <- corr_matrix["gest_weight_gain_m_preg", -which(colnames(corr_matrix) == "gest_weight_gain_m_preg")]
pvals_with_gwg <- p.mat["gest_weight_gain_m_preg", -which(colnames(p.mat) == "gest_weight_gain_m_preg")]

# Sort by absolute correlation
sorted_idx <- order(abs(cor_with_gwg), decreasing = TRUE)
cor_sorted <- cor_with_gwg[sorted_idx]
p_sorted <- pvals_with_gwg[sorted_idx]

# Create matrix for plotting
cor_matrix_plot <- matrix(cor_sorted, nrow = 1)
p_matrix_plot <- matrix(p_sorted, nrow = 1)
colnames(cor_matrix_plot) <- names(cor_sorted)
colnames(p_matrix_plot) <- names(p_sorted)
rownames(cor_matrix_plot) <- "GWG"
rownames(p_matrix_plot) <- "GWG"

# Plot with significance
corrplot(cor_matrix_plot,
         method = "color",
         is.corr = FALSE,
         addCoef.col = "black",  # Add correlation values
         tl.col = "black",
         tl.cex = 0.7,
         number.cex = 0.7,
         p.mat = p_matrix_plot,
         sig.level = c(0.001, 0.01, 0.05),
         insig = "label_sig",
         pch.cex = 1.2,
         pch.col = "black",
         cl.pos = "n")


# Create a dataframe from the correlation and p-value vectors
cor_p_df <- data.frame(
  variable = names(cor_sorted),
  correlation = as.numeric(cor_sorted),
  p_value = as.numeric(p_sorted)
)

# Optional: round values for easier viewing
cor_p_df$correlation <- round(cor_p_df$correlation, 3)
cor_p_df$p_value <- signif(cor_p_df$p_value, 3)

# View the dataframe
print(cor_p_df)
