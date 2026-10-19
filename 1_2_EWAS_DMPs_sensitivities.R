#==============================================================================#
#                                                                              #
#                    MATERNAL DIET MEALTIME AND PLACENTA DNA METHYLATION       #
#                                   EWAS Analysis                              #
#                                                                              #                                  
#                   Camille Lassale & Joana Llauradó Pont                      #         
#                               Janaury 2025                                   #
#                                                                              #
#==============================================================================#

# Based on maternal vegetarian/plant-based diet and coord blood DNAm by Peiyuan (Ran) Huang)

# - metadata.rds : RDS file with dataframe where colnames are covariables to 
#   adjust the EWAS models and rows are the samples/participants. Complete cases. 
#   See detailed description for metadata format.

# - methylation.rds : RDS file with a matrix of methylation Illumina beta values. 
#   Colnames are samples and rownames are probes on the array (450k or EPIC). 

# [OUTPUT]: 
# - exposure_model.rds : RDS file xith a dataframewhere rows are probes/CpGs and
#   columns have the statistics of the models. Colnames are :xEffect, SE, tvalue, 
#   pvalue, CpG.



# Description ------------------------------------------------------------------
# [INPUT DATA]: There are two input files required for this analysis:
#
# 1. metadata: A dataframe containing all the variables data needed for this project. 
#    Each row is an individual and each column is a different variable. 
#    Necessary variable names are: 
#
#    - Identifiers (1): "Basename" (sample id or mother or child id) 
#    - Covariates (11): "Sex", "maternal_educational_level", "maternal_smoking", 
#                      "maternal_age", "parity", "gestational_age", "maternal_bmi", 
#                      "maternal_physical_activity", "child_ethnicity", "gestational_diabetes",
#                      "dna_contamination", "pregnancy_complications"
#    
#    - Food groups in gr/day (8): 'cereals', 'fruit_nuts', 'meat', 'legumes','fish'
#                                 'legumes','olive_oil', 'vegetables', 'dairy'
#                                 
#    - Mealtime hour (4)  first and last time of eating on weekdays and weekends. 
#
#    - Number of eating occasions (1): average number of eating occasions (numeric) named as: "average_eating_occasions"
#
#    - Cell types (6): 'Trophoblasts','Stromal', 'Hofbauer', 'Endothelial','nRBC', 'Syncytiotrophoblast'
#       Please use the placenta reference panel from the 3st trimester (term) implemented in the planet R package (V Yuan et al.2021)
#
#
# If these columns are named differently in your dataset, please rename the columns accordingly
# 

# 2. meth: A matrix of methylation Illumina beta values. Each column is a sample and each row is a probe on the array (450k or EPIC). 
#    Column names must correspond to the Basename column in metadata.

# Which_model (input variable see section below) = choose model to run: 

# MAIN EWAS MODELS: 

# 1 = DNAm ~ nighttime_fasting + first_meal + num_eating  + covariables model 1 
# 2 = DNAm ~ nighttime_fasting + first_meal + num_eating  + covariables model 2 
# 3 = DNAm ~ nighttime_fasting + first_meal + num_eating  + covariables model 3 

# 4 = DNAm ~ eating jetlag + covariables model 1
# 5 = DNAm ~ eating jetlag + covariables model 2
# 6 = DNAm ~ eating jetlag + covariables model 3


# SENSEITIVTY EWAS MODELS (use only of covars model 1)

# Sens1: no gestational diabetes 
# 7 = DNAm ~ nighttime_fasting + first_meal + num_eating  +  model_1 
# 8 = DNAm ~ eating jetlag + model_1 -> no gestational diabetes

# Sens2: largest ethnicity 
# 9 = DNAm ~ nighttime_fasting +  first_meal + num_eating +  model_1
# 10 = DNAm ~ eating jetlag + model_1 ->  largest ethnicity 

# Sens4:non confined due to covid
# 13 = DNAm ~ nighttime_fasting(*) + first_meal + num_eating  + model_1 
# 14 = DNAm ~ eating jetlag + model_1 ->  non confined due to covid

# Sens7: last meal --> with all 3 models in separate code 
# 19, 20, 21 -> 

# Sens6: separates exposures --> with all 3 models 
# 22, 23, 24 
# 25, 26, 27
# 28, 29, 30 
# 31, 32, 33


# [NOTE]: Computing scores takes only a few minutes. 
# To compute one model takes around 3.30h :) 


# If you have any questions, please contact:
#   - Joana Llauradó Pont: joana.llaurado@isglobal.org
#   - Camille Lassale: camille.lasalle@isglobal.org
#-------------------------------------------------------------------------------

#==============================================================================#
#                                                                              #
#                             DEFINE PARAMETERS                                #----
#                                                                              #
#==============================================================================#

# Write directory where you have your data
metadata_path <- "db/metadata_plaEWAS_MT_complete_cases_with_scores_110012025.rds"
methylation_path <- "db/bVals_BISC_nosib_noout_nomask.RDS"
filter_450K_cpgs_path <- "db/filter_450K.RData" 
filter_EPIC_cpgs_path <- "db/filter_EPIC.RData"

# Write the name of your cohort and the type of array
cohort_name <- 'Bisc'
array <- 'EPIC'


# For sensitivty analaysis write the catgeory name of your largest ethnic grup 
# (e.g in BiSC most of the children are European and the category is called 'eur')
ethnicity <- 'eur'

# Directory where outputs/results will be stored. Write full path.
#results_path <- "/PROJECTES/BISC_OMICS/analyses/BiSC_23/EWAS_diet_CL/results"
#results_name <- "/ewas_mealtime/1_raw_results"
results_folder <- "results/1_raw_results/"

current_date <- "14012025"

# From the list above select a model to run 
#which_model <- 5

# [IGNORE and Don't uncoment - It's for personally running in bash array]
args <- commandArgs(trailingOnly = TRUE)
which_model <- args[1]
which_model <- as.numeric(which_model)
#which_model <- 1

# Specify any cohort specific variable you may have
# [NOTE]: BiSC cohort was recruited/Assessed during covid time and by 3 differents hopsitals

cohort_specific_variables <- c("study_center", "covid_confinement")

#-------------------------------------------------------------------------------

#==============================================================================#
#                                                                              #
#                       INIT  WORKING ENVIRONMENT                              #----
#                                                                              #
#==============================================================================#


print(paste0('YOUR WORKING DIRECTORY IS:', getwd()))

# Load libraries 
library(dplyr)
library(lubridate) # mealtime data 
library(tableone) # basic descriptives
library(psych) # basic descriptives
#library(PACEanalysis) # ewas
library(dmrff) # dmrs
library(doParallel)
library(foreach)
library(pryr)  # For memory usage tracking
source("script/Define_functions.R")
# Create results folder
if (!dir.exists(results_folder)) {dir.create(results_folder)}
out_folder <- file.path(results_folder, paste0(cohort_name, "_", current_date))
if (!dir.exists(out_folder)) {dir.create(out_folder)}



#==============================================================================#
#                                                                              #
#                     SECTION 1 : CHECK INPUT DATA                             #----
#                                                                              #
#==============================================================================#

cat("\n---------------------   START OF CODE  --------------------------------\n")


# Load data 
metadata <- readRDS(metadata_path)  #(389 x 67)
meth_data <- readRDS(methylation_path) #  565

#meth_data <- meth_data[1:200,]

cat("\n--------------------- CHECK INPUT DATA    ----------------------------\n")

# Check ids methylation and variables match in order and in length 

cat('Number of samples: ', nrow(meth_data), "\n")
cat('Number of covariables : ', ncol(metadata), "\n")


## MATCHING & ORDERING METADATA AND M/BVALS 
# Check number of obs/samples with exposures and bVals 

ids_metadata_in_meth <- sum(metadata$Basename %in% colnames(meth_data)) # 565 but 448 with Sex(TRUE)
ids_meth_in_metadata  <- sum(colnames(meth_data) %in% metadata$Basename) # 565 (TRUE)

if (ids_metadata_in_meth == ids_meth_in_metadata) {cat('IDs metadata and methylation matched :)\n')
} else {cat('IDs metadata and methylation do not match :(\n')}


# just keep those ones will methylation data # Use only those participants with both data 
metadata  <- metadata[(metadata$Basename %in% colnames(meth_data)),] # 565  59
meth_data <- meth_data[,(colnames(meth_data) %in% metadata$Basename)] # 5001  565



metadata$first_meal <- as.numeric(metadata$first_meal)/3600
metadata$nighttime_fasting <- as.numeric(metadata$nighttime_fasting)/3600
metadata$last_meal <- as.numeric(metadata$last_meal)/3600
metadata$eating_jetlag <- as.numeric(metadata$eating_jetlag)/3600
metadata$sleep_duration <- as.numeric(metadata$sleep_duration)/3600


#------------------------------------------------------------------------------#
#                          SET UP COMMON PARAMETERS                            #----
#------------------------------------------------------------------------------#


#  Define variables to include in the models 
ct_covars <- c('Trophoblasts','Stromal', 'Hofbauer', 'Endothelial','nRBC', 'Syncytiotrophoblast')

m1_covars <- c("maternal_age","parity","maternal_educational_level", "maternal_smoking", 
               "gestational_age", "child_ethnicity", "energy_intake","Sex",'rMEDp',
               'maternal_physical_activity',
               "dna_contamination")

m2_covars <- c('maternal_bmi')
m3_covars <- c('sleep_duration', 'psqi_score')

# Define models 
if (!is.null(cohort_specific_variables)){model_1 <- c(m1_covars, ct_covars, cohort_specific_variables)}
model_2 <- c(model_1, m2_covars)
model_3 <- c(model_2, m3_covars)



#------------------------------------------------------------------------------#
#                         SET UP SPECIFIC PARAMETERS                           #----
#------------------------------------------------------------------------------#



# Get exposure
if (which_model %in% c(1,2,3)) {
  exposure <- c('nighttime_fasting', "eating_occasions", "first_meal")}

if (which_model %in% c(4,5,6)) {exposure <- 'eating_jetlag'}


# Get analysis name 
if (which_model %in% c(1,4)) {analysis_name <- 'model_1'; covars <- model_1}
if (which_model %in% c(2,5)) {analysis_name <- 'model_2'; covars <- model_2}
if (which_model %in% c(3,6)) {analysis_name <- 'model_3'; covars <- model_3}



# All sensitivities with only covars from model 1 
# !!!!!! In previous version the seconda analysis were with nighttime_fasting,
# now we'll do them with first -aslt meal comb 
# Get exposure
if (which_model %in% c(7, 9, 11, 13, 15, 17)) {
  exposure <- c('nighttime_fasting', "eating_occasions", "first_meal")
  covars <- model_1}

if (which_model %in% c(8, 10, 12, 14, 16, 18)) {
  exposure <- 'eating_jetlag'
  covars <- model_1}

# Get analysis name  and sensitity covars
if (which_model %in% c(7,8)) {
  analysis_name <- 'model_noGDM_adj_first_meal'
  var_restrict <-  "gestational_diabetes" 
  var_indicator <-  "No"
}

if (which_model %in% c(9,10)) {
  analysis_name <- 'model_ethnicity_adj_first_meal'
  var_restrict <-  "child_ethnicity" 
  var_indicator <-  ethnicity
}


# to properly prepare  this ones 
if (which_model %in% c(13,14)) { # do not run, variables being updated 
  analysis_name <- 'model_noconf_adj_first_meal'
  var_restrict <-  "covid_confinement" 
  var_indicator <- "non_confined"
  
}

if (which_model %in% c(15,16)) {
  analysis_name <- 'early_eaters_adj_first_meal'
  var_restrict <-  "first_meal_2c" 
  var_indicator <-  "early"
  
}

if (which_model %in% c(17,18)) {
  analysis_name <- 'short_fasting_adj_first_meal'
  var_restrict <-  "nighttime_fasting_2c" 
  var_indicator <-  "short"
  
}


# if (which_model %in% c(19,20,21)) {
#   exposure <- c('nighttime_fasting', "eating_occasions", "last_meal")
# }

# run with last meal 
# if (which_model == 19) {analysis_name <- 'model_1_last_meal'; covars <- model_1}
# if (which_model == 20) {analysis_name <- 'model_2_last_meal'; covars <- model_2}
# if (which_model == 21) {analysis_name <- 'model_3_last_meal'; covars <- model_3}



# RE run models for separate exposures 

if (which_model %in% c(22,23,24)) {exposure <- c('nighttime_fasting')}
if (which_model %in% c(25,26,27)) {exposure <- c("first_meal")}
if (which_model %in% c(28,29,30)) {exposure <- c("eating_occasions")}
if (which_model %in% c(31,32,33)) {exposure <- c("last_meal")}

if (which_model %in% c(22,25,28,31)) {analysis_name <- 'model_1_sep_exp'; covars <- model_1}
if (which_model %in% c(23,26,29,32)) {analysis_name <- 'model_2_sep_exp'; covars <- model_2}
if (which_model %in% c(24,27,30,33)) {analysis_name <- 'model_3_sep_exp'; covars <- model_3}

if (which_model %in% c(34,35,36)) {
  exposure <- c("first_meal","eating_occasions", "last_meal")}


if (which_model == 34) {analysis_name <- 'model_1_firstlast'; covars <- model_1}
if (which_model == 35) {analysis_name <- 'model_2_firstlast'; covars <- model_2}
if (which_model == 36) {analysis_name <- 'model_3_firstlast'; covars <- model_3}


cat('\nMODEL SELECTED FROM LIST IS:', which_model,"\n")
cat("Exposure:", exposure,"\n")
cat("Analysis: ", analysis_name, "\n")




#-------------------------------------------------------------------------------



#==============================================================================#
#                                                                              #
#                     EXCLUSION UNRELAIBLE CPG PROBES                          #----
#                                                                              #
#==============================================================================#

cat("\n----------------- EXCLUSION CpGs/PROBES ----------------------\n")

# Dataframes and fucntions used are Copy-pasted from EASIER package
# References: https://github.com/isglobal-brge/EASIER

if(array == "450K") {load(filter_450K_cpgs_path)}
if(array == "EPIC") {load(filter_EPIC_cpgs_path)}

ethnic <- 'GMAF1p'

# Prepare inputs to run EASIER functions 
exclude <-  c( 'control_probes', 'noncpg_probes', 'Sex',
               'MASK_mapping', 'MASK_sub30_copy', 'MASK_extBase', 'MASK_typeINextBaseSwitch', 
               'MASK_snp5_common','MASK_snp5_ethnic',
               'Unrel_450_EPIC_pla_restrict')
# 'MASK_snp5_ethnic' -> remove SNPs common in European population 

cat('Number CpGs before exclusion', nrow(meth_data), "\n")
n <- ncol(meth_data) # number of samples 
meth_data <- cbind(meth_data, probeID=rownames(meth_data))

file_path1 <- file.path(out_folder, paste0("list_excluded_probes_",current_date,".txt"))
file_path2 <- file.path(out_folder, paste0("descriptives_excluded_probes_",current_date,".txt"))

meth_data <- exclude_CpGs(meth_data, "probeID", 
                          exclude, 
                          ethnic = ethnic, 
                          filename = file_path1, 
                          fileresume = file_path2,
                          artype = array )


cat('Number CpGs after exclusion', nrow(meth_data), "\n")


# Arrange output from exclude_CpGs to fit rlm_single_cpg function:
# - Metadata: Samples (rownames) x covariables (colnames)
# - Methylation data: CpGs (rownames) x Samples (colnames)
# - Remove added colanmes by exclude_cpgs
# - Convert again methylaiton values to numeric 

meth_data <- as.data.frame(lapply(meth_data, type.convert, as.is = TRUE))
colnames(meth_data) <- sub("^X", "", colnames(meth_data))

rownames(meth_data) <- meth_data$probeID

to_remove <- c("probeID", "CpG_chrm", "CpG_beg", "CpG_end", "Unrel_450_EPIC_blood",
               "MASK_mapping", "MASK_typeINextBaseSwitch", "MASK_rmsk15", 
               "MASK_sub40_copy", "MASK_sub35_copy", "MASK_sub30_copy",
               "MASK_sub25_copy", "MASK_snp5_common", "MASK_snp5_GMAF1p", "
               MASK_extBase", "MASK_general", "probeType", 
               "Unrel_450_EPIC_pla_restrict", "Unrel_450_EPIC_pla", "MASK_extBase")

meth_data <- meth_data[ , -which(names(meth_data) %in% to_remove)]




#stop("Stopping the script execution here")

#==============================================================================#
#                                                                              #
#                            RUN EWAS MODELS                                   #----
#                                                                              #
#==============================================================================#

cat("\n----------------------- DO EWAS MODELS ------------------------------\n")

#------------------------------------------------------------------------------#
#                         PREPARE DATA FOR MODEL                               #----
#------------------------------------------------------------------------------#

# Subset data for sensitivites 
if (dplyr::between(which_model, 7, 18)) {
  metadata <- metadata[metadata[[var_restrict]] == var_indicator, ]
  
  # Remove the `var_restrict` column from final_metadata
  metadata <- metadata[, !colnames(metadata) %in% var_restrict]
  
  # Remove `var_restrict` from `covars` if it exists
  covars <- setdiff(covars, var_restrict)
  
}

meth_data <- meth_data[,(colnames(meth_data) %in% metadata$Basename)] # 5001  565

cat('N sample meta : ', nrow(metadata),"\n")
cat('N meth data : ', ncol(meth_data),"\n")

#------------------------------------------------------------------------------#
#                          SPLIT DATA INTO CHUNKS                              #----
#------------------------------------------------------------------------------#


# Set up data for parallelization 
sys_info <- benchmarkme::get_ram() #236 GB or 540 in yamabuki 
#n_cores <- parallel::detectCores() - 1 #47 or 72 in yamabuki 
n_cores <- 8
chunk_size <- ceiling(nrow(meth_data) / n_cores)
cpg_chunks <- split(as.data.frame(meth_data), ceiling(seq_len(nrow(meth_data)) / chunk_size))
n_chunk <- length(cpg_chunks)
cat('Number of chunks:', length(cpg_chunks), "\n")
cat('Chunk size:', chunk_size, "\n")
cat('Numer of cores:', n_cores, "\n")


#------------------------------------------------------------------------------#
#                              RUN MODEL                                       #----
#------------------------------------------------------------------------------#

cat('Start parallelization')
cl <- makeCluster(n_cores)  # Create a cluster with the available cores
registerDoParallel(cl)      # Register the parallel backend

chunk_results_list <- foreach(i = 1:n_chunk, 
                              .export = c('cpg_chunks', 'exposure', 'covars', 
                                          'metadata', 'cohort_name', 'analysis_name', 
                                          'current_date', 'out_folder')) %dopar% {
                                            
                                            # Define a log file for each chunk
                                            # log_file <- file.path(out_folder, paste0("chunk_", i, "_log.txt"))
                                            # sink(log_file, append = TRUE)  # Redirect output to log file
                                            # 
                                            # Log the start of the chunk
                                            cat(sprintf("Chunk %d starting at %s\n", i, Sys.time()))
                                            
                                            # Record the start time and memory usage
                                            start_time <- proc.time()
                                            initial_memory <- pryr::mem_used()
                                            
                                            # Perform the computation
                                            chunk_result <- rlm_single_cpg(cpg_chunk = cpg_chunks[[i]], 
                                                                           exposures = exposure, 
                                                                           covariates = covars, 
                                                                           metadata = metadata)
                                            
                                            # Calculate runtime and memory usage
                                            run_time <- proc.time() - start_time
                                            final_memory <- pryr::mem_used()
                                            memory_diff <- final_memory - initial_memory
                                            
                                            # Log the completion of the chunk
                                            cat(sprintf("Chunk %d completed at %s. Runtime: %.2f seconds. Memory used: %.2f MB\n", 
                                                        i, Sys.time(), run_time['elapsed'], as.numeric(memory_diff) / (1024^2)))
                                            
                                            sink()  # Stop redirecting output to log file
                                            return(chunk_result)
                                          }

stopCluster(cl)

cat("\n--------------- FINISHED PARALLELIZATION -----------------------------\n")


# Get and save results by exposure 
all_results_list <- do.call('c', chunk_results_list)

# Combining sublists based on names
all_results_list <- lapply(unique(names(all_results_list)), function(x) {
  do.call(rbind, all_results_list[x == names(all_results_list)])
})


# Assign the unique names back to the new combined list
names(all_results_list) <- exposure

# Assign CpGs as rownames 
for (exp in names(all_results_list)) {
  rownames(all_results_list[[exp]]) <- all_results_list[[exp]]$CpG
}       



# Save each list/exposure in an RDS file while preserving the data frame structure
lapply(names(all_results_list), function(exp) {
  saveRDS(all_results_list[[exp]], 
          file.path(out_folder, 
                    paste0(exp, '_',cohort_name,'_', analysis_name,'_', current_date,'.rds')))
})

cat("\n-------------------- EWAS MODEL FINISHED -----------------------------\n")








