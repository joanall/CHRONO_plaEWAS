#==============================================================================#
#                                                                              #
#                    MATERNAL DIET MEALTIME AND PLACENTA DNA METHYLATION       #
#                                   EWAS Analysis                              #
#                                                                              #                                  
#                   Camille Lassale & Joana Llauradó Pont                      #         
#                               Janaury 2025                                   #
#                                                                              #
#==============================================================================#



#==============================================================================#
#                                                                              #
#                       INIT  WORKING ENVIRONMENT                              #----
#                                                                              #
#==============================================================================#


print(paste0('YOUR WORKING DIRECTORY IS:', getwd()))

# Load libraries 
library(dplyr)
library(dmrff) # dmrs

source("script/Define_functions.R")


#==============================================================================#
#                                                                              #
#                      Differential methylated region (DMRs)                   #----
#                                                                              #
#==============================================================================#

ewas_files <- c("results/1_raw_results/Bisc_05022025/nighttime_fasting_Bisc_model_3_sep_exp_05022025.rds", 
                "results/1_raw_results/Bisc_05022025/nighttime_fasting_Bisc_model_3_05022025.rds", # adjsuted first meal
                "results/1_raw_results/Bisc_05022025/first_meal_Bisc_model_3_sep_exp_05022025.rds",
                "results/1_raw_results/Bisc_05022025/first_meal_Bisc_model_3_firstlast_05022025.rds",
                "results/1_raw_results/Bisc_05022025/last_meal_Bisc_model_3_sep_exp_05022025.rds",
                "results/1_raw_results/Bisc_05022025/last_meal_Bisc_model_3_firstlast_05022025.rds",
                "results/1_raw_results/Bisc_05022025/eating_occasions_Bisc_model_3_sep_exp_05022025.rds",
                "results/1_raw_results/Bisc_05022025/eating_occasions_Bisc_model_3_firstlast_05022025.rds",
                "results/1_raw_results/Bisc_05022025/eating_jetlag_Bisc_model_3_05022025.rds")

meth_path <- "db/bVals_BISC_nosib_noout_nomask.RDS"

cohort_name <- 'Bisc'
array <- 'EPIC'
save_to  <- c("results/1_raw_results/Bisc_05022025/nighttime_fasting_Bisc_model_3_sep_exp_05022025_DMR.rds", 
              "results/1_raw_results/Bisc_05022025/nighttime_fasting_Bisc_model_3_05022025_DMR.rds", # adjsuted first meal
              "results/1_raw_results/Bisc_05022025/first_meal_Bisc_model_3_sep_exp_05022025_DMR.rds",
              "results/1_raw_results/Bisc_05022025/first_meal_Bisc_model_3_firstlast_05022025_DMR.rds",
              "results/1_raw_results/Bisc_05022025/last_meal_Bisc_model_3_sep_exp_05022025_DMR.rds",
              "results/1_raw_results/Bisc_05022025/last_meal_Bisc_model_3_firstlast_05022025_DMR.rds",
              "results/1_raw_results/Bisc_05022025/eating_occasions_Bisc_model_3_sep_exp_05022025_DMR.rds",
              "results/1_raw_results/Bisc_05022025/eating_occasions_Bisc_model_3_firstlast_05022025_DMR.rds",
              "results/1_raw_results/Bisc_05022025/eating_jetlag_Bisc_model_3_05022025_DMR.rds")



current_date <- "05022025"

# args <- commandArgs(trailingOnly = TRUE)
# which_file <- args[1]
# which_file  <- as.numeric(which_file)
which_file <- 9

results_path <- ewas_files[[which_file]]
print(results_path)


# Read files 
meth_data <- as.data.frame(readRDS(meth_path))
results <- readRDS(results_path)

cat("\n-----------------------  DMRs (meffil)   ----------------------------\n")

# Get annotation
if (array == "450K") {
  library('IlluminaHumanMethylation450kanno.ilmn12.hg19')
  annotation <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)}
if (array == "EPIC") {
  library('IlluminaHumanMethylationEPICanno.ilm10b2.hg19')
  annotation <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)}

annotation <- as.data.frame(annotation)

# Match data ewas result - annotation 
meth_data <- meth_data[na.omit(match(rownames(results), rownames(meth_data))), ]
result_match <- results[na.omit(match(rownames(meth_data), rownames(results))), ]
ifelse(all(result_match$CpG == rownames(meth_data)), 
       "Meth and results data successfully matched :)", 
       "Data not matched :(")

# Match data: meth data - annotation
annotation <- annotation[na.omit(match(rownames(meth_data), annotation$Name)), ]
meth_data <- meth_data[na.omit(match(annotation$Name, rownames(meth_data))), ]
ifelse(all(annotation$Name == rownames(meth_data)), 
       "Meth and annotation data successfully matched :)", 
       "Data not matched :(")



meth_data <- as.matrix(meth_data)

# Get DMRs 
# p.cutoff Unadjusted p-value cutoff for membership in a candidate DMR
# (Default: 0.05).

dmr <- dmrff (estimate = result_match$Effect,
              se = result_match$SE,
              p.value = result_match$pvalue,
              methylation = meth_data,
              chr = annotation$chr,
              pos = annotation$pos,
              p.cutoff = 0.05,
              maxgap = 500,
              verbose=T)

#write_rds(dmr.model, paste0("results/6_dmr.ewas.model.",ewaspcutoff,"cutoff.3.B.rds"))

saveRDS(dmr, file =  save_to[[which_file]])

  






