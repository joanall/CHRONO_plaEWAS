#==============================================================================#
#                                                                              #
#               MATERNAL DIET MEALTIME AND PLACENTA DNA METHYLATION            #
#                             Analyze DMRs                                     #
#                                                                              #                                  
#                   Camille Lassale & Joana Llauradó Pont                      #         
#                               January 2025                                   #
#                                                                              #
#==============================================================================#

# [Description]:
# [INPUT]:
# [OUTPUT]:

#for i in {1..9}; do
#Rscript script_name.R $i
#done
# for i in {1..9}; do Rscript script/6_Analyze_singleEWAS_DMRs.R $i; done


#==============================================================================#
#                           DEFINE PATHS                                       #
#==============================================================================#
ewas_files <- c("results/2_quality_control/Bisc_05022025/nighttime_fasting_model_3_sep_exp/nighttime_fasting_Bisc_model_3_sep_exp_05022025_adjusted.rds",
                "results/2_quality_control/Bisc_05022025/nighttime_fasting_model_3/nighttime_fasting_Bisc_model_3_05022025_adjusted.rds",
                "results/2_quality_control/Bisc_05022025/first_meal_model_3_sep_exp/first_meal_Bisc_model_3_sep_exp_05022025_adjusted.rds",
                "results/2_quality_control/Bisc_05022025/first_meal_model_3_firstlast/first_meal_Bisc_model_3_firstlast_05022025_adjusted.rds",
                "results/2_quality_control/Bisc_05022025/last_meal_model_3_sep_exp/last_meal_Bisc_model_3_sep_exp_05022025_adjusted.rds",
                "results/2_quality_control/Bisc_05022025/last_meal_model_3_firstlast/last_meal_Bisc_model_3_firstlast_05022025_adjusted.rds",
                "results/2_quality_control/Bisc_05022025/eating_occasions_model_3_sep_exp/eating_occasions_Bisc_model_3_sep_exp_05022025_adjusted.rds",
                "results/2_quality_control/Bisc_05022025/eating_occasions_model_3_firstlast/eating_occasions_Bisc_model_3_firstlast_05022025_adjusted.rds",
                "results/2_quality_control/Bisc_05022025/eating_jetlag_model_3/eating_jetlag_Bisc_model_3_05022025_adjusted.rds") 


results_folder <- "results/6_analysis_dmrs/Bisc_05022025_combp"

# array <- 'EPIC'
# args <- commandArgs(trailingOnly = TRUE)
#i <- as.numeric(args[1])

i <- 1


#==============================================================================#
#                       INIT  WORKING ENVIRONMENT                              #
#==============================================================================#

library(ENmix)
# library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
# library(data.table)
# 
# if (array == "450K") {
#   library('IlluminaHumanMethylation450kanno.ilmn12.hg19')
#   annotation <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)}
# if (array == "EPIC") {
#   library('IlluminaHumanMethylationEPICanno.ilm10b2.hg19')
#   annotation <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)}




print(paste0('YOUR WORKING DIRECTORY IS:', getwd()))


# Create folder (if it doesn't exist) to save outputs (main ones)
if (!dir.exists(results_folder)) {dir.create(results_folder, recursive=TRUE)}

#==============================================================================#
#                             GET SIGNIFICANT DMRs                             #----
#==============================================================================#


ewas_result <- readRDS(ewas_files [i])
# CpGId chromosome  position eaf      beta       se  beta_95L beta_95U
# 1 cg00000029         -9        -9  -9 -0.001353 0.001303 -0.003908 0.001202
# 2 cg00000109          3 172198246  -9  0.000525 0.000398 -0.000256 0.001306
# 3 cg00000165          1  90729116  -9 -0.000075 0.000576 -0.001204 0.001055
# 4 cg00000236          8  42405775  -9  0.000051 0.000839 -0.001593 0.001695
# 5 cg00000321          8  41310282  -9  0.001274 0.000937 -0.000562 0.003111


N_CpGs <- dim(ewas_result)[1]
df = data.frame(matrix(nrow = N_CpGs, ncol = 0)) 

df$probe <- ewas_result$CpGId

df$chr <- ewas_result$chromosome 
df$p <- ewas_result$p.value

df$start <- ewas_result$position
df$end <- df$start +1 


# head(df)

# probe chr        p     start       end
# 1 cg00000029  -9 0.299166        -9        -8
# 2 cg00000109   3 0.187815 172198246 172198247
# 3 cg00000165   1 0.896819  90729116  90729117
# 4 cg00000236   8 0.951527  42405775  42405776
# 5 cg00000321   8 0.173652  41310282  41310283
# 6 cg00000363   1 0.949015 230425046 230425047


# seed ->  	 FDR significance threshold for initial selection of DMR region
DMR1 <- combp(df, dist.cutoff=1000, bin.size=310,
              seed=0.05,region_plot=TRUE, 
              mht_plot=TRUE, 
              nCores=10,
              verbose=TRUE)


