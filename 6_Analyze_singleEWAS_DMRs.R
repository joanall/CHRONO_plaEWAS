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
input_files <- c("results/1_raw_results/Bisc_05022025/nighttime_fasting_Bisc_model_3_sep_exp_05022025_DMR.rds", 
                 "results/1_raw_results/Bisc_05022025/nighttime_fasting_Bisc_model_3_05022025_DMR.rds", # adjsuted first meal
                 "results/1_raw_results/Bisc_05022025/first_meal_Bisc_model_3_sep_exp_05022025_DMR.rds",
                 "results/1_raw_results/Bisc_05022025/first_meal_Bisc_model_3_firstlast_05022025_DMR.rds",
                 "results/1_raw_results/Bisc_05022025/last_meal_Bisc_model_3_sep_exp_05022025_DMR.rds",
                 "results/1_raw_results/Bisc_05022025/last_meal_Bisc_model_3_firstlast_05022025_DMR.rds",
                 "results/1_raw_results/Bisc_05022025/eating_occasions_Bisc_model_3_sep_exp_05022025_DMR.rds",
                 "results/1_raw_results/Bisc_05022025/eating_occasions_Bisc_model_3_firstlast_05022025_DMR.rds",
                 "results/1_raw_results/Bisc_05022025/eating_jetlag_Bisc_model_3_05022025_DMR.rds")


ewas_files <- c("results/2_quality_control/Bisc_05022025/nighttime_fasting_model_3_sep_exp/nighttime_fasting_Bisc_model_3_sep_exp_05022025_adjusted.rds",
                "results/2_quality_control/Bisc_05022025/nighttime_fasting_model_3/nighttime_fasting_Bisc_model_3_05022025_adjusted.rds",
                "results/2_quality_control/Bisc_05022025/first_meal_model_3_sep_exp/first_meal_Bisc_model_3_sep_exp_05022025_adjusted.rds",
                "results/2_quality_control/Bisc_05022025/first_meal_model_3_firstlast/first_meal_Bisc_model_3_firstlast_05022025_adjusted.rds",
                "results/2_quality_control/Bisc_05022025/last_meal_model_3_sep_exp/last_meal_Bisc_model_3_sep_exp_05022025_adjusted.rds",
                "results/2_quality_control/Bisc_05022025/last_meal_model_3_firstlast/last_meal_Bisc_model_3_firstlast_05022025_adjusted.rds",
                "results/2_quality_control/Bisc_05022025/eating_occasions_model_3_sep_exp/eating_occasions_Bisc_model_3_sep_exp_05022025_adjusted.rds",
                "results/2_quality_control/Bisc_05022025/eating_occasions_model_3_firstlast/eating_occasions_Bisc_model_3_firstlast_05022025_adjusted.rds",
                "results/2_quality_control/Bisc_05022025/eating_jetlag_model_3/eating_jetlag_Bisc_model_3_05022025_adjusted.rds") 


results_folder <- "results/6_analysis_dmrs/Bisc_05022025"

array <- 'EPIC'
args <- commandArgs(trailingOnly = TRUE)
#i <- as.numeric(args[1])

i <- 9
#==============================================================================#
#                       INIT  WORKING ENVIRONMENT                              #
#==============================================================================#

library(dmrff)
library(missMethyl)

print(paste0('YOUR WORKING DIRECTORY IS:', getwd()))


# Create folder (if it doesn't exist) to save outputs (main ones)
if (!dir.exists(results_folder)) {dir.create(results_folder, recursive=TRUE)}

#==============================================================================#
#                             GET SIGNIFICANT DMRs                             #----
#==============================================================================#

# Load dmr file 
cat('Reading file...', basename(input_files[i]))
dmrs <- readRDS(file = input_files[i]) 

# Colnames: chr, star, end, n , B, S, estimate, se, z, p.value, p.adjust

# Paraemeterrs p.cutoff = 0.05 Unadjusted p-value cutoff for membership in a candidate DMR
# and GAP = 500 CpG
 
# Get significance based on 1st and 2nd criteria: padjust (Bonferroni) < 0.05 & n > 1)
sig_dmrs <- dmrs[which(dmrs$p.adjust < 0.05 & dmrs$n > 1),]

# Create a sequence of numbers based on the number of rows in sig_dmrs
sig_dmrs$DMR <- paste0("DMR_", seq_len(nrow(sig_dmrs)))

# Annotate DMRs 
if (array == "450K") {
  library('IlluminaHumanMethylation450kanno.ilmn12.hg19')
  annotation <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)}
if (array == "EPIC") {
  library('IlluminaHumanMethylationEPICanno.ilm10b2.hg19')
  annotation <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)}

annotation <- as.data.frame(annotation)

# Add CpGs to each sig region
sig_dmrs_anno <- dmrff.annotate(sig_dmrs, annotation)

# Rename columns to match ewas dmps results 
colnames(sig_dmrs_anno)[colnames(sig_dmrs_anno) == 'Name'] <- 'CpG'
colnames(sig_dmrs_anno)[colnames(sig_dmrs_anno) == 'p.value'] <- 'pvalue_dmr'
colnames(sig_dmrs_anno)[colnames(sig_dmrs_anno) == 'p.adjust'] <- 'pval_adj_dmr'
colnames(sig_dmrs_anno)[colnames(sig_dmrs_anno) == 'estimate'] <- 'estimate_dmr'
colnames(sig_dmrs_anno)[colnames(sig_dmrs_anno) == 'se'] <- 'se_dmr'


# Subset to only interesting annotation columns 
sig_dmrs_anno <- sig_dmrs_anno[,c("DMR","chr", "start", "end", "n",
                                  "estimate_dmr", "se_dmr", "pvalue_dmr", 
                                  "pval_adj_dmr", "UCSC_RefGene_Name",'CpG' )]


# Add stats for each CpG 
stats <- readRDS(file = ewas_files[i]) 
sig_dmrs_anno$pvalue_dmp <- stats$pvalue[match(sig_dmrs_anno$CpG, stats$CpG)]
sig_dmrs_anno$Effect_dmp <- stats$Effect[match(sig_dmrs_anno$CpG, stats$CpG)]
sig_dmrs_anno$SE_dmp <- stats$SE[match(sig_dmrs_anno$CpG, stats$CpG)]
sig_dmrs_anno$bf_dmp <- stats$bf[match(sig_dmrs_anno$CpG, stats$CpG)]

save_to <- paste0(results_folder, '/', sub(".rds", '_anno.csv', basename(input_files[i])))
write.table(sig_dmrs_anno, save_to, sep = ";", row.names = FALSE)


# Get significance based on 3rd criteria: 2  CpGs sites within a region
# beta values in the same direction and nominal p-value < 0.05 

# Sublist df based on DMR 
dmr_list <- split(sig_dmrs_anno, sig_dmrs_anno$DMR)

# Check within each DMR 3rd criteria 

final_sig_dmrs <- lapply(dmr_list, function(df) {
  # Remove rows where pvalue_dmp or Effect_dmp is NA
  df <- df[!is.na(df$pvalue_dmp) & !is.na(df$Effect_dmp), ]
  
  # Filter CpGs with pvalue_dmp < 0.05
  df <- df[df$pvalue_dmp < 0.05, ]
  
  # Check if at least 2 CpGs exist, else exclude dmr 
  if (nrow(df) < 2) { return(NULL)}
  
  # Get the effect direction of the DMR
  effect_dmr <- unique(df$estimate_dmr)
  
  # Ensure there is only one unique estimate_dmr value (otherwise, exclude)
  if (length(effect_dmr) != 1) {return(NULL)}
  
  # Count how many CpGs have the same sign as estimate_dmr
  same_direction_count <- sum(sign(df$Effect_dmp) == sign(effect_dmr), na.rm = TRUE)
  
  # Keep the DMR only if at least 2 CpGs have the same direction
  if (same_direction_count >= 2) {return(df)
  } else {return(NULL)}
})



# Remove NULL elements (DMRs that didn't meet criteria)
final_sig_dmrs <- Filter(Negate(is.null), final_sig_dmrs)

# Convert the filtered list back into a dataframe
final_sig_dmrs_df <- do.call(rbind, final_sig_dmrs)

# View the final dataframe
head(final_sig_dmrs_df)
  
# Save if sig dmrs 
if (length(final_sig_dmrs_df) > 0) {
 
  save_to <- paste0(results_folder, '/', sub(".rds", '_sig_final.csv', basename(input_files[i])))
  write.table(final_sig_dmrs_df, save_to, sep = ";", row.names = FALSE)

}



#==============================================================================#
#                         DMRs ENRICHMENT (GO&KEGG COLLECTIONS)                #---- 
#==============================================================================#

if (!is.null(final_sig_dmrs_df)) {
  # Transform into GRanges object 
  dmrs_GR <- makeGRangesFromDataFrame(final_sig_dmrs_df)
  
  # GO enrichment 
  gst.region <- goregion(dmrs_GR, collection="GO", array.type="EPIC", plot.bias=TRUE) 
  
  
  if (sum(gst.region$FDR < 0.05) > 0) { 
    save_to <- file.path(results_folder, sub(".rds", '_DMR_GOenrich.csv', basename(input_files[i])))
    write.table(gst.region, save_to, sep = ";", row.names = FALSE)
  } 
  
  #topGSA(gst.region, n=15)
  
  # KEGG enrichment
  gst.region.kegg <- goregion(dmrs_GR, plot.bias=TRUE, 
                              collection="KEGG", array.type="EPIC")
  
  #topGSA(gst.region.kegg, n=15)
  if (sum(gst.region$FDR < 0.05) > 0) { 
    save_to <- file.path(results_folder, sub(".rds", '_DMR_KEGGenrich.csv', basename(input_files[i])))
    write.table(gst.region.kegg, save_to, sep = ";", row.names = FALSE)
  }

}








