#==============================================================================#
#                                                                              #
#               MATERNAL DIET MEALTIME AND PLACENTA DNA METHYLATION            #
#                             Functional Annotation                            #
#                                                                              #                                  
#                   Camille Lassale & Joana Llauradó Pont                      #         
#                               January 2025                                   #
#                                                                              #
#==============================================================================#

# Description:  We perform a series of Functional Analysis to gain insight about the 
# biology underlying a list of CpG/genes either by attaching biological information 
# to sequences of genes or proteins  (Functional Annotation) or by determining 
# classes of genes or proteins that are over-represented in a large group of 
# genes or proteins, and may have relations with disease phenotypes (Functional Enrichment)
# Additionally we look into specialized databases in EWAS data, to look at previous
# reported information (exposures, studies, related diseases etc) given a set of CpG. 

# - annotation hg19/37 
# - annotation EPIraction
# - annotation of eQTMs (not done yet) 
# - annotation of mQTLs
# - annotation 15 chromatin states 
# - annotation PMD
# - annotation gDMR
# snps ? imprinted regions? eForge? 


# INPUT FILE: RDS File (list of files) storing a dataframe where rows are CpGs and columns: 
# - CpGs: chr/string CpG Id
# - n : number of samples to compute the model 
# - coef: coefficient of the exposure of interest
# - se: standard error of the model
# - teststat: --- 
# - pval: p value

# INPUT DB: 
# Since this is an annotation script, some of the annotation are directly imported
# from libraries while some others have to be downloaded from different sources. 
# The list of annotation databases that have to be downloaded are: 
# - Hg38 annotation -> UCSC go to anotation EPIraction section and follow the steps
# - eQTMs -> Make it manually by saving reported eQTMS in this paper: XXX
# - mQTLs -> From BiSC cohort study. Email the responisble and ask for it.
# - Chramatin States -> FP_15_E091 and FP_18_E091 from the following github:
#   https://github.com/isglobal-brge/easier/tree/main/data
# - PMD -> Donwload file named PMD_placenta from https://github.com/isglobal-brge/EASIER/tree/main/data
# - GMD -> Download file named IR_Placenta from  https://github.com/isglobal-brge/EASIER/tree/main/data

# [OUTPUT]


#==============================================================================#
#                                                                              #
#                       SECTION 0: USER INPUT                                  #
#                                                                              #
#==============================================================================#

# Note: The example paths are relative paths as this script is part or an R.Project
# and it executes from root. If it doens't find the paths, write complete paths.

cohort_name <- 'Bisc'
array <- 'EPIC'

# Output folder 
results_path <- "results/4_functional_analysis_annotation/Bisc_13062025"
if (!dir.exists(results_path)) {dir.create(results_path, recursive=TRUE)}

# annotated_files <- c(annotated_file1, annotated_file2,...)
file <- "results/2_quality_control/Main_13062025/last_meal_model_1_sep_exp/last_meal_model_1_sep_exp_adjusted.rds"

# Create folder (if it doesn't exist) to save outputs (main ones)



#==============================================================================#
#                                                                              #
#                       INIT  WORKING ENVIRONMENT                              #
#                                                                              #
#==============================================================================#

# rm(list=ls())
#.rs.restartR()
# gc()

print(paste0('YOUR WORKING DIRECTORY IS:', getwd()))

# Root folder -> /PROJECTES/BISC_OMICS/analyses/BiSC_23/EWAS_diet_CL

# Load libraries 
library(httr)
library(rvest)
library(dplyr)
library(R.utils)
library(data.table)
library(GenomicRanges)
source("script/Define_functions.R")



#==============================================================================#
#                                                                              #
#                       Human annotation  (Hg19/GRCh37)                        #----
#                                                                              #
#==============================================================================#

# GRCh37, hg19 is a human genome reference assembly used for mapping and aligning DNA sequences.

cat("\n------------------  Hg19/GRCh37 annotation -----------------------------\n")

# Input folder 
#file <- "results/3_analysis_dmps/Bisc_14012025/nighttime_fasting_model_1/nighttime_fasting_Bisc_model_1_14012025_adjusted_1e-06.csv" # ok path now 

# Get annotation
if (array == "450K") {
  library('IlluminaHumanMethylation450kanno.ilmn12.hg19')
  annotation <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)}
if (array == "EPIC") {
  library('IlluminaHumanMethylationEPICanno.ilm10b2.hg19')
  annotation <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)}

annotation <- as.data.frame(annotation)
annotation$CpG <- row.names(annotation)

# Load the EWAS result and get df 
cat('Reading file', basename(file),'\n')
results <- readRDS(file)
# Add annotation to results by merging on row names = CpGs
results <- merge(results, annotation, by = 'CpG')

# Susbet to only BF CpGs 
sel_results <- subset(results, results$bf <= 0.05)
sel_results <- sel_results[, c('CpG', 'chr','pos', 'Relation_to_Island', 
                           'UCSC_RefGene_Name','UCSC_RefGene_Group', 'CpG_rs')]

# Split the UCSC_RefGene_Name column by ";", unlist the result, and take unique values
sel_results$gene_name <- sapply(strsplit(as.character(sel_results$UCSC_RefGene_Name), ";"),
                                function(x) paste(unique(x), collapse = ";"))

sel_results$gene_grup <- sapply(strsplit(as.character(sel_results$UCSC_RefGene_Group), ";"),
                                function(x) paste(unique(x), collapse = ";"))

# Remove 
sel_results <- sel_results[, !colnames(sel_results) %in% c('UCSC_RefGene_Name', 'UCSC_RefGene_Group')]

sel_results

#write.table(sel_results, "results/4_functional_analysis_annotation/Bisc_TEST/main_fdr_annot_pol.csv", sep=";")
save_to <- 'results/4_functional_analysis_annotation/Bisc_13062025/last_meal_sep_exp_bf_illumina.csv'
write.table(sel_results, save_to, sep = ";", row.names = FALSE)
  

#==============================================================================#
#                                                                              #
#                         Chromatin States                                     #----
#                                                                              #
#==============================================================================#

# Chromatin States are ...

# The DB has to be downloaded previously from the following github repository
# (https://github.com/isglobal-brge/easier/tree/main/data)
# Here I save it in "db" folder, and the Rscript exceuted from root, load the
# db from there. 

# Note: chroamtin states are differents for placenta than for blood
# Note: The Input  required is a df of EWAS result with the hg19 annotation.

cat("\n------------------ Chromatin States annotation  -----------------------\n")

# Load Fetal Placenta (FP) chromatin states annotation (15 states)
load('db/annotation_data/FP_15_E091.RData') # Object Name: FP_15_E091

# Convert sel_results dataframe to genomic ranges
results_GRange <- GRanges(seqnames = Rle(sel_results$chr),
                          ranges = IRanges(sel_results$pos, end = sel_results$pos),
                          name = sel_results$CpG,
                          chr = sel_results$chr,
                          pos = sel_results$pos)

names(results_GRange) <- results_GRange$name


# Search overlaps between CpGs (selected data) and Fetal Placenta (complete data)
over15 <- findoverlapvalues(selected = results_GRange,
                            complete = FP_15_E091)


# Create a dataframe with overlapping states and CpGs
data_chrstates <- c(mcols(over15$ranges), over15$values)
colnames(data_chrstates)[grep("States",colnames(data_chrstates))] <-  c("States15_FP")

# Add chromatin states with states to EWAS (CpG) data
sel_results <- merge(sel_results, data_chrstates[, c('name', 'States15_FP')], 
                     by.x = "CpG", by.y = "name")

#==============================================================================#
#                                                                              #
#                                 eQTMs                                        #----
#                                                                              #
#==============================================================================#

# Expression quantitative trait methylation (eQTMs), defined as associations 
# between DNA methylation levels and gene expression -> 
# In Children we have a database for eQTM in blood. In Placenta (no db) but there
# are 2 papers that report eQTMs. 

# Paper 1: Genetic variants influence on the placenta regulatory landscape 
# doi: 10.1371/journal.pgen.1007785
# To consider - Sample size is small. 


# Paper 2: In-depth characterization of the placental imprintome reveals novel 
# differentially methylated regions across birth weight categories
# doi: 10.1080/15592294.2019.1647945
# To consider - They look at imprinted regions 

#==============================================================================#
#                                                                              #
#                methylation quantitative trait loci (mQTLs)                   #----
#                                                                              #
#==============================================================================#

# mqTLS are SNPs associated with DNAm. The DB we are using is from BiSC cohort, 
# were they have expression data for placenta. We use the mqTLs all beta adjusted by CT and GA 
# Note: The Input required is a df of EWAS result with the hg19 annotation.

# cat("\n--------------------  mQTLs annotation  ---------------------------------\n")

# Copy file from: BISC_OMICS/data_final/resources/pla_fetal_cis_mQTL/
mQTLS_db <- fread("db/annotation_data/all_beta_a.mPCs_a.GA.txt", data.table = FALSE, dec=",")
mQTLS_db <- fread(paste0(mQTLS_path, mQTLS_file), data.table = FALSE, dec = ",")


# Add mQTLS data to EWAS results
mqtl_results <- merge(sel_results, mQTLS_db, all.x = TRUE, by.x = "CpG",
                      by.y ="phenotype_id")

# Aggregate/Group mQTLs by CpG
mqtl_results <- aggregate(mqtl ~ CpG, data = mqtl_results,
                                 FUN = function(x) paste(x, collapse = ";"))

# Add to sel_results 
sel_results <- merge(sel_results, mqtl_results, all.x = TRUE, by.x = "CpG")





#==============================================================================#
#                                                                              #
#                    Partially Methylated Domains (PMDs) PLACENTA              #----
#                                                                              #
#==============================================================================#

# PMDs are extended regions in the genome exhibiting a reduced average DNAm level

cat("\n----------  Partially Methylated Domains (PMDs) PLACENTA  ---------------\n")

# Load PMD annotation
load('db/annotation_data/PMD_placenta.RData') # Object Name: PMD_placenta

# Convert PMD_placenta df to Genom Ranges (getEnrichGenomicRanges does, I prefer wo function)
# Creates the name to assign to GRanges
pmd_name <- paste(PMD_placenta$Chr_PMD, 
                  PMD_placenta$Start_PMD, PMD_placenta$End_PMD, 
                  sep = ';')

# Create genomic ranges from PMD data
PMD_GRange <- GRanges(seqnames = Rle(PMD_placenta$Chr_PMD),
                      ranges = IRanges(PMD_placenta$Start_PMD, end = PMD_placenta$End_PMD),
                      name = pmd_name)
names(PMD_GRange) <- PMD_GRange$name

# Search overlaps between CpGs (selected data) and PMD (complete data)
overPMD <- findoverlapvalues(selected = results_GRange,
                             complete = PMD_GRange)

# Create a dataframe with overlapping states and CpGs
data_PMD  <- cbind.data.frame(DataFrame(CpG = results_GRange$name[overPMD$qhits]), 
                              DataFrame(PMD = PMD_GRange$name[overPMD$shits]))

# Add annotation PMD to results
sel_results <- merge(sel_results, data_PMD, by.x = "CpG", all=T)




#==============================================================================#
#                                                                              #
#            germ line differentially methylated regions (gDMRs)               #----
#                                                                              #
#==============================================================================#

# Germline Differentially Methylated Regions (gDMRs) are specific regions in the 
# genome where DNA methylation patterns are established in the germline 
# (sperm or egg cells) and are maintained throughout development and into the 
# somatic cells of an organism. These regions are crucial for regulating imprinted 
# genes, which are expressed in a parent-of-origin-specific manner.

# Convert sel_results dataframe to genomic ranges
results_GRange <- GRanges(seqnames = Rle(sel_results$chr),
                          ranges = IRanges(sel_results$pos, end = sel_results$pos),
                          name = sel_results$CpG,
                          chr = sel_results$chr,
                          pos = sel_results$pos)
names(results_GRange) <- results_GRange$name


cat("\n--------------------- gDMRs (imprinted regions)  ---------------------\n")

# Load gDMRs (impritned regions) annotation
load('db/annotation_data/IR_Placenta.RData') # Object Name: IR_Placenta

# Convert IR_Placenta df to Genomic Ranges 
ir_name <- paste(IR_Placenta$Chr_DMR, 
                 IR_Placenta$Start_DMR, IR_Placenta$End_DMR, 
                  sep = ';')

# Create genomic ranges from IR data
IR_GRange <- GRanges(seqnames = Rle(IR_Placenta$Chr_DMR),
                     ranges = IRanges(IR_Placenta$Start_DMR, end = IR_Placenta$End_DMR),
                     name = ir_name)
names(IR_GRange) <- IR_GRange$name

# Search overlaps between CpGs (selected data) and  DMR-IR (complete data)
overIR <- findoverlapvalues(selected = results_GRange,
                            complete = IR_GRange)

# Create a dataframe with overlapping Imprinted Regions and CpGs
data_IR <- cbind.data.frame(DataFrame(CpG = results_GRange$name[overIR$qhits]), 
                            DataFrame(gDMR = IR_GRange$name[overIR$shits]))

# Merge with re
sel_results <- merge(sel_results, data_IR, by.x = "CpG", all=T)

# Save 
save_to <- 'results/4_functional_analysis_annotation/Bisc_13062025/last_meal_sep_exp_bf_annot.csv'
write.table(sel_results, save_to, sep = ";", col.names = TRUE,  row.name = FALSE)  

#==============================================================================#
#                                                                              #
#                              EPIraction                                      #----
#                                                                              #
#==============================================================================#

# EPIRaction tool contains a catalogue of tissue specific gene-enhancer interactions 
# identified by applying two methods. 
# Method 1: For every broadly expressed gene within each tissue where this gene
# has expression above 1 TPM they predicted the candidate enhancer with conventional 
# ABC algorithm https://www.nature.com/articles/s41588-019-0538-0.
# Method 2:  For every tissue-specific gene within each tissue where this gene is overexpressed
# they applied logistic regression to predict candidate enhancers. 


cat("\n------------------  EPIraction annotation  -------------------------------\n")

# [DESCRIPTION] ----------------------------------------------------------------------------------
# EPIRation tool input expected is a TXT file with CpG positions (in hg38!) tab-separated.
# Output from EPIraction : Download Files and Save as txt into a folder, the path
# of such folder has to be specified at the section 0 under the variable annotation_data_folder. 

# Annotation --> Each tissue as a separate row	--> With QTL
# Placenta.Trophoblast

# Note: We only annotate the first 5000 CpGs when ordering the results by p value. 

# We need to first map our hg37 annotation to the hg38 annotation. To do so we'll
# use the Genome UCSS. This part has to be done MANUALLY! 

# THIS PART HAS TO BE DONE JUST ONCE FOR EACH TYPE OF ARRAY. To obtain the BED files
# UNCOMMENT the following lines. After that use the Browser as explained next and
# save the file obtained to the SAME FOLDER as input_folder. I also renamed them 
# as hg38_EPIC_annotation and hg38_450K_annotation. 

## [STEPS/CODE TO GET HG38 ANNOTATION BASED ON ILLUMINA HG19/37 ANNOTATION]
# # 1) Get annotation
# #meth_data <- readRDS(methylation_path) # 5001  565
# 
# if (array == "450K") {
#   library('IlluminaHumanMethylation450kanno.ilmn12.hg19')
#   annotation <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)}
# if (array == "EPIC") {
#   library('IlluminaHumanMethylationEPICanno.ilm10b2.hg19')
#   annotation <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)}
# 
# annotation <- as.data.frame(annotation)
# 
# # Transform annotation to BED format (e.g. "chr4 100000 100001", 0-based)
# annotation_bed <- data.frame(
#   chr = annotation$chr,
#   start = as.integer(annotation$pos) - 1,  # Convert to 0-based start
#   end = as.integer(annotation$pos),  # BED format uses 0-based start and 1-based end
#   CpG = annotation$Name
# )
# 
# # Ensure no scientific notation
# annotation_bed$start <- format(annotation_bed$start, scientific = FALSE)
# annotation_bed$end <- format(annotation_bed$end, scientific = FALSE)
# 
# # Write the BED file
# write.table(annotation_bed, paste0(annotation_data_folder, "/hg19_", array, "_annotation.bed"),
#             sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
# 
# 
# # 2) Go to genome ucsc and get DB to map  hg37 to hg38 (manually)
# # Web browser: https://genome.ucsc.edu/cgi-bin/hgLiftOver
# # Parameters specified:
# # - Original Genome: Human
# # - Original Assembly: Feb.2009 (GRCh37/hg19)
# # - New Genome: Human
# # - New Assembly: Dec.2013 (GRCh38/hg38)
# # - Minimum ratio of bases that must remap: 0.95 (default)
# # - BED 4 to BED 6 Options
#     # Allow multiple output regions:	off
#     # Minimum hit size in query:	0
#     # Minimum chain size in target:	0
# # - BED 12 Options
#     # Min ratio of alignment blocks or exons that must map:	1.00
#     # If thickStart/thickEnd is not mapped, use the closest mapped base:	off
# # - Data to transform: hg37_annotation.deb file
# 
# 
# # 3) Get DB mapped (Output) - EPIC
# # Successfully converted 866738 records: View Conversions (Click Download)
# # Conversion failed on 98 records.
# 
# # 3) Get DB mapped (Output) - 450K
# # Successfully converted 485441 records: View Conversions
# # Conversion failed on 71 records.  Display failure file Explain failure messages
# 
# # Reasons Conversions failed:
# # Deleted in new: Sequence intersects no chains
# # Partially deleted in new: Sequence insufficiently intersects one chain
# # Split in new: Sequence insufficiently intersects multiple chains
# # Duplicated in new: Sequence sufficiently intersects multiple chains
# # Boundary problem: Missing start or end base in an exon
#-----------------------------------------------------------------------------

# EPIRation tool input expected is a TXT file with CpG positions (in hg38!) tab-separated.
# frist 5000 cpgs ordered by p value 

# Load hg38 annotation 
if (array == "450K") {
  hg38_annotation <- read.table('db/annotation_data/hg38_450K_annotation.bed',
                                sep = "\t", header = FALSE,
                                stringsAsFactors = FALSE)}
if (array == "EPIC") {
  hg38_annotation <- read.table('db/annotation_data/hg38_EPIC_annotation.bed',
                                sep = "\t", header = FALSE,
                                stringsAsFactors = FALSE)}

# Rename colnames 
colnames(hg38_annotation) <- c("chr_hg38", "start_hg38", "end_hg38", "CpG")

# Add hg38 annotation to results by colnames = CpG
results <- merge(results, hg38_annotation, by = 'CpG', all.x = TRUE)  

# Rename columns as required by EPIaction (chr pos cpg).
colnames(results)[colnames(results) == "start_hg38"] <- "pos_hg38"
colnames(results)[colnames(results) == "CpG"] <- "CpGId"

# Select only those columns required by EPIraction  
# We'll annotate the first 5000 CpGs (ordered by p value)
results <- results[order(results$pvalue),]
results_epi <- results[1:5000, c("chr_hg38","pos_hg38","CpGId")]

# Save df 
save_to_txt <- 'results/4_functional_analysis_annotation/Bisc_13062025/last_meal_input_epiraction_hg38.txt' 
# ADD header fris tline staritng with #
custom_header <- "#chr\tpos\tcpg"

# Write the header manually
writeLines(custom_header, save_to_txt)
write.table(results_epi, save_to_txt, sep = "\t", col.names = FALSE,
            row.names = FALSE, quote = FALSE, append = TRUE)

# Compress to gzip 
save_to_gz <- 'results/4_functional_analysis_annotation/Bisc_13062025/last_meal_input_epiraction_hg38.gz' 
if (file.exists(save_to_gz)) {file.remove(save_to_gz)}
R.utils::gzip(save_to_txt, destname = save_to_gz, remove = FALSE, overwrite = TRUE)


# Get annotation from Epiraction 
epiraction_annot_db <- get_annot_EPIraction(gzip_file = save_to_gz,
                                            out_folder = 'results/4_functional_analysis_annotation/Bisc_13062025/')


# Load the epiraction annotation 
epiraction_annot <- fread(epiraction_annot_db, data.table = FALSE)
save_to <- "results/4_functional_analysis_annotation/Bisc_13062025/lastmeal_top5000_alltissues_EPIraction.csv"
write.table(epiraction_annot, save_to, sep = ";", col.names = TRUE, row.name = FALSE)  

head(epiraction_annot)

# Select Placenta Tissue and add to results 
sel_epiraction_annot <- epiraction_annot[which(epiraction_annot$tissue == "Placenta.Trophoblast"), ]
colnames(sel_epiraction_annot)[colnames(sel_epiraction_annot) == "annotation"] <- 'CpG'

# Add EPIraction annotation to results by colnames = CpG 
# Note: We cna have multiple rows for CpG
sel_results_epi <- merge(sel_results, sel_epiraction_annot, by = 'CpG', all.x = TRUE) 
head(sel_results_epi)

# Aggregate/Group  symbol/epiraction by CpG 
sel_results_epi_by_cpg <- aggregate(symbol ~ CpG, data = sel_results_epi, 
                                FUN = function(x) paste(x, collapse = ";"))

sel_results <- merge(sel_results, sel_results_epi_by_cpg, by = 'CpG', all.x = TRUE) 


# Save df 
save_to <- 'results/4_functional_analysis_annotation/Bisc_13062025/last_meal_sep_exp_bf_annot.csv'
write.table(sel_results, save_to, sep = ";", col.names = TRUE,  row.name = FALSE)  

#==============================================================================#-----



