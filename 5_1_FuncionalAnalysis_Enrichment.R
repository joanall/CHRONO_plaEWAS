#==============================================================================#
#                                                                              #
#               MATERNAL DIET MEALTIME AND PLACENTA DNA METHYLATION            #
#                             Functional Enrichment                            #
#                                                                              #                                  
#                   Camille Lassale & Joana Llauradó Pont                      #         
#                               January 2025                                   #
#                                                                              #
#==============================================================================#

# [Description]: Series of enrichment to see biological pathways relaetd to the 
# genes associated with mealtime and placental DNAm. 

# A) GO & KEGG collections: Gene Ontology defines concepts/classes used to describe 
# gene function, and relationships between these concepts. 
# 

# B) Consenus Path DB: Integrates interaction networks in Homo sapiens. Not EWAS specific.

# Overrepresentation analysis requires a list of gene/protein, metabolite and 
# is computed with Fisher’s exact test . (We don't do enrichment, it  requires 
# numerical data from two different states for comparison)

# Pathway Based Sets: These are collections of biological pathways representing a series of 
# actions among molecules in a cell that lead to a certain product or change in the cell.
# Included databases: Manual upload, Wikipathways, Biocarta, Smpdb, Ehmn, Kegg, Humancyc, 
# Reactome, Inoh, Pharmgkb, Netpath, Pid, Signalink

# Protein  Complex-based sets:  are sets of genes whose protein products are members 
# of the same annotated protein complex . 

# Ref: doi: 10.1093/nar/gkab1128


# C) Molecular Signatures Database enrichment - -> Collection
# The nerichment is pefromed with gsa function from missmehtyl using the colleciotn
# Hs.c2 from molecular signatu DB

# [INPUT]: RDS File (list of files) storing a dataframe where rows are CpGs and columns are statsistics 

# [OUTPUT]: 
# - 1 CSV with CpG-Genes included in the Functional Enrichment Analysis
# - 1 CSV with enriched GO terms
# - 1 CSV with enriched KEGG pathways 
# - 1 CSV With enriched pathways from Consensus Path DB -> manually 
# - 1 CSV with enriched pathways/terms from Molecular Signature DB
#

#==============================================================================#
#                                                                              #
#                             USER INPUT                                       #
#                                                                              #
#==============================================================================#


# Output folder 
results_path <- "results/5_functional_analysis_enrichment/Last_Meal_13062025"
if (!dir.exists(results_path)) {dir.create(results_path, recursive = TRUE)}

# Cutoff Enrichemnts (pval)
cutoff_enrich <- 1E-04
#input_cpgs_toenrich <- "results/5_functional_analysis_enrichment/Bisc_05022025/mealtime_cpgs_toenrich_model_3_1E06.rds" #hg19
input_cpgs_toenrich <- "results/5_functional_analysis_enrichment/Bisc_13062025/last_meal_cpgs_toenrich_model_1_1E04.rds"

#input_allcpgs <- "results/5_functional_analysis_enrichment/Bisc_05022025/mealtime_cpgs_all_model_3_1E06.rds"
input_allcpgs <- "results/5_functional_analysis_enrichment/Bisc_13062025/last_meal_cpgs_all_model_1_1E05.rds"


#==============================================================================#
#                                                                              #
#                       INIT  WORKING ENVIRONMENT                              #
#                                                                              #
#==============================================================================#

# rm(list=ls())
#.rs.restartR()
# gc()

print(paste0('YOUR WORKING DIRECTORY IS:', getwd()))

# Import libraries
library(missMethyl)
#install.packages("devtools")
#library(devtools)
#install_github("isglobal-brge/brgeEnrich")

#library(brgeEnrich) # Parse web Consesusn Path
library(msigdbr) # Load collection Molec.Sign DB
library(dplyr)




#==============================================================================#
#                         GET GENE SET FOR ENRICHMENT                          #
#==============================================================================#

cat("\n------------- Get CpG/Genes inlcuded in the enrichment -------------------\n")


results <- readRDS(input_cpgs_toenrich) 

enrich_cpgs <- results[results$pvalue <= cutoff_enrich, c('CpG','UCSC_RefGene_Name')]


save_to <- paste0(results_path, 'last_meal_genes_enrichment_model_1_1E04.csv')
write.table(enrich_cpgs, save_to, sep=";")

# Get all the CpG sites used in the analysis to form the background (gene universe)
allcpgs <-  readRDS(input_allcpgs) # Alternative: allcpgs <- rownames(annotation)

#==============================================================================#
#                        GO & KEGG collections (missMethyl)                    #----
#==============================================================================#

cat("\n---------------------  GO & KEGG collections ----------------------------\n")


# Run GO enrichment

enrich <- gometh(sig.cpg = enrich_cpgs$CpG, all.cpg = allcpgs$CpG, collection = "GO",
                 plot.bias = TRUE, prior.prob = TRUE)

save_to <- file.path(results_path,'last_meal_go_enrichment_model_1_1E04.csv')
write.table(enrich, save_to, sep=";")

# Run KEGG enrichment
enrich <- gometh(sig.cpg = enrich_cpgs$CpG, all.cpg = allcpgs$CpG, collection = "KEGG",
                 plot.bias = FALSE, prior.prob = TRUE)

save_to <- file.path(results_path,'last_meal_kegg_enrichment_model_1_1E04.csv')
write.table(enrich, save_to, sep=";")

#==============================================================================#
#                     Molecular Signatures Database enrichment                 #----
#==============================================================================#


# Collection used: MSigDB c2 - v7p1 
# 

cat("\n----------------- Molecular Signatures Database -------------------------\n")

# Load the 'C2' collection for Homo sapiens
Hs.c2 <- msigdbr(species = "Homo sapiens", category = "C2")

# Transform the data frame to a list format as required by gsameth
Hs.c2_list <- split(Hs.c2$entrez_gene, Hs.c2$gs_name)

msd_enrich <- gsameth(sig.cpg = enrich_cpgs$CpG,
                      all.cpg = allcpgs$CpG,
                      collection = Hs.c2_list,
                      array.type = 'EPIC')

msd_enrich
save_to <- file.path(results_path,'last_meal_MSigDB_enrichment_model_1_1E04.csv')
write.table(msd_enrich, save_to, sep=";")



#==============================================================================#
#                            Consensus Path DB                                 #----
#==============================================================================#

# Note: Get genes and parse manually, the original function to parse the webpage
# not working, i think the page was updated. 


cat("\n-------------------- Consensus Path DB -----------------------------\n")

# Available FSet types :
# 1 P     manually curated pathways from pathway databases
# 7 C     protein complex-based sets

# Overlaps Set
# A)  Pathway-based sets -> pathways as defined by pathway databases	
# Databases: Manual upload, Wikipathways, Biocarta, Smpdb, Ehmn, Kegg (unclick),
#                    Humancyc, Reactome, Inoh, Pharmgkb, Netpath, Pid, Signalink
# Parameters: minimum overlap with input list = 2 and p-value cutoff = 0.01

# B)  Protein complex-based gene sets -> sets of genes whose products are found together in protein complexes 
# Parameters: minimum complex size = 2 
#             minimum overlap with input list	= 2
#             p-value cutoff = 0.01


# C)  Putative enhancer target sets -> sets of genes predicted to be regulated by the same enhancer	
# - database: EnhancerAtlas v2.0
# Parameters :  minimum overlap with input list = 2 and p-value cutoff = 0.01


# Select Genes to be enriched 
enrich_genes <- results[results$pvalue <= cutoff_enrich, "UCSC_RefGene_Name"]
enrich_genes <- unique(unlist(sapply(enrich_genes, function(x) strsplit(x, ";"))))
enrich_genes <- na.omit(as.character(enrich_genes))
enrich_genes <- unlist(lapply(enrich_genes, function(gene_name) paste0(gene_name, '_HUMAN')))

save_to <- file.path(results_path, 'mealtime_genes_uniprot_consensus_model_1_1E04.txt')
writeLines(enrich_genes, save_to)


# ORA for manually curated pathways from pathway databases (P)  
# CPDB_enrich <- cpdbOverRepresentationAnalysis(entityType = 'genes', 
#                                               fsetType = 'P',
#                                               accNumbers = enrich_genes, 
#                                               accType = 'uniprot')













