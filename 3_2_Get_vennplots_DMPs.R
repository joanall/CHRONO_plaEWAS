# Author:  Camille&Joana
# Title: Summary Tables ()
# Date: July 2024


# INPUT: RDS or CSV with Processed EWAS results
# OUTPUT: Summary tables grouped by exposure and/or model 


# Note: Not prepared to run all together, go to each chunk and run the summ tables 
# or the venn diagram and change the exposure and model name 
#==============================================================================#
#                                                                              #
#                       INIT  WORKING ENVIRONMENT                              #
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
library(VennDiagram)
library(grid)
#------------------------------------------------------------------------------#
#                           Input from User                                    #----
#------------------------------------------------------------------------------#

date <- "16072024"
array <- 'EPIC'
cohort_name <- 'Bisc'
suggestive_cutoff <- 1E-06

# Input folder 
input_folder <- "results/EWAS_MEALTIME/3_QUALITY_CONTROL/"

# Output folder 
results_path <- "results/EWAS_MEALTIME/4_ANALYSIS_DMPS/"

#------------------------------------------------------------------------------#
#                            Output folder                                     #----
#------------------------------------------------------------------------------#


# Define folder to store the new processed files 
results_name <- paste0(cohort_name, "_",date, "_SumTables/")
results_folder <- paste0(results_path, results_name)
cat("OUTPUT FOLDER:\n",results_folder, "\n")

# Create folder (if it doesn't exist) to save outputs (main ones)
if (!dir.exists(results_folder)) {dir.create(results_folder)}


################################################################################
#                                                                              #
#                                 BY EXPOSURE                                  #
#                       (Summary tables + Venn Diagrams)                       #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
#                               Main Model  (1st meal)                         #----
#------------------------------------------------------------------------------#

# [CHOOSE EXPOSURE]: night_fasting , first_meal, eating_occasions, eating_jetlag
exposure_name <- "eating_occasions"

# Subfolders names main models 
input_files <- c(paste0(input_folder, "Bisc_", date, "_QC/", exposure_name, 
                        "_model_1/Bisc_", date, "_", exposure_name, "_model_1_adjusted.rds"),
                 
                 paste0(input_folder, "Bisc_", date, "_QC/", exposure_name, 
                        "_model_2/Bisc_", date, "_", exposure_name, "_model_2_adjusted.rds"),
                 
                 paste0(input_folder, "Bisc_", date, "_QC/", exposure_name, 
                        "_model_3/Bisc_", date, "_", exposure_name, "_model_3_adjusted.rds"))



# Read all files into a list
list_ewas_results <- lapply(input_files, function(file) readRDS(file))

# Extract top CpGs for each file base on suggestive cutoff
top_cpgs <- lapply(list_ewas_results, 
                   function(ewas) {subset(ewas, pval <= suggestive_cutoff)$CpG})

# Make venn diagram 
sink("/dev/null")


save_to_counts <- paste0(results_path, "/Bisc_", date, "_VennPlot/vennplot_", 
                         exposure_name, "_mainmodels_counts.png")

# Create the Venn diagram with counts only
venn_counts <- venn.diagram(
  x = top_cpgs,
  category.names = c("Model1", "Model2", "Model3"),
  filename = NULL,  # Set to NULL to add further customization before saving
  output = TRUE,
  
  # Appearance customization
  fill = c("#D9D9D9", "#B0B0B0", "#808080"),  # Grayscale colors
  alpha = 0.5,
  cex = 1.5,  # Size for the counts
  fontface = "bold",
  fontfamily = "serif",
  cat.cex = 1.5,  # Size for category labels
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-25, 25, 0),  # Adjust positions of the category names
  cat.dist = c(0.025, 0.025, 0.025),  # Bring category labels closer to the circles
  cat.fontfamily = "serif",
  lwd = 2,  # Line width of the circles
  lty = 'solid',  # Line type of the circles
  col = c("#404040", "#606060", "#808080"),  # Grayscale border colors
  
  # Display counts only
  print.mode = c("raw"),  # Show the counts (raw numbers) in the intersections
  #print.mode = c("percent", "raw"),  # Show the counts (raw numbers) in the intersections
  label.col = "black",
  label.fontface = "bold",
  label.cex = 1.2,
  label.just = "center",
  
  # Other customizations
  margin = 0.1  # Adjust margin to fit long text if necessary
)

# Save Diagram 1 (Counts Only)
png(save_to_counts, width = 800, height = 800)
grid.newpage()
grid.draw(venn_counts)
dev.off()
sink()


# Create the Venn diagram with CpGs ID --> needs improvement 
sink("/dev/null")

save_to_cpgs <- paste0(results_path, "/Bisc_", date, "_VennPlot/vennplot_", 
                       exposure_name, "_mainmodels_cpgs.png")

venn_cpgs <- venn.diagram(
  x = top_cpgs,
  category.names = c("Model1", "Model2", "Model3"),
  filename = NULL,  # Set to NULL to add further customization before saving
  output = TRUE,
  
  # Appearance customization (no counts, just circles)
  fill = c("#D9D9D9", "#B0B0B0", "#808080"),  # Grayscale colors
  alpha = 0.5,
  cex = 0,  # Set cex to 0 to hide intersection counts
  fontface = "bold",
  fontfamily = "serif",
  cat.cex = 1.5,  # Size for category labels
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-25, 25, 0),  # Adjust positions of the category names
  cat.dist = c(0.025, 0.025, 0.025),  # Bring category labels closer to the circles
  cat.fontfamily = "serif",
  lwd = 2,  # Line width of the circles
  lty = 'solid',  # Line type of the circles
  col = c("#404040", "#606060", "#808080"),  # Grayscale border colors
  
  # Other customizations
  margin = 0.1  # Adjust margin to fit long text if necessary
)

# Draw and save the Venn diagram with CpG IDs only
png(save_to_cpgs, width = 800, height = 800)
grid.newpage()
grid.draw(venn_cpgs)


# Extract CpG IDs for each intersection manually, ensuring uniqueness in display
unique_inters12 <- setdiff(intersect(top_cpgs[[1]], top_cpgs[[2]]), top_cpgs[[3]])
unique_inters13 <- setdiff(intersect(top_cpgs[[1]], top_cpgs[[3]]), top_cpgs[[2]])
unique_inters23 <- setdiff(intersect(top_cpgs[[2]], top_cpgs[[3]]), top_cpgs[[1]])
unique_inters123 <- Reduce(intersect, top_cpgs)

# Add the CpG IDs manually at appropriate locations, only once per intersection
grid.text(paste(unique_inters123, collapse = "\n"), x = 0.5, y = 0.5, gp=gpar(cex=1.2))  # Center
grid.text(paste(unique_inters12, collapse = "\n"), x = 0.35, y = 0.7, gp=gpar(cex=1.2))  # Top left
grid.text(paste(unique_inters13, collapse = "\n"), x = 0.65, y = 0.7, gp=gpar(cex=1.2))  # Top right
grid.text(paste(unique_inters23, collapse = "\n"), x = 0.5, y = 0.8, gp=gpar(cex=1.2))  # Top center


dev.off()  # Finish and save the PNG file
sink()



#------------------------------------------------------------------------------#
#                             Sens   Model                                     #----
#------------------------------------------------------------------------------#

# [CHOOSE EXPOSURE]: night_fasting , first_meal, eating_occasions, eating_jetlag
exposure_name <- "night_fasting"

input_files <- c(paste0(input_folder, "Bisc_", date, "_QC/", exposure_name, 
                        "_model_noGDM/Bisc_", date, "_", exposure_name, "_model_noGDM_adjusted.rds"),
                 
                 paste0(input_folder, "Bisc_", date, "_QC/", exposure_name, 
                        "_model_ethnicity/Bisc_", date, "_", exposure_name, "_model_ethnicity_adjusted.rds"),
                 
                 paste0(input_folder, "Bisc_", date, "_QC/", exposure_name, 
                        "_model_noPC/Bisc_", date, "_", exposure_name, "_model_noPC_adjusted.rds"),
                 
                 paste0(input_folder, "Bisc_", date, "_QC/", exposure_name, 
                        "_model_noconf/Bisc_", date, "_", exposure_name, "_model_noconf_adjusted.rds"),
                 
                 paste0(input_folder, "Bisc_", date, "_QC/", exposure_name, 
                        "_model_EarlyEat/Bisc_", date, "_", exposure_name, "_model_EarlyEat_adjusted.rds"))


# Read all files into a list
list_ewas_results <- lapply(input_files, function(file) readRDS(file))

# Extract top CpGs for each file base on suggestive cutoff
top_cpgs <- lapply(list_ewas_results, 
                   function(ewas) {subset(ewas, pval <= suggestive_cutoff)$CpG})

# PLOT 



#------------------------------------------------------------------------------#
#                             Separate Mealtime Exposure                       #----
#------------------------------------------------------------------------------#

# [CHOOSE EXPOSURE]: night_fasting , first_meal, eating_occasions, last_meal
exposure_name <- "last_meal"


input_file <- c(paste0(input_folder, "Bisc_", date, "_QC/", exposure_name, 
                       "_model_1/Bisc_", date, "_", exposure_name, "_model_1_SepExp_adjusted.rds"),
                
                paste0(input_folder, "Bisc_", date, "_QC/", exposure_name, 
                       "_model_1/Bisc_", date, "_", exposure_name, "_model_2_SepExp_adjusted.rds"),
                
                paste0(input_folder, "Bisc_", date, "_QC/", exposure_name, 
                       "_model_1/Bisc_", date, "_", exposure_name, "_model_3_SepEXp_adjusted.rds"))


# Read all files into a list
list_ewas_results <- lapply(input_files, function(file) readRDS(file))

# Extract top CpGs for each file base on suggestive cutoff
top_cpgs <- lapply(list_ewas_results, 
                   function(ewas) {subset(ewas, pval <= suggestive_cutoff)$CpG})

# PLOTS 



#------------------------------------------------------------------------------#
#                               Main Model  (last meal)                        #----
#------------------------------------------------------------------------------#


# [CHOOSE EXPOSURE]: night_fasting,  eating_occasions, last_meal
exposure_name <- "last_meal"


# Subfolders names main models 
input_files <- c(paste0(input_folder, "Bisc_", date, "_QC/", exposure_name, 
                        "_model_1/Bisc_", date, "_", exposure_name, "_model_1_lastmeal_adjusted.rds"),
                 
                 paste0(input_folder, "Bisc_", date, "_QC/", exposure_name, 
                        "_model_2/Bisc_", date, "_", exposure_name, "_model_2_lastmeal_adjusted.rds"),
                 
                 paste0(input_folder, "Bisc_", date, "_QC/", exposure_name, 
                        "_model_3/Bisc_", date, "_", exposure_name, "_model_3_lastmeal_adjusted.rds"))



# Read all files into a list
list_ewas_results <- lapply(input_files, function(file) readRDS(file))

# Extract top CpGs for each file base on suggestive cutoff
top_cpgs <- lapply(list_ewas_results, 
                   function(ewas) {subset(ewas, pval <= suggestive_cutoff)$CpG})

# PLOT 


################################################################################
#                                                                              #
#                          VENN DIAGRAMS BY EXPOSURE                           #
#                            (to do)                                           #
################################################################################
# We can have the following venn diagrams: # 


# Options: 
# - main models -> model_1, model_2, model_3
# - sensitivty model -> model_ethnicity, model_noPC, model_noconf, model_EarlyEat
# - Sep ExP -> model_1_SepExp, model_2_SepExp, model_3_SepExp
# - main models last meal -> model_1_lastmeal, model_2_lastmeal, model_3_lastmeal


if (exposure_name != "") {
  
  #------------------------------------------------------------------------------#
  #             CpGs Overlap for Between Exposures for each model                #----
  #------------------------------------------------------------------------------#
  
  if (which_model != "[]_lastmeal") {
    # Subfolders names main models 
    input_files <- c(paste0(input_folder, "Bisc_", date, "_QC/night_fasting_", which_model, "/Bisc_", 
                            date, "_night_fasting_", which_model, "_adjusted.rds"),
                     paste0(input_folder, "Bisc_", date, "_QC/first_meal_", which_model, "/Bisc_", 
                            date, "_night_fasting_", which_model_1, "_adjusted.rds"),
                     paste0(input_folder, "Bisc_", date, "_QC/eating_occasions_", which_model, "/Bisc_", 
                            date, "_night_fasting_", which_model, "_adjusted.rds"),
                     paste0(input_folder, "Bisc_", date, "_QC/eating_jetlag_", which_model, "/Bisc_", 
                            date, "_night_fasting_", which_model, "_adjusted.rds"))
    
    
    
    
    # Save to CSV file 
    save_to <- paste0(results_path, "/Bisc_", date, "_SumTables/vennplot_", 
                      exposure_name, "_allexp_mainmodels.png")
    
    # Read all files into a list
    list_ewas_results <- lapply(input_files, function(file) readRDS(file))
    
    # Extract top CpGs for each file base on suggestive cutoff
    top_cpgs <- lapply(list_ewas_results, 
                       function(ewas) {subset(ewas, pval <= suggestive_cutoff)$CpG})
    
    
    
    venn.diagram(
      xx = list(Night_Fasting = top_cpgs[1],
                Eating_Occasions = dmrs_eating_occasions,
                First_Meal = dmrs_time_first_meal,
                Eating_Jetlag = dmrs_eating_jetlag),
      category.names = c("Nighttime Fasting", "Eating Occasions", 
                         "Time First Meal", "Eating Jelag"),
      filename = file.path(output_dir, "venn_nighttime_fasting.png"),
      output = TRUE)
    
    
    #------------------------------------------------------------------------------#
    #             CpGs Overalp for Between Exposures for each model                #----
    #------------------------------------------------------------------------------#
    
    
    # Define the DMRs for each mealtime exposure and model
    nighttime_fasting_model1 <- c("P4HTM", "RNF5P1;RNF5;AGPAT1", "C3orf55", "TICAM1", "EDNRB")
    exp_model1 <- c("RNF5P1;RNF5;AGPAT1", "P4HTM")
    exp_model1 <- c()
    expo_model1 <- c() 
    
    # Nighttime Fasting
    venn.diagram(
      x = list(Model1 = nighttime_fasting_model1, Model2 = nighttime_fasting_model2),
      category.names = c("Model1", "Model2"),
      filename = file.path(output_dir, "venn_nighttime_fasting.png"),
      output = TRUE
    )
    
  }
  #------------------------------------------------------------------------------#
  #             CpGs Overalp for Between Exposures for all models                #----
  #------------------------------------------------------------------------------#
  
  exp1_cpgs <- c()
  exp2_cpgs <- c()
  exp3_cpgs <- c()
  
  venn.diagram(
    x = list(Model1 = nighttime_fasting_model1, Model2 = nighttime_fasting_model2),
    category.names = c("Nighttime Fasting", "Eating Occasions", "Time First Meal", "Eating Jetlag"),
    filename = file.path(output_dir, "venn_nighttime_fasting.png"),
    output = TRUE
  )
  
  
}


#==============================================================================#
#                                                                              #
#                                 BY MODEL                                     #
#                             (Venn diagrams)                                  #
#                                                                              #
#==============================================================================#











