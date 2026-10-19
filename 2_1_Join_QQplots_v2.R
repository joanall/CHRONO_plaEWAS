


#==============================================================================#
#                       INIT  WORKING ENVIRONMENT                              #
#==============================================================================#

library(ggplot2)
library(png)
library(grid)
library(patchwork) # For arranging multiple plots easily
library(cowplot)



#==============================================================================#
#                          DEFINE PATHS                                        #
#==============================================================================#


files <- c(
  
  # Night Fasting
  "results/2_quality_control/Bisc_14012025/nighttime_fasting_model_1_sep_exp/QQplot_lambda.rds",
  "results/2_quality_control/Bisc_14012025/nighttime_fasting_model_2_sep_exp/QQplot_lambda.rds",
  "results/2_quality_control/Bisc_05022025/nighttime_fasting_model_3_sep_exp/QQplot_lambda.rds",
  
  # Night Fasting Adjusted 
  "results/2_quality_control/Bisc_14012025/nighttime_fasting_model_1/QQplot_lambda.rds",
  "results/2_quality_control/Bisc_14012025/nighttime_fasting_model_2/QQplot_lambda.rds",
  "results/2_quality_control/Bisc_05022025/nighttime_fasting_model_3/QQplot_lambda.rds",
  
  # First Meal 
  "results/2_quality_control/Bisc_14012025/first_meal_model_1_sep_exp/QQplot_lambda.rds",
  "results/2_quality_control/Bisc_14012025/first_meal_model_2_sep_exp/QQplot_lambda.rds",
  "results/2_quality_control/Bisc_05022025/first_meal_model_3_sep_exp/QQplot_lambda.rds",
  
  # First Meal Adjusted 
  "results/2_quality_control/Bisc_14012025/first_meal_model_1_firstlast/QQplot_lambda.rds",
  "results/2_quality_control/Bisc_14012025/first_meal_model_2_firstlast/QQplot_lambda.rds",
  "results/2_quality_control/Bisc_05022025/first_meal_model_3_firstlast/QQplot_lambda.rds",
  
  # Eating occasions 
  "results/2_quality_control/Bisc_14012025/eating_occasions_model_1_sep_exp/QQplot_lambda.rds",
  "results/2_quality_control/Bisc_14012025/eating_occasions_model_2_sep_exp/QQplot_lambda.rds",
  "results/2_quality_control/Bisc_05022025/eating_occasions_model_3_sep_exp/QQplot_lambda.rds",
  
  # Eating occasions adjusted 
  "results/2_quality_control/Bisc_14012025/eating_occasions_model_1_firstlast/QQplot_lambda.rds",
  "results/2_quality_control/Bisc_14012025/eating_occasions_model_2_firstlast/QQplot_lambda.rds",
  "results/2_quality_control/Bisc_05022025/eating_occasions_model_3_firstlast/QQplot_lambda.rds",
  
  # Last Meal 
  "results/2_quality_control/Bisc_14012025/last_meal_model_1_sep_exp/QQplot_lambda.rds",
  "results/2_quality_control/Bisc_14012025/last_meal_model_2_sep_exp/QQplot_lambda.rds",
  "results/2_quality_control/Bisc_05022025/last_meal_model_3_sep_exp/QQplot_lambda.rds",
  
  # Last Meal adjusted 
  "results/2_quality_control/Bisc_14012025/last_meal_model_1_firstlast/QQplot_lambda.rds",
  "results/2_quality_control/Bisc_14012025/last_meal_model_2_firstlast/QQplot_lambda.rds",
  "results/2_quality_control/Bisc_05022025/last_meal_model_3_firstlast/QQplot_lambda.rds",
  
  # Eating Jetlag 
  "results/2_quality_control/Bisc_14012025/eating_jetlag_model_1/QQplot_lambda.rds",
  "results/2_quality_control/Bisc_14012025/eating_jetlag_model_2/QQplot_lambda.rds",
  "results/2_quality_control/Bisc_05022025/eating_jetlag_model_3/QQplot_lambda.rds"
)


# LIST QQ PLOTS 
# 3 columnes M1 / M2 / M3  
# each row one epxosure = 9 rows 
# exposures: NT / NT ad / FM / FM adj / EAT /EAT adj / LAST / LAST adj / JETLag



#==============================================================================#
#                       DEFINE LABELS AND TITLES                               #
#==============================================================================#

# Exposure row labels (9 rows, each with 3 models)
row_exposures <- c(
  'Nighttime fasting', 'Mutually adjusted',
  'First meal', 'Mutually adjusted',
  'Eating occasions', 'Mutually adjusted',
  'Last meal', 'Mutually adjusted',
  'Eating jetlag'
)

# Column headers
column_titles <- c("Model 1", "Model 2", "Model 3")

#==============================================================================#
#                      LOAD GGPLOT OBJECTS                                     #
#==============================================================================#

# Load all .rds ggplot objects
plot_list <- lapply(files, readRDS)


plot_list <- lapply(plot_list, function(p) {
  for (i in seq_along(p$layers)) {
    layer <- p$layers[[i]]
    
    if (inherits(layer$geom, "GeomText")) {
      # Shrink lambda text
      p$layers[[i]]$aes_params$size <- 3
      
      # Optionally: shift slightly to the right and down if x=0 and y=max
      if (!is.null(layer$mapping$x)) {
        if (layer$mapping$x == 0) {
          layer$mapping$x <- 0.1
        }
      }
    }
  }
  p
})


#==============================================================================#
#                      ARRANGE BY ROWS                                         #
#==============================================================================#

# Group plots into rows of 3
rows <- split(plot_list, ceiling(seq_along(plot_list) / 3))

# Add exposure labels on the left of each row
rows_with_titles <- mapply(function(row, title) {
  plot_grid(
    # Left label
    ggdraw() + 
      draw_label(title, x = 0, hjust = 0, fontface = "bold", size = 12),
    
    # Row of 3 QQ plots
    plot_grid(plotlist = row, ncol = 3),
    
    ncol = 1,
    rel_heights = c(0.12, 1)
  )
}, rows, row_exposures, SIMPLIFY = FALSE)


#==============================================================================#
#                      COMBINE FINAL PLOT                                      #
#==============================================================================#

# Column headers
column_titles <- c("Model 1", "Model 2", "Model 3")

# Add a top row with column headers (above the full grid)
top_labels <- wrap_plots(
  lapply(column_titles, function(title) {
    ggplot() +
      theme_void() +
      ggtitle(title) +
      theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))
  }),
  ncol = 3
)

# Combine column headers and plot rows
final_plot <- wrap_plots(
  c(list(top_labels), rows_with_titles),
  ncol = 1
)

# Save output
ggsave("results/2_quality_control/combined_QQplots3.png",
       plot = final_plot, width = 10, height = 18, dpi = 300)


stop()

################################################################################
################## BACON 
###############################

stop()


# JOIN QQ PLOTS BACON CORRECTED 
files_bacon <- c(
  "results/2_quality_control/Bisc_05022025/first_meal_model_3_sep_exp/QQplot_corrected_lambda.png",
  
  "results/2_quality_control/Bisc_14012025/first_meal_model_1_firstlast/QQplot_corrected_lambda.png",
  "results/2_quality_control/Bisc_14012025/first_meal_model_2_firstlast/QQplot_corrected_lambda.png",
  "results/2_quality_control/Bisc_05022025/first_meal_model_3_firstlast/QQplot_corrected_lambda.png",
  
  # Eat occasions
  "results/2_quality_control/Bisc_14012025/eating_occasions_model_1_sep_exp/QQplot_corrected_lambda.png",
  "results/2_quality_control/Bisc_14012025/eating_occasions_model_2_sep_exp/QQplot_corrected_lambda.png",
  "results/2_quality_control/Bisc_05022025/eating_occasions_model_3_sep_exp/QQplot_corrected_lambda.png",
  
  # Eating occasions adjsuted 
  "results/2_quality_control/Bisc_14012025/eating_occasions_model_1_firstlast/QQplot_corrected_lambda.png",
  "results/2_quality_control/Bisc_14012025/eating_occasions_model_2_firstlast/QQplot_corrected_lambda.png",
  "results/2_quality_control/Bisc_05022025/eating_occasions_model_3_firstlast/QQplot_corrected_lambda.png",
)

