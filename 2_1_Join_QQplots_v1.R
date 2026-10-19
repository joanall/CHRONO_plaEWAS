#==============================================================================#
#                                                                              #
#               MATERNAL DIET MEALTIME AND PLACENTA DNA METHYLATION            #
#                             Quality Control                                  #
#                                                                              #                                  
#                   Camille Lassale & Joana Llauradó Pont                      #         
#                               Janaury 2025                                   #
#                                                                              #
#==============================================================================#

# DESCRIPTION: Read PNG pltos with QQ plots and join them in one 


# [INPUT]: qqplots PNG indiviually 

# [OUTPUT]: QQ plots all together

#==============================================================================#
#                       INIT  WORKING ENVIRONMENT                              #
#==============================================================================#

library(ggplot2)
library(png)
library(grid)
library(patchwork) # For arranging multiple plots easily


#==============================================================================#
#                          DEFINE PATHS                                        #
#==============================================================================#

main_path <- c('results/2_quality_control/Bisc_05022025/nighttime_fasting_model_3_sep_exp/QQplot_lambda.png',
               'results/2_quality_control/Bisc_05022025/nighttime_fasting_model_3/QQplot_lambda.png',  # adjsuted first meal
               'results/2_quality_control/Bisc_05022025/first_meal_model_3_sep_exp/QQplot_lambda.png', 
               'results/2_quality_control/Bisc_05022025/first_meal_model_3_firstlast/QQplot_lambda.png',
               'results/2_quality_control/Bisc_05022025/last_meal_model_3_sep_exp/QQplot_lambda.png',
               'results/2_quality_control/Bisc_05022025/last_meal_model_3_firstlast/QQplot_lambda.png',
               'results/2_quality_control/Bisc_05022025/eating_occasions_model_3_sep_exp/QQplot_lambda.png',
               'results/2_quality_control/Bisc_05022025/eating_occasions_model_3_firstlast/QQplot_lambda.png',
               'results/2_quality_control/Bisc_05022025/eating_jetlag_model_3/QQplot_lambda.png')

#==============================================================================#
#                               JOIN PLOTS                                     #
#==============================================================================#

title <- 'MODEL 3'
sub_titles <- c('night time fasting', 'night time fasting adj. first meal',
                'first meal', 'first meal adj. last meal', 
                'last meal', 'last meal adj. first meal', 
                'eating occasions', 'eating occasions adj. by first and last meal',
                'eating jetlag')


# Initialize an empty list to store plots
plot_list <- list()

# Loop through file paths and create plots
for (i in seq_along(main_path)) {
  # Read image
  img <- readPNG(main_path[i])
  # Create a ggplot object with the image
  p <- ggplot() +
    annotation_raster(img, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
    ggtitle(sub_titles[i]) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 8)) #subtitle 
  # Append to the plot list
  plot_list[[i]] <- p
}

# Combine plots using patchwork
combined_plot <- wrap_plots(plot_list, ncol = 2) + 
  plot_annotation(title = title,
                  theme = theme(plot.title = element_text(size = 8, 
                                                          face = "bold", 
                                                          hjust = 0.5)))

# Save the combined plot to a file
save_to <- 'results/2_quality_control/Bisc_05022025/combined_QQplots_model3.png'
ggsave(save_to, combined_plot, width = 7, height = 15, dpi = 300)

# Display the combined plot
#print(combined_plot)





#==============================================================================#
#                          DEFINE PATHS                                        #
#==============================================================================#

main_path <- c('results/2_quality_control/Bisc_05022025/nighttime_fasting_model_3_sep_exp/QQplot_lambda.png',
               'results/2_quality_control/Bisc_05022025/nighttime_fasting_model_3/QQplot_lambda.png',  # adjsuted first meal
               'results/2_quality_control/Bisc_05022025/first_meal_model_3_sep_exp/QQplot_lambda.png', 
               'results/2_quality_control/Bisc_05022025/first_meal_model_3_firstlast/QQplot_lambda.png',
               'results/2_quality_control/Bisc_05022025/last_meal_model_3_sep_exp/QQplot_lambda.png',
               'results/2_quality_control/Bisc_05022025/last_meal_model_3_firstlast/QQplot_lambda.png',
               'results/2_quality_control/Bisc_05022025/eating_occasions_model_3_sep_exp/QQplot_lambda.png',
               'results/2_quality_control/Bisc_05022025/eating_occasions_model_3_firstlast/QQplot_lambda.png',
               'results/2_quality_control/Bisc_05022025/eating_jetlag_model_3/QQplot_lambda.png')

#==============================================================================#
#                               JOIN PLOTS                                     #
#==============================================================================#

title <- 'MODEL 3'
sub_titles <- c('night time fasting', 'night time fasting adj. first meal',
                'first meal', 'first meal adj. last meal', 
                'last meal', 'last meal adj. first meal', 
                'eating occasions', 'eating occasions adj. by first and last meal',
                'eating jetlag')


# Initialize an empty list to store plots
plot_list <- list()

# Loop through file paths and create plots
for (i in seq_along(main_path)) {
  # Read image
  img <- readPNG(main_path[i])
  # Create a ggplot object with the image
  p <- ggplot() +
    annotation_raster(img, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
    ggtitle(sub_titles[i]) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 8)) #subtitle 
  # Append to the plot list
  plot_list[[i]] <- p
}

# Combine plots using patchwork
combined_plot <- wrap_plots(plot_list, ncol = 2) + 
  plot_annotation(title = title,
                  theme = theme(plot.title = element_text(size = 8, 
                                                          face = "bold", 
                                                          hjust = 0.5)))

# Save the combined plot to a file
save_to <- 'results/2_quality_control/Bisc_05022025/combined_QQplots_model3.png'
ggsave(save_to, combined_plot, width = 7, height = 15, dpi = 300)

# Display the combined plot
#print(combined_plot)


