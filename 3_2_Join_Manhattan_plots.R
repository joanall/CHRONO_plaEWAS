#==============================================================================#
#                                                                              #
#               MATERNAL DIET MEALTIME AND PLACENTA DNA METHYLATION            #
#                             Analysis DMPsñ                                   #
#                                                                              #                                  
#                   Camille Lassale & Joana Llauradó Pont                      #         
#                               Janaury 2025                                   #
#                                                                              #
#==============================================================================#

# DESCRIPTION: Read PNG pltos with manhttan plots an join them


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

main_path <- c('results/3_analysis_dmps/Bisc_05022025/nighttime_fasting_model_noconf_sep_exp/nighttime_fasting_Bisc_model_noconf_sep_exp_05022025_adjusted_manhattan.png',
               'results/3_analysis_dmps/Bisc_05022025/nighttime_fasting_model_noconf_adj_first_meal/nighttime_fasting_Bisc_model_noconf_adj_first_meal_05022025_adjusted_manhattan.png',
               
               'results/3_analysis_dmps/Bisc_05022025/first_meal_model_noconf_sep_exp/first_meal_Bisc_model_noconf_sep_exp_05022025_adjusted_manhattan.png',
               'results/3_analysis_dmps/Bisc_05022025/first_meal_model_noconf_firstlast/first_meal_Bisc_model_noconf_firstlast_05022025_adjusted_manhattan.png',
               
               'results/3_analysis_dmps/Bisc_05022025/last_meal_model_noconf_sep_exp/last_meal_Bisc_model_noconf_sep_exp_05022025_adjusted_manhattan.png',
               'results/3_analysis_dmps/Bisc_05022025/last_meal_model_noconf_firstlast/last_meal_Bisc_model_noconf_firstlast_05022025_adjusted_manhattan.png',
               
               'results/3_analysis_dmps/Bisc_05022025/eating_occasions_model_noconf_sep_exp/eating_occasions_Bisc_model_noconf_sep_exp_05022025_adjusted_manhattan.png',
               'results/3_analysis_dmps/Bisc_05022025/eating_occasions_model_noconf_firstlast/eating_occasions_Bisc_model_noconf_firstlast_05022025_adjusted_manhattan.png') 


#==============================================================================#
#                               JOIN PLOTS                                     #
#==============================================================================#

title <- 'SENS.MODEL NO-CONF;'
sub_titles <- c('night time fasting', 'night time fasting adj. first meal',
                'first meal', 'first meal adj. last meal', 
                'last meal', 'last meal adj. first meal', 
                'eating occasions' ,'eating occasions adj. by first and last meal')


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
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 5)) #subtitle 
  # Append to the plot list
  plot_list[[i]] <- p
}

# Combine plots using patchwork
combined_plot <- wrap_plots(plot_list, ncol = 2) + 
  plot_annotation(title = title ,
                  theme = theme(plot.title = element_text(size = 8, 
                                                          face = "bold", 
                                                          hjust = 0.5)))

# Save the combined plot to a file
save_to <- paste0('results/3_analysis_dmps/Bisc_05022025/combined_manhattan_plots_', title,'.png')
ggsave(save_to, combined_plot, width = 7, height = 10, dpi = 300)

print(combined_plot)


