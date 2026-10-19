#==============================================================================#
#                                                                              #
#               MATERNAL DIET MEALTIME AND PLACENTA DNA METHYLATION            #
#                               Functions                                      #   
#                                                                              #                                  
#                   Camille Lassale & Joana Llauradó Pont                      #         
#                               Janaury 2025                                   #
#                                                                              #
#==============================================================================#

#------------------------------------------------------------------------------#
#                           DIET SCORES FUNCTIONS                              #----
#------------------------------------------------------------------------------#


# Positive scoring function
# - Default case 
# - Case 1 : If > 33% participants had no intake 
#   1) Set the first "tertile" to zero (i.e., the lowest intake in your data, which includes >33% of participants);
#   2) Then dividing the second and third "tertiles" using the median of the rest participants (i.e., those with non-zero intakes).
# - Case 2 : If 33% == 66 %  
#   1) Set the first and second "tertile" to 1, representing the lowest intake. (i.e.)
#   2) Set the third "tertile" to 3, representing the highest intake. (i.e)
# - Case 3: If 66 = 100% 
#   1) Set the third "tertile" to 3, representing the highest intake.
#   2) Divide the rest of the participants into the first and second "tertiles" using the median of their non-maximum intakes.

p_score_3 <- function(df,var_name){
  
  my_var <- df[,var_name]
  
  #### First obtain the cutoffs for tertiles
  P0 <- min(my_var, na.rm = T)
  P33.3 <- quantile(my_var, probs = 0.333, na.rm = T)
  P66.7 <- quantile(my_var, probs = 0.667, na.rm = T)
  P100 <- max(my_var, na.rm = T)
  my_var_3 <- rep(NA, length(my_var)) # Create empty vector 
  
  # CAS 1: When >33.3% with the same lowest intake (most likely to be 0)
  if (P0 == P33.3) {
    # Warning message and solution
    writeLines(paste0("Warning for ", var_name, ": >33.3% with the same lowest intake (", P0, ");
                      \nSolution: Set all of the lowest as Q1, then take the median of the rest to split Q2 and Q3."))
    ##### Modify the range of tertiles and assign positive scores (1-3)
    my_var_3[which(my_var == P0)] <- 1 # Set all of those with the lowest intake as Q1
    my_var_3[which(my_var > P0 & my_var <= median(my_var[my_var > P0], na.rm = T))] <- 2  # For the rest, set those with intake <=median as Q2
    my_var_3[which(my_var > median(my_var[my_var > P0], na.rm = T))] <- 3  # And then set those with intake >median as Q3
    return(my_var_3)
  }
  
  # CAS 2: P33.3 = P66.7 
  if (P33.3 == P66.7) {
    my_var_3[which(my_var < P33.3)] <- 1  # Set minimum intake <= P66.7 (Q1) 
    my_var_3[which(my_var == P66.7)] <- 2  # Set minimum intake <= P66.7 (Q1) 
    my_var_3[which(my_var > P66.7)] <- 3  # Set maximum with the highest intake as Q3
    return(my_var_3)
  }
  
  # CAS 3: P66.6 == P100 
  
  if (P66.7 == P100) {
    ##### Modify the range of tertiles and assign positive scores (1-3)
    my_var_3[which(my_var == P100)] <- 3  # Set all of those with the lowest intake as Q1
    my_var_3[which(my_var < P100 & my_var <= median(my_var[my_var < P100], na.rm = T))] <- 2  # For the rest, set those with intake <=median as Q2
    my_var_3[which(my_var > median(my_var[my_var < P100], na.rm = T))] <- 3  # And then set those with intake >median as Q3
    return(my_var_3)
  }
  
  #### In other situations - Take tertiles and assign positive scores
  else {
    my_var_3 <- findInterval(my_var, c(-Inf, quantile(my_var, probs = c(0.333, 0.667), na.rm = T), Inf))
    return(my_var_3)
  }
  
  
}

# Negative scoring function
r_score_3 <- function(df, var_name) {
  #Tertils variables 
  my_var_3 <- p_score_3(df,var_name)  # Same as above for setting tertiles
  
  my_var_3 <- dplyr::recode(my_var_3, "1" = 3, "2" = 2, "3" = 1)  # Recoded as reverse scores
  
  #Dichotomic variables
  return(my_var_3)
}

#------------------------------------------------------------------------------#
#                           MEALTIME  FUNCTIONS                                #----
#------------------------------------------------------------------------------#

get_nighttime_fasting <- function (last_meal, first_meal) {
  
  # If else to check whether last meal was before or after midnight
  result <- dplyr::if_else((!is.na(last_meal) & !is.na(first_meal) & 
                              last_meal >= as.duration(hm("12:00")) & last_meal <= as.duration(hm("24:00"))),
                           as.duration(hm("24:00")) - (last_meal - first_meal), # when last meal is day before
                           as.duration(first_meal - last_meal)) # otherwise do this
  return(result)
}

get_eating_window <- function (last_meal, first_meal) {
  
  # Init df to store results of fasting
  result <- as.data.frame(matrix(nrow=length(last_meal), ncol=0))
  
  # If else to check whether last meal was before or after midnight
  result <- dplyr::if_else(last_meal <= as.duration(hm("12:00")), # if last meal is the next day
                            as.duration(hm("24:00")) - first_meal + last_meal, 
                              as.duration(last_meal - first_meal)) #else
  
  return(result)
}

# We need a function to compute the difference of hours between last meal and the time
# of goign to bed. Teh iputs is 2 vector of a dataframe 
#Input variables are in duration format, we need to handle the following cases 
# 1) if either value is NA, return NA 
# 2) if go bed is after midnight (in here seen as 0s or 30 min etc ) we need to account
# for that so the hours between last meal and going to be make sense 

# 3) if last meal is bigger than go bed, either in the same day or in the next day 
# for exmaple las tmeal at 1am and go bed at 0:30 am, or last meal at 21pm and go bed at 20pm
# we set the fasting to 0 houors 


get_fasting_before_sleep <- function(last_meal, go_bed) {
  # Create an empty vector to store fasting hours
  fasting_hours <- numeric(length(last_meal))
  
  for (i in seq_along(last_meal)) {
    meal_time <- last_meal[i]
    bed_time <- go_bed[i]
    
    # Case 1: If either value is NA, return NA
    if (is.na(meal_time) | is.na(bed_time)) {
      fasting_hours[i] <- NA
      next
    }
    
    # Convert duration to numeric hours
    meal_hours <- as.numeric(meal_time, "hours")
    bed_hours <- as.numeric(bed_time, "hours")
    
    # Case 2: Adjust for bedtime after midnight
    if (meal_hours == 0 && bed_hours > 20) {  # Midnight meal but bedtime is late night
      fasting_hours[i] <- 0
    } else if (bed_hours < meal_hours || (meal_hours - bed_hours) >= 24) {
      fasting_hours[i] <- 0  # Case 3: Last meal is after bedtime or invalid fasting period
    } else {
      fasting_hours[i] <- bed_hours - meal_hours
    }
  }
  
  return(fasting_hours)
}


#------------------------------------------------------------------------------#
#                         DESCRIPTIVE  FUNCTIONS                               #----
#------------------------------------------------------------------------------#

show_values <- function(df,var_name) {
  my_var <- df[, var_name]
  x <- quantile(my_var, probs = c(0, 0.1, 0.2, 0.25, 0.3, 0.333, 0.4, 0.5, 0.6, 0.667, 0.7, 0.75, 0.8, 0.9, 1), na.rm = T)
  y <- data.frame("Mealtime" = var_name,
                  "N" = length(na.omit(my_var)),
                  `Min` = x[1], `P10` = x[2], `P20` = x[3], `P25` = x[4], `P30`  = x[5],
                  `P33.3` = x[6], `P40` = x[7], `P50` = x[8], `P60` = x[9],
                  `P66.7` = x[10], `P70` = x[11], `P75` = x[12], `P80` = x[13], `P90` = x[14], `Max` = x[15])
  return(y)
}

# Function to compute tertiles
tertile_cutoff <- function(df,var_name) {
  my_var <- df[, var_name]
  my_var_3 <- p_score_3(df,var_name)
  
  T1_lo <- round(min(my_var[my_var_3 == 1], na.rm = T), digits = 2)
  T1_up <- round(max(my_var[my_var_3 == 1], na.rm = T), digits = 2)
  T2_lo <- round(min(my_var[my_var_3 == 2], na.rm = T), digits = 2)
  T2_up <- round(max(my_var[my_var_3 == 2], na.rm = T), digits = 2)
  T3_lo <- round(min(my_var[my_var_3 == 3], na.rm = T), digits = 2)
  T3_up <- round(max(my_var[my_var_3 == 3], na.rm = T), digits = 2)
  x <- paste0("T1: ", T1_lo, "-", T1_up, " | T2: ", T2_lo, "-", T2_up, " | T3: ", T3_lo, "-", T3_up)
  return(x)
}

# Function to compute tertiles and scoring food groups
describe_food_group <- function(df,var_name) {
  x <- c(`Food group` = var_name, "Cutoffs" = tertile_cutoff(df,var_name)) 
  x <- c(x, table(factor(p_score_3(df,var_name), label = c("T1", "T2", "T3")), useNA = "always"))
  
  return(x)
  
}

#------------------------------------------------------------------------------#
#                       METHYLATION QC  FUNCTIONS                              #----
#------------------------------------------------------------------------------#

exclude_CpGs <- function(cohort, cpgid, exclude, ethnic, artype='450K', filename = NULL, fileresume = NULL) {
  # Test if filter data exists and gets it
  if( toupper(artype) == '450K' ) {
    #try(data("filter_450K") )
    filters <- filter_450K
  }else if(toupper(artype) == 'EPIC') {
    #try(data("filter_EPIC"))
    filters <- filter_EPIC
  }else{
    stop( paste0( "Unknown array type ", toupper(artype) ) )
  }
  
  
  # In order to minimize data size we merge only the selected ethnic column
  
  
  # Get first ethnic column position in filters
  firstEthnicPosition <-  min(grep(paste0("MASK_snp5_.*[^GMAF1p&^common]"), colnames(filters), perl = TRUE))-1
  
  # Select ethnic columns not related with our data from filters (ONLY IF ethnic != '' or NA)
  if(ethnic!='' && !is.na(ethnic)) {
    # Get first ethnic column position in filters
    fieldstodelete <- grep(paste0("MASK_snp5_.*[^",ethnic,"]"),  colnames(filters)[grep(paste0("MASK_snp5_.*[^GMAF1p&^common]"), colnames(filters), perl = TRUE)], perl = TRUE) + firstEthnicPosition
    
  } else {
    #.Remove from code.# fieldstodelete <- grep(paste0("MASK_snp5_.*[^EUR]"),  colnames(filters)[grep(paste0("MASK_snp5_.*[^GMAF1p&^common]"), colnames(filters), perl = TRUE)], perl = TRUE) + firstEthnicPosition
    fieldstodelete <- ""
  }
  
  
  fieldstomerge <- which(!seq(1:dim(filters)[2]) %in% fieldstodelete)
  
  # Merge cohort with CpGs filters
  cohort <- merge(cohort, filters[,fieldstomerge], by.x= cpgid, by.y = "probeID", all.x = TRUE )
  
  
  # CpGs id to remove
  if( is.null(exclude)) {
    excludeid <- NULL
  } else {
    excludeid <- cohort[eval(parse(text=getCritera(exclude, ethnic))), cpgid]
  }
  
  
  # Report descriptive exclussions to a descriptive file
  if(!is.null(fileresume)) {
    write(sprintf('\n# %s', strrep("-",16)), file = fileresume, append = TRUE)
    write(sprintf('# Remove "problematic"  CpGs : '), file = fileresume, append = TRUE)
    write(sprintf('# %s\n', strrep("-",16)), file = fileresume, append = TRUE)
    write(sprintf('# Criteria : \n\tArray type : %s \n\tEthnia : %s \n', toupper(artype), toupper(ethnic)), file = fileresume, append = TRUE)
    write(sprintf('# Mask : %s', exclude), file = fileresume, append = TRUE)
    write(sprintf('\n'), file = fileresume, append = TRUE)
    write(sprintf('# Total CpGs in data : %d', dim(cohort)[1]), file = fileresume, append = TRUE)
    write(sprintf('# Number of excluded CpGs: %d', length(excludeid)), file = fileresume, append = TRUE)
    write(sprintf('# Total CpGs after exclusion : %d\n', (dim(cohort)[1]) - length(excludeid)), file = fileresume, append = TRUE)
    write(sprintf('# Percent excluded CpGs: %f %%\n', ((length(excludeid)/dim(cohort)[1])*100 )), file = fileresume, append = TRUE)
    
    suppressWarnings(
      write.table( cohort[eval(parse(text=getCritera(exclude, ethnic))),],
                   filename, col.names = TRUE, row.names = FALSE, sep = '\t', append = FALSE, dec='.'))
  }
  
  # # Report CpG excluded and reason to a file
  # if(!is.null(filename)) {
  #    # write(sprintf('# Criteria : %s\n', toupper(artype)), file = filename)
  #    # write(sprintf('# Total CpGs in data : %d', dim(cohort)[1]), file = filename, append = TRUE)
  #    # write(sprintf('# Number of excluded CpGs: %d', length(excludeid)), file = filename, append = TRUE)
  #    # write(sprintf('# Total CpGs after exclusion : %d\n', (dim(cohort)[1]) - length(excludeid)), file = filename, append = TRUE)
  #    # write(sprintf('# Percent excluded CpGs: %f %%\n', ((length(excludeid)/dim(cohort)[1])*100 )), file = filename, append = TRUE)
  #    # Report exclusion reason
  #    suppressWarnings(
  #       write.table( cohort[eval(parse(text=getCritera(exclude, ethnic))),],
  #                    filename, col.names = TRUE, row.names = FALSE, sep = '\t', append = FALSE, dec='.'))
  # }
  
  # Rmove CpGs with exclusion parameters
  cohort <- cohort[ !cohort[,cpgid] %in% excludeid,]
  
  warning( length(excludeid), ' CpGs have been excluded')
  
  return(cohort)
  
}


getCritera <- function(exclude, ethnic){
  
  criter = ''
  possible_crit <- c( 'MASK_sub25_copy', 'MASK_sub30_copy', 'MASK_sub35_copy', 'MASK_sub40_copy',
                      'MASK_mapping', 'MASK_extBase', 'MASK_typeINextBaseSwitch', 'MASK_snp5_common', 'MASK_snp5_GMAF1p',
                      'MASK_general', 'cpg_probes', 'noncpg_probes', 'control_probes', 'Unrel_450_EPIC_blood', 'MASK_rmsk15',
                      'Sex', 'Unrel_450_EPIC_pla_restrict', 'Unrel_450_EPIC_pla', 'MASK_snp5_ethnic')
  
  if(!is.null(exclude[1]) && !is.na(exclude[1]) && exclude[1] != '') {
    # Test if all parameters are allowed
    if(length(which(! exclude %in% possible_crit))>=1)
      stop(paste0('Parameter(s) : ',paste(exclude[which(! exclude %in% possible_crit)], sep = ','),' not valid. Possible values are ',possible_crit ))
    
    # Get all exclude variables
    ##.. Works with lists -> convert list to matrix ..## exclusion.crit <- do.call(cbind, lapply( ls(patt="exclude"), get) )
    
    # Create formula with exclusion criteria
    ##.. Works with lists -> gets only parameters with value 'Exclude'..## criter <- paste(paste("cohort$",attributes(exclusion.crit)$dimnames[[1]][which(exclusion.crit=='Exclude')], sep = ""),"TRUE | ",sep=" == ", collapse = '')
    criter <- paste(paste("cohort$",exclude, sep = ""),"TRUE | ",sep=" == ", collapse = '')
    
    # If  MASK_general = 'Exclude' --> Change to :  ("MASK.sub30.copy", "MASK.mapping", "MASK.extBase", "MASK.typeINextBaseSwitch" and "MASK.snp5.GMAF1p") == TRUE
    if( length(grep("cohort$MASK_general == TRUE",criter, fixed = TRUE))>0 )
      criter <- sub("cohort$MASK_general == TRUE", " cohort$MASK_sub30_copy == TRUE | cohort$MASK_mapping == TRUE | cohort$MASK_extBase == TRUE | cohort$MASK_typeINextBaseSwitch == TRUE | cohort$MASK_snp5_GMAF1p == TRUE " , criter, fixed = TRUE)
    
    # If Sex = 'Exclude' --> Change to : ( CpG_chrm %in% "chrX" | CpG_chrm %in% "chrY" )
    if( length(grep("cohort$Sex == TRUE",criter, fixed = TRUE))>0 )
      criter <- sub("cohort$Sex == TRUE", " cohort$CpG_chrm %in% 'chrX' | cohort$CpG_chrm %in% 'chrY' " , criter, fixed = TRUE)
    
    newCpGcond <- vector()
    if( length( grep("| cohort$cpg_probes == TRUE",criter, fixed = TRUE)) >0 )
      newCpGcond <- append(newCpGcond, c('cg'))
    if( length(grep("| cohort$noncpg_probes == TRUE",criter, fixed = TRUE))>0 )
      newCpGcond <- append(newCpGcond, c('ch'))
    if( length(grep("| cohort$control_probes == TRUE",criter, fixed = TRUE))>0 )
      newCpGcond <- append(newCpGcond, c('rs'))  
    
    if(!is.null(vector())) {
      # Remove invalid criteria
      criter <- sub("| cohort$cpg_probes == TRUE ", "", criter, fixed = TRUE)
      criter <- sub("| cohort$noncpg_probes == TRUE ", "", criter, fixed = TRUE)
      criter <- sub("| cohort$control_probes == TRUE ", "", criter, fixed = TRUE)
      criter <- sub("| cohort$MASK_snp5_ethnic == TRUE ", "", criter, fixed = TRUE)
      criter <- sub("cohort$cpg_probes == TRUE |", "", criter, fixed = TRUE)
      criter <- sub("cohort$noncpg_probes == TRUE |", "", criter, fixed = TRUE)
      criter <- sub("cohort$control_probes == TRUE |", "", criter, fixed = TRUE)
      criter <- sub("cohort$MASK_snp5_ethnic == TRUE |", "", criter, fixed = TRUE)
      
      # Add new criteria
      criter <- paste0(criter, paste0(" cohort$probeType %in% '",newCpGcond,"' | ", collapse = ''))
    }
    
    # Add ethnicity exlcusion criteria
    criter <- paste0(criter, paste0(" cohort$probeType %in% '",newCpGcond,"' | ", collapse= ''))
    if(is.na(ethnic) | ethnic == ''){
      criter <- paste("which(", substr(criter,0,nchar(criter)-3),")",sep="")
    }
  }
  
  if(!is.na(ethnic) & ethnic != ''){
    criter <- paste0(criter, " cohort$MASK_snp5_", ethnic," == TRUE | ")
    criter <- paste("which(", substr(criter, 0, nchar(criter)-2), ")", sep="")
  }
  criter
}

#------------------------------------------------------------------------------#
#                         DEFINE EWAS FUNCTION                                 #----
#------------------------------------------------------------------------------#


rlm_single_cpg <- function(cpg_chunk, exposures, 
                           covariates, metadata) {
  
  # Get data for covariates and exposures
  exp_cov <- c(exposures, covariates)
  df_exp_cov <- metadata[, exp_cov]
  
  
  # Initialize a list to store results for each exposure
  exposure_results <- vector("list", length(exposures))
  names(exposure_results) <- exposures
  
  # Loop through CpG sites in the chunk
  for (cpg_index in seq_len(nrow(cpg_chunk))) {
    
    cpg <- rownames(cpg_chunk)[cpg_index]
    #cat('Processing CpG:', cpg, "\n")
    
    # Select the methylation values for the current CpG site
    meth_vals_cpg <- as.data.frame(t(cpg_chunk[cpg_index, ]))
    colnames(meth_vals_cpg) <- cpg
    
    
    # Combine the methylation values with covariates and exposure data
    df_for_model <- cbind(meth_vals_cpg, df_exp_cov)
    
    # Formula for robust linear regression using the CpG name
    f <- as.formula(paste0(cpg, " ~ ", paste(exp_cov, collapse = " + ")))
    
    # Fit the robust linear model with the combined data
    fitted_model <- MASS::rlm(formula = f, data = df_for_model)
    #cat('Model fitted for CpG:', cpg, "\n")
    
    # Loop through each exposure and store the result in the corresponding data frame
    for (exp in exposures) {
      pval <- as.numeric(sfsmisc::f.robftest(fitted_model, var = exp)["p.value"])
      
      # Extract the coefficient, standard error, etc., for the exposure
      result <- as.data.frame(summary(fitted_model)$coefficients[exp, , drop = FALSE])
      result$pvalue <- pval
      result$cpg <- cpg
      
      names(result) <- c('Effect', 'SE', 'tvalue', 'pvalue', 'CpG')
      
      # Store results in the exposure_results list
      if (is.null(exposure_results[[exp]])) { exposure_results[[exp]] <- result
      } else { exposure_results[[exp]] <- rbind(exposure_results[[exp]], result)}
    }
  }
  
  
  # Return the list of data frames, each corresponding to one exposure
  return(exposure_results)
}


#------------------------------------------------------------------------------#
#                  Quality Control-  Precision plot function                   #----
#------------------------------------------------------------------------------#

plot_precisionp <- function(x, filename, main= NULL, ...) {
  
  if (is.null(main))
    main = "Precision Plot -  1/median(SE) vs sqrt(n)"
  
  p <- ggplot2::ggplot( data = x, mapping = ggplot2::aes( x = round(as.numeric(sqrt_N),2),
                                                          y = round(as.numeric(invSE),2) ) ) +
    ggplot2::theme_bw() +
    ggplot2::geom_point( size = 3, ggplot2::aes( colour = cohort ) ) +
    #geom_line( aes( group = cohort , colour = cohort ) ) +
    ggplot2::ggtitle( main ) +
    ggplot2::theme( legend.position = "bottom",
                    legend.text = ggplot2::element_text(size=6),
                    legend.title = ggplot2::element_blank()) +
    ggplot2::labs( x = "sqrt(n)",
                   y = "inv SE") +
    ggplot2::scale_shape_manual( values = 0:7 )
  
  if(!is.null(filename)) {
    png(filename, type = "cairo")
    print(p)
    dev.off()
    #..# ggplot2::ggsave(filename,p)
  }
  
  
  return(p)
  
  
}

#------------------------------------------------------------------------------#
#                     Quality Control-Inflation function                       #----
#------------------------------------------------------------------------------#

# Function to calculate lambda & confidence intervals
# inflation <- function(pvector) {
#   chisq <- qchisq(1 - pvector, 1)
#   lambda <- median(chisq) / qchisq(0.5, 1)
#   SE.median <- qchisq(0.975, 1) * (1.253 * ( sd(chisq) / sqrt( length(chisq))))
#   lower_ci <-  lambda - (SE.median / qchisq(0.5, 1))
#   upper_ci <-  lambda + (SE.median / qchisq(0.5, 1))
#   lambda_ci <- c(lambda, lower_ci, upper_ci)
#   return(lambda_ci)
# }
inflation <- function(pvector) {
  # Remove invalid p-values
  pvector <- pvector[!is.na(pvector)]
  pvector[pvector <= 1e-300] <- 1e-300
  pvector[pvector >= 1 - 1e-16] <- 1 - 1e-16
  
  if (length(pvector) == 0) return(c(NA_real_, NA_real_, NA_real_))
  
  chisq <- qchisq(1 - pvector, 1)
  
  # Remove Inf/NaN from chisq directly
  chisq <- chisq[is.finite(chisq)]
  
  if (length(chisq) == 0) return(c(NA_real_, NA_real_, NA_real_))
  
  lambda <- median(chisq) / qchisq(0.5, 1)
  
  # Compute SD safely
  s <- suppressWarnings(sd(chisq))
  
  # Use fallback CI if sd is NA, NaN, or <= 0
  if (length(chisq) > 1 && is.finite(s) && !is.na(s) && s > 0) {
    SE.median <- qchisq(0.975, 1) * (1.253 * s / sqrt(length(chisq)))
    lower_ci <- lambda - (SE.median / qchisq(0.5, 1))
    upper_ci <- lambda + (SE.median / qchisq(0.5, 1))
  } else {
    # Fallback ±5% range if CI fails
    lower_ci <- lambda * 0.95
    upper_ci <- lambda * 1.05
  }
  
  return(c(lambda, lower_ci, upper_ci))
}





#------------------------------------------------------------------------------#
#                  Quality Control-  QQ plot function                          #----
#------------------------------------------------------------------------------#

qq_norm_plot <- function(df, var_name,title=NULL) {
  
  vals <- df[[var_name]]
  o <- -log10(sort(vals,decreasing=F))
  e <- -log10( 1:length(o)/length(o))
  
  qqdata <- data.frame(o=o,e=e)
  
  ggplot(qqdata, aes(e,o)) +
    geom_point() +
    
    geom_abline(intercept=0,slope=1,col="#d7301f") + 
    
    # axis labels 
    xlab(expression(Expected~(-log[10](pval)))) + 
    ylab(expression(Observed~(-log[10](pval)))) +
    
    # styling
    theme_bw() +
    ggtitle(title) +
    theme(axis.title=element_text(face="bold",size=12),
          plot.title=element_text(face="bold",size=18,hjust=0))
}

# Function to add lambda and CI intervals to a qqplot
add_lambda_to_plot <- function(pvector, qqplot) {
  lambda <- inflation(pvector)
  # Add lamba + CI 
  lambda_CI <- sprintf("λ = %.2f [%.2f,%.2f]", lambda[1], lambda[2], lambda[3]) 
  qqplot <- qqplot + annotate(geom = "text", x = -Inf, y = Inf,
                              hjust = -0.15, vjust = 1 + 0.15 * 3, label = lambda_CI, size = 5)
  return(qqplot)
}



qq_plot2 <- function(df, var_name) {
  pvals <- df[[var_name]]
  n <- length(pvals)
  
  observed <- -log10(sort(pvals))
  expected <- -log10(ppoints(n))
  
  qq_data <- data.frame(Expected = expected, Observed = observed)
  
  ggplot(qq_data, aes(x = Expected, y = Observed)) +
    geom_point(size = 0.6, color = "black") +
    geom_abline(slope = 1, intercept = 0, color = "#d7301f") +
    labs(
      x = expression(Expected~~-log[10](italic(p))),
      y = expression(Observed~~-log[10](italic(p)))
    ) +
    theme_classic() +
    theme(
      axis.title = element_text(size = 10, face = "bold"),
      panel.grid = element_blank()
    ) +
    coord_fixed()
}



add_lambda_to_plot2 <- function(lambda_ci, qqplot) {
  label_text <- sprintf("λ = %.2f [%.2f, %.2f]",
                        lambda_ci[1, "lambda"], lambda_ci[1, "lower_ci"], lambda_ci[1, "upper_ci"])
  
  y_pos <- max(qqplot$data$Observed)
  
  qqplot +
    annotate("text", x = 0, y = y_pos, label = label_text,
             hjust = 0, vjust = 1, size = 4, fontface = "bold")
}

#------------------------------------------------------------------------------#
#                  Quality Control- Distribution plot function                 #----
#------------------------------------------------------------------------------#

# Fucntion to do a dsitribution plot 
plot_distribution <- function(x, main='', xlab='',...) {
  
  h <- hist(x, breaks = 100, plot = FALSE)
  d <- density(x)
  plot(h, freq = FALSE, col = ggcolors(2)[2], border = 'white',
       main = main, xlab = xlab, ylim = c(0, max(d$y, h$density)))
  lines(d, col = ggcolors(2)[1], lwd = 2)
  
}

#------------------------------------------------------------------------------#
#                 Analysis DMPs- define summary table                          #----
#------------------------------------------------------------------------------#

get_summary_table <- function(list_cpgs, list_files, new_filename){
  
  if (length(list_cpgs) != 0) {
    # Create a summary table for each EWAS result with the specified columns
    summary_tables <- lapply(seq_along(list_files), function(i) {
      ewas <- list_files[[i]]
      summtab <- ewas[ewas$CpG %in% list_cpgs, c('CpG', 'Effect', 'SE', 'pvalue', 'bf')]
      colnames(summtab) <- c('CpG', paste0('coef_model', i), 
                             paste0('se_model', i), 
                             paste0('pval_model', i), 
                             paste0('bf_model', i))
      return(summtab)
    })
    
    # Merge all summary tables by CpG
    combined_summary_table <- Reduce(function(x, y) merge(x, y, by = "CpG", all = TRUE), summary_tables)
    
    # Order by CpG
    combined_summary_table <- combined_summary_table[order(combined_summary_table$CpG), ]
    
    write.table(combined_summary_table, file = new_filename, sep = ";", row.names = FALSE)
    
  }
}



#------------------------------------------------------------------------------#
#                     functions for Annotation                                 #----
#------------------------------------------------------------------------------#
get_annot_EPIraction <- function(gzip_file, out_folder) {
  
  # Go to EPIraction Tool and Get the form 
  url <- "https://epiraction.crg.eu/"
  epiraction_url <- read_html(url)
  
  forms <- html_form(epiraction_url)
  
  # We have 6 forms, if we go to the web browser we'll get to see the matching so 
  # 1 - interval --> genomic insights
  # 2 - rss_gene --> gene variants
  # 3 - overlap  --> gwas annotation (file) --> this is the one we want
  # 4 - enrichment  --> gwas annotation (gene name)
  # 5 - regions_list --> regions 
  # 6 - gene --> enhancer-gene pairs
  
  # Prepare form for file upload
  form_index <- 3  # The actual form index for "gwas annotation (file)"
  file_input_name <- "overlap"  # The name of the file input field
  file_path <- gzip_file
  form <- forms[[form_index]]
  
  # Check if the file exists
  if (!file.exists(file_path)) {
    stop("The file does not exist. Please check the file path.")
  } else { cat("File exists: ", file_path, "\n") }
  
  # Step 5: Create a session and set the form values using html_form_set
  session <- session(url)
  form <- html_form_set(form, !!file_input_name := file_path)
  print(form)
  
  # Prepare the form data for submission manually
  form_data <- list(overlap = upload_file(file_path))
  cat("Form data to be submitted:\n")
  print(form_data)
  
  # Step 6: Submit the form using POST request with multipart encoding
  response <- POST(url = form$action,
                   body = form_data,
                   encode = "multipart")
  
  # Check the response status
  if (response$status_code == 200) {
    cat("Form submission successful.\n")
    
    # Extract the HTML content
    content_text <- rawToChar(response$content)
    content_html <- read_html(content_text)
    
    # Extract tables
    tables <- html_table(content_html, fill = TRUE)
    
    if (length(tables) > 0) {
      print(tables[[1]]) # Parameters used (1st table it appears in the results web)
      
      # Extract the second table -> table with links to download 
      table2 <- tables[[2]]
      
      # Extract download links
      download_links <- content_html %>%
        html_nodes(xpath = "//table[2]//td/a[text()='Download']") %>%
        html_attr("href")
      
      # We want the first link -> all tissues as row with QTls 
      link <- download_links[1]
      download_url <- paste0("https://epiraction.crg.eu/", link)
      
      response <- GET(download_url)
      
      # Check if the request was successful
      if (response$status_code == 200) {
        destfile <- paste0(out_folder, "/", basename(link), ".txt")
        writeBin(response$content, destfile)
        cat("Downloaded and saved as: ", destfile, "\n")
      } else {
        cat("Failed to download: ", download_url, "\n")
      }
      
      return(destfile)
      
    } else {
      cat("No tables found in the response.\n")
    }
  } else {
    cat("Form submission failed with status code:", response$status_code, "\n")
    
    # Print the response content for debugging
    content_text <- rawToChar(response$content)
    cat("Response content:\n")
    cat(content_text)
  }
}

# Find Overlaps between two genomics ranges and returns the coincident values. 
# Search overlaps between CpGs (selected data) and Fetal Placenta (complete data)

findoverlapvalues <- function(selected, complete) {
  
  # Find overlaps
  ranges <- subsetByOverlaps(selected, complete)
  hits <- findOverlaps(selected, complete)
  
  # Find the names of states associated to each position
  idx <- subjectHits(hits)
  values <- DataFrame(States  = complete$name[idx])
  qh <- queryHits(hits)
  
  ans <- list("ranges" = ranges,
              "hits" = hits,
              "shits" = idx,
              "qhits" = qh,
              "values" = values)
  
  return(ans)}



