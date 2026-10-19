#==============================================================================#
#                                                                              #
#                MATERNAL DIET MEALTIME AND PLACENTA DNA METHYLATION           #
#                              EWAS Analysis                                   #
#                                                                              #                                  
#                   Camille Lassale & Joana Llauradó Pont                      #         
#                              April 2024                                      #
#                                                                              #
#==============================================================================#

# Description: Get the descriptive tables of the cohort used in this analysis.

# [INPUT]: RDS dataframe with Nrow samples and colnames covariables included in the EWAS
# [OUTPUT]: RDS dataframe with exposures needed for the analysis 



#==============================================================================#
#                             DEFINE PARAMETERS                                #
#==============================================================================#

metadata_path <- "db/metadata_plaEWAS_MT_incomplete_cases2_wgtest.rds"
#metadata_path <- "db/archive/metadata_plaEWAS_MT_complete_cases.rds"
cohort_specific_variables <- c("study_center", "covid_confinement")
current_date <- "test"

#==============================================================================#
#                       INIT  WORKING ENVIRONMENT                              #
#==============================================================================#

print(paste0('YOUR WORKING DIRECTORY IS:', getwd()))

# Load libraries 
library(dplyr)
library(lubridate) # mealtime data 
library(tableone) # basic descriptives
library(psych) # basic descriptives

source("script/Define_functions.R")


#==============================================================================#
#                                                                              #
#                       LOAD DATA & COMPUTE SCORES                             #----
#                                                                              #
#==============================================================================#

cat("\n--------------------   START OF CODE  -------------------------------\n")

# Load data 
df <- readRDS(metadata_path)  #389

#------------------------------------------------------------------------------#
#                         Compute rMEDp score                                  #----
#------------------------------------------------------------------------------#

# Step 1: Adjsut for energy desitity and compute tertiles  for food groups -

# STEP 2: Single-item mMDS score 
# - For each food group based on whether values are above or below the defined tertile (based on whether food are "beneficial" or "detrimental") 
# assign a value of 0, 1 or 3 -> the p_score funcitons is vuil on 1-2-3 (extract 1)

# Food groups classification
# food groups that fit more the Mediterranean diet  (5) -> vegetables, legumes, fruit (incl.nuts,seeds), fish and seafood, cereals
# food groups that fit less the Mediterranean diet  (2) -> meat, dairy product 
#  olive oil median

# STEP 3: Total Score 
# Ranges 0-16, with higher values reflecting higher adherence to a traditional Mediterranean diet.

# References: check that 
# Küpers, L. K.,et al. 2022. https://doi.org/10.1080/15592294.2022.2038412

cat("\n------------------- COMPUTE rMEDp SCORES ----------------------------\n")

# Adjust food groups for energy density (*1000/cals)
df$vegetables <- 1000*(df$vegetables/df$energy_intake)
df$legumes <- 1000*(df$legumes/df$energy_intake)
df$fruit_nuts <- 1000*(df$fruit_nuts/df$energy_intake)
df$fish <- 1000*(df$fish/df$energy_intake)
df$cereals <- 1000*(df$cereals/df$energy_intake)
df$meat <- 1000*(df$meat/df$energy_intake)
df$dairy <- 1000*(df$dairy/df$energy_intake)
df$olive_oil <- 1000*(df$olive_oil/df$energy_intake)


# Compute rMEDp
df$rMEDp <- p_score_3(df, 'vegetables') - 1 + p_score_3(df, 'legumes') - 1 + p_score_3(df, 'fruit_nuts') - 1 + 
  p_score_3(df, 'fish') - 1 + p_score_3(df,'cereals') - 1 + # 5 fit more meditteranean diet (+)
  
  r_score_3(df, 'meat') - 1 + r_score_3(df, 'dairy') - 1 + # 2 fit less meditteranean diet(-) 
  
  ifelse(df[['olive_oil']] < median(df[['olive_oil']], na.rm = TRUE),0,1)


#------------------------------------------------------------------------------#
#                     Subset data with mealtime variables                      #----
#------------------------------------------------------------------------------#

cat("\n------------------- GET MEALTIME EXPOSURES ----------------------------\n")

mealtime_vars <- c('first_meal_wday','last_meal_wday',
                   'first_meal_wend','last_meal_wend',
                   'eating_occasions_wday','eating_occasions_wend')


temp <- df[, c('Basename', mealtime_vars)]


#------------------------------------------------------------------------------#
#                     Compute Nighttime Fasting Duration                       #----
#------------------------------------------------------------------------------#

# Night time fasting = Timing of the last meal - Timing of the first meal 
# Average night time fasting (h) = [(5 × night time fasting weekdays ) + (2 × night time fasting on weekends)]/7.
# Alternative: Average night time fasting (h) = 24h - Average eating duration 

# References: 
# - Palomar-Cros A, et.al,  2021. doi: 10.3390/nu13082662
# - Kesztyüs D, et.al, 2019. doi: 10.3390/nu11122854

#  Get nigh-time fasting for week days & week ends
temp$nighttime_fasting_wday <- get_nighttime_fasting(temp$last_meal_wday, temp$first_meal_wday)
temp$nighttime_fasting_wend <- get_nighttime_fasting(temp$last_meal_wend, temp$first_meal_wend)

# Get average
temp$nighttime_fasting <- dplyr::if_else(!is.na(temp$nighttime_fasting_wday) & !is.na(temp$nighttime_fasting_wend),   # If we have both wday and wend
                                         ((temp$nighttime_fasting_wday * 5) + (temp$nighttime_fasting_wend * 2)) / 7, # Use both to calculate the average
                                         dplyr::if_else(!is.na(temp$nighttime_fasting_wday) & is.na(temp$nighttime_fasting_wend),  # If we have wday but not wend
                                                        temp$nighttime_fasting_wday,                                       # Use wday
                                                        dplyr::if_else(is.na(temp$nighttime_fasting_wday) & !is.na(temp$nighttime_fasting_wend), # If we have wend but not wday
                                                                       temp$nighttime_fasting_wend,                                     # Use wend
                                                                       NA)))         



#------------------------------------------------------------------------------#
#                         Compute Eating Jet Lag                               #----
#------------------------------------------------------------------------------#

# References: 
# Zerón-Rugerio MF, et.al, 2019. doi: 10.3390/nu11122980. 

# Get eating window 
temp$eating_window_wday <- get_eating_window(temp$last_meal_wday, temp$first_meal_wday)
temp$eating_window_wend <- get_eating_window(temp$last_meal_wend, temp$first_meal_wend)

# Get eating midpoint 
# Eating midpoint (local time) = ([Timing of the last meal − Timing of the first meal]/2) + Timing of the first meal.

temp$eating_midpoint_wday <- (temp$eating_window_wday)/2 + temp$first_meal_wday
temp$eating_midpoint_wend <- (temp$eating_window_wend)/2 + temp$first_meal_wend

# Step 3: Get Eating Jetlag
# Eating JetLag (h) = Eating midpoint on weekends − Eating midpoint on weekdays.

temp$eating_jetlag <- temp$eating_midpoint_wend - temp$eating_midpoint_wday


# Add mealtime variable to metadata 
df <- df[,!names(df) %in% mealtime_vars] # remove original mealtime vars


df <- merge(df, temp, by = 'Basename', all = TRUE) 

#------------------------------------------------------------------------------#
#                          Sleep Midpoint                                      #----
#------------------------------------------------------------------------------#

# Note: Variables shold be in numeric format (already are in here, but if 
# necessary check beforehand)

# Sleep midpoint = time midway between sleep onset and offset.  
# Ref: https://pmc.ncbi.nlm.nih.gov/articles/PMC7969471/ 

# Note: sleep duration is in duartion fortmat: "31114.2857142857s (~8.64 hours)" 
# We need it as numeric and in hours (/3600)
# Now apply the conditional logic to calculate the average sleep duration
df$sleep_midpoint_imp <- duration(as.numeric(df$sleep_duration)/3600/2,units = "hours")

# Sleep duration midpoint for "31114.2857142857s (~8.64 hours)"  is = 4.321429 


#------------------------------------------------------------------------------#
#                 Compute Time Between Last Meal and Going to Bed              #----
#------------------------------------------------------------------------------#

df$fasting_before_sleep_wday <- get_fasting_before_sleep(df$last_meal_wday, df$gobed_wday)

df$fasting_before_sleep_wend <- get_fasting_before_sleep(df$last_meal_wend, df$gobed_wend)


# Get average
df$fasting_before_sleep <- dplyr::if_else(!is.na(df$fasting_before_sleep_wday) & !is.na(df$fasting_before_sleep_wend),   # If we have both wday and wend
                                         ((df$fasting_before_sleep_wday * 5) + (df$fasting_before_sleep_wend * 2)) / 7, # Use both to calculate the average
                                         dplyr::if_else(!is.na(df$fasting_before_sleep_wday) & is.na(df$fasting_before_sleep_wend),  # If we have wday but not wend
                                                        df$fasting_before_sleep_wday,                                       # Use wday
                                                        dplyr::if_else(is.na(df$fasting_before_sleep_wday) & !is.na(df$fasting_before_sleep_wend), # If we have wend but not wday
                                                                       df$fasting_before_sleep_wend,                                     # Use wend
                                                                       NA)))         

df$fasting_before_sleep <- duration(df$fasting_before_sleep,units = "minutes")
df$fasting_before_sleep_wday <- duration(df$fasting_before_sleep,units = "hours")
df$fasting_before_sleep_wend <- duration(df$fasting_before_sleep,units = "hours")

df$fasting_before_sleep[1:100]

df$fasting_before_sleep_wday[1:100]


# Note: Sleep duration was imputed using the "raw" (cleaned) non-imputed wake up and go bed variables. 
# This means we have complete cases for sleep duration. Because the fasting is derived from the wake up/go bed
# variables, we do have missing data in fasting_before_sleep. 

# Option 1: Impute raw variables, wake up and go bed. With the imputed values compute sleep duration, 
# psqi_score (which needs in one of the item sleep duration), and fasting 

# Option 2 (what we're doing for now): Keep imputed sleep duration. We lose N = 50 for the fasting. 

sum(is.na(df$sleep_duration)) # 0
sum(is.na(df$sleep_midpoint)) # 0
sum(is.na(df$fasting_before_sleep)) # 48
sum(is.na(df$fasting_before_sleep_wday)) # 48
sum(is.na(df$fasting_before_sleep_wend)) # 52



#------------------------------------------------------------------------------#
#                 Compute Time Between Waking Up and First Meal                #----
#------------------------------------------------------------------------------#

get_fasting_after_sleep <- function(wake_up, first_meal) {
  # Create an empty vector to store fasting hours
  fasting_hours <- numeric(length(first_meal))
  
  for (i in seq_along(first_meal)) {
    wake_time <- wake_up[i]
    meal_time <- first_meal[i]

    # Case 1: If either value is NA, return NA
    if (is.na(meal_time) | is.na(wake_time)) {
      fasting_hours[i] <- NA
      next
    }
    
    # Convert duration to numeric hours
    meal_hours <- as.numeric(meal_time, "hours")
    wake_hours <- as.numeric(wake_time, "hours")
    
    
    
    # Case 2: Adjust for bedtime after midnight
    # if (meal_hours == 0 && bed_hours > 20) {  # Midnight meal but bedtime is late night
    #   fasting_hours[i] <- 0
    if (meal_hours < wake_hours) {
      fasting_hours[i] <- 0  # Case : First meal is before wake up time - set up to 0 fasting hours 
    } else {
    fasting_hours[i] <- meal_hours - wake_hours
    }
  }
  
  return(fasting_hours)
}


df$fasting_after_sleep_wday <- get_fasting_after_sleep(df$wakeup_wday, df$first_meal_wday)
df$fasting_after_sleep_wend <- get_fasting_after_sleep(df$wakeup_wend, df$first_meal_wend)


# Get average
df$fasting_after_sleep <- dplyr::if_else(!is.na(df$fasting_after_sleep_wday) & !is.na(df$fasting_after_sleep_wend),   # If we have both wday and wend
                                          ((df$fasting_after_sleep_wday * 5) + (df$fasting_after_sleep_wend * 2)) / 7, # Use both to calculate the average
                                          dplyr::if_else(!is.na(df$fasting_after_sleep_wday) & is.na(df$fasting_after_sleep_wend),  # If we have wday but not wend
                                                         df$fasting_after_sleep_wday,                                       # Use wday
                                                         dplyr::if_else(is.na(df$fasting_after_sleep_wday) & !is.na(df$fasting_after_sleep_wend), # If we have wend but not wday
                                                                        df$fasting_after_sleep_wend,                                     # Use wend
                                                                        NA)))         

df$fasting_after_sleep <- duration(df$fasting_after_sleep,units = "minutes")
df$fasting_after_sleep_wday <- duration(df$fasting_after_sleep,units = "hours")
df$fasting_after_sleep_wend <- duration(df$fasting_after_sleep,units = "hours")


sum(is.na(df$fasting_after_sleep)) # 48
sum(is.na(df$fasting_after_sleep_wday)) # 48
sum(is.na(df$fasting_after_sleep_wend)) # 55


#------------------------------------------------------------------------------#
#                     Get  early - non early eaters                            #----
#------------------------------------------------------------------------------#

# Get first meal into 2c : early and non early eaters  
median_first_meal <- quantile(df$first_meal, probs = 0.5, na.rm = TRUE)
# median is 8.64 (numeric) --> 8:38h 
df$first_meal_2c <- ifelse(df$first_meal < median_first_meal, "early", "non_early")

#------------------------------------------------------------------------------#
#                     Get  short - long fasting                                #----
#------------------------------------------------------------------------------#

# Get nighttime fasting into 2c : short and long fasting  
median_nighttime_fasting <- quantile(df$nighttime_fasting, probs = 0.5, na.rm = TRUE)

# median is  11.57 hours (numeric) --> 11:34h
df$nighttime_fasting_2c <- ifelse(df$nighttime_fasting < median_nighttime_fasting, "short", "long")


#==============================================================================#
#                                                                              #
#                         SAVE DATA WITH SCORES                                #----
#                                                                              #
#==============================================================================#


save_to <- paste0("db/metadata_plaEWAS_MT_incomplete_cases2_with_scores_",current_date,".rds")
saveRDS(df, save_to) 
colnames(df)
# [1] "Basename"                   "maternal_age"               "maternal_bmi"              
# [4] "maternal_educational_level" "maternal_smoking"           "parity"                    
# [7] "maternal_physical_activity" "gestational_age"            "Sex"                       
# [10] "dna_contamination"          "study_center"               "gestational_diabetes"      
# [13] "pregnancy_complications"    "covid_confinement"          "child_ethnicity"           
# [16] "sleep_duration"             "psqi_score"                 "first_meal"                
# [19] "last_meal"                  "eating_occasions"           "cereals"                   
# [22] "fruit_nuts"                 "meat"                       "legumes"                   
# [25] "fish"                       "legumes.1"                  "olive_oil"                 
# [28] "vegetables"                 "dairy"                      "energy_intake"             
# [31] "Trophoblasts"               "Stromal"                    "Hofbauer"                  
# [34] "Endothelial"                "nRBC"                       "Syncytiotrophoblast"       
# [37] "wakeup_wday"                "wakeup_wend"                "gobed_wday"                
# [40] "gobed_wend"                 "sleep_duration_wday"        "sleep_duration_wend"       
# [43] "gest_weight_gain_m_preg"    "rMEDp"                      "first_meal_wday"           
# [46] "last_meal_wday"             "first_meal_wend"            "last_meal_wend"            
# [49] "eating_occasions_wday"      "eating_occasions_wend"      "nighttime_fasting_wday"    
# [52] "nighttime_fasting_wend"     "nighttime_fasting"          "eating_window_wday"        
# [55] "eating_window_wend"         "eating_midpoint_wday"       "eating_midpoint_wend"      
# [58] "eating_jetlag"              "sleep_midpoint"             "fasting_before_sleep_wday" 
# [61] "fasting_before_sleep_wend"  "fasting_before_sleep"       "fasting_after_sleep_wday"  
# [64] "fasting_after_sleep_wend"   "fasting_after_sleep"        "first_meal_2c"             


