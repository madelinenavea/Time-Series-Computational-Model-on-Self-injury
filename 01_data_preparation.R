## ----setup------------------------------------------------------------------------------------------------------------------------
# Load Libraries
library(magrittr)
library(tidyverse)
library(lme4)
library(lmerTest) # Adds significance asterices

# Unneeded packages
# library(gridExtra)
# library(ggplot2)
# library(tibble)
# library(purrr)

# Person-centered means function from Kleiman's EMAtools package (deprecated)
pcenter<-function(ID,var){
  centered<- var-ave(var, ID,FUN=function(x) mean(x, na.rm=T))
  return(centered)
}


## ----data_prep--------------------------------------------------------------------------------------------------------------------
# Read raw df
raw_data <- read.csv("Raw data/CompleteDay3Up.csv")

# Select required variables
selected <- raw_data %>% 
  dplyr::select(
    PersonID,                                                                   # Identifier
    Male, Ethnicity, Race, Height_totalin, Weight, BMI, Age, Gen,               # Demographics
    SessionTime, studydays, EMA_instances, EMA_total,                           # EMA Metadata
    U_Binge, U_Vomit, U_Laxative, U_Exercise, U_Fast, U_Restrict, U_NSSI, U_NA, # Binary SIU
    Distress, Upset, Nervous, Scared, Afraid,                                   # Negative Affect
    Excited, Inspired, Determined, Enthusiast, Alert,                           # Positive Affect
    ends_with("UrgInt"),                                                        # SIU intensity
    SIB_Binge, SIB_Vomit, SIB_Laxative, SIB_Exercise, SIB_Fast, SIB_Restrict,   # Binary SIB ED
    SIB_NSSI) %>%                                                               # Binary SIB NSSI
    dplyr::rename(day_number = studydays) %>% 
    mutate_all(~ ifelse(. == -999, NA, .),
             ~ ifelse(. == "NaN", NA, .)) 
  
df <- selected %>% 
# Create composite items for all variables at all observations
  mutate(
    SIU_intensity_t0 = rowMeans(select(., c(BinUrgInt, VomUrgInt, LaxUrgInt, ExUrgInt, FastUrgInt, ResUrgInt, NSSIUrgInt)), na.rm=TRUE),
    SIU_count_t0 = rowSums(dplyr::select(., c(U_Binge, U_Vomit, U_Laxative, U_Exercise, U_Fast, U_Restrict, U_NSSI)), na.rm = TRUE),
    SIB_count_t0 = rowSums(dplyr::select(., c(SIB_Binge, SIB_Vomit, SIB_Laxative, SIB_Exercise, SIB_Fast, SIB_Restrict, SIB_NSSI)), na.rm = TRUE),
    NA_intensity_t0 = rowMeans(select(., c(Distress, Upset, Nervous, Scared, Afraid)), na.rm = TRUE),
    PA_intensity_t0 = rowMeans(select(., c(Excited, Inspired, Determined, Enthusiast, Alert)), na.rm = TRUE),
    
    # log NA due to zero inflaction
    NA_intensity_t0_log = log(NA_intensity_t0+1),
    
    # Create binary variables for SIBs and urges
    SIB_binary_t0 = ifelse(SIB_count_t0 > 0, 1, 0),
    SIU_binary_t0 = ifelse(SIU_count_t0 > 0, 1, 0))

# Create meta data
# Use SessionTime to count off observations
# Convert SessionTime to POSIXct
df$SessionTime <- as.POSIXct(df$SessionTime, format = "%m/%d/%y %H:%M")
  
# Group by ID and count # sessions
df %<>%
  arrange(PersonID, day_number, SessionTime) %>%          # Arrange data for counting
  dplyr::group_by(PersonID) %>%                           # Group by PersonID
  dplyr::mutate(observation_number = row_number())        # Count observations within each group

# Calculate maximum day and observation number for each PersonID
max_num <- df %>%
  dplyr::group_by(PersonID) %>%
  dplyr::summarise(day_total = max(day_number, na.rm = TRUE),
            observation_total = max(observation_number, na.rm = TRUE))

# Join day_total into df
df %<>% dplyr::left_join(max_num, by = "PersonID")

# Person-centered means for the predictor vars
df %<>%
  mutate(
    SIU_intensity_t0_c = pcenter(PersonID, SIU_intensity_t0),
    NA_intensity_t0_c  = pcenter(PersonID, NA_intensity_t0_log),
    PA_intensity_t0_c  = pcenter(PersonID, PA_intensity_t0))

# Leading the outcome variables to create later timepoint variables
df %<>%
  group_by(PersonID, day_number) %>%  
  mutate(
    SIB_binary_t1 = lead(SIB_binary_t0, n = 1),
    SIB_binary_t2 = lead(SIB_binary_t0, n = 2),
    SIB_count_t1 = lead(SIB_count_t0, n = 1),
    SIB_count_t2 = lead(SIB_count_t0, n = 2),
    
    SIU_binary_t1 = lead(SIU_binary_t0, n = 1),
    SIU_count_t1 = lead(SIU_count_t0, n = 1),
    SIU_intensity_t1_c = lead(SIU_intensity_t0_c, n = 1),
    
    NA_intensity_t1_c = lead(NA_intensity_t0_c, n = 1),
    PA_intensity_t1_c = lead(PA_intensity_t0_c, n = 1)) %>% 
  ungroup()

# Creating the prepped_data df to run the models in Mplus
prepped_data <- df %>%
  dplyr::select(
    PersonID, observation_number, observation_total, day_number, day_total, 
    SIU_count_t0, SIU_count_t1, SIU_intensity_t0_c, SIU_intensity_t1_c, SIU_binary_t0, SIU_binary_t1, 
    SIB_binary_t0, SIB_binary_t1, SIB_binary_t2, SIB_count_t0, SIB_count_t2,
    NA_intensity_t0_c, PA_intensity_t0_c, NA_intensity_t1_c, PA_intensity_t1_c)

# write.csv(prepped_data, "Outputs/prepped_data.csv", row.names = F) #write to csv for potential later usage
# prepped_data[is.na(prepped_data)] <- 999  #for Mplus to understand 'NAs'
# write.table(prepped_data, "prepped_data.dat", row.names = F, sep = ",") #for Mplus to read the file


## ---------------------------------------------------------------------------------------------------------------------------------
#### SIU_int0 -> SIB_bin1
SIU_int0_to_SIB_bin1 <- glmer(SIB_binary_t1 ~ SIU_intensity_t0_c + (1|PersonID),
data = prepped_data, family = binomial(link = "logit"))
summary(SIU_int0_to_SIB_bin1)
exp(fixef(SIU_int0_to_SIB_bin1))
exp(confint(SIU_int0_to_SIB_bin1, method = "Wald"))

#### SIU_bin0 -> SIB_bin1
SIU_bin0_to_SIB_bin1 <- glmer(SIB_binary_t1 ~ SIU_binary_t0 + (1|PersonID), # remove random slopes for low variability
data = prepped_data, family = binomial(link = "logit"))
summary(SIU_bin0_to_SIB_bin1)
exp(fixef(SIU_bin0_to_SIB_bin1))
exp(confint(SIU_bin0_to_SIB_bin1, method = "Wald"))


#? SIU_int0 -> SIB_bin2
SIU_int0_to_SIB_bin2 <- glmer(SIB_binary_t2 ~ SIU_intensity_t0_c + (1|PersonID), # remove random slopes for low variability
data = prepped_data, family = binomial(link = "logit"))
summary(SIU_int0_to_SIB_bin2)
exp(fixef(SIU_int0_to_SIB_bin2))
exp(confint(SIU_int0_to_SIB_bin2, method = "Wald"))

## SIU_bin0 -> SIB_bin2
SIU_bin0_to_SIB_bin2 <- glmer(SIB_binary_t2 ~ SIU_binary_t0 + (1|PersonID), #**better w random slopes
data = prepped_data, family = binomial(link = "logit"))
summary(SIU_bin0_to_SIB_bin2)
exp(fixef(SIU_bin0_to_SIB_bin2))
exp(confint(SIU_bin0_to_SIB_bin2, method = "Wald"))


### SIU_bin1 -> SIB_bin2
SIU_bin1_to_SIB_bin2 <- glmer(SIB_binary_t2 ~ SIU_binary_t1 + (1|PersonID),
data = prepped_data, family = binomial(link = "logit"))
summary(SIU_bin1_to_SIB_bin2)
exp(fixef(SIU_bin1_to_SIB_bin2))
exp(confint(SIU_bin1_to_SIB_bin2, method = "Wald"))

### SIU_int1 -> SIB_bin2
SIU_int1_to_SIB_bin2 <- glmer(SIB_binary_t2 ~ SIU_intensity_t1_c + (1|PersonID), # remove random slopes for low variability
data = prepped_data, family = binomial(link = "logit"))
summary(SIU_int1_to_SIB_bin2)
exp(fixef(SIU_int1_to_SIB_bin2))
exp(confint(SIU_int1_to_SIB_bin2, method = "Wald"))


## ---------------------------------------------------------------------------------------------------------------------------------
#X NA_int0 -> SIB_bin1
NA_to_SIB_bin1 <- glmer(SIB_binary_t1 ~ NA_intensity_t0_c + (1|PersonID),
data = prepped_data, family = binomial(link = "logit"))
summary(NA_to_SIB_bin1)
exp(fixef(NA_to_SIB_bin1))
exp(confint(NA_to_SIB_bin1, method = "Wald"))

### NA_int0 -> SIB_bin2
NA_to_SIB_bin2 <- glmer(SIB_binary_t2 ~ NA_intensity_t0_c + (1|PersonID), #***with random slopes
data = prepped_data, family = binomial(link = "logit"))
summary(NA_to_SIB_bin2)
exp(fixef(NA_to_SIB_bin2))
exp(confint(NA_to_SIB_bin2, method = "Wald"))


# NA AND URGES
#### NA_int0 -> SIU_int1
NA_to_SIU_int1 <- lmer(SIU_intensity_t1 ~ NA_intensity_t0_c + (1|PersonID),
data = prepped_data, control = lmerControl(optimizer = "bobyqa"))
summary(NA_to_SIU_int1)
exp(fixef(NA_to_SIU_int1))
exp(confint(NA_to_SIU_int1, method = "Wald"))

#X NA_int0 -> SIU_binary_t1
NA_to_SIU_bin1 <- glmer(SIU_binary_t1 ~ NA_intensity_t0_c + (1|PersonID),
data = prepped_data, family = binomial(link = "logit"))
summary(NA_to_SIU_bin1)
exp(fixef(NA_to_SIU_bin1))
exp(confint(NA_to_SIU_bin1, method = "Wald"))



## ---------------------------------------------------------------------------------------------------------------------------------
#X PA_int0 -> SIB_bin1
PA_to_SIB_bin1 <- glmer(SIB_binary_t1 ~ PA_intensity_t0_c + (1|PersonID),
data = prepped_data, family = binomial(link = "logit"))
summary(PA_to_SIB_bin1)
exp(fixef(PA_to_SIB_bin1))
exp(confint(PA_to_SIB_bin1, method = "Wald"))

#X PA_int0 -> SIB_bin2
PA_to_SIB_bin2 <- glmer(SIB_binary_t2 ~ PA_intensity_t0_c + (1|PersonID),
data = prepped_data, family = binomial(link = "logit"))
summary(PA_to_SIB_bin2)
exp(fixef(PA_to_SIB_bin2))
exp(confint(PA_to_SIB_bin2, method = "Wald"))


# NA AND URGES
#X PA_int0 -> SIU_int1
PA_to_SIU_int1 <- lmer(SIU_intensity_t1 ~ PA_intensity_t0_c + (1|PersonID), # remove random slopes for low variability
data = prepped_data, control = lmerControl(optimizer = "bobyqa"))
summary(PA_to_SIU_int1)
exp(fixef(PA_to_SIU_int1))
exp(confint(PA_to_SIU_int1, method = "Wald"))

#X PA_int0 -> SIU_binary_t1
PA_to_SIU_bin1 <- glmer(SIU_binary_t1 ~ PA_intensity_t0_c + (1|PersonID),
data = prepped_data, family = binomial(link = "logit"))
summary(PA_to_SIU_bin1)
exp(fixef(PA_to_SIU_bin1))
exp(confint(PA_to_SIU_bin1, method = "Wald"))

# table_summary_MLMsforDSEM <- modelsummary(list(
#  "PA_to_SIU_int1" = PA_to_SIU_int1), output = "Outputs/table_summary_MLMsforDSEM.csv")


## ---------------------------------------------------------------------------------------------------------------------------------
#### SIU_binary_t0 -> NA_intensity_t1
SIU_bin0_to_NA1 <- lmer(NA_intensity_t1 ~ SIU_binary_t0 + (1|PersonID),
data = prepped_data)
summary(SIU_bin0_to_NA1)
exp(fixef(SIU_bin0_to_NA1))
exp(confint(SIU_bin0_to_NA1, method = "Wald"))

#### SIU_int_t0 -> NA_intensity_t1
SIU_int0_to_NA1 <- lmer(NA_intensity_t1 ~ SIU_intensity_t0_c + (1|PersonID),
data = prepped_data)
summary(SIU_int0_to_NA1)
exp(fixef(SIU_int0_to_NA1))
exp(confint(SIU_int0_to_NA1, method = "Wald"))

## SIU_binary_t0 -> PA_intensity_t1
SIU_bin0_to_PA1 <- lmer(PA_intensity_t1 ~ SIU_binary_t0 + (1|PersonID),
data = prepped_data)
summary(SIU_bin0_to_PA1)
exp(fixef(SIU_bin0_to_PA1))
exp(confint(SIU_bin0_to_PA1, method = "Wald"))

### SIU_int_t0 -> PA_intensity_t1
SIU_int0_to_PA1 <- lmer(PA_intensity_t1 ~ SIU_intensity_t1_c + (1|PersonID),
data = prepped_data)
summary(SIU_int0_to_PA1)
exp(fixef(SIU_int0_to_PA1))
exp(confint(SIU_int0_to_PA1, method = "Wald"))

