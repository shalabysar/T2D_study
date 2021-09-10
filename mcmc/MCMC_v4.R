library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(ggdendro)
library(cluster)
library(purrr)
library(tibble)
library(ggradar)
library(naniar)
library(lubridate)
library(fossil)
library(matrixStats)
library(kableExtra)
library(pastecs)

###############################################################################
################################ TABLES LOADING ###############################
###############################################################################

#load tables
setwd("~/ukbb_diabetes/msm/obs_selection")

df_main <- fread("~/ukbb_diabetes/data/obs_selected/df_main_msm_sel.csv")
df_bio <- fread("~/ukbb_diabetes/data/obs_selected/df_bio_msm_sel.csv")
df_genes <- fread("~/ukbb_diabetes/data/obs_selected/df_genes_msm_sel.csv")

df_main <- as.data.frame(df_main)
df_bio <- as.data.frame(df_bio)
df_genes <- as.data.frame(df_genes)

###############################################################################
################################# GET FUNCTIONS ###############################
###############################################################################
source("~/ukbb_diabetes/functions.R")

###############################################################################
############################### PATIENTS SELECTION ############################
###############################################################################

##################### 1. SELECT PRE-DIABETIC PATIENTS AT INITIAL
#####################################################

#eid of pre-diabetes at initial
eid_diabetes_initial <- df_main %>% 
  filter(status_init == "prediabetes") %>% 
  select(eid) %>% 
  unlist()

#dataframe of pre-diabetes at initial
df_main_sel <- df_main %>% 
  filter(eid %in% eid_diabetes_initial)

df_bio_sel <- df_bio %>% 
  filter(eid %in% eid_diabetes_initial)

df_genes_sel <- df_genes %>% 
  filter(eid %in% eid_diabetes_initial)


##################### 2. REMOVE INCONSISTENT PEOPLE
#####################################################

#eid of age_diabetes_diagnosis < age at initial assessment
eid_inconsistent <- df_main_sel %>% 
  filter(`21003-0.0` > diabetes_age_diag_fu ) %>% 
  select(eid) %>% 
  unlist()

#dataframe cleaned
df_main_sel <- df_main_sel %>% 
  filter(!(eid %in% eid_inconsistent))

df_bio_sel <- df_bio_sel %>% 
  filter(!(eid %in% eid_inconsistent))

df_genes_sel <- df_genes_sel %>% 
  filter(!(eid %in% eid_inconsistent))


###############################################################################
################################# MATRICES INIT ###############################
###############################################################################

##################### 1. PREPARE TABLES
#####################################################

#variables
var_main_init <- c("21001-0.0", "4079-0.0", "4080-0.0", "diabetes_age_diag_init", "status_init", 
                   "diabetes_class_source_init")
var_bio_init <- c("30750-0.0", "30760-0.0", "30780-0.0", "30870-0.0")
var_main_fu <- c("21001-1.0", "4079-1.0", "4080-1.0", "diabetes_age_diag_fu", "status_fu", 
                 "diabetes_class_source_fu")
var_bio_fu <- c("30750-1.0", "30760-1.0", "30780-1.0", "30870-1.0")
var_genes <- c("score_secr_weighted", "score_res_weighted")

#create table joined
df_all_sel_init <- prepare_tables_3df(
  df_main_sel, c("eid", var_main_init), 
  df_bio_sel, c("eid", var_bio_init),
  df_genes_sel,c("eid", var_genes),
  "eid")

df_all_sel_fu <- prepare_tables_3df(
  df_main_sel, c("eid", var_main_fu), 
  df_bio_sel, c("eid", var_bio_fu),
  df_genes_sel,c("eid", var_genes),
  "eid")

#rename columns
colnames(df_all_sel_init) <- 
  c("eid", "bmi", "dbp", "sbp", "diabetes_age_diag", "status", "diabetes_class_source",
    "hba1c", "hdl_chol", "ldl_direct","triglycerides", 
    "pgs_secr_w", "pgs_res_w")
colnames(df_all_sel_fu) <- 
  c("eid", "bmi", "dbp", "sbp", "diabetes_age_diag", "status", "diabetes_class_source",
    "hba1c", "hdl_chol", "ldl_direct","triglycerides", 
    "pgs_secr_w", "pgs_res_w")

#complete_cases
var_compl_cases <- 
  c("bmi", "dbp", "sbp", "hba1c", "hdl_chol", "ldl_direct", "triglycerides", 
    "pgs_secr_w", "pgs_res_w")
compl_cases_init <- df_all_sel_init[complete.cases(df_all_sel_init[ ,var_compl_cases]),
                                    "eid"]
compl_cases_fu <- df_all_sel_fu[complete.cases(df_all_sel_fu[ ,var_compl_cases]),
                                "eid"]
compl_cases <-  intersect(compl_cases_init, compl_cases_fu)

#table with complete cases
df_all_sel_init <- df_all_sel_init %>% 
  filter(eid %in% compl_cases)
df_all_sel_fu <- df_all_sel_fu %>% 
  filter(eid %in% compl_cases)

##################### 2. ASSESSMENTS DATES
#####################################################

#Year and Date of Initial assessment
DOA1 <- df_main_sel %>% 
  filter(eid %in% compl_cases) %>% 
  select(`53-0.0`)

YOA1 <- year(ymd(DOA1$`53-0.0`))

#Year and Date of FU assessment
DOA2 <- df_main_sel %>% 
  filter(eid %in% compl_cases) %>% 
  select(`53-1.0`)

YOA2 <- year(ymd(DOA2$`53-1.0`))

#First and last year
first_year <- min(YOA1)
last_year <- max(YOA2)

#Year of follow-up
n_year_fu <- last_year - first_year + 1


##################### 3. DATE OF BIRTH
#####################################################

#Year of birth
YOB <- df_main_sel %>% 
  filter(eid %in% compl_cases) %>% 
  select(`34-0.0`) %>% 
  unlist()

#Month of birth
MOB <- df_main_sel %>% 
  filter(eid %in% compl_cases) %>% 
  select(`52-0.0`) %>% 
  unlist()


##################### 4. INITIALIZE MATRICES - ALL BUT AGE
#####################################################

#Variables selected
var_selected <- c("bmi", "dbp", "sbp", "hba1c", "hdl_chol", "ldl_direct","triglycerides", 
                  "pgs_secr_w", "pgs_res_w")

#create matrices 
matrices <- lapply(1:length(var_selected), matrix, data= NA, 
                   nrow = nrow(df_all_sel_init), ncol = n_year_fu)

names(matrices) <- var_selected

for (ind in 1:nrow(df_all_sel_init)){
  #obtain YOA
  YOA1_ind <- YOA1[ind]
  YOA2_ind <- YOA2[ind]
  
  #obtain position of YOA
  YOA1_pos <- YOA1_ind - first_year + 1
  YOA2_pos <- YOA2_ind - first_year + 1
  
  #fill matrices
  for (col in 1:length(var_selected)){
    
    #value at A1
    matrices[[col]][ind, YOA1_pos] <- df_all_sel_init[ind, var_selected[col]]
    #value at A2
    matrices[[col]][ind, YOA2_pos] <- df_all_sel_fu[ind, var_selected[col]]
  }
  
}

##################### 5. OBTAIN AGE
#####################################################
mat_age <- matrix(NA, nrow = nrow(df_all_sel_init), ncol = n_year_fu)

for (ind in 1:nrow(df_all_sel_init)){
  #obtain YOA
  YOA1_ind <- YOA1[ind]
  YOA2_ind <- YOA2[ind]
  
  #obtain position of YOA
  YOA1_pos <- YOA1_ind - first_year + 1
  YOA2_pos <- YOA2_ind - first_year + 1
  
  #fill matrix
  for (year in YOA1_pos:YOA2_pos){
    mat_age[ind, year] <- first_year - YOB[ind] + year - 1
  }
}

###############################################################################
################################# INTERPOLATIONS ##############################
###############################################################################

##################### 1. INTERPOLATE ALL BUT AGES
#####################################################

mat_interpolated <- lapply(matrices, interpolate_linear)


###############################################################################
################################# DATE OUT OF PD ##############################
###############################################################################

##################### 1. DATE OF DIABETES DIAGNOSIS
#####################################################

#initialize
YOD <- rep(NA, nrow(df_all_sel_fu))

#Age of diagnosis from main matrix
AOD <- df_all_sel_fu$diabetes_age_diag

#for patients identified with quest use age of diag
YOD[which(df_all_sel_fu$diabetes_class_source == "quest")] <-
  YOB[which(df_all_sel_fu$diabetes_class_source == "quest")] + 
  AOD[which(df_all_sel_fu$diabetes_class_source == "quest")]

#for patients identified with verbal use age of diag
YOD[which(df_all_sel_fu$diabetes_class_source == "verbal")] <-
  YOB[which(df_all_sel_fu$diabetes_class_source == "verbal")] + 
  AOD[which(df_all_sel_fu$diabetes_class_source == "verbal")]

#for patients born in H2 add +1 to YOD
for (ind in 1:length(MOB)){
  if (MOB[ind] > 6){
    YOD[ind] <- YOD[ind] + 1
  }
}


#for patients identified with hba1c use first year hba1c become > 48
for (ind in 1:nrow(df_all_sel_fu)){
  if(is.na(df_all_sel_fu$diabetes_class_source[ind])){
    
  } else if (df_all_sel_fu$diabetes_class_source[ind] == "hba1c"){
    pos_diabetes <- min(which(mat_interpolated$hba1c[ind, ] > 48))
    YOD[ind] <- first_year + pos_diabetes - 1
  }
}

#check no patients with YOD<YOA1 - ok
sum(YOD<YOA1, na.rm = TRUE)
for (ind in 1:length(YOD)){
  if (!is.na(YOD[ind])){
    if (YOD[ind]<YOA1[ind]){
      YOD[ind] <- YOD[ind] + 1
    }
  }
}

#check no patients with YOD>YOA2 - 1 obs.
sum(YOD>YOA2, na.rm = TRUE)
for (ind in 1:length(YOD)){
  if (!is.na(YOD[ind])){
    if (YOD[ind]>YOA2[ind]){
      YOD[ind] <- YOD[ind] - 1
    }
  }
}


##################### 2. DATE OF HEALTHY DIAGNOSIS
#####################################################

#initialize
YOR <- rep(NA, nrow(df_all_sel_fu))

#use first year hba1c become <39
for (ind in 1:nrow(df_all_sel_fu)){
  if (df_all_sel_fu$status[ind] == "healthy"){
    pos_healthy <- min(which(mat_interpolated$hba1c[ind, ] < 39))
    YOR[ind] <- first_year + pos_healthy - 1
  }
}

###############################################################################
############################### PREPARE MATRICES ##############################
###############################################################################

##################### 1. CREATE LOGS
#####################################################

#age
mat_age_log <- log(mat_age)

#all others
mat_interpolated_log <- lapply(mat_interpolated, log)


##################### 2. SCALE
#####################################################

#age
mat_age_scaled <- matrix(scale(as.vector(mat_age)),
                            ncol = ncol(mat_age))

#all others
mat_interpolated_scaled <- lapply(mat_interpolated, scale_matrix)

##################### 2. CLEAN MATRICES OF VARIABLES
#####################################################

#age
mat_age[is.na(mat_age)] <- 0
mat_age_log[is.na(mat_age_log)] <- 0
mat_age_scaled[is.na(mat_age_scaled)] <- 0

#all others
for(i in 1:length(var_selected)){
  mat_interpolated[[i]][is.na(mat_interpolated[[i]])] <- 0
}

for(i in 1:length(var_selected)){
  mat_interpolated_log[[i]][is.na(mat_interpolated_log[[i]])] <- 0
}

for(i in 1:length(var_selected)){
  mat_interpolated_scaled[[i]][is.na(mat_interpolated_scaled[[i]])] <- 0
}


##################### 3. CLEAN DATES
#####################################################
#deal with missing values
YOD[is.na(YOD)] <- 1900
YOR[is.na(YOR)] <- 1900

#create matrix
dates <- data.frame(YOA1, YOA2, YOD, YOR)

##################### 3. INFO ON PATIENTS
#####################################################

df_all_sel_fu %>% 
  filter(status == "diabetes") %>% 
  group_by(diabetes_class_source) %>% 
  summarise(tot=n())

df_all_sel_fu %>% 
  group_by(status) %>% 
  summarise(tot=n())

generate_sum_stat(df_all_sel_fu %>% 
                    filter(status == "healthy") %>% 
                    as.data.frame(), 
                  colnames(df_all_sel_fu)[-c(1, 3, 5, 6, 7, 8, 10)], 
                  colnames(df_all_sel_fu)[-c(1, 3, 5, 6, 7, 8, 10)])

generate_sum_stat(df_all_sel_fu %>% 
                    filter(status == "prediabetes") %>% 
                    as.data.frame(), 
                  colnames(df_all_sel_fu)[-c(1, 3, 5, 6, 7, 8, 10)], 
                  colnames(df_all_sel_fu)[-c(1, 3, 5, 6, 7, 8, 10)])

generate_sum_stat(df_all_sel_fu %>% 
                    filter(status == "diabetes") %>% 
                    as.data.frame(), 
                  colnames(df_all_sel_fu)[-c(1, 3, 5, 6, 7, 8, 10)], 
                  colnames(df_all_sel_fu)[-c(1, 3, 5, 6, 7, 8, 10)])


##################### 4. STORE
#####################################################
out_dir <- "~/ukbb_diabetes/msm/MCMC_V4/Data/"
simul_name <- "_simul_V4"

save_for_mcmc(dates, out_dir, "dates", simul_name)
save_for_mcmc(mat_interpolated_scaled$bmi, out_dir, "bmi", simul_name)
save_for_mcmc(mat_interpolated_scaled$sbp, out_dir, "sbp", simul_name)
save_for_mcmc(mat_interpolated_scaled$hdl_chol, out_dir, "hdl", simul_name)
save_for_mcmc(mat_interpolated_scaled$triglycerides, out_dir, "trigl", simul_name)
save_for_mcmc(mat_interpolated_scaled$pgs_secr_w, out_dir, "pgs_is", simul_name)
save_for_mcmc(mat_interpolated_scaled$pgs_res_w, out_dir, "pgs_ir", simul_name)
save_for_mcmc(mat_age_scaled, out_dir, "age", simul_name)
