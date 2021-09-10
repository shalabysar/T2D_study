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
library(flexclust)

###############################################################################
################################ TABLES LOADING ###############################
###############################################################################

#load tables - XXX to be filled as required
setwd("~/XXX")

df_main <- fread("~/XXX/df_main_msm_sel.csv")
df_bio <- fread("~/XXX/df_bio_msm_sel.csv")
df_genes <- fread("~/XXX/df_genes_msm_sel.csv")
df_hes <- fread("~/XXX/df_hes_msm_sel.csv")

df_main <- as.data.frame(df_main)
df_bio <- as.data.frame(df_bio)
df_genes <- as.data.frame(df_genes)
df_hes <- as.data.frame(df_hes)

#load model
load("~/XXX/kmeans_3.RData")

kcca_model <- kcca_kmeans_3
cluster_names_ordered <- cluster_names_ordered_kmeans_3
cluster_ordered <- cluster_ordered_kmeans_3


###############################################################################
################################# GET FUNCTIONS ###############################
###############################################################################
source("~/functions.R")

###############################################################################
############################### PATIENTS SELECTION ############################
###############################################################################

##################### 1. SELECT DIABETIC PATIENTS AT INITIAL
#####################################################

#eid of diabetes at initial
eid_diabetes_initial <- df_main %>% 
  filter(status_init == "diabetes") %>% 
  filter(diabetes_class_source_init != "hba1c") %>% 
  filter(status_fu == "diabetes") %>% 
  dplyr::select(eid) %>% 
  unlist()

#dataframe of diabetes at initial
df_main_sel <- df_main %>% 
  filter(eid %in% eid_diabetes_initial)

df_bio_sel <- df_bio %>% 
  filter(eid %in% eid_diabetes_initial)

df_genes_sel <- df_genes %>% 
  filter(eid %in% eid_diabetes_initial)

df_hes_sel <- df_hes %>% 
  filter(eid %in% eid_diabetes_initial)

###############################################################################
################################# PREPARE DATA ################################
###############################################################################

##################### 1. PREPARE TABLES
#####################################################

#variables
var_main_init <- c("21001-0.0", "4080-0.0", "diabetes_age_diag_init", "status_init", 
                   "diabetes_class_source_init")
var_bio_init <- c("30750-0.0", "30760-0.0", "30870-0.0")
var_main_fu <- c("21001-1.0", "4080-1.0", "diabetes_age_diag_fu", "status_fu", 
                 "diabetes_class_source_fu")
var_bio_fu <- c("30750-1.0", "30760-1.0", "30870-1.0")
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
  c("eid", "bmi", "sbp", "diabetes_age_diag", "status", "diabetes_class_source",
    "hba1c", "hdl_chol","triglycerides", 
    "pgs_secr_w", "pgs_res_w")
colnames(df_all_sel_fu) <- 
  c("eid", "bmi", "sbp", "diabetes_age_diag", "status", "diabetes_class_source",
    "hba1c", "hdl_chol", "triglycerides", 
    "pgs_secr_w", "pgs_res_w")

#complete_cases
var_compl_cases <- 
  c("bmi", "sbp", "diabetes_age_diag", "hba1c", "hdl_chol", "triglycerides", 
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
  dplyr::select(`53-0.0`)

YOA1 <- year(ymd(DOA1$`53-0.0`))

#Year and Date of FU assessment
DOA2 <- df_main_sel %>% 
  filter(eid %in% compl_cases) %>% 
  dplyr::select(`53-1.0`)

YOA2 <- year(ymd(DOA2$`53-1.0`))

#First and last year
first_year_a <- min(YOA1)
last_year_a <- max(YOA2)

#Year of follow-up
n_year_fu_a <- last_year_a - first_year_a + 1


##################### 3. DATE OF BIRTH
#####################################################

#Year of birth
YOB <- df_main_sel %>% 
  filter(eid %in% compl_cases) %>% 
  dplyr::select(`34-0.0`) %>% 
  unlist()


#Month of birth
MOB <- df_main_sel %>% 
  filter(eid %in% compl_cases) %>% 
  dplyr::select(`52-0.0`) %>% 
  unlist()

##################### 4. DATE OF DIABETES DIAGNOSIS
#####################################################

#initialize
YOD <- rep(NA, nrow(df_all_sel_init))

#Age of diagnosis from main matrix
AOD <- df_all_sel_init$diabetes_age_diag

df_all_sel_init %>% 
  filter(eid %in% unlist(df_all_sel_init$eid[which(is.na(AOD))]))

YOD <- as.integer(unlist(YOB + AOD))

#for patients born in H2 add +1 to YOD
for (ind in 1:length(MOB)){
  if (MOB[ind] > 6){
    YOD[ind] <- YOD[ind] + 1
  }
}

#check no patients with YOD>YOA1 - ok
sum(YOD>YOA1, na.rm = TRUE)
for (ind in 1:length(YOD)){
  if (!is.na(YOD[ind])){
    if (YOD[ind]<YOA1[ind]){
      YOD[ind] <- YOD[ind] + 1
    }
  }
}

#check no patients with YOD>YOA2 - ok
sum(YOD>YOA2, na.rm = TRUE)
for (ind in 1:length(YOD)){
  if (!is.na(YOD[ind])){
    if (YOD[ind]>YOA2[ind]){
      YOD[ind] <- YOD[ind] - 1
    }
  }
}


##################### 5. DATE OF CVD DIAGNOSIS
#####################################################

#cvd diagnostics code
icd10_cardiovascular_list <- c(paste0("I0",0:2), # acute rheumatic fever
                               paste0("I0",5:9), # chronic rheumatic
                               paste0("I",20:25), # ischemic heart
                               paste0("I",26:28), # pulmonary
                               paste0("I",30:52), # other heart disease
                               paste0("I",60:69), # cerebrovascular
                               paste0("I",70:76), # arteries -- excluding I77.6 (in autoimmune)
                               paste0("I77", 0:5),
                               paste0("I77", 7:9),
                               paste0("I",78:79), 
                               paste0("I",81:82), # vein thrombosis
                               paste0("I",95:99)) # other circulatory

icd10_cardiovascular <- c()
for(i in 1:length(icd10_cardiovascular_list)){
  tmp <- grep(paste0("^", icd10_cardiovascular_list[i]), 
              df_hes$diag_icd10, 
              value=TRUE)
  icd10_cardiovascular <- c(icd10_cardiovascular, tmp)
}
icd10_cardiovascular <- unique(icd10_cardiovascular)

#select patients with cvd diagnosis at anytime
case_cvd <- df_hes_sel %>% 
  filter(eid %in% unlist(df_all_sel_init$eid)) %>% 
  filter(diag_icd10 %in% icd10_cardiovascular) %>% 
  dplyr::select(eid, epistart) %>% 
  arrange(eid, epistart)

#patients with cvd episode after diag.
case_cvd$year_epistart <- year(dmy(case_cvd$epistart))

case_cvd <- inner_join(data.frame(eid = df_all_sel_init$eid, 
                                  yod = YOD), 
                       case_cvd, by="eid") %>% 
  filter(year_epistart >= yod) %>% 
  dplyr::select(eid, year_epistart, yod)

#date of first episode of cvd
eid_cvd <- unique(unlist(case_cvd$eid))
YOC <- rep(NA, nrow(df_all_sel_init))

for(i in 1:nrow(df_all_sel_init)){
  if(df_all_sel_init$eid[i] %in% eid_cvd){
    pos_sel <- min(which(case_cvd$eid == df_all_sel_init$eid[i]))
    YOC[i] <- case_cvd$year_epistart[pos_sel] 
  }else{
    YOC[i] <- NA
  }
}

###############################################################################
################################# CREATE MATRIX ###############################
###############################################################################

##################### 1. GET BOUND DATES
#####################################################

#First year
##check min is diagnostic - ok
min(YOC, na.rm =TRUE)
min(YOD)
min(YOA1)
##assign first year
first_year <- min(YOD, na.rm = TRUE)

#Last year
##check max is compl.or YOA2
max(YOC, na.rm =TRUE)
max(YOD)
max(YOA2)
##assign last year
last_year <- max(YOC, na.rm = TRUE)

#Years of follow-up
n_year_fu <- last_year - first_year + 1

##################### 2. INITIALIZE MATRICES - ALL BUT AGE AT DIAGNOSIS
#####################################################

#Variables selected
var_selected <- c("bmi", "sbp", "hba1c", "hdl_chol","triglycerides", 
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




##################### 5. OBTAIN AGE AT DIAGNOSIS
#####################################################
mat_age_diag <- matrix(NA, nrow = nrow(df_all_sel_init), ncol = n_year_fu)

for (ind in 1:nrow(df_all_sel_init)){
  #obtain first year
  first_year_ind <- YOD[ind]
  
  #obtain last year
  last_year_ind <- YOC[ind]
  
  #obtain position of first year
  first_pos <- first_year_ind - first_year + 1
  
  #obtain position of last year
  if(is.na(last_year_ind)){
    last_pos <- n_year_fu
  }else{
    last_pos <- last_year_ind - first_year + 1
  }
  
  #fill matrix
  for (year in first_pos:last_pos){
    mat_age_diag[ind, year] <- df_all_sel_init[ind, "diabetes_age_diag"]
  }
}



###############################################################################
################################# INTERPOLATIONS ##############################
###############################################################################

##################### 1. INTERPOLATE ALL BUT AGES
#####################################################

mat_interpolated <- lapply(matrices, interpolate_linear)

##################### 2. EXTRAPOLATE ALL BUT AGES
#####################################################

mat_extrapolated <- lapply(mat_interpolated, function(mat){
  extrapolate(mat, YOD, YOC, first_year)})

###############################################################################
############################### CLUSTER DATA ##################################
###############################################################################

##################### 1. SCALE
#####################################################

#age
mat_age_diag_scaled <- matrix(scale(as.vector(mat_age_diag)),
                            ncol = ncol(mat_age_diag))

#all others
mat_extrapolated_scaled <- lapply(mat_extrapolated, scale_matrix)


##################### 2. PREDICT
#####################################################
#prepare matrix
mat_clustered <- mat_age_diag_scaled


#loop over individuals
for(ind in 1:nrow(mat_clustered)){
  
  #find pos with first non NA value
  pos_start <- 1
  while(is.na(mat_clustered[ind, pos_start])){
    pos_start <- pos_start + 1
  }
  
  #find pos with last non NA value
  if(!is.na(mat_clustered[ind, ncol(mat_clustered)])){
    pos_end <- ncol(mat_clustered)
  } else {
    pos_end <- pos_start + 1
    while(!is.na(mat_clustered[ind, pos_end])){
      pos_end <- pos_end + 1
    }
    pos_end <- pos_end - 1
  }
  
  #loop over positions
  for(pos in pos_start:pos_end){
    mat_clustered[ind, pos] <-
      predict(kcca_model,
            data.frame(bmi = mat_extrapolated_scaled$bmi[ind, pos],
                       sbp = mat_extrapolated_scaled$sbp[ind, pos],
                       diabetes_age_diag = mat_age_diag_scaled[ind, pos],
                       hba1c = mat_extrapolated_scaled$hba1c[ind, pos],
                       hdl_chol = mat_extrapolated_scaled$hdl_chol[ind, pos],
                       triglycerides = mat_extrapolated_scaled$triglycerides[ind, pos]))
    
  }
}

##################### 3. REORDER CLUSTERING
#####################################################

#clean matrix
mat_clustered_ordered <- matrix(NA, nrow = nrow(mat_clustered),
                                ncol = ncol(mat_clustered))

#reorder
for(i in 1:length(cluster_ordered)){
  mat_clustered_ordered[
    which(mat_clustered == cluster_ordered[i])] <- cluster_names_ordered[i]
}
apply(mat_clustered_ordered, 1, function(x)length(unique(x)))


##################### 2. CLEAN CLUSTERING AND CODE FACTOR
#####################################################

#create matrices of variables
mat_variables <- lapply(1:(length(cluster_ordered)-1), matrix, data= 0, 
                        nrow = nrow(mat_clustered_ordered), 
                        ncol = ncol(mat_clustered_ordered))

names(mat_variables) <- cluster_names_ordered[-1]



#fill matrices
for(mat in 1:length(mat_variables)){
  for(col in 1:ncol(mat_clustered_ordered)){
    #replace NA with -1
    mat_variables[[mat]][ ,col][
      which(is.na(mat_clustered_ordered[ ,col]))] <- -1
    
    #replace value of cluster 
    mat_variables[[mat]][ ,col][
      which(mat_clustered_ordered[ ,col] == 
              cluster_names_ordered[mat+1])] <- 1
  }
}



##################### 3. CLEAN DATES
#####################################################
#deal with missing values
YOC[is.na(YOC)] <- 1900

#create matrix
dates <- data.frame(YOA1, YOA2, YOD, YOC)


##################### 4. GET INITIAL AND LAST CLUSTER CODED
#####################################################
cluster_init <- rep(NA, nrow(mat_clustered_ordered))
cluster_fu <- rep(NA, nrow(mat_clustered_ordered))

for(ind in 1:nrow(mat_clustered_ordered)){
  
  #find pos with first non NA value
  pos_start <- 1
  while(is.na(mat_clustered_ordered[ind, pos_start])){
    pos_start <- pos_start + 1
  }
  
  #find pos with last non NA value
  if(!is.na(mat_clustered_ordered[ind, ncol(mat_clustered_ordered)])){
    pos_end <- ncol(mat_clustered_ordered)
  } else {
    pos_end <- pos_start + 1
    while(!is.na(mat_clustered_ordered[ind, pos_end])){
      pos_end <- pos_end + 1
    }
    pos_end <- pos_end - 1
  }
  
  
  #get values
  cluster_init[ind] <- mat_clustered_ordered[ind, pos_start]
  cluster_fu[ind] <- mat_clustered_ordered[ind, pos_end]
  
  #code
  for(i in 1:length(cluster_names_ordered)){
    cluster_init[which(cluster_init == cluster_names_ordered[i])] <- i
    cluster_fu[which(cluster_fu == cluster_names_ordered[i])] <- i
  }
}

clusters <- data.frame(cluster_init, cluster_fu)

##################### 5. STORE
#####################################################
out_dir <- "~/XXX/"
simul_name <- "_simul_V6"

save_for_mcmc(dates, out_dir, "dates", simul_name)
save_for_mcmc(mat_variables$MARD1, out_dir, "mard1", simul_name)
save_for_mcmc(mat_variables$MARD2, out_dir, "mard2", simul_name)
save_for_mcmc(mat_variables$MOD, out_dir, "mod", simul_name)
save_for_mcmc(mat_variables$SIRD, out_dir, "sird", simul_name)
save_for_mcmc(clusters, out_dir, "clusters", simul_name)


