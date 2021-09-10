library(data.table)
library(dplyr)
library(tidyr)
library(NbClust)
library(ggplot2)
library(gridExtra)
library(ggdendro)
library(cluster)
library(purrr)
library(tibble)
library(ggradar)
library(naniar)
library(lubridate)

###############################################################################
################################ TABLES LOADING ###############################
###############################################################################

#load tables
setwd("~/ukbb_diabetes/data/obs_selection")

df_main <- fread("~/ukbb_diabetes/data/obs_selected/df_main_msm_presel.csv")
df_bio <- fread("~/ukbb_diabetes/data/obs_selected/df_bio_msm_presel.csv")
df_genes <-fread("~/ukbb_diabetes/data/obs_selected/df_genes_msm_presel.csv")
df_hes <-fread("~/ukbb_diabetes/data/obs_selected/df_hes_msm_presel.csv")
df_genetic_dict <- fread("~/ukbb_diabetes/data/obs_selected/genes_dict.csv")
df_diag_dict <- fread("~/ukbb_diabetes/data/obs_selected/diag_dict.csv")

df_main <- as.data.frame(df_main)
df_bio <- as.data.frame(df_bio)
df_genes <- as.data.frame(df_genes)
df_hes <- as.data.frame(df_hes)
df_genetic_dict <- as.data.frame(df_genetic_dict)
df_diag_dict <- as.data.frame(df_diag_dict)

###############################################################################
################################# GET FUNCTIONS ###############################
###############################################################################
source("~/ukbb_diabetes/functions.R")


###############################################################################
############################### PATIENTS SELECTION ############################
###############################################################################

##################### 1. REMOVE WOMEN PREGNANT OR UNSURE ABOUT PREGANCY
#####################################################
eid_preg_init <- df_main %>% 
  filter(`31-0.0` == 0 ) %>% #select women
  filter(`3140-0.0` == 1 | `3140-0.0` == 2 ) %>% #remove women preg or unsure
  dplyr::select(eid)

eid_preg_fu <- df_main %>% 
  filter(`31-0.0` == 0 ) %>% #select women
  filter(`3140-1.0` == 1 | `3140-1.0` == 2 ) %>% #remove women preg or unsure
  dplyr::select(eid)

#combine eid to remove
eid_to_remove <- unique(unlist(c(eid_preg_init,eid_preg_fu)))

#remove from dataframes
df_main <- df_main %>% 
  filter (!(eid %in% eid_to_remove))

df_bio <- df_bio %>% 
  filter (!(eid %in% eid_to_remove))

df_genes <- df_genes %>% 
  filter (!(eid %in% eid_to_remove))

df_hes <- df_hes %>% 
  filter (!(eid %in% eid_to_remove))


###############################################################################
############################# DATAFRAMES PREPARATION ##########################
###############################################################################

#get info on patients at init and fu
info_init <- split_by_status(df_main, df_bio, df_genes, df_hes, 0)
info_fu <- split_by_status(df_main, df_bio, df_genes, df_hes, 1)

#add info columns in df_main
df_main <- df_main %>% 
  inner_join(info_init, by = "eid")

df_main <- df_main %>% 
  inner_join(info_fu, by = "eid")

###############################################################################
########################### ANALYSIS OF TRANSITIONS ###########################
###############################################################################

##################### 1. TRANSITIONS INITIAL - FU
#####################################################
df_main %>% 
  filter(status_init == "diabetes") %>% 
  group_by(status_fu) %>% 
  dplyr::summarise(tot = n())

df_main %>% 
  filter(status_init == "prediabetes") %>% 
  group_by(status_fu) %>% 
  dplyr::summarise(tot = n())

df_main %>% 
  filter(status_init == "healthy") %>% 
  group_by(status_fu) %>% 
  dplyr::summarise(tot = n())

##################### 2. TRANSITIONS FU - HES
#####################################################

###### a. Convert to date

df_hesl$epistart <- dmy(df_hes$epistart)
df_hes$epiend <- dmy(df_hes$epiend)

df_main$`53-0.0` <- ymd(df_main$`53-0.0`)
df_main$`53-1.0` <- ymd(df_main$`53-1.0`)

###### b. Pre-diabetes at fu
eid_diabetes_hes_pred <- df_hes %>% 
  left_join(df_main[ ,c("eid", "53-1.0", "status_fu")], by = "eid") %>% 
  filter(status_fu == "prediabetes") %>% 
  filter(diag_icd10 %in% c(codes_compl_t1d, codes_compl_t2d)) %>% 
  filter(epistart > `53-1.0`) %>% 
  select(eid) %>% 
  unique() %>% 
  unlist()
length(eid_diabetes_hes_pred)

###### c. Healthy at fu
eid_diabetes_hes_healthy <- df_hes_sel %>% 
  left_join(df_main_sel[ ,c("eid", "53-1.0", "status_fu")], by = "eid") %>% 
  filter(status_fu == "healthy") %>% 
  filter(diag_icd10 %in% c(codes_compl_t1d, codes_compl_t2d)) %>% 
  filter(epistart > `53-1.0`) %>% 
  select(eid) %>% 
  unique() %>% 
  unlist()
length(eid_diabetes_hes_healthy)


###############################################################################
############################ DATAFRAMES SAVING ################################
###############################################################################

#save data frames
write.csv(df_main,"~/ukbb_diabetes/data/obs_selected/df_main_msm_sel.csv", row.names = FALSE)
write.csv(df_bio,"~/ukbb_diabetes/data/obs_selected/df_bio_msm_sel.csv", row.names = FALSE)
write.csv(df_genes,"~/ukbb_diabetes/data/obs_selected/df_genes_msm_sel.csv", row.names = FALSE)
write.csv(df_hes,"~/ukbb_diabetes/data/obs_selected/df_hes_msm_sel.csv", row.names = FALSE)
