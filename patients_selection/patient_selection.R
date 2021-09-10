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


###############################################################################
################################ TABLES LOADING ###############################
###############################################################################

#load tables
setwd("~/ukbb_diabetes/data/obs_selection")

df <- fread("~/ukbb_diabetes/data/raw/ukb26390.csv")
df_diag <- fread("~/ukbb_diabetes/data/raw/hesin_diag.txt")
df_bio <- fread("~/ukbb_diabetes/data/raw/ukb27725.csv")
df_genetic <- readRDS("~/ukbb_diabetes/data/raw/genetic_data_extracted.rds")
df_genetic_dict <- fread("~/ukbb_diabetes/data/obs_selected/genes_dict.csv")
df_hes_main <- fread("~/ukbb_diabetes/data/raw/hesin.txt")
df_hes_diag <- fread("~/ukbb_diabetes/data/raw/hesin_diag.txt")
df_death_main <- fread("~/ukbb_diabetes/data/raw/death.txt")
df_death_cause <- fread("~/ukbb_diabetes/data/raw/death_cause.txt")


df <- as.data.frame(df)
df_diag <- as.data.frame(df_diag)
df_bio <- as.data.frame(df_bio)
df_genetic <- as.data.frame(df_genetic)
df_hes_main <- as.data.frame(df_hes_main)
df_hes_diag <- as.data.frame(df_hes_diag)
df_death_main <- as.data.frame(df_death_main)
df_death_cause <- as.data.frame(df_death_cause)
df_genetic_dict <- as.data.frame(df_genetic_dict)

df_genetic <- rownames_to_column(df_genetic, var = "eid")
df_genetic[ ,-1] <- sapply(df_genetic[ ,-1], as.integer)
df_genetic <- df_genetic %>% 
  na_if(0)


###############################################################################
################################# GET FUNCTIONS ###############################
###############################################################################
source("~/ukbb_diabetes/functions.R")

###############################################################################
############################# GENES DATA PROCESSING ###########################
###############################################################################

#remove snps not wanted
snps_to_remove <- c("rs5215")

#create dataframe with number of risk allele
df_genes <- df_genetic %>% 
  dplyr::select(!snps_to_remove)

#populate dataframe  with number of risk alleles
for (j in 2:ncol(df_genes)){
  #save info of dict
  dict <- df_genetic_dict %>% 
    filter(snp == names(df_genes)[j])
  
  risk_allele <- dict[5]
  minor_allele <- dict[8]
  
  #select column
  col <- df_genes[, j]
  
  #update values with number of risk alleles
  if (minor_allele == risk_allele){
    col[col == 1] = 2
    col[col == 2] = 1
    col[col == 3] = 0
  } else {
    col[col == 1] = 0
    col[col == 2] = 1
    col[col == 3] = 2
  }
  
  df_genes[, j] <- col
}

#calculate scores

##select snps
df_genetic_dict_secr <- df_genetic_dict %>% 
  filter(type == "Secretion") %>% 
  filter(snp != snps_to_remove) %>% 
  dplyr::select(c("snp", "insulin_effect"))

df_genetic_dict_res <- df_genetic_dict %>% 
  filter(type == "Resistance") %>% 
  filter(snp != snps_to_remove) %>% 
  dplyr::select(c("snp", "insulin_effect"))

##weighted scores
score_secr_weighted <- as.matrix(df_genes[ ,df_genetic_dict_secr$snp]) %*% 
  as.matrix(df_genetic_dict_secr$insulin_effect)
score_secr_weighted <- 
  score_secr_weighted/(sum(df_genetic_dict_secr$insulin_effect)) * nrow(df_genetic_dict_secr)

score_res_weighted <- as.matrix(df_genes[ ,df_genetic_dict_res$snp]) %*% 
  as.matrix(df_genetic_dict_res$insulin_effect)
score_res_weighted <- 
  score_res_weighted/(sum(df_genetic_dict_res$insulin_effect)) * nrow(df_genetic_dict_res)

##non-weighted scores
score_secr_nweighted <- rowSums(df_genes[ ,df_genetic_dict_secr$snp])
score_res_nweighted <- rowSums(df_genes[ ,df_genetic_dict_res$snp])

##save in dataframe
df_genes$score_secr_weighted  <- score_secr_weighted
df_genes$score_res_weighted  <- score_res_weighted
df_genes$score_secr_nweighted  <- score_secr_nweighted
df_genes$score_res_nweighted  <- score_res_nweighted


###############################################################################
############################### HES DATA PROCESSING ###########################
###############################################################################

#join the two tables
df_hes <- left_join(df_hes_diag[ ,c("eid", "ins_index", "arr_index", 
                                    "level", "diag_icd10")], 
                    df_hes_main[ ,c("eid", "ins_index", "epistart",
                                     "epiend", "mainspef_uni", "tretspef_uni")], 
                    by = c("eid", "ins_index"))


###############################################################################
############################## MULTI STATES MODEL #############################
###############################################################################

##################### 1.SELECT PATIENTS FOR STUDY
#####################################################

#filter out rows with no bio data
df_msm <- df %>% 
  filter(eid %in% df_bio$eid)

nrow(df_msm)

#filter out rows with no 2 assessment centers
df_msm <- df_msm%>% 
  filter((!(is.na(`53-0.0`))) & (!(is.na(`53-1.0`))))

nrow(df_msm)

#filter out rows with missing glycated haemoglobin
eid_msm <- df_msm$eid

df_bio_msm <- df_bio %>% 
  filter(eid %in% eid_msm)

nrow(df_bio_msm) #to check

df_bio_msm <- df_bio_msm %>% 
  filter((!(is.na(`30750-0.0`))) & (!(is.na(`30750-1.0`))))

nrow(df_bio_msm)

df_msm <- df_msm %>% 
  filter(eid %in% df_bio_msm$eid)

nrow(df_msm) #to check


##################### 2. CREATE GENES AND HES DATAFRAMES
#####################################################
#create genes dataset
df_genes_msm <- df_genes %>% 
  filter(eid %in% unlist(df_bio_msm$eid))

nrow(df_genes_msm) #to check

#filter out rows with no genes data
df_bio_msm <- df_bio_msm %>% 
  filter(eid %in% unlist(df_genes_msm$eid))
nrow(df_bio_msm)

df_msm <- df_msm %>% 
  filter(eid %in% unlist(df_genes_msm$eid))
nrow(df_msm)

#create hes dataset
df_hes_msm <- df_hes %>% 
  filter(eid %in% unlist(df_bio_msm$eid))

##################### 3. DATA FRAMES SAVING
#####################################################

#save dataframes
write.csv(df_msm,"~/ukbb_diabetes/data/obs_selected/df_main_msm_presel.csv", row.names = FALSE)
write.csv(df_bio_msm,"~/ukbb_diabetes/data/obs_selected/df_bio_msm_presel.csv", row.names = FALSE)
write.csv(df_genes_msm,"~/ukbb_diabetes/data/obs_selected/df_genes_msm_presel.csv", row.names = FALSE)
write.csv(df_hes_msm,"~/ukbb_diabetes/data/obs_selected/df_hes_msm_presel.csv", row.names = FALSE)

###############################################################################
############################### DIABETES CLUSTERING ###########################
###############################################################################

##################### 1. BASIC PATIENTS SELECTION
#####################################################

#filter out rows with no biochemestry data
df_main_dc <- df %>% 
  filter(eid %in% df_bio$eid)

nrow(df_main_dc)

df_bio_dc <- df_bio %>% 
  filter(eid %in% df_main_dc$eid)

nrow(df_bio_dc)

#filter out rows with missing glycated haemoglobin at initial
df_bio_dc <- df_bio_dc %>% 
  filter((!(is.na(`30750-0.0`))))

nrow(df_bio_dc)

df_main_dc <- df_main_dc %>% 
  filter(eid %in% df_bio_dc$eid)

nrow(df_main_dc) #to check

##################### 2. CREATE GENES AND HES DATAFRAMES
#####################################################

#create genes dataset
df_genes_dc <- df_genes %>% 
  filter(eid %in% unlist(df_bio_dc$eid))

nrow(df_genes_dc) #to check

#filter out rows with no genes data
df_bio_dc <- df_bio_dc %>% 
  filter(eid %in% unlist(df_genes_dc$eid))
nrow(df_bio_dc)

df_main_dc <- df_main_dc %>% 
  filter(eid %in% unlist(df_genes_dc$eid))
nrow(df_main_dc)

#create hes dataset
df_hes_dc <- df_hes %>% 
  filter(eid %in% unlist(df_bio_dc$eid))

##################### 3. GET STATUS
#####################################################

#get info on patients at init and fu
info_init_dc <- split_by_status(df_main_dc, df_bio_dc, df_genes_dc, df_hes_dc, 0)
info_fu_dc <- split_by_status(df_main_dc, df_bio_dc, df_genes_dc, df_hes_dc, 1)

#add info columns in df_main
df_main_dc <- df_main_dc %>% 
  inner_join(info_init_dc, by = "eid")

df_main_dc <- df_main_dc %>% 
  inner_join(info_fu_dc, by = "eid")


##################### 4. SELECT PATIENTS WITH DIABETES
#####################################################
df_main_dc_sel <- df_main_dc %>% 
  filter(status_init == "diabetes" | status_fu == "diabetes")

df_bio_dc_sel <- df_bio_dc %>% 
  filter(eid %in% unlist(df_main_dc_sel$eid))

df_genes_dc_sel <- df_genes_dc %>% 
  filter(eid %in% unlist(df_main_dc_sel$eid))

df_hes_dc_sel <- df_hes_dc %>% 
  filter(eid %in% unlist(df_main_dc_sel$eid))

##################### 5. DATA FRAMES SAVING
#####################################################

write.csv(df_main_dc_sel,"~/ukbb_diabetes/data/obs_selected/df_main_dc_sel.csv", row.names = FALSE)
write.csv(df_bio_dc_sel,"~/ukbb_diabetes/data/obs_selected/df_bio_dc_sel.csv", row.names = FALSE)
write.csv(df_genes_dc_sel,"~/ukbb_diabetes/data/obs_selected/df_genes_dc_sel.csv", row.names = FALSE)
write.csv(df_hes_dc_sel,"~/ukbb_diabetes/data/obs_selected/df_hes_dc_sel.csv", row.names = FALSE)

