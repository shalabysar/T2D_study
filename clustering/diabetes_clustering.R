library(data.table)
library(dplyr)
library(tidyr)
library(NbClust)
library(ggplot2)
library(gridExtra)
library(cluster)
library(purrr)
library(tibble)
library(kableExtra)
library(fpc)
library(fossil)
library(matrixStats)
library(clusterCrit)
library(logr)
library(clusterSim)
library(ggdendro)
library(dendextend)
library(shipunov)
library(clv)
library(ggbiplot)
library(ggalluvial)
library(GPArotation)
library(psych)
library (pastecs)
library(lubridate)
library(SwarmSVM)
library(flexclust)

###############################################################################
################################### LOAD TABLES ###############################
###############################################################################

#load tables - XXX to be filled as required
setwd("~/XXX")

df_main_sel <- fread("~/XXX/df_main_dc_sel.csv")
df_bio_sel <- fread("~/XXX/df_bio_dc_sel.csv")
df_genes_sel <- fread("~/XXX/df_genes_dc_sel.csv")
df_hes_sel <- fread("~/XXX/df_hes_dc_sel.csv")
df_genetic_dict <- fread("~/XXX.csv")
df_diag_dict <- fread("~/XXX.csv")


df_main_sel <- as.data.frame(df_main_sel)
df_bio_sel <- as.data.frame(df_bio_sel)
df_genes_sel <- as.data.frame(df_genes_sel)
df_hes_sel <- as.data.frame(df_hes_sel)
df_genetic_dict <- as.data.frame(df_genetic_dict)
df_diag_dict <- as.data.frame(df_diag_dict)

###############################################################################
################################# GET FUNCTIONS ###############################
###############################################################################
source("~/functions.R")

###############################################################################
############################### BUILD DATAFRAMES ##############################
###############################################################################

##################### 1. SELECT PATIENTS IN MAIN
#####################################################

###### a. Diabetes at initial 
df_main <- df_main_sel %>% 
  filter((status_init == "diabetes" & diabetes_class_source_init != "hba1c"))

###### b. Remove pregnant women or unsure about pregnancy 

#get eid to remove
eid_preg <- df_main %>% 
  filter(`31-0.0` == 0 ) %>% #select women
  filter(`3140-0.0` == 1 | `3140-0.0` == 2 ) %>% #remove women preg or unsure
  dplyr::select(eid)

###### c. Remove non-adults
#transform all to integer
df_main$diabetes_age_diag_init <- as.integer(df_main$diabetes_age_diag_init)

#select eid of children - to remove
eid_children <- df_main %>% 
  filter(!(eid %in% unlist(eid_preg))) %>% 
  filter(diabetes_age_diag_init < 18) %>% 
  dplyr::select(eid) 

###### d. Remove obese and diagnosed below 30 
eid_young_notobese <- df_main %>% 
  filter(!(eid %in% unlist(eid_preg))) %>% 
  filter(!(eid %in% unlist(eid_children))) %>% 
  filter(diabetes_age_diag_init <=30) %>% 
  filter(`21001-0.0`<30) %>% 
  dplyr::select(eid)
  
  
###### e. Remove from dataframe
eid_to_remove <- unique(unlist(c(eid_preg,
                                 eid_children,
                                 eid_young_notobese)))
  
df_main <- df_main %>% 
  filter (!(eid %in% unlist(eid_to_remove)))

df_bio <- df_bio_sel %>% 
  filter (eid %in% unlist(df_main$eid))

df_genes <- df_genes_sel %>% 
  filter (eid %in% unlist(df_main$eid))

df_hes <- df_hes_sel %>% 
  filter (eid %in% unlist(df_main$eid))

##################### 2. DEFINE VARIABLES AND TABLE
#####################################################
var_main <- c("21001-0.0", "4079-0.0", "4080-0.0", "21003-0.0","diabetes_age_diag_init", "74-0.0")
var_bio <- c("30750-0.0", "30760-0.0", "30780-0.0", "30870-0.0", "30740-0.0")
var_genes <- c("score_secr_weighted", "score_res_weighted", "score_secr_nweighted", "score_res_nweighted")

df_all <- prepare_tables_3df(
  df_main, c("eid", var_main), 
  df_bio, c("eid", var_bio),
  df_genes,c("eid", var_genes),
  "eid")

var_all <- c("eid", 
             "bmi", "dbp", "sbp", "age_at_assessment", "diabetes_age_diag", "fasting_time", 
             "hba1c", "hdl_chol", "ldl_direct","triglycerides", "glucose", 
             "pgs_secr_w", "pgs_res_w", "pgs_secr_nw", "pgs_res_nw")

colnames(df_all) <- var_all

var_names_display <- c("Eid", 
                       "BMI (Kg/m2)", "Diastolic BP (mmHg)", "Systolic BP (mmHg)",
                       "Age at assessment (years)", "Age at diagnosis (years)", "Fasting time (hours)",
                       "HbA1c (mmol/mol)", "HDL Cholesterol (mmol/L)", 
                       "LDL Cholesterol (mmol/L)", "Triglycerides (mmol/L)", "Glucose (mmol/L)",
                       "PGS Secretion (weighted)", "PGS Resistance (weighted)",
                       "PGS Secretion (non-weighted)", "PGS Resistance (non-weighted)")


##################### 3. REMOVE NON-COMPLETE CASES
#####################################################

###### a. Select relevant columns
var_all_compl <- c("bmi", "dbp", "sbp", "age_at_assessment", "diabetes_age_diag",
                   "hba1c", "hdl_chol", "ldl_direct", "triglycerides", 
                   "pgs_secr_w", "pgs_res_w")

var_names_display_compl <- c("BMI (Kg/m2)", "DBP (mmHg)", "SBP (mmHg)",
                             "Age (years)", "Age at diagnosis (years)",
                             "HbA1c (mmol/mol)", "HDL Cholesterol (mmol/L)", 
                             "LDL Cholesterol (mmol/L)", "Triglycerides (mmol/L)",
                             "PGS IS", "PGS IR")

###### b. Analyze missingness
sum(complete.cases(df_all[ , var_all_compl]))
sapply(df_all[ , var_all_compl] , function(x) sum(is.na(x)))

###### c. Remove non-complete cases
df_all_compl <- df_all[complete.cases(df_all[ , var_all_compl]), ]

###############################################################################
############################# DISTRIBUTION ANALYSIS ###########################
###############################################################################

##################### 1. GENERAL OVERVIEW
#####################################################
#get summary statistics
sum_compl <- generate_sum_stat(df_all_compl,
                               var_all_compl, var_names_display_compl)

#print summary statistics
kbl(t(round(sum_compl, 3))) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), font_size = 12) %>%
  scroll_box(height = "300px")

#check if all age of diagnosis are =< age at assessment - OK
sum((df_all_compl$diabetes_age_diag - df_all_compl$age_at_assessment) > 0)

##################### 3. PAIRS PLOTS
#####################################################
plt_pairs_comp <- ggpairs(df_all_compl[ ,var_all_compl])

png("plots/plt_pairs_comp.png", width = 1700, height = 850)
plt_pairs_comp
dev.off()

###############################################################################
################################## OUTLIERS ###################################
###############################################################################

##################### 1. GET BOUNDS
#####################################################

###### a. Get bounds - automated
#get upper bound - 99% quantile
upper_bounds <- sapply(c(2:16), 
                       function(i){quantile(df_all_compl[ ,i], 0.99, na.rm = TRUE)})
names(upper_bounds) <- colnames(df_all_compl)[2:16]

#get lower bound - 1% quantile
lower_bounds <- sapply(c(2:16), 
                       function(i){quantile(df_all_compl[ ,i], 0.01, na.rm = TRUE)})
names(lower_bounds) <- colnames(df_all_compl[2:16])

#merge
bounds_outliers <- cbind(lower_bounds, upper_bounds)

kbl(round(bounds_outliers,2)) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), font_size = 12) %>%
  scroll_box(height = "300px")


##################### 2. REMOVE OUTLIERS
#####################################################
df_all_compl_no <- remove_outliers(df_all_compl, 
                                   var_all_compl,
                                   bounds_outliers)


##################### 3. GET SUMMARY STAT
#####################################################

#get summary statistics
sum_compl_no <- generate_sum_stat(df_all_compl_no,
                                  var_all_compl, var_names_display_compl)

#print summary statistics
kbl(t(round(sum_compl_no, 3))) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), font_size = 12) %>%
  scroll_box(height = "300px")


##################### 4. PAIRS PLOTS
#####################################################
# pairs plot
plt_pairs_comp_no <- ggpairs(df_all_compl_no[ ,var_all_compl], 
                             columnLabels = c("BMI", "DBP", "SBP",
                                              "Age", "Age at diagnosis",
                                              "HbA1c", "HDL Chol", 
                                              "LDL Chol", "Triglycerides",
                                              "PGS IS", "PGS IR"))

png("plots/plt_pairs_comp_no.png", width = 1700, height = 850)
plt_pairs_comp_no
dev.off()

###############################################################################
############################## FEATURES SELECTION #############################
###############################################################################

##################### 0. SAVE CURRENT IMAGE TO USE IN STEPWISE ALGO
#####################################################
save.image("~/XXX/data_for_stepwise_210725.RData")

##################### 1. FORWARD STEPWISE CLUSTERING
#####################################################

#p1 Run for k_range=c(2,6) on df_all_compl_no for var_all_compl variables
#p2 Run for k_range=c(7,10) on df_all_compl_no for var_all_compl variables

###### a. Load results from cluster run
#load("~/XXX/diabetes_clusters_forw_step_p1_210725.RData")
#load("~/XXX/diabetes_clusters_forw_step_p2_210725.RData")


#format results
cluster_step_kmeans_sil <- cluster_step_kmeans_sil_26
cluster_step_kmeans_ch <- cluster_step_kmeans_ch_26
cluster_step_clara_sil <- cluster_step_clara_sil_26
cluster_step_clara_ch <- cluster_step_clara_ch_26
cluster_step_kmeans_sil

###### b. Analyze results
sum_step_kmeans_sil <- create_summary_stepwise(cluster_step_kmeans_sil, var_all_compl, 3)
sum_step_kmeans_ch <- create_summary_stepwise(cluster_step_kmeans_ch, var_all_compl, 3)
sum_step_clara_sil <- create_summary_stepwise(cluster_step_clara_sil, var_all_compl, 3)
sum_step_clara_ch <- create_summary_stepwise(cluster_step_clara_ch, var_all_compl, 3)

###### c. Print results
round(cluster_step_kmeans_sil[[1]][[1]],3) %>% 
  as.data.frame() %>% 
  kbl() %>% 
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), font_size = 12) %>%
  scroll_box(height = "300px")

 sum_step_kmeans_sil %>% 
  as.data.frame() %>% 
  kbl() %>% 
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), font_size = 12) %>%
  scroll_box(height = "300px")

sum_step_kmeans_ch %>% 
  as.data.frame() %>% 
  kbl() %>% 
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), font_size = 12) %>%
  scroll_box(height = "300px")


sum_step_clara_sil %>% 
  as.data.frame() %>% 
  kbl() %>% 
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), font_size = 12) %>%
  scroll_box(height = "300px")

sum_step_clara_ch %>% 
  as.data.frame() %>% 
  kbl() %>% 
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), font_size = 12) %>%
  scroll_box(height = "300px")


##################### 2. BACKWARD STEPWISE CLUSTERING
#####################################################

# Run for k_range=c(2,10) on df_all_compl_no for var_all_compl variables

###### a. Load results from cluster run
#load("~/XXX/diabetes_clusters_back_step_210725.RData")

###### b. Format results

#create matrix
sum_step_back_kmeans <- matrix(NA, 
                               ncol = ncol(cluster_one_step_backward_kmeans_210),
                               nrow = nrow(cluster_one_step_backward_kmeans_210))
colnames(sum_step_back_kmeans) <- colnames(cluster_one_step_backward_kmeans_210)

#fill scores
for (j in 1:ncol(sum_step_back_kmeans)){
  sum_step_back_kmeans[-1, j] <- paste0(var_all_compl[order(cluster_one_step_backward_kmeans_210[-1, j])],
                                        " (",
                                        round(cluster_one_step_backward_kmeans_210[-1, j]
                                              [order(cluster_one_step_backward_kmeans_210[-1, j])], 3),
                                        ")")
}

#fill score for all
sum_step_back_kmeans[1, ] <- round(cluster_one_step_backward_kmeans_210[1, ], 3)

#Separate silhouette and ch
sum_step_back_kmeans_sil <- sum_step_back_kmeans[ ,c(TRUE,FALSE)]
sum_step_back_kmeans_sil[1, ] <- paste0(c(2:10), " centers (", sum_step_back_kmeans_sil[1, ], ")")

sum_step_back_kmeans_ch <- sum_step_back_kmeans[ ,c(FALSE,TRUE)]
sum_step_back_kmeans_ch[1, ] <- paste0(c(2:10), " centers (", sum_step_back_kmeans_ch[1, ], ")")

###### c. Print results
sum_step_back_kmeans_sil %>% 
  as.data.frame() %>% 
  kbl() %>% 
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), font_size = 12) %>%
  scroll_box(height = "300px")

sum_step_back_kmeans_ch %>% 
  as.data.frame() %>% 
  kbl() %>% 
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), font_size = 12) %>%
  scroll_box(height = "300px")


##################### 3. PCA Analysis
#####################################################
data_for_pca <- process_for_cluster(df_all_compl_no, var_sel_1)
colnames(data_for_pca) <- c("BMI", "DBP", "SBP",
                                     "Age", "Age at diagnosis",
                                     "HbA1c", "HDL Chol", 
                                     "LDL Chol", "Triglycerides",
                                     "PGS-IS", "PGS-IR")
res.pca <- prcomp(data_for_pca,  scale = FALSE)
biplot_pca <- ggbiplot(res.pca, alpha = 0) + xlim(-1, 2) + ylim(-1, 1.75) +
  xlab("Standardised PC1 (17.6% explained variance)") +
  ylab("Standardised PC2 (14.6% explained variance)")

pdf("plots/biplot_pca.pdf")
biplot_pca
dev.off()


##################### 4. FINAL VARIABLES SELECTED
#####################################################

###### a. All variables
var_sel_1 <- var_all_compl

###### b. Remove age at assessment, DBP and LDL
var_sel_2 <- var_all_compl[-c(2, 4, 8)]

###### c. Remove age at assessment, DBP, LDL and PGS
var_sel_3 <- var_all_compl[-c(2, 4, 8, 10, 11)]

###### d. Original Ahlqvist variables
var_sel_4 <- var_all_compl[-c(2:4, 7:9)]


###############################################################################
################################ BEST K SELECTION #############################
###############################################################################

##################### 1. CLUSTERING
#####################################################

###### a. K-means
cluster_kmeans_1 <- cluster_kmeans(process_for_cluster(df_all_compl_no, var_sel_1),
                                   k_range = c(2:10),
                                   n_runs = 50,
                                   iter_max = 500,
                                   algo = "MacQueen")

cluster_kmeans_2 <- cluster_kmeans(process_for_cluster(df_all_compl_no, var_sel_2),
                                   k_range = c(2:10),
                                   n_runs = 50,
                                   iter_max = 500,
                                   algo = "MacQueen")

cluster_kmeans_3 <- cluster_kmeans(process_for_cluster(df_all_compl_no, var_sel_3),
                                   k_range = c(2:10),
                                   n_runs = 50,
                                   iter_max = 500,
                                   algo= "MacQueen")

cluster_kmeans_4 <- cluster_kmeans(process_for_cluster(df_all_compl_no, var_sel_4),
                                   k_range = c(2:10),
                                   n_runs = 50,
                                   iter_max = 500,
                                   algo = "MacQueen")
###### b. K-medoids euclidean

cluster_clara_1 <- cluster_clara(process_for_cluster(df_all_compl_no, var_sel_1),
                                 k_range = c(2:10),
                                 n_samples = 500, 
                                 n_sampsize = 1000)
cluster_clara_2 <- cluster_clara(process_for_cluster(df_all_compl_no, var_sel_2),
                                 k_range = c(2:10),
                                 n_samples = 500, 
                                 n_sampsize = 1000)
cluster_clara_3 <- cluster_clara(process_for_cluster(df_all_compl_no, var_sel_3),
                                 k_range = c(2:10),
                                 n_samples = 500, 
                                 n_sampsize = 1000)
cluster_clara_4 <- cluster_clara(process_for_cluster(df_all_compl_no, var_sel_4),
                                 k_range = c(2:10),
                                 n_samples = 500, 
                                 n_sampsize = 1000)

###### c. K-medoids Manhattan

cluster_clara_m_1 <- cluster_clara(process_for_cluster(df_all_compl_no, var_sel_1),
                                   k_range = c(2:10),
                                   n_samples = 500, 
                                   n_sampsize = 1000,
                                   distance_metric = "manhattan")
cluster_clara_m_2 <- cluster_clara(process_for_cluster(df_all_compl_no, var_sel_2),
                                   k_range = c(2:10),
                                   n_samples = 500, 
                                   n_sampsize = 1000,
                                   distance_metric = "manhattan")
cluster_clara_m_3 <- cluster_clara(process_for_cluster(df_all_compl_no, var_sel_3),
                                   k_range = c(2:10),
                                   n_samples = 500, 
                                   n_sampsize = 1000,
                                   distance_metric = "manhattan")
cluster_clara_m_4 <- cluster_clara(process_for_cluster(df_all_compl_no, var_sel_4),
                                   k_range = c(2:10),
                                   n_samples = 500, 
                                   n_sampsize = 1000,
                                   distance_metric = "manhattan")

###### d. K-medoids Jaccard

cluster_clara_j_1 <- cluster_clara(process_for_cluster(df_all_compl_no, var_sel_1),
                                   k_range = c(2:10),
                                   n_samples = 500, 
                                   n_sampsize = 1000,
                                   distance_metric = "jaccard")
cluster_clara_j_2 <- cluster_clara(process_for_cluster(df_all_compl_no, var_sel_2),
                                   k_range = c(2:10),
                                   n_samples = 500, 
                                   n_sampsize = 1000,
                                   distance_metric = "jaccard")
cluster_clara_j_3 <- cluster_clara(process_for_cluster(df_all_compl_no, var_sel_3),
                                   k_range = c(2:10),
                                   n_samples = 500, 
                                   n_sampsize = 1000,
                                   distance_metric = "jaccard")
cluster_clara_j_4 <- cluster_clara(process_for_cluster(df_all_compl_no, var_sel_4),
                                   k_range = c(2:10),
                                   n_samples = 500, 
                                   n_sampsize = 1000,
                                   distance_metric = "jaccard")

##################### 2. INDICES FOR BEST K SELECTION
#####################################################

best_k_kmeans_1 <- select_best_k(df_all_compl_no,
                                 var_sel_1,
                                 cluster_kmeans_1,
                                 center_type = "centroids")
best_k_kmeans_2 <- select_best_k(df_all_compl_no,
                                 var_sel_2,
                                 cluster_kmeans_2,
                                 center_type = "centroids")
best_k_kmeans_3 <- select_best_k(df_all_compl_no,
                                 var_sel_3,
                                 cluster_kmeans_3,
                                 center_type = "centroids")
best_k_kmeans_4 <- select_best_k(df_all_compl_no,
                                 var_sel_4,
                                 cluster_kmeans_4,
                                 center_type = "centroids")

best_k_clara_1 <- select_best_k(df_all_compl_no,
                                var_sel_1,
                                cluster_clara_1,
                                center_type = "medoids")
best_k_clara_2 <- select_best_k(df_all_compl_no,
                                var_sel_2,
                                cluster_clara_2,
                                center_type = "medoids")
best_k_clara_3 <- select_best_k(df_all_compl_no,
                                var_sel_3,
                                cluster_clara_3,
                                center_type = "medoids")
best_k_clara_4 <- select_best_k(df_all_compl_no,
                                var_sel_4,
                                cluster_clara_4,
                                center_type = "medoids")

best_k_clara_m_1 <- select_best_k(df_all_compl_no,
                                  var_sel_1,
                                  cluster_clara_m_1,
                                  center_type = "medoids")
best_k_clara_m_2 <- select_best_k(df_all_compl_no,
                                  var_sel_2,
                                  cluster_clara_m_2,
                                  center_type = "medoids")
best_k_clara_m_3 <- select_best_k(df_all_compl_no,
                                  var_sel_3,
                                  cluster_clara_m_3,
                                  center_type = "medoids")
best_k_clara_m_4 <- select_best_k(df_all_compl_no,
                                  var_sel_4,
                                  cluster_clara_m_4,
                                  center_type = "medoids")

best_k_clara_j_1 <- select_best_k(df_all_compl_no,
                                  var_sel_1,
                                  cluster_clara_j_1,
                                  center_type = "medoids")
best_k_clara_j_2 <- select_best_k(df_all_compl_no,
                                  var_sel_2,
                                  cluster_clara_j_2,
                                  center_type = "medoids")
best_k_clara_j_3 <- select_best_k(df_all_compl_no,
                                  var_sel_3,
                                  cluster_clara_j_3,
                                  center_type = "medoids")
best_k_clara_j_4 <- select_best_k(df_all_compl_no,
                                  var_sel_4,
                                  cluster_clara_j_4,
                                  center_type = "medoids")


##################### 3. HIERARCHICAL CLUSTERING
#####################################################

###### a. Define distance measures
hclust_args <- expand.grid(dist_method=c("euclidean", "manhattan"), 
                           linkage_method=c("complete", "average"))

###### b. Get clustering
hclust_list_1 <- get_hclust_list(hclust_args, df_all_compl_no, var_sel_1)
hclust_list_2 <- get_hclust_list(hclust_args, df_all_compl_no, var_sel_2)
hclust_list_3 <- get_hclust_list(hclust_args, df_all_compl_no, var_sel_3)
hclust_list_4 <- get_hclust_list(hclust_args, df_all_compl_no, var_sel_4)

###### c. Get plots
#get plots
plot_hclust_1 <- plot_hclust_list(hclust_list_1, hclust_args)
plot_hclust_2 <- plot_hclust_list(hclust_list_2, hclust_args)
plot_hclust_3 <- plot_hclust_list(hclust_list_3, hclust_args)
plot_hclust_4 <- plot_hclust_list(hclust_list_4, hclust_args)

#arrange plots
plot_hclust1_all <- grid.arrange(plot_hclust_1[[1]], plot_hclust_1[[2]],
                                 plot_hclust_1[[3]], plot_hclust_1[[4]])

plot_hclust2_all <- grid.arrange(plot_hclust_2[[1]], plot_hclust_2[[2]], 
                                 plot_hclust_2[[3]], plot_hclust_2[[4]], 
                                 ncol=2)

plot_hclust3_all <- grid.arrange(plot_hclust_3[[1]], plot_hclust_3[[2]],
                                 plot_hclust_3[[3]], plot_hclust_3[[4]])

plot_hclust4_all <- grid.arrange(plot_hclust_4[[1]], plot_hclust_4[[2]],
                                 plot_hclust_4[[3]], plot_hclust_4[[4]])

#save plots
png("plots/plt_hclust1.png", width = 1700, height = 850)
plot(plot_hclust1_all)
dev.off()
png("plots/plt_hclust2.png", width = 1700, height = 850)
plot(plot_hclust2_all)
dev.off()
png("plots/plt_hclust3.png", width = 1700, height = 850)
plot(plot_hclust3_all)
dev.off()
png("plots/plt_hclust4.png", width = 1700, height = 850)
plot(plot_hclust4_all)
dev.off()


##################### 4. IBM SPSS CRITERIA
#####################################################

hclust_list_cut_eucl_compl_1 <- list()
for (i in 1:9){
  hclust_list_cut_eucl_compl_1[[i]] <- cutree(hclust_list_1[[1]], k = i+1)
}
best_k_bic_eucl_compl_1 <- calculate_IC(df_all_compl_no, 
                                        var_sel_1, 
                                        hclust_list_cut_eucl_compl_1)

hclust_list_cut_eucl_compl_2 <- list()
for (i in 1:9){
  hclust_list_cut_eucl_compl_2[[i]] <- cutree(hclust_list_2[[1]], k = i+1)
}
best_k_bic_eucl_compl_2 <- calculate_IC(df_all_compl_no, 
                                        var_sel_2, 
                                        hclust_list_cut_eucl_compl_2)

hclust_list_cut_eucl_compl_3 <- list()
for (i in 1:9){
  hclust_list_cut_eucl_compl_3[[i]] <- cutree(hclust_list_3[[1]], k = i+1)
}
best_k_bic_eucl_compl_3 <- calculate_IC(df_all_compl_no, 
                                        var_sel_3, 
                                        hclust_list_cut_eucl_compl_3)

hclust_list_cut_eucl_compl_4 <- list()
for (i in 1:9){
  hclust_list_cut_eucl_compl_4[[i]] <- cutree(hclust_list_4[[1]], k = i+1)
}
best_k_bic_eucl_compl_4 <- calculate_IC(df_all_compl_no, 
                                        var_sel_4, 
                                        hclust_list_cut_eucl_compl_4)


###############################################################################
############################ CLUSTERING COMPARISON ############################
###############################################################################

##################### 1. CLUSTER DATA
#####################################################

cluster_sel_kmeans_1 <- cluster_kmeans(
  data_processed = process_for_cluster(df_all_compl_no, var_sel_1),
  k_range = 6, 
  n_runs = 100, 
  iter_max = 100)

cluster_sel_kmeans_2 <- cluster_kmeans(
  data_processed = process_for_cluster(df_all_compl_no, var_sel_2),
  k_range = 4, 
  n_runs = 100, 
  iter_max = 1000,
  algo = "MacQueen")

cluster_sel_kmeans_3 <- cluster_kmeans(
  data_processed = process_for_cluster(df_all_compl_no, var_sel_3),
  k_range = 5, 
  n_runs = 100, 
  iter_max = 1000,
  algo = "MacQueen")

##################### 2. CREATE PROFILES
#####################################################

var_names_display_compl_prof <- c("BMI (Kg/m2)", "DBP (mmHg)", "SBP (mmHg)",
                                    "Age (years)", "Age at diagnosis (years)",
                                    "HbA1c (mmol/mol)", "HDL Chol (mmol/L)", 
                                    "LDL Chol (mmol/L)", "Triglycerides (mmol/L)",
                                    "PGS IS", "PGS IR")

var_names_display_compl_ordered <- var_names_display_compl_prof[c(5, 1, 6, 10, 11,
                                                           3, 9, 7,
                                                           4, 2, 8)]
#create
profile_sel_kmeans_1 <- create_profile(
  data = df_all_compl_no, 
  variables = var_sel_1, 
  variables_names = var_names_display_compl_prof, 
  variables_names_ordered = var_names_display_compl_ordered ,
  clustering = cluster_sel_kmeans_1[[1]]$partition, 
  cluster_order = c(1:6),
  cluster_name_order = c(1:6),
  n_col = 4, 
  add_title = FALSE)

profile_sel_kmeans_2 <- create_profile(
  data = df_all_compl_no, 
  variables = var_sel_1, 
  variables_names = var_names_display_compl_prof, 
  variables_names_ordered = var_names_display_compl_ordered,
  clustering = cluster_sel_kmeans_2[[1]]$partition,
  cluster_order = c(4, 3, 1, 2),
  cluster_name_order = c("SIDD", "SIRD", "MOD", "MARD"),
  n_col = 4, 
  add_title = FALSE)

profile_sel_kmeans_3 <- create_profile(
  data = df_all_compl_no, 
  variables = var_sel_1, 
  variables_names = var_names_display_compl_prof,
  variables_names_ordered = var_names_display_compl_ordered ,
  clustering = cluster_sel_kmeans_3[[1]]$partition, 
  cluster_order = c(1, 5, 3, 2, 4),
  cluster_name_order = c("SIDD", "SIRD", "MOD", "MARD1", "MARD2"),
  n_col = 4, 
  add_title = FALSE)

#save
pdf("plots/plot_prof_km_1.pdf")
profile_sel_kmeans_1
dev.off()

pdf("plots/plot_prof_km_2.pdf")
profile_sel_kmeans_2
dev.off()

pdf("plots/plot_prof_km_3.pdf")
profile_sel_kmeans_3
dev.off()

#obtain centres
xtable(t(round(cluster_sel_kmeans_1[[1]]$result$centers, 2)))
xtable(t(round(cluster_sel_kmeans_2[[1]]$result$centers [c(4, 3, 1, 2), ], 2)))
xtable(t(round(cluster_sel_kmeans_3[[1]]$result$centers [c(1, 5, 3, 2, 4), ], 2)))


#table characteristics
profile_sel_kmeans_1_table <- create_profile_char(df_all_compl_no, 
                                                  cluster_sel_kmeans_1[[1]]$partition, c(1:6),
                                                  c(1:6),
                                                  var_sel_1, var_names_display_compl)


profile_sel_kmeans_2_table <- create_profile_char(df_all_compl_no, 
                    cluster_sel_kmeans_2[[1]]$partition, c(4, 3, 1, 2),
                    c("SIDD", "SIRD", "MOD", "MARD"),
                    var_sel_1, var_names_display_compl)

profile_sel_kmeans_3_table <- create_profile_char(df_all_compl_no, 
                                                  cluster_sel_kmeans_3[[1]]$partition, c(1, 5, 3, 2, 4),
                                                  c("SIDD", "SIRD", "MOD", "MARD1", "MARD2"),
                                                  var_sel_1, var_names_display_compl)

##################### 3. BOOTSTRAPS
#####################################################

boot_sel_kmeans_1 <- clusterboot(process_for_cluster(df_all_compl_no, var_sel_1),
                                 k = 6,
                                 B = 200, 
                                 bootmethod = "boot",
                                 clustermethod = kmeansCBI, 
                                 scaling = FALSE,
                                 runs = 100,
                                 iter.max = 1000, 
                                 algorithm = "MacQueen",
                                 criterion="ch",
                                 seed = 1234567)

boot_sel_kmeans_2 <- clusterboot(process_for_cluster(df_all_compl_no, var_sel_2),
                                 k = 4,
                                 B = 300, 
                                 bootmethod = "boot",
                                 clustermethod = kmeansCBI, 
                                 scaling = FALSE,
                                 runs = 100,
                                 iter.max = 1000, 
                                 algorithm = "MacQueen",
                                 criterion="ch",
                                 seed = 1234567)

boot_sel_kmeans_3 <- clusterboot(process_for_cluster(df_all_compl_no, var_sel_3),
                                 k = 5,
                                 B = 200, 
                                 bootmethod = "boot",
                                 clustermethod = kmeansCBI, 
                                 scaling = FALSE,
                                 runs = 100,
                                 iter.max = 1000, 
                                 algorithm = "MacQueen",
                                 criterion="ch",
                                 seed = 1234567)


##################### 4. PCA VISUALISATION
#####################################################

#PCA
res.pca_1 <- prcomp(process_for_cluster(df_all_compl_no, var_sel_1),  scale = FALSE)
res.pca_2 <- prcomp(process_for_cluster(df_all_compl_no, var_sel_2),  scale = FALSE)
res.pca_3 <- prcomp(process_for_cluster(df_all_compl_no, var_sel_3),  scale = FALSE)

#plots
plot_pca_sel_kmeans_1 <-create_cluster_pca_plot(res.pca = res.pca_1, 
                                                clustering = cluster_sel_kmeans_1[[1]]$partition, 
                                                cluster_order = c(1:6),
                                                cluster_name_order = c(1:6),
                                                x_title = "Stand. PC1 (17.6% explained var.)",
                                                y_title = "Stand. PC2 (14.6% explained var.)", 
                                                add_title = TRUE, 
                                                plot_title = "Model 1")

plot_pca_sel_kmeans_2 <- create_cluster_pca_plot(res.pca = res.pca_2, 
                                                 clustering = cluster_sel_kmeans_2[[1]]$partition, 
                                                 cluster_order = c(4, 3, 1, 2),
                                                 cluster_name_order = c(1:4),
                                                 x_title = "Stand. PC1 (17.6% explained var.)",
                                                 y_title = "Stand. PC2 (14.6% explained var.)", 
                                                 add_title = TRUE, 
                                                 plot_title = "Model 2")

plot_pca_sel_kmeans_3 <- create_cluster_pca_plot(res.pca = res.pca_3, 
                                                 clustering = cluster_sel_kmeans_3[[1]]$partition, 
                                                 cluster_order = c(1, 5, 3, 2, 4),
                                                 cluster_name_order = c(1:5),
                                                 x_title = "Stand. PC1 (17.6% explained var.)",
                                                 y_title = "Stand. PC2 (14.6% explained var.)", 
                                                 add_title = TRUE, 
                                                 plot_title = "Model 3")

pdf("plots/pca_sel_kmeans_123.pdf")
grid.arrange(plot_pca_sel_kmeans_1, plot_pca_sel_kmeans_2, plot_pca_sel_kmeans_3, 
             ncol=2)
dev.off()



##################### 5. ALLUVIAL PLOT
#####################################################

###### a. 2 clusters, variables group 2
cluster_sel_kmeans_2_2cl <- cluster_kmeans(
  data_processed = process_for_cluster(df_all_compl_no, var_sel_2),
  k_range = 2,
  n_runs = 100,
  iter_max = 1000,
  algo = "MacQueen")

#alluvial plot
group_colors_2 <- c(SIDD = "#F8766D", SIRD = "#7CAE00", MOD ="#00BFC4", MARD = "#C77CFF")

plot_alluvial_2_2clto2 <- 
  plot_alluvial(cluster_sel_kmeans_2_2cl[[1]]$partition, 
                cluster_sel_kmeans_2[[1]]$partition,
                c(1, 2), c(4, 3, 1, 2),
                "2 clusters", "4 clusters",
                paste0("Cluster", c(1:2)), c("SIDD", "SIRD", "MOD", "MARD"),
                "right", group_colors_2,
                add_title = FALSE)

pdf("plots/plot_alluvial_2_2clto2.pdf")
plot_alluvial_2_2clto2
dev.off()


###### b. 2 clusters, variables group 3

#cluster
cluster_sel_kmeans_3_2cl <- cluster_kmeans(
  data_processed = process_for_cluster(df_all_compl_no, var_sel_3),
  k_range = 2, 
  n_runs = 100, 
  iter_max = 1000,
  algo = "MacQueen")

#alluvial plot
group_colors_3 <- c(SIDD = "#F8766D", SIRD = "#93AA00", 
                    MOD ="#00C19F", MARD1 = "#00B9E3", MARD2 = "#DB72FB")

plot_alluvial_3_2clto3 <- 
  plot_alluvial(cluster_sel_kmeans_3_2cl[[1]]$partition, 
                cluster_sel_kmeans_3[[1]]$partition,
                c(1, 2), c(1, 5, 3, 2, 4),
                "2 clusters", "5 clusters",
                paste0("Cluster", c(1:2)), c("SIDD", "SIRD", 
                                             "MOD", "MARD1", "MARD2"),
                "right", group_colors_3,
                add_title = FALSE)

pdf("plots/plot_alluvial_3_2clto3.pdf")
plot_alluvial_3_2clto3
dev.off()

###### c. Model 2 to 3
plot_alluvial_2to3 <- plot_alluvial(cluster_sel_kmeans_2[[1]]$partition, 
                                    cluster_sel_kmeans_3[[1]]$partition,
                                    c(4, 3, 1, 2),c(1, 5, 3, 2, 4),
                                    "Model 2", "Model 3",
                                    c("SIDD", "SIRD", "MOD", "MARD"), c("SIDD", "SIRD", "MOD", "MARD1", "MARD2"),
                                    "left", group_colors_2,
                                    add_title = FALSE)


pdf("plots/plot_alluvial_2to3.pdf")
plot_alluvial_2to3
dev.off()


###############################################################################
######################### ASSOCIATIONS WITH COMPLICATIONS #####################
###############################################################################

##################### 1. ALL ICD10 CODES
#####################################################

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

##################### 2. MODEL 2
#####################################################

#excl. episodes of CVD before diagnosis of T2D
assoc_kmeans_2_all <- calculate_associations(df_main,
                                   df_all_compl_no,
                                   df_hes,
                                   icd10_cardiovascular,
                                   cluster_sel_kmeans_2[[1]]$partition,
                                   c(4, 3, 1, 2),
                                   c("SIDD", "SIRD", "MOD", "MARD"),
                                   var_sel_2,
                                   FALSE,
                                   FALSE)

#excl. patients with episodes of CVD before diagnosis of T2D 
assoc_kmeans_2_new_only <- calculate_associations(df_main,
                                             df_all_compl_no,
                                             df_hes,
                                             icd10_cardiovascular,
                                             cluster_sel_kmeans_2[[1]]$partition,
                                             c(4, 3, 1, 2),
                                             c("SIDD", "SIRD", "MOD", "MARD"),
                                             var_sel_2,
                                             TRUE,
                                             FALSE)

##################### 3. MODEL 3
#####################################################

#excl. episodes of CVD before diagnosis of T2D
assoc_kmeans_3_all <- calculate_associations(df_main,
                                         df_all_compl_no,
                                         df_hes,
                                         icd10_cardiovascular,
                                         cluster_sel_kmeans_3[[1]]$partition,
                                         c(1, 5, 3, 2, 4),
                                         c("SIDD", "SIRD", "MOD", "MARD1", "MARD2"),
                                         var_sel_3,
                                         FALSE,
                                         FALSE)

#excl. patients with episodes of CVD before diagnosis of T2D 
assoc_kmeans_3_new_only <- calculate_associations(df_main,
                                             df_all_compl_no,
                                             df_hes,
                                             icd10_cardiovascular,
                                             cluster_sel_kmeans_3[[1]]$partition,
                                             c(1, 5, 3, 2, 4),
                                             c("SIDD", "SIRD", "MOD", "MARD1", "MARD2"),
                                             var_sel_3,
                                             TRUE,
                                             FALSE)

##################### 4. ANALYSIS
#####################################################

assoc_kmeans_2_all$cont
assoc_kmeans_2_all$glm_res
assoc_kmeans_2_new_only$cont
assoc_kmeans_2_new_only$glm_res
chisq.test(assoc_kmeans_2$cluster, assoc_kmeans_2$cvd)
sum(assoc_kmeans_2_all$cvd)/length(assoc_kmeans_2_all$cvd)


assoc_kmeans_3_all$cont
assoc_kmeans_3_all$glm_res
assoc_kmeans_3_new_only$cont
assoc_kmeans_3_new_only$glm_res
chisq.test(assoc_kmeans_3$cluster, assoc_kmeans_3$cvd)
sum(assoc_kmeans_3_all$cvd)/length(assoc_kmeans_3_all$cvd)



###############################################################################
############################## EVOLUTION OF CLUSTERS ##########################
###############################################################################

##################### 1. DATA PREPARATION
#####################################################

#eid with diabetes at both assessments
eid_evol <- df_main %>% 
  filter(status_init == "diabetes") %>% 
  filter(status_fu == "diabetes") %>% 
  dplyr::select(eid) %>% 
  unlist()
 
#variables  
var_main_fu <- c("21001-1.0", "4079-1.0", "4080-1.0", "21003-1.0","diabetes_age_diag_fu", "74-1.0")
var_bio_fu <- c("30750-1.0", "30760-1.0", "30780-1.0", "30870-1.0", "30740-1.0")
var_genes_fu <- var_genes

#fu dataframe
df_all_fu <- prepare_tables_3df(
  df_main, c("eid", var_main_fu), 
  df_bio, c("eid", var_bio_fu),
  df_genes,c("eid", var_genes_fu),
  "eid")

colnames(df_all_fu) <- var_all

df_all_compl_fu <- df_all_fu[complete.cases(df_all_fu[ , var_sel_2]), ]

#eid with complete case at follow-up and init, and diabetes at initial and fu
eid_compl_init_fu <- intersect(
  intersect(unlist(df_all_compl_fu$eid), unlist(df_all_compl_no$eid)),
  eid_evol)

length(eid_compl_init_fu)

#time between assessments
dates_ass <- 
  df_main %>% 
  filter(eid %in% eid_compl_init_fu) %>% 
  dplyr::select(eid, `53-0.0`, `53-1.0`)

dates_ass$diff_ass <- time_length(
  interval(ymd(dates_ass$`53-0.0`), ymd(dates_ass$`53-1.0`)), "year")


mean_dates_ass <- mean(dates_ass$diff_ass)

##################### 2. MODEL 2
#####################################################

#prediction
kcca_kmeans_2 <- as.kcca(cluster_sel_kmeans_2[[1]]$result,
                process_for_cluster(df_all_compl_no, var_sel_2))


pred_kmeans_2_init <- predict(kcca_kmeans_2, 
                              process_for_cluster(df_all_compl_no %>% 
                                                          filter(eid %in% eid_compl_init_fu), 
                                                        var_sel_2))

pred_kmeans_2_fu <- predict(kcca_kmeans_2, 
                            process_for_cluster(df_all_compl_fu %>% 
                                                 filter(eid %in% eid_compl_init_fu), 
                                               var_sel_2))
#transition analysis
kmeans_2_inittofu <- 
  get_evol_cluster(pred_init = pred_kmeans_2_init,
                         pred_fu = pred_kmeans_2_fu,
                         cluster_order = c(4, 3, 1, 2),
                         cluster_name_order = c("SIDD", "SIRD", "MOD", "MARD"),
                         color_group = group_colors_2)

length(which(pred_kmeans_2_init == pred_kmeans_2_fu))/length(pred_kmeans_2_init)

xtable(kmeans_2_inittofu$trans_percent)

#save
pdf("plots/plot_redistr_km2.pdf")
kmeans_2_inittofu$trans_plot
dev.off()


##################### 3. MODEL 3
#####################################################
#predictions
kcca_kmeans_3 <- as.kcca(cluster_sel_kmeans_3[[1]]$result,
                         process_for_cluster(df_all_compl_no, var_sel_3))

pred_kmeans_3_init <- predict(kcca_kmeans_3, 
                              process_for_cluster(df_all_compl_no %>% 
                                                          filter(eid %in% eid_evol) %>% 
                                                          filter(eid %in% unlist(df_all_compl_fu$eid)), 
                                                        var_sel_3))

pred_kmeans_3_fu <- predict(kcca_kmeans_3,
                            process_for_cluster(df_all_compl_fu %>% 
                                                        filter(eid %in% eid_evol) %>% 
                                                        filter(eid %in% unlist(df_all_compl_no$eid)), 
                                                      var_sel_3))

#transition analysis
kmeans_3_inittofu <- 
  get_evol_cluster(pred_init = pred_kmeans_3_init,
                   pred_fu = pred_kmeans_3_fu,
                   cluster_order = c(1, 5, 3, 2, 4),
                   cluster_name_order = c("SIDD", "SIRD", "MOD", "MARD1", "MARD2"),
                   color_group = group_colors_3)

length(which(pred_kmeans_3_init == pred_kmeans_3_fu))/length(pred_kmeans_3_init)

xtable(kmeans_3_inittofu$trans_percent)
#save
pdf("plots/plot_redistr_km3.pdf")
kmeans_3_inittofu$trans_plot
dev.off()


###############################################################################
########################## SAVE SELECTED MODELS ###############################
###############################################################################

cluster_ordered_kmeans_2 <- c(4, 3, 1, 2)
cluster_names_ordered_kmeans_2 <- c("SIDD", "SIRD", "MOD", "MARD")

cluster_ordered_kmeans_3 <- c(1, 5, 3, 2, 4)
cluster_names_ordered_kmeans_3 <- c("SIDD", "SIRD", "MOD", "MARD1", "MARD2")

save(kcca_kmeans_2, cluster_ordered_kmeans_2,
     cluster_names_ordered_kmeans_2,
     file = "~/XXX/kmeans_2.RData")

save(kcca_kmeans_3, cluster_ordered_kmeans_3,
     cluster_names_ordered_kmeans_3,
     file = "~/XXX/kmeans_3.RData")
