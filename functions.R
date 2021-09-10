###############################################################################
############################# PATIENTS CLASSIFICATION #########################
###############################################################################

##################### 1. SPLIT PATIENT BY STATUS
#####################################################

split_by_status <- function(df_main, df_bio, df_genes, df_hes, assessment){
  #df_main is the main df to be splitted
  #df_bio is the biochemestry df to be splitted
  #df_genes is the genes df to be splitted
  #df_hes is the HES df to be splitted
  #assessment is the assessment # (0 for initial, 1 for follow-up)
  
  ##### 1. DECLARE PARAMETERS
  ##############################
  
  ###### a. Columns of diabetes diagnosis and age
  
  #col diag diabetes - quest
  if(assessment == 0){
    col_diabetes_quest <- expr(`2443-0.0`)
  } else {
    col_diabetes_quest <- expr(`2443-1.0`) 
  }
  
  #cols age diagnosis - quest
  if(assessment == 0){
    col_age_diabetes_quest <- "2976-0.0"
  } else {
    col_age_diabetes_quest <- "2976-1.0"
  }
  
  #cols diag diabetes - verbal
  if(assessment == 0){
    col_diabetes_verbal <- grep("^20002-0.", names(df_main), value=TRUE) 
  } else {
    col_diabetes_verbal <- grep("^20002-1.", names(df_main), value=TRUE)
  }
  
  #cols age diagnosis - verbal
  if(assessment == 0){
    col_age_diabetes_verbal <- grep("^20009-0.", names(df_main), value=TRUE) 
  } else {
    col_age_diabetes_verbal <- grep("^20009-1.", names(df_main), value=TRUE) 
  }
  
  #col hba1c
  if(assessment == 0){
    col_hba1c <- expr(`30750-0.0`)
  } else {
    col_hba1c <- expr(`30750-1.0`) 
  }
  
  ###### b. Columns of diabetes type
  
  #col gestational diabetes
  if(assessment == 0){
    col_gest <- expr(`4041-0.0`)
  } else {
    col_gest <- expr(`4041-1.0`) 
  }
  
  ###### c. Columns of Medications
  
  #medications
  if(assessment == 0){
    col_medication_ts <- grep("^6177-0.", names(df_main), value=TRUE)
    col_medication_ts2 <- grep("^6153-0.", names(df_main), value=TRUE)
    col_medication_verbal <- grep("^20003-0.", names(df_main), value=TRUE)
  } else {
    col_medication_ts <- grep("^6177-1.", names(df_main), value=TRUE)
    col_medication_ts2 <- grep("^6153-1.", names(df_main), value=TRUE)
    col_medication_verbal <- grep("^20003-1.", names(df_main), value=TRUE) 
  }
  
  #started ins. within one year of diag
  if(assessment == 0){
    col_ins_within_year <- expr(`2986-0.0`) 
  } else {
    col_ins_within_year <- expr(`2986-1.0`) 
  }
  
  ###### d. Codes
  
  #codes condition
  codes_diabetes <- c(1220, 1221, 1222, 1223, 1521)
  codes_gest <- c(1221)
  codes_t1d <- c(1222)
  codes_insipidus <- c(1521)
  
  #codes medication
  codes_insulin_ts <- c(3)
  codes_insulin_verbal <- c(1140883066)
  
  #codes diabetes related complications
  codes_compl_t1d <- c("E100","E101","E102","E103","E104","E105","E106","E107","E108","E109")
  codes_compl_t2d <- c("E110","E111","E112","E113","E114","E115","E116","E117","E118","E119",
                       "E120","E121","E122","E123","E124","E125","E126","E127","E128","E129",
                       "E130","E131","E132","E133","E134","E135","E136","E137","E138","E139",
                       "E140","E141","E142","E143","E144","E145","E146","E147","E148","E149")
  
  #### 2. SPLIT HEALTHY, PRE-DIABETES, DIABETES
  ##############################
  
  ###### a. Diabetic patients
  
  #rows with diagnosis of diabetes - questionnaire
  eid_diabetes_quest <- df_main %>% 
    filter({{col_diabetes_quest}} == 1) %>% 
    dplyr::select(eid)
  
  #rows with diagnosis of diabetes - verbal interview
  eid_diabetes_verbal <- df_main %>% 
    filter_at(col_diabetes_verbal, any_vars(. %in% codes_diabetes)) %>% 
    dplyr::select(eid)
  
  #rows with diabetic values of glycated haemogoblyn
  eid_diabetes_hba1c <- df_bio %>% 
    filter({{col_hba1c}} > 48) %>% 
    dplyr::select(eid)
  
  #rows with diabetic patients (quest, verbal, glycated haemoglobyn)
  eid_diabetes <-
    unique(unlist(c(eid_diabetes_quest, 
                    eid_diabetes_verbal, 
                    eid_diabetes_hba1c)))
  
  ###### b. Pre-diabetic patients
  eid_prediabetes <- df_main %>% 
    filter(!(eid %in% eid_diabetes)) %>% 
    dplyr::select(eid)
  
  eid_prediabetes <- df_bio %>% 
    filter(eid %in% unlist(eid_prediabetes)) %>% 
    filter({{col_hba1c}} >= 39) %>% 
    dplyr::select(eid) %>% 
    unlist()
  
  ###### c. Healthy patients
  eid_healthy <- df_main %>% 
    filter(!(eid %in% eid_diabetes)) %>% 
    filter(!(eid %in% unlist(eid_prediabetes))) %>%
    dplyr::select(eid) %>% 
    unlist()
  
  #### 3. ADD SOURCES FOR DIABETES
  ##############################
  
  #create dataframe for source
  source <- data.frame(eid = df_main$eid)
  if(assessment == 0){
    source$diabetes_class_source_init <- rep(NA, nrow(source))
  } else {
    source$diabetes_class_source_fu <- rep(NA, nrow(source))
  }
  
  # patients identified with questionnaire
  indx_quest <- which(source$eid %in% unlist(eid_diabetes_quest))
  source[indx_quest, 2] <- "quest"
  
  # patients identified with verbal
  indx_verbal <- which(source$eid %in% unlist(eid_diabetes_verbal))
  indx_verbal <- indx_verbal[which(!(indx_verbal %in% indx_quest))]
  source[indx_verbal, 2] <- "verbal"
  
  # patients identified with HbA1c level
  indx_hba1c <- which(source$eid %in% unlist(eid_diabetes_hba1c))
  indx_hba1c <- indx_hba1c[which(!(indx_hba1c %in% indx_quest))]
  indx_hba1c <- indx_hba1c[which(!(indx_hba1c %in% indx_verbal))]
  source[indx_hba1c, 2] <- "hba1c"
  
  #### 4. AGE OF DIAGNOSIS
  ##############################
  
  ###### 0. Create dataframe for age
  age <- data.frame(eid = df_main$eid)
  if(assessment == 0){
    age$diabetes_age_diag_init <- rep(NA, nrow(age))
  } else {
    age$diabetes_age_diag_fu <- rep(NA, nrow(age))
  }
  
  ###### a. With questionnaire
  indx_quest <- which(source[ ,2] == "quest")
  age[indx_quest, 2] <- 
    df_main[indx_quest, col_age_diabetes_quest]
  
  ##remove -3 and -1 values
  age[indx_quest, 2][which(age[indx_quest, 2] == -3)] <- NA
  age[indx_quest, 2][which(age[indx_quest, 2] == -1)] <- NA
  
  ###### b. With verbal interview
  indx_verbal <- which(source[ , 2] == "verbal")
  
  ##matrix with 1 if diabetes diag 0 otherwise
  verbal_all_diag <- df_main[indx_verbal, col_diabetes_verbal]
  
  for (col in col_diabetes_verbal){
    for (i in 1:nrow(verbal_all_diag)){
      if(verbal_all_diag[i, col] %in% codes_diabetes){
        verbal_all_diag[i, col] <- 1
      } else{
        verbal_all_diag[i, col] <- 0
      }
    }
  }
  
  ##age at diagnosis of for each condition
  verbal_all_diag_ages <- df_main[indx_verbal, col_age_diabetes_verbal]
  
  verbal_diag_ages <- rep(NA, nrow(verbal_all_diag))
  
  for (i in 1:nrow(verbal_all_diag)) {
    col <- 1
    while (is.na(verbal_diag_ages[i])) {
      if (verbal_all_diag[i, col] == 1) {
        verbal_diag_ages[i] <- verbal_all_diag_ages[i, col]
      } else {
        col = col + 1
      }
    }
  }
  
  age[indx_verbal, 2] <- verbal_diag_ages
  
  ##remove -3 and -1 values
  age[indx_verbal, 2][which(age[indx_verbal, 2] == -3)] <- NA
  age[indx_verbal, 2][which(age[indx_verbal, 2] == -1)] <- NA
  
  #### 4. SPLIT DIABETES PATIENTS FURTHER
  ##############################
  
  ###### a. Gestational diabetes
  
  #based on touchscreen
  eid_gest_diabetes_ts <- df_main %>% 
    filter(eid %in% eid_diabetes) %>% 
    filter({{col_gest}} == 1) %>% 
    dplyr::select(eid)
  
  #based on verbal interview
  eid_gest_diabetes_verbal <- df_main %>% 
    filter(eid %in% eid_diabetes) %>% 
    filter_at(col_diabetes_verbal, any_vars(. %in% codes_gest)) %>% 
    dplyr::select(eid)
  
  #all
  eid_gest_diabetes <- unique(unlist(c(eid_gest_diabetes_ts,
                                       eid_gest_diabetes_verbal)))
  
  ###### b. Insipidus diabetes
  eid_insipidus_diabetes <- df_main %>% 
    filter(eid %in% eid_diabetes) %>% 
    filter_at(col_diabetes_verbal, any_vars(. %in% codes_insipidus)) %>% 
    dplyr::select(eid)
  
  ###### c. Type I Diabetes main
  
  #based on verbal interview
  eid_t1d_verbal <- df_main%>%
    filter(eid %in% eid_diabetes) %>% 
    filter_at(col_diabetes_verbal, any_vars(. %in% codes_t1d)) %>%
    dplyr::select(eid)
  
  #age <=30, insulin within 1 year of diag. and currently on insulin according to ts
  eid_t1d_age_med_ts1 <- df_main %>%
    filter(eid %in% eid_diabetes) %>% 
    filter({{col_ins_within_year}} == 1) %>%
    filter_at(col_medication_ts, any_vars(. %in% codes_insulin_ts)) %>%
    dplyr::select(eid)
  
  if(assessment == 0){
    eid_t1d_age_med_ts1 <- age %>% 
      filter(eid %in% unlist(eid_t1d_age_med_ts1)) %>% 
      filter(diabetes_age_diag_init <= 30) %>% 
      dplyr::select(eid)
  } else {
    eid_t1d_age_med_ts1 <- age %>% 
      filter(eid %in% unlist(eid_t1d_age_med_ts1)) %>% 
      filter(diabetes_age_diag_fu <= 30) %>% 
      dplyr::select(eid)
  }
  
  #age <=30, insulin within 1 year of diag. and currently on insulin according to ts2
  eid_t1d_age_med_ts2 <- df_main %>%
    filter(eid %in% eid_diabetes) %>% 
    filter({{col_ins_within_year}} == 1) %>% 
    filter_at(col_medication_ts2, any_vars(. %in% codes_insulin_ts)) %>%
    dplyr::select(eid)
  
  if(assessment == 0){
    eid_t1d_age_med_ts2 <- age %>% 
      filter(eid %in% unlist(eid_t1d_age_med_ts2)) %>% 
      filter(diabetes_age_diag_init <= 30) %>% 
      dplyr::select(eid)
  } else {
    eid_t1d_age_med_ts2 <- age %>% 
      filter(eid %in% unlist(eid_t1d_age_med_ts2)) %>% 
      filter(diabetes_age_diag_fu <= 30) %>% 
      dplyr::select(eid)
  }
  
  #age <30, insulin within 1 year of diag. and currently on insulin according to vi
  eid_t1d_age_med_verbal <- df_main %>%
    filter(eid %in% eid_diabetes) %>% 
    filter({{col_ins_within_year}} == 1) %>% 
    filter_at(col_medication_verbal, any_vars(. %in% codes_insulin_verbal)) %>%
    dplyr::select(eid)
  
  if(assessment == 0){
    eid_t1d_age_med_verbal <- age %>% 
      filter(eid %in% unlist(eid_t1d_age_med_verbal)) %>% 
      filter(diabetes_age_diag_init <= 30) %>% 
      dplyr::select(eid)
  } else {
    eid_t1d_age_med_verbal <- age %>% 
      filter(eid %in% unlist(eid_t1d_age_med_verbal)) %>% 
      filter(diabetes_age_diag_fu <= 30) %>% 
      dplyr::select(eid)
  }
  
  #all eid with t1d
  eid_t1d <-
    unique(unlist(c(eid_t1d_verbal, eid_t1d_age_med_ts1, 
                    eid_t1d_age_med_ts2, 
                    eid_t1d_age_med_verbal)))
  
  ###### c. Type I Diabetes - hes
  
  #eid with codes of t2d related compl.
  eid_t2d_compl <- df_hes %>% 
    filter(eid %in% eid_diabetes) %>% 
    filter(diag_icd10 %in% codes_compl_t2d) %>% 
    dplyr::select(eid) %>% 
    unique()
  
  #eid with codes of t1d related compl. only (ie no compl related to other type of diabetes)
  eid_t1d_compl_only <- df_hes %>% 
    filter(eid %in% eid_diabetes) %>% 
    filter(diag_icd10 %in% codes_compl_t1d) %>% 
    filter(!(eid %in% unlist(eid_t2d_compl))) %>% 
    dplyr::select(eid) %>% 
    unique()
  
  ###### d. eid of other diabetes types
  
  #combine all
  eid_other_diabetes <- unique(unlist(c(eid_gest_diabetes,
                                        eid_insipidus_diabetes,
                                        eid_t1d,
                                        eid_t1d_compl_only)))
  
  #### 7. CREATE STATUS FLAG
  ##############################
  
  ###### 0. Create dataframe for status
  status <- data.frame(eid = df_main$eid)
  
  ###### 1. Fill status
  if(assessment == 0){
    status <- status %>% 
      mutate(status_init = ifelse(eid %in% eid_other_diabetes, "other diabetes", 
                                  ifelse(eid %in% eid_diabetes, "diabetes", 
                                         ifelse(eid %in% eid_prediabetes, "prediabetes", "healthy"))))
  } else {
    status <- status %>% 
      mutate(status_fu = ifelse(eid %in% eid_other_diabetes, "other diabetes", 
                                ifelse(eid %in% eid_diabetes, "diabetes", 
                                       ifelse(eid %in% eid_prediabetes, "prediabetes", "healthy"))))
  }
  
  
  #### 8. CREATE DATAFRAME
  ##############################
  result <- inner_join(status, source, by="eid")
  result <- inner_join(result, age, by="eid")
  
  return(result)
}


###############################################################################
############################### DATA PROCESSING ###############################
###############################################################################

##################### 2. PREPARE TABLES - 2 TABLES
#####################################################

prepare_tables_2df <- function(table1, variables1, table2, variables2, join_var){
  #table 1 is the first table to use
  #variables 1 are the names of variables to select for table 1
  #table 2 is the second table to use
  #variables 2 are the names of variables to select for table 2
  #join_var is the variable to use for the inner join
  
  t1 <- table1 %>% 
    dplyr::select(variables1)
  t2<- table2 %>% 
    dplyr::select(variables2)
  
  res <- inner_join(t1, t2, by = join_var)
  return(res)
}

##################### 3. PREPARE TABLES - 3 TABLES
#####################################################

prepare_tables_3df <- function(table1, variables1, table2, variables2, 
                               table3, variables3, join_var){
  #table 1 is the first table to use
  #variables 1 are the names of variables to select for table 1
  #table 2 is the second table to use
  #variables 2 are the names of variables to select for table 2
  #table 3 is the second table to use
  #variables 3 are the names of variables to select for table 3
  #join_var is the variable to use for the inner join
  
  t1 <- table1 %>% 
    dplyr::select(variables1)
  t2 <- table2 %>% 
    dplyr::select(variables2)
  t3 <- table3 %>% 
    dplyr::select(variables3)
  
  res <- inner_join(t1, t2, by = join_var)
  res <- inner_join(res, t3, by = join_var)
  return(res)
}

##################### 4. REMOVE OUTLIERS
#####################################################

remove_outliers <- function(data, variables, bounds){
  #data is the dataframe
  #variables is a vector with the name of variables for outliers removing
  #bounds is a matrix with lower bounds and upper bounds (ie remove obs out of bounds)
  
  X <- data
  
  #filter for each variable
  for (var in variables){
    selected_indx <- ((X[ ,var] > bounds[var, 1]) & (X[ ,var] < bounds[var, 2]))
    X <- X [selected_indx, ]
  }
  
  return(X)
}

##################### 5. PROCESS DATA FOR CLUSTER
#####################################################

process_for_cluster <- function(data, variables){
  #data is the dataframe
  #variables is a vector with the name of variables to use for clustering
  
  X <- data[, variables]
  X <- scale(X)
  
  return (X)
}

#################### 6. CLEAN FOR LATEX
#####################################################

clean_for_latex <- function(matrix){
  #matrix is the matrix to be cleaned
  
  res <- matrix
  
  for(j in 1:ncol(res)){
    #select col
    col <- res[ ,j]
    
    #apply changes
    col[which(col == "age_at_assessment")] <- "Age"
    col[which(col == "diabetes_age_diag")] <- "Age at diagnosis"
    col[which(col == "dbp")] <- "DBP"
    col[which(col == "sbp")] <- "SBP"
    col[which(col == "bmi")] <- "BMI"
    col[which(col == "hba1c")] <- "HbA1c"
    col[which(col == "triglycerides")] <- "Triglycerides"
    col[which(col == "hdl_chol")] <- "HDL Chol"
    col[which(col == "ldl_direct")] <- "LDL Chol"
    col[which(col == "pgs_secr_w")] <- "PGS-IS"
    col[which(col == "pgs_res_w")] <- "PGS-IR"
    
    #save changes
    res[ ,j] <- col
  }
  
  
  return(res) 
  
}


###############################################################################
################################### CLUSTERING ################################
###############################################################################

##################### 7. CLUSTER KMEANS
#####################################################

cluster_kmeans <- function(data_processed, k_range, n_runs, 
                           iter_max = 10, algo = "Hartigan-Wong"){
  #data_processed is the dataframe with only variables to be clustered scaled
  #k_range is the range to be tested
  #n_runs is the nb of random initializations
  #iter_max is the max number of iterations
  #algo is the algo to use
  
  #create result list
  res <- list()
  
  for (i in 1:length(k_range)){
    res[[i]] <- kmeansCBI(data_processed,
                          k = k_range[i],
                          scaling = FALSE,
                          runs = n_runs,
                          criterion="ch",
                          iter.max = iter_max, 
                          algorithm = algo)
  }
  
  names(res) <- k_range
  return (res)
}

##################### 8. CLUSTER PAM
#####################################################

cluster_pam <- function(data_processed, k_range, n_starts){
  #data_processed is the dataframe with only variables to be clustered scaled
  #k_range is the range to be tested
  #n_start is the number of different start
  
  #create result list
  res <- list()
  
  for (i in 1:length(k_range)){
    res[[i]] <-  pamkCBI(data_processed,
                         k = k_range[i],
                         scaling = FALSE, 
                         usepam = TRUE, 
                         diss = FALSE,
                         criterion="ch",
                         metric = "euclidean",
                         nstart = n_starts)
  }
  
  names(res) <- k_range
  return (res)
}

##################### 9. CLUSTER CLARA
#####################################################

cluster_clara <- function(data_processed, k_range, n_samples, 
                          n_sampsize, distance_metric = "euclidean"){
  #data_processed is the dataframe with only variables to be clustered scaled
  #k_range is the range to be tested
  #n_samples is the number of samples
  #n_sampsize is the size of each sample
  #distance_metric is the distance to use
  
  #create result list
  res <- list()
  
  for (i in 1:length(k_range)){
    res[[i]] <-  claraCBI(data_processed,
                          k = k_range[i],
                          usepam = FALSE, 
                          samples = n_samples,
                          sampsize = n_sampsize,
                          diss = FALSE, 
                          metric = distance_metric)
  }
  
  names(res) <- k_range
  return (res)
}

##################### 10. GET HIERARCHICAL CLUSTERING LIST
#####################################################

get_hclust_list <- function(hclust_args, data, variables){
  #hclust_args is the list of distance x linkage to use
  #data is the dataset to use
  #variables is the variables set to use
  
  X <- process_for_cluster(data, variables)
  
  get_clust_assignment <- function(dist_method, linkage_method) {
    hclust(dist(X, method=dist_method), method=linkage_method)
  }
  
  res <- pmap(hclust_args, get_clust_assignment)
  
  return(res)
}

#################### 11. PLOT HIERACRHCICAL CLUSTERING
#####################################################

plot_hclust_list <- function(hclust_list, hclust_args){
  #hclust list is a list of hclustering as obtained in get_hclust_list
  #hclust_args is the list of distance x linkage to use
  
  res <- list()
  
  for(i in 1:length(hclust_list)) {
    #print status
    log_print(i)
    
    res[[i]] <- hclust_list[[i]] %>% 
      as.dendrogram() %>% 
      dendro_data(type = "rectangle") %>% 
      segment() %>% 
      ggplot() + 
      geom_segment(aes(x = x, y = y, xend = xend, yend = yend))+
      ggtitle(paste(hclust_args[i, 1], "distance and",  hclust_args[i, 2], "linkage" ))
  }
  
  return(res)
}

#################### 12. CHANGE ORDER OF CLUSTERING
#####################################################

reorder_clustering <- function(clustering, cluster_order, cluster_name_order){
  #clustering is the partition
  #cluster_order is the reorder clustering
  #cluster_name_order is the names to be given to the clusters
  
  clustering_ordered <- c(NA, length(clustering))
  
  for(k in 1:max(clustering)){
    clustering_ordered[which(clustering == k)] <- cluster_name_order[
      which(cluster_order == k)]
  }
  
  return(clustering_ordered)
}


###############################################################################
####################### FEATURES AND CLUSTERING SELECTION #####################
###############################################################################

#################### 13. CREATE BIPLOT
#####################################################
create_cluster_pca_plot <- function(res.pca, 
                                    clustering, cluster_order, cluster_name_order,
                                    x_title, y_title, add_title , plot_title = NULL){
  #res.pca is the pca of the variable
  #clustering is the clustering
  #cluster_order is the order of the clusters
  #cluster_name_order is the name of the clusters in order
  #x_title is the name of the x axis
  #y_title is the name of the y axis
  #add_title is wether yes or no a title shoule be added
  #plot_title is the global name of the title
  
  clustering_ordered <- factor(
    reorder_clustering(clustering, cluster_order, cluster_name_order), 
    levels = cluster_name_order)
  
  
  
  res <- ggbiplot(res.pca, ellipse = TRUE,
                  groups = clustering_ordered,
                  var.axes=FALSE) + 
    xlab(x_title) + 
    ylab(y_title) + 
    theme(axis.title.x = element_text(size = 9),
          axis.title.y = element_text(size = 9),
          legend.position = "bottom") +
    guides(color = guide_legend(title = "Cluster"))
  
  if(add_title == TRUE){
    res <- res + ggtitle(plot_title)
  }
  
  return(res)
}

##################### 14. GET VALIDATION STAT
#####################################################

get_validation_stat <- function(data_processed, clustering, 
                                index_measures = c("Silhouette","Calinski_Harabasz", "Dunn")){
  #data_processed is the dataframe with only variables to be clustered scaled
  #clustering is the clustering result
  
  #create result vector
  res <- rep(NA, 3)
  
  res <- intCriteria(data_processed, 
                     clustering, 
                     index_measures)
  return (res)
}

##################### 15. FORWARD STEPWISE CLUSTERING
#####################################################
cluster_stepwise <- function(data, variables, k_range, crit, algo, n_runs = NULL,
                             n_starts = NULL, n_samples = NULL, n_sampsize = NULL){
  #data is the dataframe
  #variables is a vector with the name of variables to use for clustering
  #k_range is the range to be tested
  #crit is the criteria to be used for selection of best variable
  #algo is the algo to use (kmeans, pam, clara)
  #n_runs is the nb of random initializations
  #n_starts is the number of different starts
  #n_samples is the number of samples
  #n_sampsize is the size of each sample
  
  ###### a. Create dataframe
  X <- process_for_cluster(data, variables)
  
  
  ###### b. Create parameters
  var_indx <- rep(NA, length(variables))
  var_out_all <- list()
  var_in_all <- list()
  
  ###### c. Create results
  res_var_all <- list()
  var_sel_all <- list()
  
  ###### d. Initialize loop
  #print log
  log_print(paste0("INIT"))
  
  for(i_k in 1:length(k_range)){
    
    #print log
    log_print(paste0("k: ",i_k))
    
    #initialize values
    var_out_all[[i_k]] <- variables
    res_var_all[[i_k]] <- matrix(rep(rep(NA, length(variables)), 
                                     length(variables)), 
                                 nrow=length(variables))
    var_sel_all[[i_k]] <- matrix(rep(rep(0, length(variables)), 
                                     length(variables)), 
                                 nrow=length(variables))
    
    #loop over variables
    for (i_var in 1:length(var_out_all[[i_k]])){
      
      #print log
      log_print(paste0("var: ",i_var))
      
      #cluster based on variable selected
      var_current <- var_out_all[[i_k]][i_var]
      X_current <- X[ ,var_current]
      
      if (algo=="kmeans"){
        cl_current <- cluster_kmeans(X_current, k_range[i_k], n_runs)
      } else if(algo == "pam"){
        cl_current <- cluster_pam(X_current, k_range[i_k], n_starts)
      } else {
        cl_current <- cluster_clara(X_current, k_range[i_k], n_samples, n_sampsize)
      }
      
      #fill index
      var_indx[i_var] <- intCriteria(as.matrix(X_current), 
                                     cl_current[[1]]$partition, crit)[[1]]
    }
    
    #remove unwanted values
    var_indx[which(var_indx=="NaN")] <- NA
    
    #fill result matrix
    res_var_all[[i_k]][ ,1] <- var_indx
    indx_best_var <- which(var_indx == max(var_indx, na.rm=TRUE))
    var_sel_all[[i_k]][indx_best_var, 1] <- 1
    
    #update variables out and in
    var_out_all[[i_k]] <- var_out_all[[i_k]][-indx_best_var]
    var_in_all[[i_k]] <- variables[indx_best_var]
  }
  
  #loop over all steps
  for(step in 2:length(variables)){
    #print log
    log_print(paste0("STEP: ",step))
    
    for(i_k in 1:length(k_range)){
      #print log
      log_print(paste0("k: ",i_k))
      
      #clean var_indx
      var_indx <- res_var_all[[i_k]][, step]
      
      for (i_var in 1:length(var_out_all[[i_k]])){
        
        #print log
        log_print(paste0("var: ",i_var))
        
        #cluster based on variable selected
        var_current <- c(var_in_all[[i_k]], var_out_all[[i_k]][i_var])
        X_current <- X[ ,var_current]
        
        if (algo=="kmeans"){
          cl_current <- cluster_kmeans(X_current, k_range[i_k], n_runs)
        } else if(algo == "pam"){
          cl_current <- cluster_pam(X_current, k_range[i_k], n_starts)
        } else {
          cl_current <- cluster_clara(X_current, k_range[i_k], n_samples, n_sampsize)
        }
        
        #fill index
        indx_var <- which(variables == var_out_all[[i_k]][i_var])
        var_indx[indx_var] <- intCriteria(as.matrix(X_current), 
                                          cl_current[[1]]$partition, crit)[[1]]
      }
      
      #remove unwanted values
      var_indx[which(var_indx=="NaN")] <- NA
      
      #fill result matrix
      res_var_all[[i_k]][ ,step] <- var_indx
      indx_best_var <- which(var_indx == max(var_indx, na.rm=TRUE))
      var_sel_all[[i_k]][indx_best_var, step] <- 1
      
      #update variables out and in
      indx_to_remove <- which(var_out_all[[i_k]] == variables[indx_best_var])
      var_out_all[[i_k]] <- var_out_all[[i_k]][-indx_to_remove]
      var_in_all[[i_k]] <- c(var_in_all[[i_k]], variables[indx_best_var])
      
    }
  }
  
  #name list
  names(res_var_all) <- paste0(k_range,"_k")
  names(var_sel_all) <- paste0(k_range,"_k")
  
  for(i_k in 1:length(k_range)){
    rownames(res_var_all[[i_k]]) <- variables
    rownames(var_sel_all[[i_k]]) <- variables
  }
  
  #return results
  res <- list(res_var_all, var_sel_all)
  return(res)
}

##################### 16. COMBINE FORWARD STEPWISE CLUSTERING RESULTS
#####################################################
combine_cluster_stepwise <- function(list1, list2){
  #list1 is the first list of results of cluster_stepwise
  #list2 is the second list of results of cluster_stepwise
  
  #get size of lists
  n_el1 <- lengths(list1)[1]
  n_el2 <- lengths(list2)[1]
  
  #create results
  res <- list1
  
  #add elements to list1
  for(i in 1:n_el2){
    res[[1]][[n_el1 + i]] <- list2[[1]][[i]]
    res[[2]][[n_el1 + i]] <- list2[[2]][[i]]
  }
  
  return(res)
}

##################### 17. CREATE SUMMARY FORWARD STEPWISE
#####################################################
create_summary_stepwise <- function(data, variables, n_round){
  #data is the result from cluster_stepwise function
  #variables is the vector of variables used in cluster_stepwise
  #n_round is the number of digits for rounding
  
  #get size of list
  n_el <- lengths(data)[1]
  n_var <- length(variables)
  
  #create results
  res <- matrix(NA, nrow = n_el * 2, ncol = n_var)
  
  #fill matrix
  for (i in 1:n_el){
    for (j in 1:n_var){
      indx_max <- which(data[[2]][[i]][ ,j] == 1)
      res[2*i - 1, j] <- variables[indx_max]
      res[2*i, j] <- round(data[[1]][[i]][indx_max, j], n_round)
    }
  }
  return(res)
}

##################### 18. BACKWARD ONE STEP STEPWISE CLUSTERING
#####################################################

cluster_one_step_backward <- function(data, variables, k_range, algo, n_runs = NULL,
                                      n_starts = NULL, n_samples = NULL, n_sampsize = NULL){
  #data is the datatframe
  #variables are the columns to use
  #k_range is the range to be tested
  #algo is the algo to use (kmeans, pam, clara)
  #n_runs is the nb of random initializations
  #n_starts is the number of different starts
  #n_samples is the number of samples
  #n_sampsize is the size of each sample
  
  #create dataframe
  X <- process_for_cluster(data, variables)
  
  
  #matrix to save result
  res <- matrix(NA, nrow = length(variables) + 1, ncol = 2* length(k_range))
  rownames(res) <- c("all", paste0(variables, "_out"))
  grid_param <- expand.grid(c("sil", "ch"), paste0(k_range,"_k"))
  colnames(res) <- paste0(grid_param [ ,2], ".", grid_param [ ,1])
  
  #initialize
  #print log
  log_print(paste0("INIT"))
  
  X_current <- X
  ## loop over k
  for(i_k in 1:length(k_range)){
    
    #print log
    log_print(i_k)
    
    #cluster data
    if (algo=="kmeans"){
      cl_current <- cluster_kmeans(X_current, k_range[i_k], n_runs)
    } else if(algo == "pam"){
      cl_current <- cluster_pam(X_current, k_range[i_k], n_starts)
    } else {
      cl_current <- cluster_clara(X_current, k_range[i_k], n_samples, n_sampsize)
    }
    
    #get validation stat
    val_stat_current <- get_validation_stat(as.matrix(X_current), 
                                            cl_current[[1]]$partition)
    res[1, 2*i_k - 1] <- val_stat_current[[1]]
    res[1, 2*i_k] <- val_stat_current[[2]]
  }
  
  
  #when one variable removed
  for (i_var in 1:length(variables)){
    
    #print log
    log_print(i_var)
    
    X_current <- X[ , -i_var]
    
    ## loop over k
    for(i_k in 1:length(k_range)){
      
      #print log
      log_print(i_k)
      
      #cluster data
      if (algo=="kmeans"){
        cl_current <- cluster_kmeans(X_current, k_range[i_k], n_runs)
      } else if(algo == "pam"){
        cl_current <- cluster_pam(X_current, k_range[i_k], n_starts)
      } else {
        cl_current <- cluster_clara(X_current, k_range[i_k], n_samples, n_sampsize)
      }
      
      #get validation stat
      val_stat_current <- get_validation_stat(as.matrix(X_current), 
                                              cl_current[[1]]$partition)
      res[i_var + 1, 2*i_k - 1] <- val_stat_current[[1]]
      res[i_var + 1, 2*i_k] <- val_stat_current[[2]]
    }
    
  }
  
  return(res)
}

##################### 19. CALCULATE GD33
#####################################################

calculate_GD33 <- function(data_processed, clustering, centers){
  #data_processed is the datatframe
  #clustering is a vector with the membership
  #centers is the vector of the centers for each cluster (either medoids or centroids)
  
  #cluster_number
  n_clusters <- max(clustering)
  
  #Compactness with centroids
  upper_delta_3_centroids <- rep(NA, n_clusters)
  
  for (cluster in 1:n_clusters){
    selected_obs <- data_processed[which(clustering == cluster), ]
    centroid <- colMeans(selected_obs)
    
    current_dist <- apply(selected_obs, 1, function(x)sqrt(sum((x - centroid)^2)))
    upper_delta_3_centroids[cluster] <- 2/nrow(selected_obs) * sum(current_dist)
  }
  
  
  
  #Compactness
  upper_delta_3 <- rep(NA, n_clusters)
  
  for (cluster in 1:n_clusters){
    selected_obs <- data_processed[which(clustering == cluster), ]
    current_dist <- apply(selected_obs, 1, function(x)sqrt(sum((x - centers[cluster,])^2)))
    upper_delta_3[cluster] <- 2/nrow(selected_obs) * sum(current_dist)
  }
  
  #GD33 from intCriteria function
  indx_no_corr <- intCriteria(as.matrix(data_processed),
                              clustering,
                              "GDI33")
  
  res <- indx_no_corr[[1]] * max(upper_delta_3_centroids)/max(upper_delta_3)
  
  return(res)
}

##################### 20. CALCULATE GD43
#####################################################

calculate_GD43 <- function(data_processed, clustering, centers){
  #data_processed is the datatframe
  #clustering is a vector with the membership
  #centers is the vector of the centers for each cluster (either medoids or centroids)
  
  #cluster_number
  n_clusters <- max(clustering)
  
  #Separation
  lower_delta_4 <- dist(centers, method="euclidean")
  
  #Compactness
  upper_delta_3 <- rep(NA, n_clusters)
  
  for (cluster in 1:n_clusters){
    selected_obs <- data_processed[which(clustering == cluster), ]
    current_dist <- apply(selected_obs, 1, function(x)sqrt(sum((x - centers[cluster,])^2)))
    upper_delta_3[cluster] <- 2/nrow(selected_obs) * sum(current_dist)
  }
  
  res <- min(lower_delta_4) / max(upper_delta_3)
  
  return(res)
}

##################### 21. CALCULATE GD53
#####################################################

calculate_GD53 <- function(data_processed, clustering, centers){
  #data_processed is the datatframe
  #clustering is a vector with the membership
  #centers is the vector of the centers for each cluster (either medoids or centroids)
  
  #cluster_number
  n_clusters <- max(clustering)
  
  #distance list
  dist_list <- list()
  
  
  #calculate distances to centers
  for (cluster in 1:n_clusters){
    selected_obs <- data_processed[which(clustering == cluster), ]
    current_dist <- apply(selected_obs, 1, function(x)sqrt(sum((x - centers[cluster,])^2)))
    dist_list[[cluster]] <- current_dist
  }
  
  #Separation
  lower_delta_5 <- matrix(NA, nrow = n_clusters, ncol = n_clusters)
  
  for(cluster1 in 1:n_clusters){
    for(cluster2 in 1:n_clusters){
      if(cluster1 != cluster2){
        lower_delta_5[cluster1, cluster2] <- 1/(length(dist_list[[cluster1]]) + length(dist_list[[cluster2]]))*
          (sum(dist_list[[cluster1]]) + sum(dist_list[[cluster2]]))
      }
    }
  }
  
  
  #Compactness
  upper_delta_3 <- rep(NA, n_clusters)
  
  for (cluster in 1:n_clusters){
    upper_delta_3[cluster] <- 2/length(dist_list[[cluster]]) * sum(dist_list[[cluster]])
  }
  
  res <- min(lower_delta_5, na.rm = TRUE) / max(upper_delta_3)
  
  return(res)
}

##################### 22. CALCULATE VARIANT OF DB
#####################################################

calculate_DB_var <- function(data_processed, clustering, centers){
  #data_processed is the datatframe
  #clustering is a vector with the membership
  #centers is the vector of the centers for each cluster (either medoids or centroids)
  
  #cluster_number
  n_clusters <- max(clustering)
  
  #distances to centers
  dist_to_centers <- rep(NA, n_clusters)
  
  for (cluster in 1:n_clusters){
    selected_obs <- data_processed[which(clustering == cluster), ]
    current_dist <- apply(selected_obs, 1, function(x)sqrt(sum((x - centers[cluster,])^2)))
    dist_to_centers[cluster] <- sum(current_dist)/nrow(selected_obs)
  }
  
  ratio_vector <- rep(NA, n_clusters)
  sum_vector <- rep(NA, n_clusters-1)
  
  for(cluster1 in 1:n_clusters){
    
    #numerator
    sum_vector <- dist_to_centers + dist_to_centers[cluster1]
    sum_vector <- sum_vector[-cluster1]
    
    #denominator
    selected_centers <- centers[-cluster1, ]
    if(n_clusters == 2){
      all_distances_to_centers <- sqrt(sum((selected_centers - centers[cluster1,])^2))
    }else {all_distances_to_centers <- apply(selected_centers, 1, 
                                             function(x)sqrt(sum((x - centers[cluster1,])^2)))
    }
    
    ratio_vector[cluster1] <- max(sum_vector)/min(all_distances_to_centers)
  }
  
  res <- 1/n_clusters * sum(ratio_vector)
  
  return(res)
}

##################### 23. SELECT BEST K BASED ON MEASURES
#####################################################

select_best_k <- function(data, variables, cluster_object, center_type){
  #data is the datatframe
  #variables are the columns to use
  #cluster_object is the clustering object obtained from cluster_kmeans, cluster_pam or cluster_clara
  #center_type is the type of centers (centroids or medoids)
  
  #create dataframe
  X <- process_for_cluster(data, variables)
  
  #index_measures
  index_measures = c("Silhouette",
                     "Davies_Bouldin",
                     "Calinski_Harabasz",
                     "GD33",
                     "GD43",
                     "GD53")
  
  #number of k
  n_k <- length(cluster_object)
  
  #matrix to save result
  res <- matrix(NA, nrow = n_k, ncol = length(index_measures))
  rownames(res) <- names(cluster_object)
  colnames(res) <- index_measures
  
  ## loop over k
  for(i_k in 1:n_k){
    
    #print log
    log_print(i_k)
    
    #get silhouette
    log_print("Silhouette")
    current_sil <- intCriteria(as.matrix(X), 
                               cluster_object[[i_k]]$partition,
                               "Silhouette")[[1]]
    
    #get DB
    log_print("DB")
    if(center_type == "medoids"){
      current_db <- calculate_DB_var(X,
                                     cluster_object[[i_k]]$partition,
                                     cluster_object[[i_k]]$result$medoids)
    } else {
      current_db <- calculate_DB_var(X,
                                     cluster_object[[i_k]]$partition,
                                     cluster_object[[i_k]]$result$centers)
    }
    
    
    #get CH
    log_print("CH")
    if(center_type == "medoids"){
      current_ch <- index.G1(X, 
                             cluster_object[[i_k]]$partition, 
                             d = dist(X, method="euclidean"),
                             centrotypes = center_type)
    } else {
      current_ch <- index.G1(X, 
                             cluster_object[[i_k]]$partition,
                             centrotypes = center_type)
    }
    
    #get GDI33
    log_print("GD33")
    if(center_type == "medoids"){
      current_gd33 <- calculate_GD33(X,
                                     cluster_object[[i_k]]$partition,
                                     cluster_object[[i_k]]$result$medoids)
    } else {
      current_gd33 <- intCriteria(as.matrix(X), 
                                  cluster_object[[i_k]]$partition,
                                  "GDI33")[[1]]
    }
    
    #get GDI43
    log_print("GD43")
    if(center_type == "medoids"){
      current_gd43 <- calculate_GD43(X,
                                     cluster_object[[i_k]]$partition,
                                     cluster_object[[i_k]]$result$medoids)
    } else {
      current_gd43 <- intCriteria(as.matrix(X), 
                                  cluster_object[[i_k]]$partition,
                                  "GDI43")[[1]]
    }
    
    #get GDI53
    log_print("GD53")
    if(center_type == "medoids"){
      current_gd53 <- calculate_GD53(X,
                                     cluster_object[[i_k]]$partition,
                                     cluster_object[[i_k]]$result$medoids)
    } else {
      current_gd53 <- intCriteria(as.matrix(X), 
                                  cluster_object[[i_k]]$partition,
                                  "GDI53")[[1]]
    }
    
    
    
    #fill result matrix
    log_print("Fill matrix")
    res[i_k, 1] <- current_sil
    res[i_k, 2] <- current_db
    res[i_k, 3] <- current_ch
    res[i_k, 4] <- current_gd33
    res[i_k, 5] <- current_gd43
    res[i_k, 6] <- current_gd53
    
  }
  
  
  return(res)
}

#################### 24. SELECT BEST K BASED ON BIC CRITERIA (SPSS)
#####################################################

calculate_IC <- function(data, variables, clustering_list){
  #data is the data object
  #variables is the variables
  #clustering_list is the list of clusters
  
  #1. CALCULATE BIC
  X <- process_for_cluster(data, variables)
  
  #parameters
  n_clustering <- length(clustering_list)
  
  #result
  res_ic <- matrix(NA, ncol = 2, nrow = n_clustering + 1)
  colnames(res_ic) <- c("BIC", "AIC")
  
  #loop over all clustering
  for(clust in 1:n_clustering){
    
    #obtain clustering
    current_clustering <- clustering_list[[clust]]
    #obtain number of different clusters
    current_n_clust <- max(clustering_list[[clust]])
    #obtain m 
    m <- 2 * current_n_clust * ncol(X) 
    
    #obtain zeta
    zeta <- rep(NA, current_n_clust)
    for (k in 1:current_n_clust){
      n_k <- length(which(current_clustering == k))
      selected_obs <- X[which(current_clustering == k), ]
      sd_k <- colSds(selected_obs)
      zeta[k] <- -n_k/2 * sum(log(1 + sd_k^2))
    }
    
    #obtain bic
    res_ic[clust + 1, 1] <- -2* sum(zeta[1:current_n_clust]) + m * log(nrow(X))
    
    #obtain aic
    res_ic[clust + 1, 2] <- -2* sum(zeta[1:current_n_clust]) + m * 2
    
  }
  
  #2. SELECT BEST MODEL
  
  #calculate bic 1
  BIC_1 <- 2 * nrow(X)/2 * ncol(X) * log(2) + 2 * 1 * ncol(X) * log(nrow(X))
  
  #diff vector
  diff_bic <- rep(NA, n_clustering + 1)
  diff_bic[1] <- BIC_1 - res_ic[2, 1]
  
  #define R2
  R2 <- c()
  
  #if diffbic1 <0 best clust is 1
  if(diff_bic[1]<0){
    best_clust <- 1
  }else{
    
    #fill diff vector
    for (i in 2:(n_clustering)){
      diff_bic[i] <- res_ic[i, 1] - res_ic[i+1, 1] 
    }
    
    #R1 vector
    R1 <- rep(NA, n_clustering + 1)
    for (i in 1:n_clustering){
      R1[i] <- diff_bic[i]/diff_bic[1]
    }
    
    #initialisation
    if(length(which(R1 <0.04)) == 0){
      init_best_clust <- length(R1)
    } else{
      init_best_clust <- min(which(R1 <0.04))
    }
    
    
    #refining of best clust
    
    ##calculate distances
    dist_min <- rep(NA, init_best_clust)
    for(i in 2:init_best_clust){
      dist_centers <- cls.scatt.data(X, 
                                     clustering_list[[i-1]], 
                                     dist="euclidean")$intercls.centroid
      dist_centers[which(dist_centers == 0)] <- NA
      dist_min[i] <- min(dist_centers, na.rm = TRUE)
      
    }
    
    ##fill R2
    R2 <- rep(NA, init_best_clust)
    for(i in 2:init_best_clust - 1){
      R2[i] <- dist_min[i]/dist_min[i+1]
    }
    
    #compare the two largest
    max1 <- sort(R2, decreasing = TRUE)[1]
    max2 <- sort(R2, decreasing = TRUE)[2]
    
    if(max1>1.15*max2){
      best_clust <- max(which(R2 == max1))
    } else{
      best_clust <- max(max(which(R2 == max1)), max(which(R2 == max2)))
    }
    
  }
  
  res <- list(ic = res_ic, 
              best_clust = best_clust, 
              bic_1 = BIC_1,
              diff_bic = diff_bic, 
              R2 = R2)
  return(res)
  
}


###############################################################################
########################## CLUSTERS CHARACTERISTICS ###########################
###############################################################################

##################### 25. CREATE PROFILES
#####################################################

create_profile <- function(data, variables, variables_names, variables_names_ordered,
                           clustering, cluster_order, cluster_name_order, 
                           n_col, model_name = NULL, add_title){
  #data is the dataframe
  #variables are the actual names of variables in data to use in the profile
  #variables_names are the names of variables to show on the plot
  #variables_names_ordered are the names of the variables but ordered in the way to be shown
  #clustering is the vector of clustering
  #cluster_order is the order to be used for the clusters
  #cluster_name_order is the order of the clusters ordered
  #n_col is the number of columns for the facet wrap graph
  #model_name is the name of the model to be displayed as subtitle
  #add_title whether yes or no to add the title
  
  X <- data[, variables]
  
  #reorder clustering
    clustering_ordered <- reorder_clustering(clustering,
                                             cluster_order,
                                             cluster_name_order)
  
  #create data frame
  X$cluster <- factor(clustering_ordered, levels = cluster_name_order)
  colnames(X) <- c(variables_names, "Cluster")
  
  plt <- X %>%
    pivot_longer(variables_names, names_to  = "variable",
                 values_to = "Value") %>% 
    ggplot(aes(x = `variable`, y = `Value`, fill = `Cluster`)) + 
    geom_boxplot() +
    facet_wrap(~factor(`variable`,
                       levels = variables_names_ordered), 
               scale ="free", ncol = n_col) +
    scale_x_discrete(name ="",
                     breaks = c(1:max(as.numeric(X$Cluster))))+
    theme(legend.position="bottom", 
          plot.title = element_text(face = "bold")) 
  
  if(add_title == TRUE){
    plt <- plt + ggtitle(label = "Profile by cluster",
                         subtitle = model_name)
  }
  return(plt)
}

##################### 26. GENERATE SUMMARY STATISTICS
#####################################################
generate_sum_stat <- function(data, variables, variables_names){
  #data is the datatframe
  #variables are the actual names of variables in dataframe
  #variables_names are the names of variables to show in res table
  
  #calculate summary stats
  stat <- stat.desc(data[ ,variables])
  
  #create result matrix
  res <- matrix(NA, ncol = length(variables), nrow = 5)
  res <- as.data.frame(res)
  colnames(res) <- variables_names
  rownames(res) <- c("Mean", "Median", "Min", "Max", "SD")
  
  #fill result matrix
  res[1, ] <- stat["mean", ]
  res[2, ] <- stat["median", ]
  res[3, ] <- stat["min", ]
  res[4, ] <- stat["max", ]
  res[5, ] <- stat["std.dev", ]
  
  return(res)
}

#################### 27. CREATE COHOROT CARACHTERISTICS
#####################################################
create_profile_char <- function(data, clustering, cluster_order, cluster_name_order,
                                variables, variables_names){
  #data is the dataframe
  #clustering is the partition
  #cluster_order is the reorder clustering
  #cluster_name_order is the names to be given to the clusters
  #variables are the variables to be used for the charcterization
  #variables_names are the names of the variables for display
  
  #create dataframe
  df <- data.frame(data, 
                   cluster = factor(reorder_clustering(clustering,
                                                       cluster_order,
                                                       cluster_name_order),
                                    levels = cluster_name_order))
  
  #calculate means and sd
  res_mean <- df %>%
    dplyr::select(variables, cluster) %>%
    group_by(cluster) %>%
    summarise_all(mean) %>%
    as.data.frame()
  
  res_sd <- df %>%
    dplyr::select(variables, cluster) %>%
    group_by(cluster) %>%
    summarise_all(sd) %>%
    as.data.frame()
  
  #create res matrix
  res <- as.data.frame(matrix(NA, nrow = nrow(res_mean), ncol = ncol(res_mean) - 1))
  rownames(res) <- res_mean[ ,1]
  colnames(res) <- variables_names
  
  for(j in 1:ncol(res)){
    res[ ,j] <- paste0(round(res_mean[ ,j+1],2), 
                       " (",
                       round(res_sd[ ,j+1], 2),
                       ")")
  }
  
  return(res)
}

#################### 28. GET SIZE BY CLUSTER
#####################################################

get_size_by_cluster <- function(clustering,
                                cluster_order,
                                cluster_name_order){
  #clustering is the partition
  #cluster_order is the reorder clustering
  #cluster_name_order is the names to be given to the clusters
  
  cluster_ordered <- reorder_clustering(clustering,
                                        cluster_order,
                                        cluster_name_order)
  
  
  res <- data.frame(cluster = factor(cluster_ordered, 
                                     levels = cluster_name_order)) %>% 
    group_by(cluster) %>% 
    dplyr::summarise(tot = n(), 
                     freq = round(n()/length(clustering) * 100, 1)) %>% 
    as.data.frame()
  
  return(res)
}

#################### 29. CREATE ALLUVIAL PLOT
#####################################################

plot_alluvial <- function(left_clustering, right_clustering,
                          left_cluster_order, right_cluster_order,
                          left_model_name, right_model_name,
                          left_clusters_names_order, right_clusters_names_order,
                          fill_clusters, group_colors = NULL, 
                          add_title,
                          plot_name = NULL){
  #left_clustering is the clustering to use on the left of the plot
  #right_clustering is the clustering to use on the right of the plot
  #left_cluster_order is the order of the left clusters
  #right_cluster_order is the order of the right clusters
  #left_model_name is the name of the left model
  #right_model_name is the name of the right model
  #left_clusters_names_order are the names of the clusters of the left model ordered
  #right_clusters_names_order are the names of the clusters of the right model ordered
  #fill_clusters is the clusters to be used for fill of stratum and alluvium
  #group_colors is the colors for each group
  #add_title is whether yes or no to add a title
  #plot_name is the title of the plot
  
  #clustering ordered
  left_clustering_ordered <- c(NA, length(left_clustering))
  
  for(k in 1:max(left_clustering)){
    left_clustering_ordered[which(left_clustering == k)] <- which(
      left_cluster_order == k)
  }
  
  right_clustering_ordered <- c(NA, length(right_clustering))
  
  for(k in 1:max(right_clustering)){
    right_clustering_ordered[which(right_clustering == k)] <- which(
      right_cluster_order == k)
  }
  
  
  #create dataframes of clustering
  df_alluvial <- data.frame(
    left_model = as.factor(left_clustering_ordered), 
    right_mode = as.factor(right_clustering_ordered))
  
  #create dataframe of all combinations
  df_alluvial_sum <- as.data.frame(expand.grid(
    left_model = as.factor(c(1:max(left_clustering_ordered))),
    right_model = as.factor(c(1:max(right_clustering_ordered)))))
  
  #create dataframes of cluster names
  clusters_names <- expand.grid(
    left_names = left_clusters_names_order,
    right_names = right_clusters_names_order)
  
  #count number of occurence of combinations
  df_alluvial_sum$tot <- rep(NA, nrow(df_alluvial_sum))
  
  for (i in 1:nrow(df_alluvial_sum)){
    df_alluvial_sum[i, 3] <- length(which((df_alluvial[ ,1] == df_alluvial_sum[i, 1]) & 
                                            (df_alluvial[ ,2] == df_alluvial_sum[i, 2])))
  }
  
  #add cluster names
  df_alluvial_sum$left_names <- clusters_names[ ,1]
  df_alluvial_sum$right_names <- clusters_names[ ,2]
  
  #add fill names
  if(fill_clusters == "left"){
    df_alluvial_sum$fill_names <- df_alluvial_sum$left_names
  } else{
    df_alluvial_sum$fill_names <- df_alluvial_sum$right_names
  }
  
  
  #plot
  res <- ggplot(df_alluvial_sum,
                aes(y = tot, axis1 = left_names, axis2 = right_names)) +
    geom_alluvium(aes(fill = fill_names), curve_type = "sigmoid", alpha = 0.7) +
    geom_stratum(aes(fill = fill_names), color = "white") +
    geom_text(aes(label = paste0(..stratum.., "\n", 
                                 scales::percent(..count.., accuracy = 1,
                                                 scale = 100 / sum(df_alluvial_sum$tot)))), 
              stat = "stratum", size = 3.5) +
    scale_x_discrete(limits = c(left_model_name, right_model_name), expand = c(.15, .15)) +
    scale_y_continuous(label = scales::percent_format(scale = 100 / sum(df_alluvial_sum$tot)),
                       breaks = c(0, sum(df_alluvial_sum$tot/4), 
                                  sum(df_alluvial_sum$tot)/2, sum(df_alluvial_sum$tot)/4*3,
                                  sum(df_alluvial_sum$tot))) +
    theme(legend.position="bottom", 
          plot.title = element_text(face = "bold")) +
    guides(fill = guide_legend(title=paste0("Clusters"))) +
    labs(y = "Frequency") 
  
  if(!is.null(group_colors)){
    res <-  res + 
      scale_fill_manual(values = group_colors)}
  
  if(add_title == TRUE){
    res <- res + ggtitle(plot_name)
  }
  
  #return plot
  return(res)
  
}

#################### 30. GET EVOLUTION OF CLUSTER TABLE
#####################################################

get_evol_cluster <- function(pred_init, pred_fu,
                             cluster_order, cluster_name_order, color_group){
  
  #pred_init are the prediction at initial assessment
  #pred_fu are the prediction at fu assessment
  #cluster_order is the order to be used for the clusters
  #cluster_name_order is the order of the clusters ordered
  #color_group is the color scheme to be used in the plot
  
  #create alluvial plot
  trans_plot <- 
    plot_alluvial(pred_init, 
                  pred_fu,
                  cluster_order, cluster_order,
                  "Baseline", "Repeated",
                  cluster_name_order, cluster_name_order,
                  "left", color_group,
                  add_title = FALSE)
  
  #create tale
  trans_count <- data.frame(
    cluster_init = factor(reorder_clustering(
      pred_init, cluster_order, cluster_name_order),
      levels = cluster_name_order),
    cluster_fu = factor(reorder_clustering(
      pred_fu, cluster_order, cluster_name_order),
      levels = cluster_name_order)) %>% 
    group_by(cluster_init, cluster_fu) %>% 
    dplyr::summarise(tot = n()) %>% 
    pivot_wider(names_from = cluster_fu, values_from = tot) %>% 
    as.data.frame()
  
  rownames(trans_count) <- cluster_name_order
  trans_count <- trans_count[, -1]
  
  trans_percent <- round(trans_count / rowSums(trans_count, na.rm = TRUE) *100, 1)
  
  #create res
  res <- list(
    trans_plot = trans_plot,
    trans_count = trans_count,
    trans_percent = trans_percent)
  
  return(res)
  
}


###############################################################################
############################# ASSOCIATION STUDY ###############################
###############################################################################

#################### 31. CALCULATE ODD RATIO
#####################################################
calculate_odd_ratio <- function(exposure, outcome){
  #exposure is the exposure or group
  #outcome is the response
  
  tab <- table(exposure, outcome)
  
  #create exposure / outcome table
  tab2 <- matrix(NA, 2, 2)
  tab2[1, 1] <- tab[which(rownames(tab) == TRUE), which(colnames(tab) == 1)]
  tab2[2, 1] <- tab[which(rownames(tab) == TRUE), which(colnames(tab) == 0)]
  tab2[1, 2] <- tab[which(rownames(tab) == FALSE), which(colnames(tab) == 1)]
  tab2[2, 2] <- tab[which(rownames(tab) == FALSE), which(colnames(tab) == 0)]
  rownames(tab2) <- c("case", "control")
  colnames(tab2) <- c("exposed", "not_exposed")
  
  #calculate odd_ratio
  odd_ratio <- (tab2[1,1] / tab2[1, 2]) / (tab2[2,1] / tab2[2, 2])
  
  #calculate log odd
  log_odd_ratio <- log(odd_ratio)
  
  #calculate ci
  log_sd_ci <- qnorm(0.975)*sqrt(sum(1/tab2))
  
  
  #print results
  res <- list(tab2, odd_ratio, log_odd_ratio, log_sd_ci)
  
  return(res)
}

#################### 32. CALCULATE ASSOCIATIONS
#####################################################
calculate_associations <- function(data_main,
                                   data_sel,
                                   data_hes,
                                   icd10_cardiovascular,
                                   clustering,
                                   cluster_order,
                                   cluster_names_order,
                                   var_selected,
                                   exclude_before,
                                   include_freq){
  #data_main is the main dataframe
  #data_sel is the selected dataframe with selected eid patients
  #data_hes is the hes data
  #icd10 is the list of cardiovasuclar disease
  #clustering is the clustering
  #cluster_order is the clustering ordered
  #cluster_names_order are the names of the clustered ordered
  #var_selected is the selected variables for logistic regression
  #exclude_before is whether yes or no to exclude patients who had a previous diagnosis of CVD
  #include_freq is whether yes or no to include freq
  
  
  ##################### 1. PREPARE TABLE
  #####################################################
  
  ###### a. Year of diagnosis
  
  #year of assessment
  df_dates <- data_main %>% 
    filter(eid %in% unlist(data_sel$eid)) %>% 
    dplyr::select(eid, `53-0.0`)
  
  df_dates$yoa <- year(ymd(df_dates$`53-0.0`))
  
  #year of diagnosis
  df_dates$yod <- df_dates$yoa - (data_sel$age_at_assessment - data_sel$diabetes_age_diag)
  
  ###### b. HES data
  
  #eid with cvd any date
  case_cvd <- data_hes %>%
    filter(eid %in% unlist(data_sel$eid)) %>% 
    filter(diag_icd10 %in% icd10_cardiovascular) %>% 
    dplyr::select(eid, epistart) %>% 
    arrange(eid, epistart)
  
  case_cvd$year_epistart <- year(dmy(case_cvd$epistart))
  
  #eid with cvd before diagnosis
  eid_cvd_before <- inner_join(df_dates, case_cvd, by="eid") %>% 
    filter(year_epistart < yod) %>% 
    dplyr::select(eid) %>% 
    unlist() %>% 
    unique() 
  
  #eid with cvd after diagnosis
  if(exclude_before == TRUE){
    eid_cvd <- inner_join(df_dates, case_cvd, by="eid") %>% 
      filter(!(eid %in% eid_cvd_before)) %>% 
      dplyr::select(eid) %>% 
      unlist() %>% 
      unique()
  } else {
    eid_cvd <- inner_join(df_dates, case_cvd, by="eid") %>% 
      filter(year_epistart >= yod) %>% 
      dplyr::select(eid) %>% 
      unlist() %>% 
      unique()
  }
  
  
  #numb of cvd in cvd after diagnosis
  if(exclude_before == TRUE){
    cvd_freq_case <- inner_join(df_dates, case_cvd, by="eid") %>% 
      filter(!(eid %in% eid_cvd_before)) %>% 
      group_by(eid) %>% 
      dplyr::summarise(tot=n()) %>% 
      as.data.frame()
  } else{
    cvd_freq_case <- inner_join(df_dates, case_cvd, by="eid") %>% 
      filter(year_epistart >= yod) %>%
      group_by(eid) %>% 
      dplyr::summarise(tot=n()) %>% 
      as.data.frame()
  }
  
  #flag patients with cvd
  flag_cvd <- rep(0, nrow(data_sel))
  
  for(i in 1:nrow(data_sel)){
    if(data_sel$eid[i] %in% eid_cvd){
      flag_cvd[i] <- 1
    }
    if(exclude_before == TRUE){
      if(data_sel$eid[i] %in% eid_cvd_before){
        flag_cvd[i] <- -1
      }}
  }
  
  #vector of frequency of cvd
  cvd_freq <- rep(0, nrow(data_sel))
  
  for(i in 1:nrow(data_sel)){
    if(data_sel$eid[i] %in% unlist(cvd_freq_case$eid)){
      cvd_freq[i] <- cvd_freq_case$tot[which(cvd_freq_case$eid == data_sel$eid[i])]
    }
  }
  
  
  ###### b. Cluster data
  cluster_ordered <- factor(reorder_clustering(clustering,
                                               cluster_order,
                                               cluster_names_order),
                            levels = cluster_names_order)
  
  ###### c. Data frame
  df_cvd <- data.frame(cluster = cluster_ordered,
                       cvd = flag_cvd,
                       freq = cvd_freq)
  df_cvd <- df_cvd %>% 
    filter(cvd != -1)
  
  ##################### 2. CONTINGENCY TABLE
  #####################################################
  
  freq_cvd <- table(df_cvd$cluster, df_cvd$cvd)/
    rowSums(table(df_cvd$cluster, df_cvd$cvd)) * 100
  
  avg_case_cvd <- df_cvd  %>% 
    filter(cvd == 1) %>% 
    group_by(cluster) %>% 
    dplyr::summarise(mean(freq)) %>% 
    as.data.frame()
  
  if(include_freq == TRUE){
    cont <- data.frame(tot = round(rowSums(table(df_cvd$cluster, df_cvd$cvd)), 0),
                       freq_cvd = round(freq_cvd[ ,2], 1),
                       avg_case_cvd = round(avg_case_cvd[ ,2], 1))
  } else{
    cont <- data.frame(tot = round(rowSums(table(df_cvd$cluster, df_cvd$cvd)), 0),
                       freq_cvd = round(freq_cvd[ ,2], 1))
  }
  
  
  
  
  ##################### 3. LOG ODD RATIO
  #####################################################
  log_odd <- c()
  for (i in 1:length(cluster_names_order)){
    tmp <- calculate_odd_ratio(df_cvd$cluster == cluster_names_order[i], 
                               df_cvd$cvd)[[3]]
    log_odd <- c(log_odd, tmp)
  }
  
  sd_log_odd <- c()
  for (i in 1:length(cluster_names_order)){
    tmp <- calculate_odd_ratio(df_cvd$cluster == cluster_names_order[i], 
                               df_cvd$cvd)[[4]]
    sd_log_odd <- c(sd_log_odd, tmp)
  }
  
  
  ##################### 4. LOGISTIC REGRESSION
  #####################################################
  glm_unadjusted <- list()
  glm_adjusted <- list()
  
  for (i in 1:length(cluster_names_order)){
    if(exclude_before == TRUE){
      data_glm <- data.frame(data_sel %>% 
                               filter(!(eid %in% eid_cvd_before)), 
                             cluster = (df_cvd$cluster == cluster_names_order[i]),
                             outcome = df_cvd$cvd)
    } else {
      data_glm <- data.frame(data_sel, 
                             cluster = (df_cvd$cluster == cluster_names_order[i]),
                             outcome = df_cvd$cvd)
    }
    
    glm_unadj <- glm(as.formula(paste("outcome", "~", "cluster")),
                     family = "binomial", data = data_glm)
    
    glm_adj <- glm(as.formula(paste("outcome", "~", "cluster", "+", 
                                    paste0(var_selected, collapse = "+"))),
                   family = "binomial", data = data_glm)
    
    glm_unadjusted[[i]] <- glm_unadj
    glm_adjusted[[i]] <- glm_adj
  }
  
  #create res matrix
  glm_res <- matrix(NA, nrow = length(cluster_names_order), ncol = 2)
  rownames(glm_res) <- cluster_names_order
  colnames(glm_res) <- c("pval_unadj", "pval_adj")
  
  for(i in 1:length(cluster_names_order)){
    glm_res[i, 1] <- coef(summary(glm_unadjusted[[i]]))[2,4]
    glm_res[i, 2] <- coef(summary(glm_adjusted[[i]]))[2,4]
  }
  
  ##################### 5. RECAP TABLES
  #####################################################
  cont$log_odd <- paste0(round(log_odd, 2),
                         " (",
                         round(sd_log_odd, 2),
                         ")")
  
  
  res <- list(cluster = df_cvd$cluster,
              eid_sel = eid_cvd,
              cvd = df_cvd$cvd,
              freq = df_cvd$freq,
              glm_unadjusted = glm_unadjusted,
              glm_adjusted = glm_adjusted,
              cont = cont,
              glm_res = glm_res)
  
  
  return(res)
  
}


###############################################################################
############################# MCMC DATA PREPARATION ###########################
###############################################################################

##################### 33. INTERPOLATE LINEARALY
#####################################################

interpolate_linear <- function(data_mat){
  #data_mat is the matrix to interpolate (each row)
  
  res <- data_mat
  
  #loop over ind
  for (ind in 1:nrow(res)){
    
    #find pos with first non NA value
    pos_start <- 1
    while(is.na(res[ind, pos_start])){
      pos_start <- pos_start + 1
    }
    
    #find pos with last non NA value
    pos_end <- pos_start + 1
    while(is.na(res[ind, pos_end])){
      pos_end <- pos_end + 1
    }
    
    #calculate increment
    val_start <- res[ind, pos_start]
    val_end <- res[ind, pos_end]
    incr <- (val_end-val_start)/(pos_end-pos_start)
    
    #fill values interpolated
    for (pos in (pos_start+1):(pos_end-1)){
      res[ind, pos] <- res[ind, pos-1] + incr
    }
    
  }
  return(res)
}

##################### 34. EXTRAPOLATION
#####################################################

extrapolate <- function(data_mat, first_year_extrapol, 
                        last_year_extrapol, first_year){
  #data_mat is the interpolated matrix to extrapolate (each row)
  #first_year_extrapol is the first year to start the extrapolation
  #last_year is the last_year to extrapolate
  #first_year is the year of start of study
  #first_pos is the pos of the first point to extrapolate
  #last_pos is the pos of the last point to extrapolate
  
  res <- data_mat
  
  
  #loop over ind
  for (ind in 1:nrow(res)){
    
    #find pos with first non NA value
    pos_start <- 1
    while(is.na(res[ind, pos_start])){
      pos_start <- pos_start + 1
    }
    
    #find pos with last non NA value
    pos_end <- pos_start + 1
    while(!is.na(res[ind, pos_end])){
      pos_end <- pos_end + 1
    }
    
    #get position of first_year_extrapol
    first_pos <- first_year_extrapol[ind] - first_year + 1
    
    #get position of last_year_extrapol
    if(is.na(last_year_extrapol[ind])){
      last_pos <- ncol(res)
    } else {
      last_pos <- first_pos + (last_year_extrapol[ind] - first_year_extrapol[ind]) 
    }
    
    
    #fill left values
    for(pos in first_pos:(pos_start-1)){
      res[ind, pos] <- res[ind, pos_start]
    }
    
    #fill right values
    if(last_pos < (pos_end-1)){
      for(pos in (last_pos+1):(pos_end-1)){
        res[ind, pos] <- NA
      }
    } else {
      for(pos in pos_end:(last_pos)){
        res[ind, pos] <- res[ind, pos_end-1]
      }
    }
    
    
    
  }
  return(res)
}

#################### 35. SAVE DATA FOR MCMC
#####################################################

save_for_mcmc <- function(data, out_dir, name_to_save, simul_name){
  #data is the dataframe to be save
  #out_dir is the directory in which to save the dataframe
  #name_to_save is the name under which to save the dataframe
  #simul name is the name of the simul to be added as suffix
  
  #save numb of rows and of columns
  write.table(matrix(c(nrow(data), ncol(data)), ncol=1), 
              paste0(out_dir, "predetails.txt"),
              quote=FALSE, row.names=FALSE, col.names=FALSE)
  
  #save content
  write.table(data,
              paste0(out_dir, "details.txt"),
              quote=FALSE, row.names=FALSE, col.names=FALSE)
  
  #combine dim of matrix and content
  system(paste0("cat ", out_dir, "predetails.txt ", 
                out_dir, "details.txt", " > ", 
                out_dir, paste0(name_to_save, simul_name, ".txt")))
  
  #remove transition matrix
  system(paste0("rm ", out_dir, "predetails.txt"))
  system(paste0("rm ", out_dir, "details.txt"))
}

#################### 36. SCALE MATRIX
#####################################################

scale_matrix <- function(data_matrix){
  #data_matrix is the matrix of data
  
  res <- matrix(scale(as.vector(data_matrix)),
                ncol = ncol(data_matrix))
  
  return(res)
}


###############################################################################
############################# MCMC RESULTS SUMMARY ############################
###############################################################################

#################### 37. GET MCMC MODEL SUMMARY
#####################################################

get_sum_mcmc <- function(simul_name, simul_seed, seed, n_iter, BI){
  #simul_name is the name of the simul
  #simul_seed is the seed of the simul
  #seed is the seed of the run
  #n_iter is the number of iterations for the model
  #BI is the number of BI
  
  #get file name
  file_name <- paste0("~/ukbb_diabetes/msm/MCMC_", simul_name,
                      "/Results_simulation_", simul_name,
                      "_", simul_seed ,"/", 
                      simul_name, "_", simul_seed, "_",  seed ,
                      "_History_iter_", n_iter ,".txt")
  
  #read content
  data <- fread(file_name, skip = 1)
  data <- as.data.frame(data)
  data <- data[ , -ncol(data)]
  
  #read header
  colnames(data) <- strsplit(readLines(file_name, n = 1), 
                             split ="\t")[[1]]
  
  #get index of first sigma
  indx_first_sig <- min(which(grepl("sigma", 
                                    colnames(data)) == TRUE))
  
  #create result matrix
  res <- matrix(NA, nrow = 4, ncol = indx_first_sig - 2)
  colnames(res) <- colnames(data)[2: (indx_first_sig - 1)]
  rownames(res) <- c("mean", "median", "2.5%", "97.5%")
  
  #fill matrix
  for (i in 1:ncol(res)){
    res[1, i] <- mean(data[BI:nrow(data), i + 1])
    res[2, i] <- quantile(data[BI:nrow(data), i + 1],
                          probs = 0.5)
    res[3, i] <- quantile(data[BI:nrow(data), i + 1],
                          probs = 0.025)
    res[4, i] <- quantile(data[BI:nrow(data), i + 1],
                          probs = 0.975)
  }
  
  return(res)
  
}

#################### 38. GET BIC MCMC
#####################################################

get_bic_mcmc <- function(sum_model, stat_to_use, n_obs){
  #sum_model is an output from function get_sum_mcmc
  #stat_to_use is the statistic to use (mean or median)
  #n_obs is the number of observation
  
  #get numb of params
  n_params = ncol(sum_model) - 1
  
  #select mean or median
  if(stat_to_use == "mean"){
    selected_row <- 1
  } else {
    selected_row <- 2
  }
  
  #calculate bic
  sum_model <- as.data.frame(sum_model)
  res <- n_params * log(n_obs) - 2 * sum_model[selected_row, ncol(sum_model)]
  
  return(res)
}

#################### 39. GET MCMC CHAINS
#####################################################

get_mcmc_chains <- function(simul_name, simul_seed, seed, n_iter, BI){
  #simul_name is the name of the simul
  #simul_seed is the seed of the simul
  #seed is the seed of the run
  #n_iter is the number of iterations for the model
  #BI is the number of BI
  
  #get file name
  file_name <- paste0("~/ukbb_diabetes/msm/MCMC_", simul_name,
                      "/Results_simulation_", simul_name,
                      "_", simul_seed ,"/", 
                      simul_name, "_", simul_seed, "_",  seed ,
                      "_History_iter_", n_iter ,".txt")
  
  #read content
  data <- fread(file_name, skip = 1)
  data <- as.data.frame(data)
  data <- data[ , -ncol(data)]
  
  #read header
  colnames(data) <- strsplit(readLines(file_name, n = 1), 
                             split ="\t")[[1]]
  
  #get index of first sigma
  indx_first_sig <- min(which(grepl("sigma", 
                                    colnames(data)) == TRUE))
  
  #res
  res <- data[(as.integer(BI)+1):nrow(data), 2: (indx_first_sig - 2)]
  colnames(res) <- colnames(data)[2: (indx_first_sig - 2)]
  
  
  return(res)
  
}

#################### 40. GET MCMC SUMMARY CLUSTER
#####################################################

get_mcmc_sum_cluster <- function(simul_name, simul_seed, seed, n_iter, BI,
                                   cluster_names_ordered){
  #simul_name is the name of the simul
  #simul_seed is the seed of the simul
  #seed is the seed of the run
  #n_iter is the number of iterations for the model
  #BI is the number of BI
  #cluster_names_ordered is the name of the clusteres ordered
  
  #get file name
  file_name <- paste0("~/ukbb_diabetes/msm/MCMC_", simul_name,
                      "/Results_simulation_", simul_name,
                      "_", simul_seed ,"/", 
                      simul_name, "_", simul_seed, "_",  seed ,
                      "_History_iter_", n_iter ,".txt")
  
  #read content
  data <- fread(file_name, skip = 1)
  data <- as.data.frame(data)
  data <- data[ , -ncol(data)]
  
  #read header
  colnames(data) <- strsplit(readLines(file_name, n = 1), 
                             split ="\t")[[1]]
  
  #get index of first sigma
  indx_first_sig <- min(which(grepl("sigma", 
                                    colnames(data)) == TRUE))
  
  #chains
  chains <- data[(as.integer(BI)+1):nrow(data), 2: (indx_first_sig - 2)]
  colnames(chains) <- colnames(data)[2: (indx_first_sig - 2)]
  
  
  #proba matrix
  proba <- chains
  
  #initialize
  proba[ ,1] <- exp(chains[, 1])/(1 + exp(chains[, 1]))
  
  #fill other cluster
  for(i in 2:ncol(proba)){
    proba[ ,i] <- exp(chains[, 1] + chains[, i])/(1 + exp(chains[, 1] + chains[, i]))
    
  }
  
  
  #create result matrix
  res <- matrix(NA, nrow = 4, ncol = ncol(proba))
  colnames(res) <- cluster_names_ordered
  rownames(res) <- c("mean", "median", "2.5%", "97.5%")
  
  
  for (i in 1:ncol(res)){
    res[1, i] <- mean(proba[ ,i])
    res[2, i] <- quantile(proba[ ,i],
                          probs = 0.5)
    res[3, i] <- quantile(proba[ ,i],
                          probs = 0.025)
    res[4, i] <- quantile(proba[ ,i],
                          probs = 0.975)
  }
  
  return(res)
  
}

#################### 41. GET MCMC SIMULATON SUMMARY
#####################################################

get_sum_mcmc_simul <- function(simul_name, simul_seed, seed, n_simul,
                               col_proba, x_title, groups_names){
  #simul_name is the name of the simul
  #simul_seed is the seed of the simul
  #seed is the seed of the run
  #n_simul is the number of simulation
  #col_proba is the column of the proba to plot in the dataset
  #x_title is the title for the plot
  #groups in the names of the groups to use in the legend
  
  #get file name
  file_name <- paste0("~/ukbb_diabetes/msm/MCMC_", simul_name,
                      "/Results_simulation_", simul_name,
                      "_", simul_seed ,"/", 
                      simul_name, "_", simul_seed, "_",  seed ,
                      "_Summary_simul_", n_simul ,".txt")
  
  #read content
  data <- fread(file_name, skip = 1)
  data <- as.data.frame(data)
  data <- data[ , -ncol(data)]
  
  #read header
  colnames(data) <- strsplit(readLines(file_name, n = 1), 
                             split ="\t")[[1]]
  
  #create df
  df <- data.frame(status = as.factor(data[ ,2]), 
                   proba = data[ ,col_proba])
  
  #create plot
  plot_sum <- 
    ggplot(df, aes(x=proba, color=status)) +
    geom_density() +
    scale_color_discrete(name = "State", 
                         labels = groups_names) +
    xlab(x_title) +
    ylab("Density") +
    theme(legend.position="bottom")
  
  #create mean and median tables
  res_sum <- matrix(NA, ncol = 2, nrow = length(groups_names))
  rownames(res_sum) <- groups_names
  colnames(res_sum) <- c("Mean", "Median")
  
  for(i in 1:nrow(res_sum)){
    res_sum[i, 1] <- mean(df$proba[which(df$status == i-1)])
    res_sum[i, 2] <- median(df$proba[which(df$status == i-1)])
  }
  
  #return res
  res <- list(plot_sum = plot_sum,
              res_sum = res_sum)
  return(res)
  
}

#################### 42. GET MCMC SIMULATON CLUSTER
#####################################################

get_sum_mcmc_cluster <- function(simul_name, simul_seed, seed, n_simul,
                                 clusters_names){
  #simul_name is the name of the simul
  #simul_seed is the seed of the simul
  #seed is the seed of the run
  #n_simul is the number of simulation
  #clusters_names is the names of the clusters
  
  #get file name - cluster
  file_name <- paste0("~/ukbb_diabetes/msm/MCMC_", simul_name,
                      "/Results_simulation_", simul_name,
                      "_", simul_seed ,"/", 
                      simul_name, "_", simul_seed, "_",  seed ,
                      "_Cluster_D_C_", n_simul ,".txt")
  
  #read content
  data_cluster <- fread(file_name)
  data_cluster <- as.data.frame(data_cluster)

  #get file name - summary
  file_name <- paste0("~/ukbb_diabetes/msm/MCMC_", simul_name,
                      "/Results_simulation_", simul_name,
                      "_", simul_seed ,"/", 
                      simul_name, "_", simul_seed, "_",  seed ,
                      "_Summary_simul_", n_simul ,".txt")
  
  #read content
  data_sum <- fread(file_name, skip = 1)
  data_sum <- as.data.frame(data_sum)
  data_sum <- data_sum[ , -ncol(data_sum)]
  
  #read header
  colnames(data_sum) <- strsplit(readLines(file_name, n = 1), 
                             split ="\t")[[1]]
  
  #calculate ratio
  res <- rep(NA, length(clusters_names))
  names(res) <- clusters_names
  
  for (i in 1:length(res)){
    res[i] <- length(which(data_cluster == i)) / (sum(
      length(which(data_sum[, 5] == i)),
      length(which(data_sum[, 6] == i)))/2*
        ncol(data_cluster))*100
  }
  
  return(res)
}
