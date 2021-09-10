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
library(coda)


###############################################################################
################################# GET FUNCTIONS ###############################
###############################################################################
source("~/ukbb_diabetes/functions.R")

###############################################################################
############################ GET SUMMARY BY MODEL #############################
###############################################################################

sum_v1 <- get_sum_mcmc(simul_name = "V1", 
              simul_seed = 1, seed = 1, 
              n_iter = "100000", BI = "10000")
sum_v2 <- get_sum_mcmc(simul_name = "V2", 
                        simul_seed = 1, seed = 1, 
                        n_iter = "100000", BI = "10000")
sum_v3 <- get_sum_mcmc(simul_name = "V3", 
                        simul_seed = 1, seed = 1, 
                        n_iter = "100000", BI = "10000")
sum_v4 <- get_sum_mcmc(simul_name = "V4", 
                        simul_seed = 1, seed = 1, 
                        n_iter = "100000", BI = "10000")
sum_v5 <- get_sum_mcmc(simul_name = "V5", 
                       simul_seed = 1, seed = 1, 
                       n_iter = "100000", BI = "10000")
sum_v6 <- get_sum_mcmc(simul_name = "V6", 
                       simul_seed = 1, seed = 1, 
                       n_iter = "100000", BI = "10000")

###############################################################################
################################### GET BIC ###################################
###############################################################################

#based on means
get_bic_mcmc(sum_v1, "mean", 922)
get_bic_mcmc(sum_v2, "mean", 922)
get_bic_mcmc(sum_v3, "mean", 922)
get_bic_mcmc(sum_v4, "mean", 922)
get_bic_mcmc(sum_v5, "mean", 284)
get_bic_mcmc(sum_v6, "mean", 284)

#based on medians
get_bic_mcmc(sum_v1, "median", 922)
get_bic_mcmc(sum_v2, "median", 922)
get_bic_mcmc(sum_v3, "median", 922)
get_bic_mcmc(sum_v4, "median", 922)
get_bic_mcmc(sum_v5, "median", 284)
get_bic_mcmc(sum_v6, "median", 284)


###############################################################################
######################## CONVERGENCE DIAGNOSTIC V3 ############################
###############################################################################

chains_list_v3 <- list(mcmc(get_mcmc_chains(simul_name = "V3", 
                                       simul_seed = 1, seed = 1, 
                                       n_iter = "100000", BI = "10000")), 
                       mcmc(get_mcmc_chains(simul_name = "V3", 
                                                   simul_seed = 2, seed = 1, 
                                                   n_iter = "100000", BI = "10000")),
                       mcmc(get_mcmc_chains(simul_name = "V3", 
                                            simul_seed = 3, seed = 1, 
                                            n_iter = "100000", BI = "10000")),
                       mcmc(get_mcmc_chains(simul_name = "V3", 
                                            simul_seed = 4, seed = 1, 
                                            n_iter = "100000", BI = "10000")),
                       mcmc(get_mcmc_chains(simul_name = "V3", 
                                            simul_seed = 5, seed = 1, 
                                            n_iter = "100000", BI = "10000")))
gelman.diag(mcmc.list(chains_list_v3))

###############################################################################
######################## CONVERGENCE DIAGNOSTIC V5 ############################
###############################################################################

chains_list_v5 <- list(mcmc(get_mcmc_chains(simul_name = "V5", 
                                            simul_seed = 1, seed = 1, 
                                            n_iter = "100000", BI = "10000")), 
                       mcmc(get_mcmc_chains(simul_name = "V5", 
                                            simul_seed = 2, seed = 1, 
                                            n_iter = "100000", BI = "10000")),
                       mcmc(get_mcmc_chains(simul_name = "V5", 
                                            simul_seed = 3, seed = 1, 
                                            n_iter = "100000", BI = "10000")),
                       mcmc(get_mcmc_chains(simul_name = "V5", 
                                            simul_seed = 4, seed = 1, 
                                            n_iter = "100000", BI = "10000")),
                       mcmc(get_mcmc_chains(simul_name = "V5", 
                                            simul_seed = 5, seed = 1, 
                                            n_iter = "100000", BI = "10000")))
gelman.diag(mcmc.list(chains_list_v5))



###############################################################################
######################## CONVERGENCE DIAGNOSTIC V6 ############################
###############################################################################

chains_list_v6 <- list(mcmc(get_mcmc_chains(simul_name = "V6", 
                                            simul_seed = 1, seed = 1, 
                                            n_iter = "100000", BI = "10000")), 
                       mcmc(get_mcmc_chains(simul_name = "V6", 
                                            simul_seed = 2, seed = 1, 
                                            n_iter = "100000", BI = "10000")),
                       mcmc(get_mcmc_chains(simul_name = "V6", 
                                            simul_seed = 3, seed = 1, 
                                            n_iter = "100000", BI = "10000")),
                       mcmc(get_mcmc_chains(simul_name = "V6", 
                                            simul_seed = 4, seed = 1, 
                                            n_iter = "100000", BI = "10000")),
                       mcmc(get_mcmc_chains(simul_name = "V6", 
                                            simul_seed = 5, seed = 1, 
                                            n_iter = "100000", BI = "10000")))
gelman.diag(mcmc.list(chains_list_v6))


###############################################################################
######################## GET SUMMARY OF PARAMETERS PD-D #######################
###############################################################################
#create summary table
summary_mcmc_simple <- matrix(NA, nrow = 16 + 1, ncol = 3)
colnames(summary_mcmc_simple) <- c("Model 1", "Model 2", "Model 3")
rownames(summary_mcmc_simple) <- colnames(sum_v2)

#fill info
summary_mcmc_simple[, 1] <- paste0(round(sum_v2[2, ],2 ), 
                                   " [", round(sum_v2[3, ], 2), ", ",
                                   round(sum_v2[4, ], 2), "]")

summary_mcmc_simple[, 3] <- paste0(round(sum_v4[2, ],2 ), 
                                   " [", round(sum_v4[3, ], 2), ", ",
                                   round(sum_v4[4, ], 2), "]")

summary_mcmc_simple[8:17, 2] <- paste0(round(sum_v3[2, ],2 ), 
                                       " [", round(sum_v3[3, ], 2), ", ",
                                       round(sum_v3[4, ], 2), "]")
summary_mcmc_simple[1, 2] <- summary_mcmc_simple[8, 2]
summary_mcmc_simple[8, 2] <- NA

kbl(summary_mcmc_simple) %>% 
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), font_size = 12) %>%
  scroll_box(height = "300px")

###############################################################################
######################## GET SUMMARY OF PARAMETERS D-C ########################
###############################################################################

sum_v5 <- get_mcmc_sum_cluster(simul_name = "V5", 
                     simul_seed = 3, seed = 1, 
                     n_iter = "100000", BI = "10000",
                     cluster_names_ordered = c("SIDD", "SIRD", "MOD", "MARD"))

sum_v6 <- get_mcmc_sum_cluster(simul_name = "V6", 
                               simul_seed = 3, seed = 1, 
                               n_iter = "100000", BI = "10000",
                               cluster_names_ordered = c("SIDD", "SIRD", "MOD", "MARD1", "MARD2"))
get_sum_mcmc()
#fill info
sum_v5 <- rbind(sum_v5,
paste0(round(sum_v5[2, ],3 ), 
       " [", round(sum_v5[3, ], 3), ", ",
       round(sum_v5[4, ], 3), "]"))

sum_v6 <- rbind(sum_v6,
                paste0(round(sum_v6[2, ],3 ), 
                       " [", round(sum_v6[3, ], 3), ", ",
                       round(sum_v6[4, ], 3), "]"))

###############################################################################
######################## ABILITY TO RECOVER TRAJ PD-D #########################
###############################################################################

#simul pdd
sum_simul_pdd_v3 <- get_sum_mcmc_simul(simul_name = "V3", 
             simul_seed = 1, seed = 1, 
             n_simul = "10000",
             col_proba = 4,
             x_title  = "Simulated PD-D transition probability",
             groups_names = c("Pre-Diabetes", "T2D", "Healthy"))

pdf("~/ukbb_diabetes/msm/Figures/sum_simul_pdd_v3.pdf")
sum_simul_pdd_v3$plot_sum
dev.off()

#simul pdh
sum_simul_pdh_v3 <- get_sum_mcmc_simul(simul_name = "V3", 
                                       simul_seed = 1, seed = 1, 
                                       n_simul = "10000",
                                       col_proba = 6,
                                       x_title  = "Simulated PD-H transition probability",
                                       groups_names = c("Pre-Diabetes", "T2D", "Healthy"))

pdf("~/ukbb_diabetes/msm/Figures/sum_simul_pdh_v3.pdf")
sum_simul_pdh_v3$plot_sum
dev.off()

#recap table
xtable(data.frame(pdd = sum_simul_pdd_v3$res_sum[,1],
                  pdh = sum_simul_pdh_v3$res_sum[,2]))


###############################################################################
######################## ABILITY TO RECOVER TRAJ D-C ##########################
###############################################################################
sum_simul_dc_v5 <- get_sum_mcmc_cluster(simul_name = "V5", 
                           simul_seed = 1, seed = 1, 
                           n_simul = "10000",
                           clusters_names = c("SIDD", "SIRD", "MOD", "MARD"))

sum_simul_dc_v6 <- get_sum_mcmc_cluster(simul_name = "V6", 
                              simul_seed = 1, seed = 1, 
                              n_simul = "10000",
                              clusters_names = c("SIDD", "SIRD", "MOD", "MARD1", "MARD2"))
sum_simul_dc_v6 <- as.data.frame(sum_simul_dc_v6)

xtable(sum_simul_dc_v6)
