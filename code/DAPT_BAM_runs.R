#Merritt Harlan
#BAM Sensitivity for TOWNS dataset
library(bamr)
library(ggplot2)
library(viridis)
library(tidyr)
library(ggsci)
library(scales)
library(ggpubr)

#load a river as a single bam data object
Sac_rch1_data <- readRDS(file = "Sac_rch1_data.rds")
Sac_rch2_data <- readRDS(file = "Sac_rch2_data.rds")
Sac_rch3_data <- readRDS(file = "Sac_rch3_data.rds")
Olen_data <- readRDS(file = "Olen_data.rds")
Will_data <- readRDS(file = "Will_data.rds")
Tan_data <- readRDS(file = "Tan_data.rds")
North_Sask_data <- readRDS(file = "North_Sask_data.rds")

#lapply generate bam priors and run bam onto the giant prior matrix
#Sac_rch1----
Sac_rch1_prior_df <-
  generate_prior_df(bamdata = Sac_rch1_data,
                    n_levels = 10,
                    true = 0, var_scale = 1)
Sac_rch1_prior_df = Sac_rch1_prior_df[Sac_rch1_prior_df$variant == "manning",]
Sac_rch1_sensitivity =  lapply(as.list(1:nrow(Sac_rch1_prior_df)), function(x) {
  run_and_validate_bam(
    bamdata = Sac_rch1_data,
    new_df = Sac_rch1_prior_df[x, ],
    n_levels = 10,
    qval = exp(Sac_rch1_data$logQ_hat)
  )
})

#Sac_rch2----
Sac_rch2_prior_df <-
  generate_prior_df(bamdata = Sac_rch2_data,
                    n_levels = 10,
                    true = 0, var_scale = 1)
Sac_rch2_prior_df = Sac_rch2_prior_df[Sac_rch2_prior_df$variant == "manning",]
Sac_rch2_sensitivity =  lapply(as.list(1:nrow(Sac_rch1_prior_df)), function(x) {
  run_and_validate_bam(
    bamdata = Sac_rch2_data,
    new_df = Sac_rch2_prior_df[x, ],
    n_levels = 10,
    qval = exp(Sac_rch2_data$logQ_hat)
  )
  
})

#Sac_rch3----
Sac_rch3_prior_df <-
  generate_prior_df(bamdata = Sac_rch3_data,
                    n_levels = 10,
                    true = 0, var_scale = 1)
Sac_rch3_sensitivity =  lapply(as.list(1:nrow(Sac_rch3_prior_df)), function(x) {
  run_and_validate_bam(
    bamdata = Sac_rch3_data,
    new_df = Sac_rch3_prior_df[x, ],
    n_levels = 10,
    qval = exp(Sac_rch3_data$logQ_hat)
  )
  print(x)
})

#Olentangy----
Olen_prior_df <-
  generate_prior_df(bamdata = Olen_data,
                    n_levels = 10,
                    true = 0, var_scale = 1)
Olen_sensitivity =  lapply(as.list(1:nrow(Olen_prior_df)), function(x) {
  run_and_validate_bam(
    bamdata = Olen_data,
    new_df = Olen_prior_df[x, ],
    n_levels = 10,
    qval = exp(Olen_data$logQ_hat)
  )
})

#Willamette----
Will_prior_df <-
  generate_prior_df(bamdata = Will_data,
                    n_levels = 10,
                    true = 0, var_scale = 1)
Will_sensitivity =  lapply(as.list(1:nrow(Will_prior_df)), function(x) {
  run_and_validate_bam(
    bamdata = Will_data,
    new_df = Will_prior_df[x, ],
    n_levels = 10,
    qval = exp(Will_data$logQ_hat)
  )
  print(x)
})

#Tanana----
Tan_prior_df <-
  generate_prior_df(bamdata = Tan_data,
                    n_levels = 10,
                    true = 0, var_scale = 1)
Tan_sensitivity =  lapply(as.list(1:nrow(Tan_prior_df)), function(x) {
  run_and_validate_bam(
    bamdata = Tan_data,
    new_df = Tan_prior_df[x, ],
    n_levels = 10,
    qval = exp(Tan_data$logQ_hat)
  )
})

#North Saskatchewan----
North_Sask_prior_df <-
  generate_prior_df(bamdata = North_Sask_data,
                    n_levels = 10,
                    true = 0,var_scale = 1)
North_Sask_sensitivity =  lapply(as.list(1:nrow(North_Sask_prior_df)), function(x) {
  run_and_validate_bam(
    bamdata = North_Sask_data,
    new_df = North_Sask_prior_df[x, ],
    n_levels = 10,
    qval = exp(North_Sask_data$logQ_hat)
  )
})
North_Sask_time <- proc.time() - ptm

#Plot TOWNSS-----
Sac_rch1_plot <-
  plot_bam_priors(
    bamdata = Sac_rch1_data,
    bam_priors_df = Sac_rch1_prior_df,
    bam_sensitivity = Sac_rch1_sensitivity,
    rivname = "Sac_rch1",
    n_levels = 10,
    color = viridis(1)
  )
Sac_rch2_plot <-
  plot_bam_priors(
    bamdata = Sac_rch2_data,
    bam_priors_df = Sac_rch2_prior_df,
    bam_sensitivity = Sac_rch2_sensitivity,
    rivname = "Sac_rch2",
    n_levels = 10,
    color = viridis(1)
  )

Sac_rch3_plot <-
  plot_bam_priors(
    bamdata = Sac_rch3_data,
    bam_priors_df = Sac_rch3_prior_df,
    bam_sensitivity = Sac_rch3_sensitivity,
    rivname = "Sac_rch3",
    n_levels = 10,
    color = viridis(1)
  )

Olen_plot <-
  plot_bam_priors(
    bamdata = Olen_data,
    bam_priors_df = Olen_prior_df,
    bam_sensitivity = Olen_sensitivity,
    rivname = "Olen",
    n_levels = 10,
    color = viridis(1)
  )

Will_plot <-
  plot_bam_priors(
    bamdata = Will_data,
    bam_priors_df = Will_prior_df,
    bam_sensitivity = Will_sensitivity,
    rivname = "Will",
    n_levels = 10,
    color = viridis(1)
  )

Tan_plot <-
  plot_bam_priors(
    bamdata = Tan_data,
    bam_priors_df = Tan_prior_df,
    bam_sensitivity = Tan_sensitivity,
    rivname = "Tan",
    n_levels = 10,
    color = viridis(1)
  )

North_Sask_plot <-
  plot_bam_priors(
    bamdata = North_Sask_data,
    bam_priors_df = North_Sask_prior_df,
    bam_sensitivity = North_Sask_sensitivity,
    rivname = "North_Sask",
    n_levels = 10,
    color = viridis(1)
  )

#Test # of PTs-----
Sac_rch1_num_pt = test_num_pt(bamdata = Sac_rch1_data)
Sac_rch2_num_pt = test_num_pt(bamdata = Sac_rch2_data)
Sac_rch3_num_pt = test_num_pt(bamdata = Sac_rch3_data)
Olen_num_pt = test_num_pt(bamdata = Olen_data)
Will_num_pt = test_num_pt(bamdata = Will_data)
Tan_num_pt = test_num_pt(bamdata = Tan_data)
North_Sask_num_pt = test_num_pt(bamdata = North_Sask_data)

#Turn into df
Sac_rch1_pt_test = number_df(num_pt_list = Sac_rch1_num_pt, rivname = "Sac_rch1")
Sac_rch2_pt_test = number_df(num_pt_list = Sac_rch2_num_pt, rivname = "Sac_rch2")
Sac_rch3_pt_test = number_df(num_pt_list = Sac_rch3_num_pt, rivname = "Sac_rch3")
Olen_pt_test = number_df(num_pt_list = Olen_num_pt, rivname = "Olen")
Tan_pt_test = number_df(num_pt_list = Tan_num_pt, rivname = "Tan")
Will_pt_test = number_df(num_pt_list = Will_num_pt, rivname = "Will")
North_Sask_pt_test = number_df(num_pt_list = North_Sask_num_pt, rivname = "North_Sask")

combined_pt_test = rbind(
  Sac_rch1_pt_test,
  Sac_rch2_pt_test,
  Sac_rch3_pt_test,
  Olen_pt_test,
  Tan_pt_test,
  Will_pt_test,
  North_Sask_pt_test
)

#Test Spacing of PTs------
Sac_rch1_spacing_pt = test_spacing_pt(bamdata = Sac_rch1_data)
Sac_rch2_spacing_pt = test_spacing_pt(bamdata = Sac_rch2_data)
Sac_rch3_spacing_pt = test_spacing_pt(bamdata = Sac_rch3_data)
Olen_spacing_pt = test_spacing_pt(bamdata = Olen_data)
Will_spacing_pt = test_spacing_pt(bamdata = Will_data)
Tan_spacing_pt = test_spacing_pt(bamdata = Tan_data)
North_Sask_spacing_pt = test_spacing_pt(bamdata = North_Sask_data)

Sac_rch1_xvec = readRDS("Sac_rch1_xvec.rds")
Sac_rch2_xvec = readRDS("Sac_rch2_xvec.rds")
Sac_rch3_xvec = readRDS("Sac_rch3_xvec.rds")
Olen_xvec = readRDS("Olen_xvec.rds")
Will_xvec = readRDS("Will_xvec.rds")
North_Sask_xvec = readRDS("North_Sask_xvec.rds")
Tan_xvec = readRDS("Tan_xvec.rds")

Sac_rch1_spacing = spacing_df(spacing_pt = Sac_rch1_spacing_pt,
                              rivname = "Sac_rch1",
                              riv_dist = Sac_rch1_xvec)
Sac_rch2_spacing = spacing_df(spacing_pt = Sac_rch2_spacing_pt,
                              rivname = "Sac_rch2",
                              riv_dist = Sac_rch2_xvec)
Sac_rch3_spacing = spacing_df(spacing_pt = Sac_rch3_spacing_pt,
                              rivname = "Sac_rch3",
                              riv_dist = Sac_rch3_xvec)
Olen_spacing = spacing_df(
  spacing_pt = Olen_spacing_pt,
  rivname = "Olen",
  riv_dist = rev(Olen_xvec)
)
Tan_spacing = spacing_df(spacing_pt = Tan_spacing_pt,
                         rivname = "Tan",
                         riv_dist = Tan_xvec)
Will_spacing = spacing_df(
  spacing_pt = Will_spacing_pt,
  rivname = "Will",
  riv_dist = rev(Will_xvec)
)
North_Sask_spacing = spacing_df(spacing_pt = North_Sask_spacing_pt,
                                rivname = "North_Sask",
                                riv_dist = North_Sask_xvec)

combined_spacing = rbind(
  Sac_rch1_spacing,
  Sac_rch2_spacing,
  Sac_rch3_spacing,
  Olen_spacing,
  Tan_spacing,
  Will_spacing,
  North_Sask_spacing
)

#Test Timing--------
Sac_rch1_data_6hr <- readRDS(file = "Sac_rch1_data_6hr.rds")
Sac_rch2_data_6hr <- readRDS(file = "Sac_rch2_data_6hr.rds")
Sac_rch3_data_6hr <- readRDS(file = "Sac_rch3_data_6hr.rds")
Olen_data_6hr <- readRDS(file = "Olen_data_6hr.rds")
Will_data_6hr <- readRDS(file = "Will_data_6hr.rds")
Tan_data_6hr <- readRDS(file = "Tan_data_6hr.rds")
North_Sask_data_6hr <- readRDS(file = "North_Sask_data_6hr.rds")

Sac_rch1_data_hourly <- readRDS(file = "Sac_rch1_data_hourly.rds")
Sac_rch2_data_hourly <- readRDS(file = "Sac_rch2_data_hourly.rds")
Sac_rch3_data_hourly <- readRDS(file = "Sac_rch3_data_hourly.rds")
Olen_data_hourly <- readRDS(file = "Olen_data_hourly.rds")
Will_data_hourly <- readRDS(file = "Will_data_hourly.rds")
Tan_data_hourly <- readRDS(file = "Tan_data_hourly.rds")
North_Sask_data_hourly <-readRDS(file = "North_Sask_data_hourly.rds")

Sac_rch1_timing <-
  test_timing(bamdata_list = list(Sac_rch1_data, Sac_rch1_data_6hr, Sac_rch1_data_hourly))

Sac_rch2_timing <-
  test_timing(bamdata_list = list(Sac_rch2_data, Sac_rch2_data_6hr, Sac_rch2_data_hourly))

Sac_rch3_timing <-
  test_timing(bamdata_list = list(Sac_rch3_data, Sac_rch3_data_6hr, Sac_rch3_data_hourly))

Olen_timing <-
  test_timing(bamdata_list = list(Olen_data, Olen_data_6hr, Olen_data_hourly))

Tan_timing <-
  test_timing(bamdata_list = list(Tan_data, Tan_data_6hr, Tan_data_hourly))

Will_timing <-
  test_timing(bamdata_list = list(Will_data, Will_data_6hr, Will_data_hourly))

Sac_rch1_timing_df = timing_df(Sac_rch1_timing, "Sac_rch1")
Sac_rch2_timing_df = timing_df(Sac_rch2_timing, "Sac_rch2")
Sac_rch3_timing_df = timing_df(Sac_rch3_timing, "Sac_rch3")
Olen_timing_df = timing_df(Olen_timing, "Olen")
Tan_timing_df = timing_df(Tan_timing, "Tan")
Will_timing_df = timing_df(Will_timing, "Will")

combined_timing_df = rbind.data.frame(Sac_rch1_timing_df, Sac_rch2_timing_df, Sac_rch3_timing_df,
                                           Olen_timing_df, Will_timing_df, Tan_timing_df)
#Test Q for TOWNSSS-------
Sac_rch1_Q_test = test_Q_num(bamdata = Sac_rch1_data)
Sac_rch2_Q_test = test_Q_num(bamdata = Sac_rch2_data)
Sac_rch3_Q_test = test_Q_num(bamdata = Sac_rch3_data)
Olen_Q_test = test_Q_num(bamdata = Olen_data)
Will_Q_test = test_Q_num(bamdata = Will_data)
Tan_Q_test = test_Q_num(bamdata = Tan_data)
North_Sask_Q_test = test_Q_num(bamdata = North_Sask_data)

Sac_rch1_Q_num = num_Q_df(num_Q_list = Sac_rch1_Q_test, rivname = "Sac_rch1")
Sac_rch2_Q_num = num_Q_df(num_Q_list = Sac_rch2_Q_test, rivname = "Sac_rch2")
Sac_rch3_Q_num = num_Q_df(num_Q_list = Sac_rch3_Q_test, rivname = "Sac_rch3")
Olen_Q_num = num_Q_df(num_Q_list = Olen_Q_test, rivname = "Olen")
Will_Q_num = num_Q_df(num_Q_list = Will_Q_test, rivname = "Will")
Tan_Q_num = num_Q_df(num_Q_list = Tan_Q_test, rivname = "Tan")
North_Sask_Q_num = num_Q_df(num_Q_list = North_Sask_Q_test, rivname = "North_Sask")

combined_Q_num= rbind(Sac_rch1_Q_num, Sac_rch2_Q_num, Sac_rch3_Q_num, Olen_Q_num, Will_Q_num, Tan_Q_num, North_Sask_Q_num)

Sac_rch1_TG = TG_all(num_days = ncol(Sac_rch1_data$Wobs),
                     pt_count = nrow(Sac_rch1_data$Wobs),
                     WSE_mat = as.matrix(Sac_WSE_rch1_daily_df),
                     Qobs = exp(Sac_rch1_data$logQ_hat),
                     Q_id = Q_ids(Sac_rch1_Q_test))

Sac_rch2_TG = TG_all(num_days = ncol(Sac_rch2_data$Wobs),
                     pt_count = nrow(Sac_rch2_data$Wobs),
                     WSE_mat = as.matrix(Sac_WSE_rch2_daily_df),
                     Qobs = exp(Sac_rch2_data$logQ_hat),
                     Q_id = Q_ids(Sac_rch2_Q_test))

Sac_rch3_TG = TG_all(num_days = ncol(Sac_rch3_data$Wobs),
                     pt_count = nrow(Sac_rch3_data$Wobs),
                     WSE_mat = as.matrix(Sac_WSE_rch3_daily_df),
                     Qobs = exp(Sac_rch3_data$logQ_hat),
                     Q_id = Q_ids(Sac_rch3_Q_test))

Olen_TG = TG_all(num_days = ncol(Olen_data$Wobs),
                 pt_count = nrow(Olen_data$Wobs),
                 WSE_mat = as.matrix(Olen_WSE_daily_df)[,c(seq(1,38,3))],
                 Qobs = exp(Olen_data$logQ_hat),
                 Q_id = Q_ids(Olen_Q_test))

Will_TG = TG_all(num_days = ncol(Will_data$Wobs),
                 pt_count = nrow(Will_data$Wobs),
                 WSE_mat = as.matrix(Will_WSE_daily_df),
                 Qobs = exp(Will_data$logQ_hat),
                 Q_id = Q_ids(Will_Q_test))

Tan_TG = TG_all(num_days = ncol(Tan_data$Wobs),
                pt_count = nrow(Tan_data$Wobs),
                WSE_mat = as.matrix(Tan_WSE_daily_df),
                Qobs = exp(Tan_data$logQ_hat),
                Q_id = Q_ids(Tan_Q_test))

North_Sask_TG = TG_all(num_days = ncol(North_Sask_data$Wobs),
                       pt_count = nrow(North_Sask_data$Wobs),
                       WSE_mat = as.matrix(Sask_WSE_daily_df),
                       Qobs = exp(North_Sask_data$logQ_hat),
                       Q_id = Q_ids(North_Sask_Q_test))

Sac_rch1_TG_nse = TG_NSE(tg_all = Sac_rch1_TG, rivname = "Sacramento R1")
Sac_rch2_TG_nse = TG_NSE(tg_all = Sac_rch2_TG, rivname = "Sacramento R2")
Sac_rch3_TG_nse = TG_NSE(tg_all = Sac_rch3_TG, rivname = "Sacramento R3")
Olen_TG_nse = TG_NSE(tg_all = Olen_TG, rivname = "Olentangy")
Will_TG_nse = TG_NSE(tg_all = Will_TG, rivname = "Willamette")
Tan_TG_nse = TG_NSE(tg_all = Tan_TG, rivname = "Tanana")
North_Sask_TG_nse = TG_NSE(tg_all = North_Sask_TG, rivname = "N. Saskatchewan")

Sac_rch1_comb = rbind(Sac_rch1_Q_num[c(1:3,8)], Sac_rch1_TG_nse)
Sac_rch2_comb = rbind(Sac_rch2_Q_num[c(1:3,8)], Sac_rch2_TG_nse)
Sac_rch3_comb = rbind(Sac_rch3_Q_num[c(1:3,8)], Sac_rch3_TG_nse)
Olen_comb = rbind(Olen_Q_num[c(1:3,8)], Olen_TG_nse)
Tan_comb = rbind(Tan_Q_num[c(1:3,8)], Tan_TG_nse)
Will_comb = rbind(Will_Q_num[c(1:3,8)], Will_TG_nse)
North_Sask_comb = rbind(North_Sask_Q_num[c(1:3,8)], North_Sask_TG_nse)

combined_TOWNS = rbind(
  Sac_rch1_comb,
  Sac_rch2_comb,
  Sac_rch3_comb, 
  Olen_comb,
  Tan_comb,
  Will_comb,
  North_Sask_comb
)
Sac_rch1_stats = aggregate(nse ~ variant + Q, Sac_rch1_comb, quantile)
Sac_rch2_stats = aggregate(nse ~ variant + Q, Sac_rch2_comb, quantile)
Sac_rch3_stats = aggregate(nse ~ variant + Q, Sac_rch3_comb, quantile)
Olen_stats = aggregate(nse ~ variant + Q, Olen_comb, quantile)
Tan_stats = aggregate(nse ~ variant + Q, Tan_comb, quantile)
Will_stats = aggregate(nse ~ variant + Q, Will_comb, quantile)
North_Sask_stats = aggregate(nse ~ variant + Q, North_Sask_comb, quantile)

Sac_rch1_hydro = plot_tradeoff(Sac_rch1_stats, "Sacramento R1")

#Plot Hydrographs-----
Sac_rch1_num_hydro = number_hydro(Sac_rch1_num_pt, "Sac_rch1")
Sac_rch2_num_hydro = number_hydro(Sac_rch2_num_pt, "Sac_rch2")
Sac_rch3_num_hydro = number_hydro(Sac_rch3_num_pt, "Sac_rch3")
Olen_num_hydro = number_hydro(Olen_num_pt, "Olen")
Will_num_hydro = number_hydro(Will_num_pt, "Will")
Tan_num_hydro = number_hydro(Tan_num_pt, "Tan")
North_Sask_num_hydro = number_hydro(North_Sask_num_pt, "North_Sask")

Sac_rch1_spacing_hydro = spacing_hydro(Sac_rch1_spacing_pt, Sac_rch1_xvec, "Sac_rch1")
Sac_rch2_spacing_hydro = spacing_hydro(Sac_rch2_spacing_pt, Sac_rch2_xvec, "Sac_rch2")
Sac_rch3_spacing_hydro = spacing_hydro(Sac_rch3_spacing_pt, Sac_rch3_xvec, "Sac_rch3")
Olen_spacing_hydro = spacing_hydro(Olen_spacing_pt, Olen_xvec, "Olen")
Will_spacing_hydro = spacing_hydro(Will_spacing_pt, Will_xvec, "Will")
Tan_spacing_hydro = spacing_hydro(Tan_spacing_pt, Tan_xvec, "Tan")
North_Sask_spacing_hydro = spacing_hydro(North_Sask_spacing_pt, North_Sask_xvec, "North_Sask")

Sac_rch1_Q_hydro = Q_hydro(Sac_rch1_Q_test, "Sac_rch1")
Sac_rch2_Q_hydro = Q_hydro(Sac_rch2_Q_test, "Sac_rch2")
Sac_rch3_Q_hydro = Q_hydro(Sac_rch3_Q_test, "Sac_rch3")
Olen_Q_hydro = Q_hydro(Olen_Q_test, "Olen")
Will_Q_hydro = Q_hydro(Will_Q_test, "Will")
Tan_Q_hydro = Q_hydro(Tan_Q_test, "Tan")
North_Sask_Q_hydro = Q_hydro(North_Sask_Q_test, "North_Sask")

Sac_rch1_df = rbind(Sac_rch1_plot[[2]][,2:4], Sac_rch1_num_hydro[,3:5], Sac_rch1_spacing_hydro[,3:5], Sac_rch1_Q_hydro[,3:5])
Sac_rch2_df = rbind(Sac_rch2_plot[[2]][,2:4], Sac_rch2_num_hydro[,3:5], Sac_rch2_spacing_hydro[,3:5], Sac_rch2_Q_hydro[,3:5])
Sac_rch3_df = rbind(Sac_rch3_plot[[2]][,2:4], Sac_rch3_num_hydro[,3:5], Sac_rch3_spacing_hydro[,3:5])#, Sac_rch3_Q_hydro[,3:5])
Olen_df = rbind(Olen_plot[[2]][,2:4], Olen_num_hydro[,3:5], Olen_spacing_hydro[,3:5], Olen_Q_hydro[,3:5])
Will_df = rbind(Will_plot[[2]][,2:4], Will_num_hydro[,3:5], Will_spacing_hydro[,3:5], Will_Q_hydro[,3:5])
Tan_df = rbind(Tan_plot[[2]][,2:4], Tan_num_hydro[,3:5], Tan_spacing_hydro[,3:5], Tan_Q_hydro[,3:5])
North_Sask_df = rbind(North_Sask_plot[[2]][,2:4], North_Sask_num_hydro[,3:5], North_Sask_spacing_hydro[,3:5], North_Sask_Q_hydro[,3:5])

Tan_hydro = plot_hydro(Tan_df, Tan_data, "Tanana")+annotate("text", x = 24, y = 865, label = "NSE = 0.97", size = 3)
Olen_hydro = plot_hydro(Olen_df, Olen_data, "Olentangy")+annotate("text", x = 25, y = 11, label = "NSE = 0.88", size = 3)
Will_hydro = plot_hydro(Will_df, Will_data, "Willamette")+annotate("text", x = 35, y = 280, label = "NSE = 0.88", size = 3)
North_Sask_hydro = plot_hydro(North_Sask_df, North_Sask_data, "N Saskatchewan")+annotate("text", x = 30, y = 580, label = "NSE = 0.81", size = 3)
Sac_rch1_hydro = plot_hydro(Sac_rch1_df, Sac_rch1_data, "Sacramento R1")+annotate("text", x = 9, y = 153, label = "NSE = 0.99", size = 3)
Sac_rch2_hydro = plot_hydro(Sac_rch2_df, Sac_rch2_data, "Sacramento R2")+annotate("text", x = 20, y = 135, label = "NSE = 0.84", size = 3)
Sac_rch3_hydro = plot_hydro(Sac_rch3_df, Sac_rch3_data, "Sacramento R3")+annotate("text", x = 35, y = 170, label = "NSE = 0.81", size = 3)


png(
  "TOWNSS_Hydrograph.png",
  width = 6.5,
  height = 3.5,
  units = 'in',
  res = 300
)
TOWNS_hydrographs = 
  ggarrange(
    Tan_hydro,
    Olen_hydro,
    Will_hydro,
    Sac_rch1_hydro,
    Sac_rch2_hydro,
    Sac_rch3_hydro,
    ncol = 3,
    nrow = 2,
    common.legend = TRUE,
    legend = "right",
    labels = c("A", "B", "C", "D", "E", "F"))

annotate_figure(TOWNS_hydrographs,
                bottom = text_grob("Number of Days", size = 10),
                left = text_grob(expression(paste(Discharge~(m ^ 3 / s))), size = 10, rot = 90))
dev.off()

#Full distribution boxplots
bam_sensitivity_comb = rbind.data.frame(
  Sac_rch1_plot[[1]],
  Sac_rch2_plot[[1]],
  Sac_rch3_plot[[1]],
  Olen_plot[[1]],
  North_Sask_plot[[1]],
  Tan_plot[[1]],
  Will_plot[[1]]
)
bam_sensitivity_comb = bam_sensitivity_comb[order(bam_sensitivity_comb$prior_names),]
combined_timing_df$variant = "manning"
combined_Q_num$variant = "manning"
full_nash = data.frame(
  test = c(
    rep("PT_Spacing", nrow(combined_spacing)),
    rep("PT_Number", nrow(combined_pt_test)),
    rep("PT_Timing", nrow(combined_timing_df)),
    rep("Q_obs", nrow(combined_Q_num)),
    as.character(bam_sensitivity_comb$prior_names)
  ),
  variant = c(
    as.character(combined_spacing$variant),
    as.character(combined_pt_test$variant),
    as.character(combined_timing_df$variant),
    as.character(combined_Q_num$variant),
    as.character(bam_sensitivity_comb$variant)
  ),
  rivname = c(
    as.character(combined_spacing$rivname),
    as.character(combined_pt_test$rivname),
    as.character(combined_timing_df$rivname),
    as.character(combined_Q_num$rivname),
    as.character(bam_sensitivity_comb$rivname)
  ),
  nse = c(
    combined_spacing$nse,
    combined_pt_test$nse,
    combined_timing_df$nse,
    combined_Q_num$nse,
    bam_sensitivity_comb$nash
  ),
  kge = c(
    combined_spacing$kge,
    combined_pt_test$kge,
    combined_timing_df$kge,
    combined_Q_num$kge,
    bam_sensitivity_comb$kge
  ),
  nrmse = c(
    combined_spacing$nrmse,
    combined_pt_test$nrmse,
    combined_timing_df$nrmse,
    combined_Q_num$nrmse,
    bam_sensitivity_comb$n_rmse
  ),
  rrmse = c(
    combined_spacing$rrmse,
    combined_pt_test$rrmse,
    combined_timing_df$rrmse,
    combined_Q_num$rrmse,
    bam_sensitivity_comb$nash
  ),
  pbias = c(
    combined_spacing$pbias,
    combined_pt_test$pbias,
    combined_timing_df$pbias,
    combined_Q_num$pbias,
    bam_sensitivity_comb$p_bias
  ),stringsAsFactors = FALSE
)

full_nash = full_nash[full_nash$variant=="manning",]
full_nash = full_nash[full_nash$test != "b_sd",]
full_nash = full_nash[full_nash$test != "b_hat",]

full_nash_E1 = full_nash[1:1553,]
full_nash_E2 = full_nash[1554:nrow(full_nash),]

saveRDS(full_nash_E1, "full_nash_E1.rds")
saveRDS(full_nash_E2, "full_nash_E2.rds")

full_nash_plot_E1 <-
  ggplot(data = full_nash_E1, aes(x = test, y = nse, fill = variant)) +
  geom_violin() + xlab("Experiment 1") +
  ylab("NSE") +
  ggtitle("BAM Sensitivity to PT Resolution \n (TOWNS dataset)") +
  scale_fill_manual(
    values = "dodgerblue") +theme_bw() + theme(text = element_text(size = 10),  plot.title = element_text(size = 11)) + 
  coord_flip(ylim = c(0.5, 1))+theme(legend.position = "none")+
  scale_x_discrete(labels = c("PT Timing", "PT Spacing", "PT Number"))

full_nash_plot_E1

full_nash_plot_E2 <-
  ggplot(data = full_nash_E2, aes(x = test, y = nse, fill = variant)) +
  geom_violin() + xlab("Experiment 2") +
  ylab("NSE") +
  ggtitle("BAM Sensitivity to BAM Parameters \n (TOWNS dataset)") +
  scale_fill_manual(
    values = "dodgerblue") +theme_bw() + theme(text = element_text(size = 10),  plot.title = element_text(size = 11)) + 
  coord_flip(ylim = c(0.5, 1))+theme(legend.position = "none")+
  scale_x_discrete(labels = c("sigmas", expression(sigma[italic(Werr)]~(9)), expression(sigma[italic(Serr)]~(2)),
                              expression(sigma[italic(dAerr)]~(8)), expression(log~sigma[italic(Ao)]~(2)),
                              expression(log~hat(italic(Ao))~(23)), expression(log~sigma[italic(n)]~(11) ),
                              expression(log~hat(italic(n))~(9)), expression(log~sigma[italic(Q)]~(7) ),
                              expression(log~hat(italic(Q))~(19)), expression(italic(numQ[obs])~(111))))

full_nash_plot_E2
png(
  "TOWNSS_NSE.png",
  width = 6.5,
  height = 6,
  units = 'in',
  res = 300
)
ggarrange(
  full_nash_plot_E1,
  full_nash_plot_E2,
  labels = c("A", "B"),
  nrow = 2,
  heights = c(5, 11)
)
dev.off()

#Plot four largest sensitivity results----
saveRDS(bam_sensitivity_comb, "bam_sensitivity_comb.rds")
bam_sensitivity_comb = readRDS("bam_sensitivity_comb.rds")
npg_pal = pal_npg("nrc")(9)
TOWNSS_colors = c(npg_pal[c(2, 3, 9, 6, 1, 5, 8)])
TOWNSS_colors = c(npg_pal[c(1, 5, 8, 3, 6, 2, 9)])
Q_hat = data.frame(bam_sensitivity_comb[bam_sensitivity_comb$prior_names == "logQ_hat",])
Q_hat_plot = ggplot(Q_hat[Q_hat$variant == "manning",], 
                    aes(x = prior_vec, y = nash, col = rivname))+
  geom_point()+geom_line()+xlab(expression(log~italic(hat(Q))))+
  ylab("NSE")+ggtitle(expression(Mean~Discharge~(italic(hat(Q)))))+theme_bw()+
  #theme(legend.position = "none")+
  theme(axis.title.y = element_blank(), text = element_text(size = 10), plot.title = element_text(size=11))+
  scale_color_manual(
    values = TOWNSS_colors[c(1:3, 5:7)],
    name = "River",
    breaks = c("Tan", "Olen", "Will", "Sac_rch1", "Sac_rch2", "Sac_rch3"),
    labels = c("Tanana", "Olentangy", "Willamette", "Sacramento R1", "Sacramento R2", "Sacramento R3"))
Q_hat_plot
Q_sd = data.frame(bam_sensitivity_comb[bam_sensitivity_comb$prior_names == "logQ_sd",])
Q_sd_plot = ggplot(Q_sd[Q_sd$variant == "manning",], aes(x = prior_vec, y = nash, col = rivname))+
  geom_point()+geom_line()+xlab(expression(log~sigma[italic(Q)]))+
  ylab("NSE")+ggtitle(expression(Discharge~Std.~Dev.~(sigma[italic(Q)])))+theme_bw()+
  #theme(legend.position = "none")+
  theme(axis.title.y = element_blank(), text = element_text(size = 10), plot.title = element_text(size=11))+
  coord_cartesian(ylim = c(-0.5,1), xlim = c(0, 0.5))+
  scale_color_manual(
    values = TOWNSS_colors[c(1:3, 5:7)],
    name = "River",
    breaks = c("Tan", "Olen", "Will","Sac_rch1", "Sac_rch2", "Sac_rch3"),
    labels = c("Tanana", "Olentangy", "Willamette","Sacramento R1", "Sacramento R2", "Sacramento R3"))
Q_sd_plot
pt_num_agg = aggregate(nse ~ pt_count+rivname, combined_pt_test[combined_pt_test$variant == "manning",], quantile)
pt_num_agg$rivname = factor(pt_num_agg$rivname, levels = unique(pt_num_agg$rivname)[c(5,4,6,7,1,2,3)])
PT_num_plot = ggplot(pt_num_agg, aes(x = pt_count, col = rivname,fill = rivname))+
  geom_ribbon(aes(ymin = pt_num_agg$nse[, 2], ymax = pt_num_agg$nse[, 4])) +
  geom_line(aes(y = pt_num_agg$nse[, 3]), size = 1.25)+ xlab("PT count")+
  ylab("NSE")+ggtitle("Number of PTs")+theme_bw()+
  #theme(legend.position = "none")+
  theme(axis.title.y = element_blank(), text = element_text(size = 10), plot.title = element_text(size=11))+
  scale_color_manual(
    values = TOWNSS_colors[c(6,7,5,1,2,3)],#6,4,7
    name = "River",
    breaks = c("Tan", "Olen", "Will", "Sac_rch1", "Sac_rch2", "Sac_rch3"),
    labels = c("Tanana", "Olentangy", "Willamette", "Sacramento R1", "Sacramento R2", "Sacramento R3"))+
  scale_fill_manual(values = alpha(TOWNSS_colors[c(6,7,5,1,2,3)], 0.2),
                    name = "River",
                    breaks = c("Tan", "Olen", "Will", "North_Sask", "Sac_rch1", "Sac_rch2", "Sac_rch3"),
                    labels = c("Tanana", "Olentangy", "Willamette", "N. Saskatchewan", "Sacramento R1", "Sacramento R2", "Sacramento R3"))
PT_num_plot
pt_spacing_agg = aggregate(nse ~ rounded_dist+rivname, combined_spacing[combined_spacing$variant == "manning",], quantile)
pt_spacing_agg$rivname = factor(pt_spacing_agg$rivname, levels = unique(pt_spacing_agg$rivname)[c(5,4,6,7,1,2,3)])
pt_spacing_agg$rounded_dist = pt_spacing_agg$rounded_dist/1000
PT_spacing_plot = ggplot(pt_spacing_agg, aes(x = rounded_dist, col = rivname, fill = rivname))+
  geom_ribbon(aes(ymin = pt_spacing_agg$nse[, 2], ymax = pt_spacing_agg$nse[, 4])) +
  geom_line(aes(y = pt_spacing_agg$nse[, 3]), size = 1.25)+ xlab("PT distance (km)")+
  ylab("NSE")+ggtitle("Spacing of PTs")+theme_bw()+
  #theme(legend.position = "none")+
  theme(axis.title.y = element_blank(), text = element_text(size = 10), plot.title = element_text(size=11))+
  scale_color_manual(
    values = TOWNSS_colors[c(6,7,5,1,2,3)],
    name = "River",
    breaks = c("Tan", "Olen", "Will", "Sac_rch1", "Sac_rch2", "Sac_rch3"),
    labels = c("Tanana", "Olentangy", "Willamette","Sacramento R1", "Sacramento R2", "Sacramento R3"))+
  scale_fill_manual(values = alpha(TOWNSS_colors[c(6,7,5,1,2,3)], 0.2),
                    name = "River",
                    breaks = c("Tan", "Olen", "Will", "Sac_rch1", "Sac_rch2", "Sac_rch3"),
                    labels = c("Tanana", "Olentangy", "Willamette",  "Sacramento R1", "Sacramento R2", "Sacramento R3"))
PT_spacing_plot

combined_timing_df$Time = factor(combined_timing_df$Time, levels = unique(combined_timing_df$Time)[c(1,2,3)])
timing_plot = ggplot(combined_timing_df, aes(x = Time, y = nse, col = rivname, group = rivname))+
  geom_point()+geom_line()+xlab("PT timestep")+
  ylab("NSE")+ggtitle("PT Temporal Res.")+theme_bw()+
  ylim(c(0.8, 1))+
  #theme(legend.position = "none")+
  theme(axis.title.y = element_blank(),text = element_text(size = 10), plot.title = element_text(size=11))+
  scale_color_manual(
    values = TOWNSS_colors[c(1:4,7,6)],
    name = "River",
    breaks = c("Tan", "Olen", "Will", "North_Sask", "Sac_rch1", "Sac_rch2", "Sac_rch3"),
    labels = c("Tanana", "Olentangy", "Willamette", "N. Saskatchewan", "Sacramento R1", "Sacramento R2", "Sacramento R3"))
timing_plot
Qobs_agg = aggregate(nse ~ Q+rivname, combined_Q_num, quantile)
Qobs_agg$rivname = factor(Qobs_agg$rivname, levels = unique(Qobs_agg$rivname)[c(5,3,4,6,1,2)])
Qobs_plot = ggplot(Qobs_agg, aes(x = Q, col = rivname, fill = rivname))+
  geom_ribbon(aes(ymin = Qobs_agg$nse[, 2], ymax = Qobs_agg$nse[, 4])) +
  geom_line(aes(y = Qobs_agg$nse[, 3]), size = 1.25)+ xlab("Discharge count")+
  ylab("NSE")+ggtitle(expression(Discharge~(italic(numQ[obs]))))+theme_bw()+
  theme(legend.position = "none")+
  theme(axis.title.y = element_blank(),text = element_text(size = 10), plot.title = element_text(size=11))+
  coord_cartesian(ylim = c(0,1))+
  scale_color_manual(
    values = TOWNSS_colors[c(6,7,5,1,2,3)],
    name = "River",
    breaks = c("Tan", "Olen", "Will", "Sac_rch1", "Sac_rch2", "Sac_rch3"),
    labels = c("Tanana", "Olentangy", "Willamette", "Sacramento R1", "Sacramento R2", "Sacramento R3"))+
  scale_fill_manual(values = alpha(TOWNSS_colors[c(6,7,5,1,2,3)], 0.2),
                    name = "River",
                    breaks = c("Tan", "Olen", "Will", "Sac_rch1", "Sac_rch2", "Sac_rch3"),
                    labels = c("Tanana", "Olentangy", "Willamette", "Sacramento R1", "Sacramento R2", "Sacramento R3"))
Qobs_plot


png(
  "Sensitivity.png",
  width = 6.5,
  height = 5,
  units = 'in',
  res = 300
)
Sensitivity = 
  ggarrange(
    PT_num_plot,
    PT_spacing_plot,
    timing_plot,
    Q_hat_plot,
    Q_sd_plot,
    Qobs_plot,
    labels = c("A", "B", "C", "D", "E", "F"),
    ncol = 3,
    nrow = 2,
    common.legend = TRUE,
    legend = "bottom"
  )
annotate_figure(Sensitivity,
                left = text_grob("                     Experiment 2                               Experiment 1 \n               (NSE)                                           (NSE)", size = 10, rot = 90))
dev.off()

#All error metrics combined-----
bam_sensitivity_metrics <- full_nash %>% gather(Metric, Value, 4:8)
bam_sensitivity_metrics$rivname <- factor(bam_sensitivity_metrics$rivname, 
                                          levels = unique(bam_sensitivity_metrics$rivname)[c(5, 4, 6, 7, 1, 2, 3)])
bam_sensitivity_metrics$Metric <- factor(bam_sensitivity_metrics$Metric, 
                                         levels = unique(bam_sensitivity_metrics$Metric)[c(1,2,3,4,5)])
npg_pal = pal_npg("nrc")(9)
TOWNSS_colors = c(npg_pal[c(2, 3, 9, 6, 1, 5, 8)])
TOWNSS_metrics_plot <-
  ggplot(data = bam_sensitivity_metrics, aes(x = Metric, y = Value, fill = rivname)) +
  geom_boxplot(outlier.size = 0.25) + xlab("Error Metrics") +
  ggtitle("BAM Sensitivity Error Metrics by River\n (TOWNSS dataset)") +
  scale_x_discrete(labels = c("nse" = "NSE (47)", "kge" = "KGE (4)", "nrmse" = "NRMSE", 
                              "rrmse" = "RRMSE (47)","pbias" = "rBIAS"))+
  scale_fill_manual(
    values = TOWNSS_colors[c(1,2,3,5:7)],
    name = "River",
    breaks = c("Tan", "Olen", "Will","Sac_rch1", "Sac_rch2", "Sac_rch3"),
    labels = c("Tanana", "Olentangy", "Willamette", "Sacramento R1", "Sacramento R2", "Sacramento R3")
  )+
  theme_bw() + theme(text = element_text(size = 10), plot.title = element_text(size = 11)) +
  coord_cartesian(ylim = c(-1, 1))#+theme(legend.position = "right")
TOWNSS_metrics_plot

png("TOWNS_Metrics.png",
    width = 6.5,
    height = 4,
    units = 'in',
    res = 300
)
TOWNSS_metrics_plot
dev.off()






