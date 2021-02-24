#Merritt Harlan
#Discharge from Arrays of Pressure Transducers (DAPT)
#Plotting Functions

library(ggplot2)
library(viridis)
library(tidyr)
library(ggsci)
library(scales)
library(ggpubr)

#Functions --------
plot_bam_priors <-
  function(bamdata,
           bam_priors_df,
           bam_sensitivity,
           n_levels,
           rivname,
           color) {
    #add bam variant to dataframe
    add_variant <- function(bam_sensitivity_df) {
      variant = bam_sensitivity_df[[1]]
      new_df = data.frame(cbind(variant, bam_sensitivity_df[[3]]))
      return(new_df)
    }
    add_prior <- function(bam_priors_col, n_levels = 10) {
      prior_values = unique(bam_priors_col)
      if (length(prior_values) > n_levels) {
        prior_values = prior_values[!(prior_values %in% prior_values[which.max(tabulate(match(bam_priors_col, prior_values)))])]
      }
      return(prior_values)
    }
    add_hydrograph <- function(bam_sensitivity_df) {
      hydro = data.frame(bam_sensitivity_df[[2]])
      return(hydro)
    }
    
    sensitivity_df <- data.frame()
    prior_values = apply(bam_priors_df[1:ncol(bam_priors_df) - 1], 2, add_prior)
    error_df = do.call(rbind.data.frame, (lapply(bam_sensitivity, add_variant)))
    prior_names = rep(colnames(prior_values), each = 3 * n_levels)
    prior_vec = rep(as.vector(apply(bam_priors_df[1:ncol(bam_priors_df) -
                                                    1], 2, add_prior)), each = 3)
    riv_df = cbind(rivname, prior_names, prior_vec, error_df)
    sensitivity_df = rbind(riv_df, sensitivity_df)
    hydro_length = length(bam_sensitivity[[1]][[2]])
    hydrograph_df = data.frame(
      prior_names = rep(sensitivity_df$prior_names, each = hydro_length),
      variant = rep(sensitivity_df$variant, each = hydro_length),
      Q = do.call("rbind", lapply(bam_sensitivity, add_hydrograph)),
      time = c(1:hydro_length),
      id = rep(1:(n_levels * 36), each = hydro_length)
    )
    colnames(hydrograph_df)[3] = "Q"
    
    sensitivity_plots = list()
    hydrograph_plots = list()
    
    for (p in 1:length(unique(sensitivity_df$prior_names))) {
      prior = unique(sensitivity_df$prior_names)[p]
      true_prior = unique(bam_priors_df[, p])[which.max(tabulate(match(
        bam_priors_df[, p], unique(bam_priors_df[, p])
      )))]
      subset_df = sensitivity_df[sensitivity_df$prior_names == prior, ]
      subset_hydro = hydrograph_df[hydrograph_df$prior_names == prior, ]
      require(ggplot2)
      sensitivity_plots[[p]] = ggplot(subset_df, aes(prior_vec, nash)) + theme_classic(base_size = 16) +
        geom_line(aes(linetype = variant), size = 1.2, color = color) + geom_point(size = 3, color = color, aes(shape = variant)) +
        ylab("NSE") + xlab(prior) + labs(fill = "Variant") + ggtitle(rivname) +
        geom_vline(xintercept = true_prior)
      hydrograph_plots[[p]] = ggplot() + theme_classic(base_size = 16) +
        geom_line(
          data = subset_hydro,
          aes(
            time,
            Q,
            group = id,
            linetype = variant,
            color = variant
          ),
          size = 1.2
        ) +
        geom_point(aes(x = c(1:length(
          bamdata$logQ_hat
        )), y = exp(bamdata$logQ_hat))) +
        ylab("Discharge") + xlab("Days") + labs(fill = "Variant") + ggtitle(paste0(rivname, "\n", prior))
    }
    pdf(paste0(rivname, "_Sensitivity.pdf"))
    print(sensitivity_plots)
    print(hydrograph_plots)
    dev.off()
    return(list(sensitivity_df, hydrograph_df))
  }

number_df <- function(num_pt_list, rivname) {
  manning_nse = unlist(x = lapply(num_pt_list, function(num_pt) {
    lapply(num_pt, function(x) {
      return(x[[2]][[3]]$nash)
    })
  }))
  manning_kge = unlist(x = lapply(num_pt_list, function(num_pt) {
    lapply(num_pt, function(x) {
      return(x[[2]][[3]]$kge)
    })
  }))
  manning_nrmse = unlist(x = lapply(num_pt_list, function(num_pt) {
    lapply(num_pt, function(x) {
      return(x[[2]][[3]]$n_rmse)
    })
  }))
  manning_rrmse = unlist(x = lapply(num_pt_list, function(num_pt) {
    lapply(num_pt, function(x) {
      return(x[[2]][[3]]$rrmse)
    })
  }))
  manning_pbias = unlist(x = lapply(num_pt_list, function(num_pt) {
    lapply(num_pt, function(x) {
      return(x[[2]][[3]]$p_bias)
    })
  }))
  pt_count = as.numeric(unlist(lapply(num_pt_list, function(num_pt) {
    lapply(num_pt, function(x) {
      return(length(x[[1]]))
    })
  })))
  combined_df = data.frame(
    rivname = rivname,
    pt_count = pt_count,
    nse = manning_nse,
    kge = manning_kge,
    nrmse = manning_nrmse,
    rrmse = manning_rrmse,
    pbias = manning_pbias,
    variant = rep("manning",length(pt_count))
  )
  return(combined_df)
}

number_hydro <- function(num_pt_list, rivname){
  manning_hydro = unlist(lapply(num_pt_list, function(num_pt) {
    lapply(num_pt, function(x) {
      return(x[[2]][[2]])
    })
  }))
  pt_count = as.numeric(unlist(lapply(num_pt_list, function(num_pt) {
    lapply(num_pt, function(x) {
      return(length(x[[1]]))
    })
  })))
  combined_df = data.frame(
    rivname = rivname,
    pt_count = rep(pt_count, each = length(num_pt_list[[1]][[1]][[2]][[2]])),
    time = c(1:length(num_pt_list[[1]][[1]][[2]][[2]])),
    Q = c(manning_hydro),
    variant = "manning"
  )
  return(combined_df)
}

spacing_df <- function(spacing_pt, rivname, riv_dist) {
  manning_nse = unlist(lapply(spacing_pt, function(x) {
    return(x[[2]][[3]]$nash)
  }))
  manning_kge = unlist(lapply(spacing_pt, function(x) {
    return(x[[2]][[3]]$kge)
  }))
  manning_nrmse = unlist(lapply(spacing_pt, function(x) {
    return(x[[2]][[3]]$n_rmse)
  }))
  manning_rrmse = unlist(lapply(spacing_pt, function(x) {
    return(x[[2]][[3]]$rrmse)
  }))
  manning_pbias = unlist(lapply(spacing_pt, function(x) {
    return(x[[2]][[3]]$p_bias)
  }))
  pt_space = as.numeric(lapply(spacing_pt, function(x) {
    return(riv_dist[x[[1]][2]] - riv_dist[x[[1]][1]])
  }))
  rounded_dist = 5000 * round(round(pt_space,-3) / 5000)
  combined_df = data.frame(
    rivname = rivname,
    pt_space = pt_space,
    rounded_dist = rounded_dist,
    nse = manning_nse,
    kge = manning_kge,
    nrmse = manning_nrmse,
    rrmse = manning_rrmse,
    pbias = manning_pbias,
    variant = rep("manning", length(pt_space))
  )
  return(combined_df)
}

spacing_hydro <- function(spacing_pt, riv_dist, rivname){
  manning_hydro = unlist(lapply(spacing_pt, function(x) {
    return(x[[2]][[2]])
  }))
  pt_space = as.numeric(lapply(spacing_pt, function(x) {
    return(riv_dist[x[[1]][2]] - riv_dist[x[[1]][1]])
  }))
  rounded_dist = 5000 * round(round(pt_space,-3) / 5000)
  combined_df = data.frame(
    rivname = rivname,
    rounded_dist = rep(rounded_dist, each = length(spacing_pt[[1]][[2]][[2]])),
    time = c(1:length(spacing_pt[[1]][[2]][[2]])),
    Q = c(manning_hydro),
    variant = rep(c("manning"), each = length(manning_hydro))
  )
  return(combined_df)
}

timing_df <- function(timing, rivname) {
  manning_amhg_nse = unlist(lapply(timing, function(x) {
    return(x[[1]][[3]]$nash)
  }))
  manning_amhg_kge = unlist(lapply(timing, function(x) {
    return(x[[1]][[3]]$kge)
  }))
  manning_amhg_nrmse = unlist(lapply(timing, function(x) {
    return(x[[1]][[3]]$n_rmse)
  }))
  manning_amhg_rrmse = unlist(lapply(timing, function(x) {
    return(x[[1]][[3]]$rrmse)
  }))
  manning_amhg_pbias = unlist(lapply(timing, function(x) {
    return(x[[1]][[3]]$p_bias)
  }))
  combined_df = data.frame(
    rivname = rivname,
    nse = manning_amhg_nse,
    kge = manning_amhg_kge,
    nrmse = manning_amhg_nrmse,
    rrmse = manning_amhg_rrmse,
    pbias = manning_amhg_pbias,
    variant = c("manning_amhg"),
    Time = rep(c("Daily", "6-hr", "Hourly"))
  )
  return(combined_df)
}

num_Q_df <- function(num_Q_list, rivname) {
  manning_amhg_nse = unlist(lapply(num_Q_list, function(num_Q) {
    lapply(num_Q, function(x) {
      return(x[[2]][[3]]$nash)
    })
  }))
  manning_amhg_kge = unlist(x = lapply(num_Q_list, function(num_Q) {
    lapply(num_Q, function(x) {
      return(x[[2]][[3]]$kge)
    })
  }))
  manning_amhg_nrmse = unlist(x = lapply(num_Q_list, function(num_Q) {
    lapply(num_Q, function(x) {
      return(x[[2]][[3]]$n_rmse)
    })
  }))
  manning_amhg_rrmse = unlist(x = lapply(num_Q_list, function(num_Q) {
    lapply(num_Q, function(x) {
      return(x[[2]][[3]]$rrmse)
    })
  }))
  manning_amhg_pbias = unlist(x = lapply(num_Q_list, function(num_Q) {
    lapply(num_Q, function(x) {
      return(x[[2]][[3]]$p_bias)
    })
  }))
  Q = as.numeric(unlist(lapply(num_Q_list, function(num_Q) {
    lapply(num_Q, function(x) {
      return(length(x[[1]]))
    })
  })))
  combined_df = data.frame(
    rivname = rivname,
    Q = Q,
    nse = manning_amhg_nse,
    kge = manning_amhg_kge,
    nrmse = manning_amhg_nrmse,
    rrmse = manning_amhg_rrmse,
    pbias = manning_amhg_pbias,
    variant = rep("manning_amhg",length(Q))
  )
  return(combined_df)
}

Q_hydro <- function(num_Q_list, rivname){
  manning_amhg_hydro = unlist(lapply(num_Q_list, function(num_Q) {
    lapply(num_Q, function(x) {
      return(x[[2]][[2]])
    })
  }))
  Q_count = as.numeric(unlist(lapply(num_Q_list, function(num_Q) {
    lapply(num_Q, function(x) {
      return(length(x[[1]]))
    })
  })))
  combined_df = data.frame(
    rivname = rivname,
    pt_count = rep(Q_count, each = length(num_Q_list[[1]][[1]][[2]][[2]])),
    time = c(1:length(num_Q_list[[1]][[1]][[2]][[2]])),
    Q = c(manning_amhg_hydro),
    variant = rep("manning_amhg",length(manning_amhg_hydro))
  )
  return(combined_df)
}

Q_ids <- function(num_Q_list) {
  id = lapply(num_Q_list, function(num_Q) {
    lapply(num_Q, function(x) {
      return(x[[1]])
    })
  })
}

TG_NSE <- function(tg_all, rivname) {
  nse_df = data.frame()
  for (Q in 1:length(tg_all)) {
    for (i in 1:length(tg_all[[Q]])) {
      for (pt in 1:length(tg_all[[Q]][[i]])) {
        nse = tg_all[[Q]][[i]][[pt]][[2]]
        df = data.frame(
          nse = nse,
          Q = Q + 1,
          variant = "TG",
          rivname = rivname
        )
        nse_df = rbind(df, nse_df)
      }
    }
  }
  return(nse_df)
}

plot_hydro <- function(df, bamdata, rivname) {
  df = df[df$variant != "manning_amhg",]
  df$variant = "DAPT Method"
  hydro = aggregate(Q ~ time + variant, df, quantile)
  obs_df = data.frame(variant = "Observed", Q = exp(bamdata$logQ_hat), time = c(1:length(bamdata$logQ_hat)))
  hydro = rbind(hydro, obs_df)
  colnames(hydro)[2] = "Legend"
  background_col = alpha("dodgerblue", 0.2)
  ggplot(hydro, aes(x = time, col = Legend, fill = Legend, shape = Legend)) +
    geom_ribbon(aes(ymin = hydro$Q[, 2], ymax = hydro$Q[, 4])) +
    geom_line(aes(y = hydro$Q[, 3])) + 
    geom_point(aes(y = hydro$Q[,3]), color = "black")+
    ggtitle(paste0(rivname)) + 
    theme_bw() + theme(axis.title.x=element_blank(), axis.title.y = element_blank()) +
    theme(legend.position = "none")+
    theme(text = element_text(size = 10), plot.title = element_text(size = 11)) +
    scale_colour_manual(name = "Legend", values = c("dodgerblue", NA))+
    scale_fill_manual(values = c(alpha("dodgerblue", 0.2),NA))+
    scale_shape_manual(values = c(NA,19))+
    scale_x_continuous(breaks = 2*round(c(seq(0, (length(bamdata$logQ_hat)), length.out = 4))/2), 
                       labels = 2*round(c(seq(0, (length(bamdata$logQ_hat)/2), length.out = 4))/3))+
    scale_y_continuous(breaks = 5*round(c(seq(min(exp(bamdata$logQ_hat)), max(exp(bamdata$logQ_hat)), length.out = 4))/5))
}

plot_pad_hydro <- function(num_Q_list, bamdata, Q_df, Q_id, rivname){
  PAD_hydro_df = data.frame(time = c(1:length(bamdata$logQ_hat)),
                            best = c(exp(bamdata$logQ_hat), num_Q_list[[length(num_Q_list)]][[1]][[2]][[2]]),
                            Q1 = c(exp(bamdata$logQ_hat),aggregate(Q~time, Q_df, quantile)$Q[,2]),
                            Q3 = c(exp(bamdata$logQ_hat), aggregate(Q~time, Q_df, quantile)$Q[,4]),
                            Legend = c(rep(c("Observed", "DAPT Method"), each = length(bamdata$logQ_hat)))
  )
  PAD_hydro_df$Legend = factor(PAD_hydro_df$Legend, levels = unique(PAD_hydro_df$Legend)[c(2,3,1)])
  ggplot(PAD_hydro_df, aes(x = time, col = Legend, fill = Legend, shape = Legend))+
    ggtitle(rivname)+
    geom_line(aes(y = best), size = 1.25)+
    theme_bw()+
    theme(text = element_text(size = 10), plot.title = element_text(size = 11)) +
    theme(axis.title.x=element_blank(),axis.title.y = element_blank())+
    theme(legend.position = "none")+
    geom_ribbon(aes(ymin = Q1, ymax = Q3))+
    geom_point(aes(y = best), col = "black", size = 1.25)+
    scale_shape_manual(values = c("Observed" = 19, "DAPT Method" = NA))+
    scale_color_manual(values = c("Observed" = NA, "DAPT Method" = "dodgerblue"))+
    scale_fill_manual(values = c("Observed" = NA, "DAPT Method" = alpha("dodgerblue", 0.2)))+
    scale_x_continuous(breaks = 2*round(c(seq(0, length(bamdata$logQ_hat), length.out = 4))/2))
}

PAD_hydrographs <- function(WSE_all, WSE_Q, Q, bamdata, xvec, rivname, QMOD){
  #Calculate dA
  calcdA_mat <- function(w, h) {
    stopifnot(all(dim(w) == dim(h)))
    dA <- w
    for (i in 1:nrow(dA)) {
      dA[i, ] <- calcdA_vec(w[i, ], h[i, ])
    }
    dA
  }
  calcdA_vec <- function(w, h) {
    words <- order(w)
    warr <- w[words]
    harr <- h[words]
    delh <- c(0, diff(harr))
    delA <- cumsum(warr * delh)
    dA <- 1:length(w)
    dA[words] <- delA
    dA
  }
  
  #Calculate slope
  calcslope <- function(xvec, hmat) {
    dH <- apply(hmat, 2, function(x) c(NA, -diff(x)))
    dX_vec <- c(NA, diff(xvec))
    dX <- matrix(dX_vec, ncol = ncol(hmat), nrow = nrow(hmat),
                 byrow = FALSE)
    out <- dH / dX
  }
  if(rivname == "Rochers"){
    WSE_all$date = as.Date(as.character(WSE_all$new_time), format = "%m/%d/%Y")
  }else if(rivname == "Fletcher"){
    WSE_all$date = as.Date(as.character(WSE_all$new_time), format = "%m/%d/%Y")
  }else{
    WSE_all$date = as.Date(WSE_all$new_time)
  }
  
  WSE_all = aggregate(.~date, data = WSE_all, FUN = mean)
  height = t(WSE_all[, 4:ncol(WSE_all)])
  width = matrix(bamdata$Wobs[1,1], nrow = (ncol(WSE_all)-3), ncol = nrow(WSE_all))
  dA = calcdA_mat(width, height)
  slope = calcslope(xvec = xvec, hmat = height)
  slope = rbind(colMeans(slope, na.rm = TRUE), slope)
  slope = slope[c(1, 3:nrow(slope)),]
  slope[slope<=0] = min(slope[slope>0], na.rm = TRUE)
  slope[is.na(slope)] = mean(slope, na.rm = TRUE)
  data = bam_data(w = width, s = slope, dA = dA, Qhat = exp(bamdata$logQ_hat))
  prior_df = generate_prior_df(var_scale = 1,
                               bamdata = data,
                               n_levels = 1,
                               true = 1
  )
  data_mod = bam_data(w = width, s = slope, dA = dA, Qhat = QMOD)
  prior_df_mod = generate_prior_df(var_scale = 1,
                                   bamdata = data_mod,
                                   n_levels = 1,
                                   true = 1
  )
  manning_amhg = run_and_validate_bam(bamdata = data, new_df = prior_df[2,],n_levels = 1, qval = c(1:ncol(width)))
  manning_amhg_mod = run_and_validate_bam(bamdata = data_mod, new_df = prior_df_mod[2,],n_levels = 1, qval = c(1:ncol(width)))
  dat = data.frame(stage = as.numeric(rowMeans(WSE_Q)), discharge = exp(bamdata$logQ_hat))
  full_dat = data.frame(stage = as.numeric(colMeans(height)))
  model = lm(log(discharge) ~ log(stage), data = dat)
  pred = data.frame(p = exp(predict(model, full_dat)))
  Q$new_date = as.POSIXct(Q$new_date, format = "%Y-%m-%d %H:%M")
  WSE_all$date = as.POSIXct(WSE_all$date, format = "%Y-%m-%d %H:%M")
  combined_Q = data.frame(time = rep(WSE_all$date, 3), Q = c(pred$p,manning_amhg[[2]], manning_amhg_mod[[2]]),
                          model = c(rep("Temporary Gauge", length(WSE_all$date)),
                                    rep("Manning-AMHG", length(WSE_all$date)),
                                    rep("Model", length(WSE_all$date))))
  #png("hydro_legend.png", res = 300, width = 6.5, height = 4, units = "in")
  ggplot(combined_Q, aes(x = time, y = Q))+geom_line(lwd = 1, aes(col = model, linetype = model))+theme_bw()+
    ylab("")+#ylab("Discharge (m3/s)")+xlab("Date")+
    scale_color_manual(values = c(viridis(6)[c(3)],viridis(6)[c(3)], "red"),name = "Legend") + 
    theme(text = element_text(size = 10)) + ggtitle(paste0(rivname)) + 
    theme(legend.position = "none", axis.title.x=element_blank())+
    #coord_cartesian(ylim = c(0, 200))+
    geom_point(data = Q, aes(x = new_date, y = Q))+scale_x_datetime(labels = date_format("%m/%d"), date_breaks = "10 days")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))#+theme(legend.position = "bottom")+xlab("ADCP")
}

plot_full_hydro <- function(full_BAM,rivname,WSE_all, Q){
  full_BAM_Q = Q_hydro(full_BAM, rivname)
  combined_stats = aggregate(Q ~ variant + time, full_BAM_Q, quantile)
  if(rivname == "Rochers R1"){
    WSE_all$date = as.Date(as.character(WSE_all$new_time), format = "%m/%d/%Y")
  }else if(rivname == "Fletcher"){
    WSE_all$date = as.Date(WSE_all$new_time, format = "%m/%d/%Y")
  }else if(rivname == "Rochers R2"){
    WSE_all$date = as.Date(WSE_all$new_time, format = "%m/%d/%Y")
  }else{
    WSE_all$date = as.Date(WSE_all$new_time)
  }
  
  WSE_all = aggregate(.~date, data = WSE_all, FUN = mean)
  WSE_all$date = as.POSIXct(WSE_all$date, format = "%Y-%m-%d %H:%M")
  combined_stats$date = as.POSIXct(WSE_all$date, format = "%Y-%m-%d %H:%M")[1:nrow(combined_stats)]
  
  Q$new_date = as.POSIXct(Q$new_date, format = "%Y-%m-%d %H:%M")
  Q_df = data.frame(variant = "Observed", date = Q$new_date, Q = Q$Q)
  
  #if(rivname == "Rochers R2"){
  # combined_stats = combined_stats[1:36,]
  #}
  combined_stats = rbind(combined_stats[,c(1,3,4)], Q_df)
  if(rivname == "Mamawi"){
    combined_stats = rbind(combined_stats, data.frame(variant = "Gauge", date = WSE_all$date, Q = c(188, 174, 164, 160, 155, 152, 154, 157, 157, 150, 142, 135, 131)))
    ggplot(combined_stats, aes(x = date, col = variant, fill = variant, shape = variant, linetype = variant))+
      geom_ribbon(aes(ymin = combined_stats$Q[, 2], ymax = combined_stats$Q[, 4])) +
      geom_line(aes(y = combined_stats$Q[, 3]), size = 1.25) +
      ggtitle(paste0(rivname)) + ylab("")+
      theme_bw() + theme(legend.position = "none", axis.title.x=element_blank(), axis.title.y = element_blank()) + 
      theme(text = element_text(size = 10), plot.title = element_text(size = 11)) +
      geom_point(aes(y = Q[,3]), col = "black")+
      scale_fill_manual(values = c(alpha("dodgerblue",0.2), NA, NA),
                        name = "Legend",
                        breaks = c("manning_amhg","Observed", "Gauge"),
                        labels = c("DAPT Method", "ADCP", "Nearby Gauge (not mass conserved)")
      ) +
      #ylim(0,200)+
      scale_colour_manual(
        values = c("dodgerblue", NA, "gray60"),
        name = "Legend",
        breaks = c("manning_amhg", "Observed", "Gauge"),
        labels = c("DAPT Method", "ADCP", "Nearby Gauge (not mass conserved)")
      )+
      scale_shape_manual(
        values = c(NA, 19, NA),
        name = "Legend",
        breaks = c("manning_amhg", "Observed", "Gauge"),
        labels = c("DAPT Method", "ADCP", "Nearby Gauge (not mass conserved)")
      )+#scale_x_datetime(labels = date_format("%m/%d"), date_breaks = "10 days")+
      scale_linetype_manual(
        values = c(1,1,2),
        name = "Legend",
        breaks = c("manning_amhg","Observed", "Gauge"),
        labels = c("DAPT Method", "ADCP", "Nearby Gauge (not mass conserved)") 
      )+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }else if(rivname == "Richardson"){
    combined_stats = rbind(combined_stats, data.frame(variant = "Gauge", date = WSE_all$date, 
                                                      Q = c(13.6, 13.4, 13.3, 13.1, 13, 12.9, 12.8, 12.7, 12.5, 
                                                            12.3, 12.2, 12, 11.9, 11.8, 12, 12.5, 12.4, 12.3, 12.8, 12.9,
                                                            12.8, 12.6, 12.7, 13, 13.2, 13.6, 13.9, 13.8, 14.4, 14.6,
                                                            14.9, 15.2, 15.3, 15.1, 14.9, 18.3, 21.5)))
    ggplot(combined_stats, aes(x = date, col = variant, fill = variant, shape = variant, linetype = variant))+
      geom_ribbon(aes(ymin = combined_stats$Q[, 2], ymax = combined_stats$Q[, 4])) +
      geom_line(aes(y = combined_stats$Q[, 3]), size = 1.25) +
      ggtitle(paste0(rivname)) + ylab("")+
      theme_bw() + theme(legend.position = "none", axis.title.x=element_blank(), axis.title.y = element_blank()) + 
      theme(text = element_text(size = 10), plot.title = element_text(size = 11)) +
      geom_point(aes(y = Q[,3]), col = "black")+
      scale_fill_manual(values = c(alpha("dodgerblue",0.2), NA, NA),
                        name = "Legend",
                        breaks = c("manning_amhg", "Observed", "Gauge"),
                        labels = c("DAPT Method", "ADCP", "Nearby Gauge (not mass conserved)")
      ) +
      #ylim(0,200)+
      scale_colour_manual(
        values = c("dodgerblue", NA, "gray60"),
        name = "Legend",
        breaks = c("manning_amhg", "Observed", "Gauge"),
        labels = c("DAPT Method", "ADCP", "Nearby Gauge (not mass conserved)")
      )+
      scale_shape_manual(
        values = c(NA, 19, NA),
        name = "Legend",
        breaks = c("manning_amhg", "Observed", "Gauge"),
        labels = c("DAPT Method",  "ADCP", "Nearby Gauge (not mass conserved)")
      )+#scale_x_datetime(labels = date_format("%m/%d"), date_breaks = "10 days")+
      scale_linetype_manual(
        values = c(1,1,2),
        name = "Legend",
        breaks = c("manning_amhg", "Observed", "Gauge"),
        labels = c("DAPT Method", "ADCP", "Nearby Gauge (not mass conserved)") 
      )+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }else if(rivname == "Athabasca R1"){
    combined_stats = rbind(combined_stats, data.frame(variant = "Gauge", date = WSE_all$date, 
                                                      Q = c(1190, 1200, 1200, 1220, 1200, 1140, 1070, 1010, 959, 924, 905,
                                                            882, 876, 882, 894, 889, 852, 806, 775, 758, 744, 722, 696, 674,663,
                                                            674, 687, 701, 693, 680, 678, 671, 667, 676, 670, 663, 661, 653, 654, 
                                                            655, 647, 638, 636, 644, 650, 646, 639, 637, 631, 618, 610, 606, 608,
                                                            612, 613, 615)))
    ggplot(combined_stats, aes(x = date, col = variant, fill = variant, shape = variant, linetype = variant))+
      geom_ribbon(aes(ymin = combined_stats$Q[, 2], ymax = combined_stats$Q[, 4])) +
      geom_line(aes(y = combined_stats$Q[, 3]), size = 1.25) +
      ggtitle(paste0(rivname)) + ylab("")+
      theme_bw() + theme(legend.position = "none", axis.title.x=element_blank(), axis.title.y = element_blank()) + 
      theme(text = element_text(size = 10), plot.title = element_text(size = 11)) +
      geom_point(aes(y = Q[,3]), col = "black")+
      scale_fill_manual(values = c(alpha("dodgerblue",0.2), NA, NA),
                        name = "Legend",
                        breaks = c("manning_amhg","Observed", "Gauge"),
                        labels = c("DAPT Method", "ADCP", "Nearby Gauge (not mass conserved)")
      ) +
      #ylim(0,200)+
      scale_colour_manual(
        values = c("dodgerblue", NA, "gray60"),
        name = "Legend",
        breaks = c("manning_amhg", "Observed", "Gauge"),
        labels = c("DAPT Method", "ADCP", "Nearby Gauge (not mass conserved)")
      )+
      scale_shape_manual(
        values = c(NA, 19, NA),
        name = "Legend",
        breaks = c("manning_amhg", "Observed", "Gauge"),
        labels = c("DAPT Method", "ADCP", "Nearby Gauge (not mass conserved)")
      )+#scale_x_datetime(labels = date_format("%m/%d"), date_breaks = "10 days")+
      scale_linetype_manual(
        values = c(1,1,2),
        name = "Legend",
        breaks = c("manning_amhg", "Observed", "Gauge"),
        labels = c("DAPT Method", "ADCP", "Nearby Gauge (not mass conserved)") 
      )+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }else if(rivname == "Athabasca R2"){
    combined_stats = rbind(combined_stats, data.frame(variant = "Gauge", date = WSE_all$date, 
                                                      Q = c(1190, 1200, 1200, 1220, 1200, 1140, 1070, 1010, 959, 924, 905,
                                                            882, 876, 882, 894, 889, 852, 806, 775, 758, 744, 722, 696, 674,663,
                                                            674, 687, 701, 693)))
    ggplot(combined_stats, aes(x = date, col = variant, fill = variant, shape = variant, linetype = variant))+
      geom_ribbon(aes(ymin = combined_stats$Q[, 2], ymax = combined_stats$Q[, 4])) +
      geom_line(aes(y = combined_stats$Q[, 3]), size = 1.25) +
      ggtitle(paste0(rivname)) + ylab("")+
      theme_bw() + theme(legend.position = "none", axis.title.x=element_blank(), axis.title.y = element_blank()) + 
      theme(text = element_text(size = 10), plot.title = element_text(size = 11)) +
      geom_point(aes(y = Q[,3]), col = "black")+
      scale_fill_manual(values = c(alpha("dodgerblue",0.2), NA, NA),
                        name = "Legend",
                        breaks = c("manning_amhg", "Observed", "Gauge"),
                        labels = c("DAPT Method", "ADCP", "Gauge")
      ) +
      #ylim(0,200)+
      scale_colour_manual(
        values = c("dodgerblue", NA, "gray60"),
        name = "Legend",
        breaks = c("manning_amhg", "Observed", "Gauge"),
        labels = c("DAPT Method", "ADCP", "Gauge")
      )+
      scale_shape_manual(
        values = c(NA, 19, NA),
        name = "Legend",
        breaks = c("manning_amhg", "Observed", "Gauge"),
        labels = c("DAPT Method", "ADCP", "Gauge")
      )+#scale_x_datetime(labels = date_format("%m/%d"), date_breaks = "10 days")+
      scale_linetype_manual(
        values = c(1,1,2),
        name = "Legend",
        breaks = c("manning_amhg", "Observed", "Gauge"),
        labels = c("DAPT Method", "ADCP", "Nearby Gauge (not mass conserved)") 
      )+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }else{
    ggplot(combined_stats, aes(x = date, col = variant, fill = variant, shape = variant))+
      geom_ribbon(aes(ymin = combined_stats$Q[, 2], ymax = combined_stats$Q[, 4])) +
      geom_line(aes(y = combined_stats$Q[, 3]), size = 1.25) +
      ggtitle(paste0(rivname)) + ylab("")+
      theme_bw() + theme(legend.position = "none", axis.title.x=element_blank(), axis.title.y = element_blank()) + 
      theme(text = element_text(size = 10), plot.title = element_text(size = 11)) +
      geom_point(aes(y = Q[,3]), col = "black")+
      scale_fill_manual(values = c(alpha("dodgerblue",0.2), NA),
                        name = "Legend",
                        breaks = c("manning_amhg", "Observed"),
                        labels = c("DAPT Method", "ADCP")
      ) +
      #ylim(0,200)+
      scale_colour_manual(
        values = c("dodgerblue",  NA),
        name = "Legend",
        breaks = c("manning_amhg", "Observed"),
        labels = c("DAPT Method", "ADCP")
      )+
      scale_shape_manual(
        values = c(NA, 19),
        name = "Legend",
        breaks = c("manning_amhg", "Observed"),
        labels = c("DAPT Method", "ADCP")
      )+#scale_x_datetime(labels = date_format("%m/%d"), date_breaks = "10 days")+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
}