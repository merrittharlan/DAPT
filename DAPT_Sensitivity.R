#Merritt Harlan
#Discharge from Pressure Transducers (DAPT)

#Sensitivity Analysis + Temporary Gauge Comparison
library(bamr)

#Functions --------
#Run BAM
generate_prior_df <- function(bamdata, n_levels, true, var_scale) {
  #define physical bounds-----------
  n_phys_bounds = c(0.02, 0.15)
  logn_phys_bounds = log(c(0.02, 0.15))
  logQ_phys_bounds = c(max(mean(bamdata$logQ_hat) - 4 * sd(bamdata$logQ_hat), 0),
                       mean(bamdata$logQ_hat) + 4 * sd(bamdata$logQ_hat))
  Q_phys_bounds = exp(logQ_phys_bounds)
  A0_phys_bounds = c((
    min(n_phys_bounds) * min(Q_phys_bounds) * min(bamdata$Wobs) ^ (2 / 3) * max(bamdata$Sobs) ^
      (-1 / 2)
  ) ^ (3 / 5)   ,
  (
    max(n_phys_bounds) * max(Q_phys_bounds) * max(bamdata$Wobs) ^ (2 / 3) * min(bamdata$Sobs[bamdata$Sobs>0]) ^
      (-1 / 2)
  ) ^ (3 / 5))
  logA0_phys_bounds = log(A0_phys_bounds)
  b_phys_bounds = c(0.01, 0.8)
  #------------------------------
  
  #define hats-------------------
  generate_hats = function(phys_bounds, n_levels) {
    hats = seq(min(phys_bounds), max(phys_bounds), length.out = n_levels)
  }
  
  logn_hat = generate_hats(logn_phys_bounds, n_levels)
  b_hat = generate_hats(b_phys_bounds, n_levels)
  logA0_hat = generate_hats(logA0_phys_bounds, n_levels)
  logQ_hat = generate_hats(logQ_phys_bounds, n_levels)
  #------------------------------
  
  #define range of sigmas that define bounds-------
  range_sds = seq(0.1, 10, length.out = n_levels)
  
  #get true values for all parameters------
  true_df = data.frame(
    logn_hat = log(0.035),
    b_hat = 0.01,
    logA0_hat = log((
      0.035 * mean(exp(bamdata$logQ_hat)) * mean(bamdata$Wobs) ^ (2 / 3) * mean(bamdata$Sobs) ^
        (-1 / 2)
    ) ^ (3 / 5)),
    logQ_hat = median(bamdata$logQ_hat),
    logn_sd = 0.5,
    b_sd = 0.005,
    logA0_sd = sd(log(abs(bamdata$dAobs[bamdata$dAobs!=0]))),
    logQ_sd = max(var_scale*sd(bamdata$logQ_hat),0.0025),
    Werr_sd = 30,
    dAerr_sd = 30,
    Serr_sd = 0.1 * mean(bamdata$Sobs),
    sigmas = 3
  )
  if (true == 1) {
    repped_3x = true_df[rep(seq_len(nrow(true_df)), each = 3),]
    variant = data.frame(variant = rep(c("amhg", 'manning_amhg', 'manning'), nrow(true_df)))
    true_df = cbind(repped_3x, variant)
    return(true_df)
  } else{
    #define sigmas-------------------
    logn_sd = range_sds * (log(0.035) - log(0.02)) / 2
    b_sd = range_sds * (0.25 - 0.01) / 9
    logA0_sd = abs(range_sds /5 * (true_df$logA0_hat - min(logA0_phys_bounds)))
    logQ_sd = range_sds / 5 * sd(bamdata$logQ_hat)
    Werr_sd = range_sds * 3
    dAerr_sd = range_sds * 3
    Serr_sd = range_sds * 0.1 * mean(bamdata$Sobs)
    #--------------------------------
    
    #generate a prior combo bamdataframe--------
    #60 tests, using true_df value as xxxhat value
    giant_df = true_df[rep(seq_len(nrow(true_df)), each = n_levels * 12),]
    giant_df$logn_hat[1:n_levels] = logn_hat
    giant_df$b_hat[(n_levels + 1):(2 * n_levels)] = b_hat
    giant_df$logA0_hat[(2 * n_levels + 1):(3 * n_levels)] = logA0_hat
    giant_df$logQ_hat[(3 * n_levels + 1):(4 * n_levels)] = logQ_hat
    giant_df$logn_sd[(4 * n_levels + 1):(5 * n_levels)] = logn_sd
    giant_df$b_sd[(5 * n_levels + 1):(6 * n_levels)] = b_sd
    giant_df$logA0_sd[(6 * n_levels + 1):(7 * n_levels)] = logA0_sd
    giant_df$logQ_sd[(7 * n_levels + 1):(8 * n_levels)] = logQ_sd
    giant_df$Werr_sd[(8 * n_levels + 1):(9 * n_levels)] = Werr_sd
    giant_df$dAerr_sd[(9 * n_levels + 1):(10 * n_levels)] = dAerr_sd
    giant_df$Serr_sd[(10 * n_levels + 1):(11 * n_levels)] = Serr_sd
    giant_df$sigmas[(11 * n_levels + 1):(12 * n_levels)] = range_sds
    
    #now, triple it to 180 tests.
    repped_3x = giant_df[rep(seq_len(nrow(giant_df)), each = 3),]
    variant = data.frame(variant = rep(c("amhg", 'manning_amhg', 'manning'), nrow(giant_df)))
    new_df = cbind(repped_3x, variant)
    
    return(new_df)
  }
}

run_and_validate_bam <- function(bamdata, new_df, n_levels, qval) {
  #-----------------------------------------
  
  generate_bam_priors <- function(bamdata, new_df) {
    #we don't set sigma Man, sigma AMHG. these are theoretical
    #we don't set Wc params, these are calc from data
    
    #we are deriving A0hat, bhat, nhat from other priors
    
    #each_lowerbound=max(physbound,mean- x*sd)
    #each_upperbound=min(physbound,mean+ x*sd)
    
    #all hats passed to function
    priors_out = bam_priors(
      bamdata = bamdata,
      #these are calculated from input hats and sds
      lowerbound_logQ = max(new_df$logQ_hat - (new_df$sigmas *
                                                 new_df$logQ_sd), -3),
      upperbound_logQ = new_df$logQ_hat + (new_df$sigmas *
                                             new_df$logQ_sd),
      lowerbound_A0 = max(exp(new_df$logA0_hat) - (exp(
        new_df$sigmas * new_df$logA0_sd
      )), 0),
      upperbound_A0 = exp(new_df$logA0_hat) + exp(new_df$sigmas *
                                                    new_df$logA0_sd),
      lowerbound_logn = max(
        log(0.02),
        new_df$logn_hat - (new_df$sigmas * new_df$logn_sd)
      ),
      upperbound_logn  = min(
        log(0.15),
        new_df$logn_hat + (new_df$sigmas * new_df$logn_sd)
      ),
      lowerbound_b = max(0.01, new_df$b_hat - (new_df$sigmas *
                                                 new_df$b_sd)),
      upperbound_b = min(0.8, new_df$b_hat + (new_df$sigmas *
                                                new_df$b_sd)),
      
      #-------------------------------------------
      
      #this is passed to the function
      logQc_hat = new_df$logQ_hat ,
      b_hat = rep(new_df$b_hat, nrow(bamdata$Wobs)),
      logA0_hat = rep(new_df$logA0_hat, nrow(bamdata$Wobs)),
      logn_hat = new_df$logn_hat,
      #-------------------------------
      
      #this is passed to the function
      logQ_sd = new_df$logQ_sd,
      b_sd = new_df$b_sd,
      logA0_sd = new_df$logA0_sd,
      logn_sd = new_df$logn_sd,
      Werr_sd = new_df$Werr_sd,
      Serr_sd = new_df$Serr_sd,
      dAerr_sd = new_df$dAerr_sd,
      
      #calc from other priors-------------------
      lowerbound_logQc = max(new_df$logQ_hat - (new_df$sigmas *
                                                  new_df$logQ_sd), -3),
      upperbound_logQc = new_df$logQ_hat + (new_df$sigmas *
                                              new_df$logQ_sd),
      logQc_sd = new_df$logQ_sd
    )
    
    return(priors_out)
  }#end function generate_bam_priors
  
  bam_val <- function(pred, qobs, variant) {
    bam_stats = data.frame()
    rrmse = sqrt(mean(((pred - qobs) / qobs) ^ 2))
    nash = 1 - (sum((qobs - pred) ^ 2) / sum((qobs - mean(qobs)) ^ 2))
    n_rmse = sqrt(mean(((pred - qobs) / mean(qobs)) ^ 2))
    p_bias = mean((pred - qobs) / mean(qobs))
    kge = KGE(sim = pred, obs = qobs)
    df = data.frame(
      nash = nash,
      rrmse = rrmse,
      n_rmse = n_rmse,
      p_bias = p_bias,
      kge = kge
    )
    bam_stats = rbind(df, bam_stats)
    return(list(variant, pred, bam_stats))
  }
  priors = generate_bam_priors(bamdata, new_df)
  variant = toString(new_df$variant)
  if (variant == "amhg") {
    bam_qpred = bam_validate(
      bam_estimate(
        bamdata = bamdata,
        variant = variant,
        bampriors = priors,
        chains = 3
      ),
      exp(bamdata$logQ_hat)
    )$valdata$qpred
  } else if (variant == "manning_amhg") {
    bam_qpred = bam_validate(
      bam_estimate(
        bamdata = bamdata,
        variant = variant,
        bampriors = priors,
        chains = 3
      ),
      exp(bamdata$logQ_hat)
    )$valdata$qpred
  } else{
    bam_qpred = bam_validate(
      bam_estimate(
        bamdata = bamdata,
        variant = variant,
        bampriors = priors,
        chains = 3
      ),
      exp(bamdata$logQ_hat)
    )$valdata$qpred
  }
  library(hydroGOF)
  bam_metrics = bam_val(bam_qpred, qval, variant = variant)
  detach("package:hydroGOF", unload = TRUE)
  return(bam_metrics)
}

#Test BAM sensitivity to number of pressure transducers (PTs), spacing, timing, and discharge (Q) constraints
test_num_pt <- function(bamdata) {
  pt_num <- lapply(as.list(2:nrow(bamdata$Wobs)), function(pt) {
    pt_num_list = list()
    if (choose(nrow(bamdata$Wobs), pt) > 10) {
      for (i in c(1:10)) {
        pt_id = sample(c(1:nrow(bamdata$Wobs)), pt)
        data = bam_data(
          w = matrix(bamdata$Wobs[pt_id, ], nrow = length(pt_id)),
          s = matrix(bamdata$Sobs[pt_id, ], nrow = length(pt_id)),
          dA = matrix(bamdata$dAobs[pt_id, ], nrow = length(pt_id)),
          Qhat = exp(bamdata$logQ_hat)
        )
        prior_df = generate_prior_df(bamdata = data,
                                     n_levels = 1,
                                     true = 1,var_scale = 1)
        amhg = run_and_validate_bam(
          bamdata = data,
          new_df = prior_df[1, ],
          n_levels = 1,
          qval = exp(bamdata$logQ_hat)
        )
        manning_amhg = run_and_validate_bam(
          bamdata = data,
          new_df = prior_df[2, ],
          n_levels = 1,
          qval = exp(bamdata$logQ_hat)
        )
        manning = run_and_validate_bam(
          bamdata = data,
          new_df = prior_df[3, ],
          n_levels = 1,
          qval = exp(bamdata$logQ_hat)
        )
        pt_num_list[[i]] = list(pt_id, amhg, manning_amhg, manning)
      }
    } else{
      for (i in 1:choose(nrow(bamdata$Wobs), pt)) {
        pt_id = combn(c(1:nrow(bamdata$Wobs)), pt)[, i]
        data = bam_data(
          w = matrix(bamdata$Wobs[pt_id, ], nrow = length(pt_id)),
          s = matrix(bamdata$Sobs[pt_id, ], nrow = length(pt_id)),
          dA = matrix(bamdata$dAobs[pt_id, ], nrow = length(pt_id)),
          Qhat = exp(bamdata$logQ_hat)
        )
        prior_df = generate_prior_df(bamdata = data,
                                     n_levels = 1,
                                     true = 1, var_scale = 1)
        amhg = run_and_validate_bam(
          bamdata = data,
          new_df = prior_df[1, ],
          n_levels = 1,
          qval = exp(bamdata$logQ_hat)
        )
        manning_amhg = run_and_validate_bam(
          bamdata = data,
          new_df = prior_df[2, ],
          n_levels = 1,
          qval = exp(bamdata$logQ_hat)
        )
        manning = run_and_validate_bam(
          bamdata = data,
          new_df = prior_df[3, ],
          n_levels = 1,
          qval = exp(bamdata$logQ_hat)
        )
        pt_num_list[[i]] = list(pt_id, amhg, manning_amhg, manning)
      }
    }
    return(pt_num_list)
  })
  return(pt_num)
}

test_spacing_pt <- function(bamdata) {
  pt_spacing <-
    lapply(as.list(1:choose(nrow(bamdata$Wobs), 2)), function(i) {
      
      pt_id = combn(c(1:nrow(bamdata$Wobs)), 2)[, i]
      data = bam_data(
        w = matrix(bamdata$Wobs[pt_id, ], nrow = length(pt_id)),
        s = matrix(bamdata$Sobs[pt_id, ], nrow = length(pt_id)),
        dA = matrix(bamdata$dAobs[pt_id, ], nrow = length(pt_id)),
        Qhat = exp(bamdata$logQ_hat)
      )
      prior_df = generate_prior_df(bamdata = data,
                                   n_levels = 1,
                                   true = 1, var_scale = 1)
      amhg = run_and_validate_bam(
        bamdata = data,
        new_df = prior_df[1, ],
        n_levels = 1,
        qval = exp(bamdata$logQ_hat)
      )
      manning_amhg = run_and_validate_bam(
        bamdata = data,
        new_df = prior_df[2, ],
        n_levels = 1,
        qval = exp(bamdata$logQ_hat)
      )
      manning = run_and_validate_bam(
        bamdata = data,
        new_df = prior_df[3, ],
        n_levels = 1,
        qval = exp(bamdata$logQ_hat)
      )
      return(list(pt_id, amhg, manning_amhg, manning))
    })
  return(pt_spacing)
}


test_timing <- function(bamdata_list) {
  timing <- lapply(as.list(1:length(bamdata_list)), function(i) {
    data = bamdata_list[[i]]
    prior_df = generate_prior_df(bamdata = data,
                                 n_levels = 1,
                                 true = 1, var_scale = 1)
    manning_amhg = run_and_validate_bam(
      bamdata = data,
      new_df = prior_df[2, ],
      n_levels = 1,
      qval = exp(data$logQ_hat)
    )
    return(list(manning_amhg))
  })
  return(timing)
}

test_Q_num <- function(bamdata){
  Q_num <- lapply(as.list(2:ncol(bamdata$Wobs)), function(Q){
    Q_num_list = list()
    if(choose(ncol(bamdata$Wobs), Q)>10){
      for(i in c(1:10)){
        time_id = sort(sample(c(1:ncol(bamdata$Wobs)), Q))
        Q_sam = approx(x = time_id,y = exp(bamdata$logQ_hat)[time_id], xout = c(1:ncol(bamdata$Wobs)), rule = 2)$y
        data = bam_data(w = bamdata$Wobs, 
                        s = bamdata$Sobs, 
                        dA = bamdata$dAobs,
                        Qhat = Q_sam)
        prior_df = generate_prior_df(bamdata = data, n_levels = 1, true = 1, var_scale = 1)
        manning = run_and_validate_bam(bamdata = data, new_df = prior_df[3,], n_levels = 1, qval = exp(bamdata$logQ_hat))
        Q_num_list[[i]] = list(time_id, manning)
        print(Q)
      }
    }else{
      for(i in c(1:choose(ncol(bamdata$Wobs), Q))){
        time_id = sort(combn(c(1:ncol(bamdata$Wobs)), Q)[,i])
        Q_sam = approx(x = time_id,y = exp(bamdata$logQ_hat)[time_id], xout = c(1:ncol(bamdata$Wobs)), rule = 2)$y
        data = bam_data(w = bamdata$Wobs, 
                        s = bamdata$Sobs, 
                        dA = bamdata$dAobs,
                        Qhat = Q_sam)
        
        prior_df = generate_prior_df(bamdata = data, n_levels = 1, true = 1, var_scale = 1)
        manning = run_and_validate_bam(bamdata = data, new_df = prior_df[3,], n_levels = 1, qval = exp(bamdata$logQ_hat))
        Q_num_list[[i]] = list(time_id, manning)
      }
    }
    return(Q_num_list)
  })
  return(Q_num)
}

Q_ids <- function(num_Q_list) {
  id = lapply(num_Q_list, function(num_Q) {
    lapply(num_Q, function(x) {
      return(x[[1]])
    })
  })
}

#Compare BAM with a temporary gauge
TG_all <- function(num_days, pt_count, WSE_mat, Qobs, Q_id) {
  rating_curve <-
    function(stage,
             discharge,
             full_stage,
             full_discharge) {
      dat = data.frame(stage = stage, discharge = discharge)
      full_dat = data.frame(stage = full_stage, discharge = full_discharge)
      model = lm(log(discharge) ~ log(stage), data = dat)
      pred = data.frame(p = exp(predict(model, full_dat)))
      qobs = full_discharge
      nash = 1 - (sum((qobs - pred) ^ 2) / sum((qobs - mean(qobs, na.rm = TRUE)) ^
                                                 2))
      return(list(pred, nash))
    }
  tg = lapply(Q_id, function(s) {
    lapply(s, function(id) {
      q_sample = id
      all = list()
      for (pt in 1:pt_count) {
        rc = rating_curve(WSE_mat[pt, q_sample], Qobs[q_sample],
                          WSE_mat[pt, ], Qobs)
        all[[pt]] = rc
      }
      return(all)
    })
  })
  return(tg)
}

#Calculate BAM input data for the PAD
PAD_bamdata <- function(WSE_all, WSE_Q, Q, bamdata, xvec, rivname){
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
  return(data)
}


