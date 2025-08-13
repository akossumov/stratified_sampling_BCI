SRSWOR <- function(population_data, n, seed){
  
  #population_data = BCI
  
  population_size <- nrow(population_data)
  
  ###### Sampling of grid cells
  #set.seed(271199)
  set.seed(seed)
  
  units_srswor <- sample(population_size, size = n, replace = FALSE)
  mysample_srswor <- population_data[units_srswor, ]
  
  
  ###### Estimation of population mean by simple random sampling without replacement
  
  mz_srswor <- mean(mysample_srswor$Ca)
  
  est_sigma2z_srswor <- var(mysample_srswor$Ca)
  
  fpc_srswor <- (population_size - n)/population_size
  est_VAR_mz_srswor <- fpc_srswor*est_sigma2z_srswor/n
  
  sigma2z <- var(population_data$Ca)
  VAR_mz_srswor <- fpc_srswor*sigma2z/n
  
  # Alternative computation of mz_srswor and sqrt(est_VAR_mz_srswor), confirming that the manual computation above is correct
  
  design_si <- svydesign(id = ~ 1, probs = ~ pi_srswor, fpc = ~ fpc, data = transform(mysample_srswor, pi_srswor = n / population_size, fpc = population_size))
  res_srswor <- svymean(~ Ca, design_si)  # is the same to print(c(mz_srswor, sqrt(est_VAR_mz_srswor)))
  
  # control manual computation with automatic computation
  if (check_not_equal(as.numeric(mean(res_srswor)),  mz_srswor) || check_not_equal(as.numeric(SE(res_srswor)), sqrt(est_VAR_mz_srswor))){
    stop("In simple random sampling without replacement, the manual and automatic estimates of the mean/std. error are not equal!")
  }
  
  ### Confidence interval for the estimate of the population mean using simple random sampling without replacement
  conf_int_srswor <- confint(svymean(~ Ca, design_si), df = degf(design_si), level = 0.95)
  
  return(list(
    
    sigma2z = sigma2z,
    
    #sample_srswor = mysample_srswor,
    
    conf_int_srswor = conf_int_srswor,
    
    fpc_srswor = fpc_srswor,
    
    mz_srswor = mz_srswor,
    
    est_VAR_mz_srswor = est_VAR_mz_srswor,
    
    VAR_mz_srswor = VAR_mz_srswor, 
    
    est_SE_mz_srswor = sqrt(est_VAR_mz_srswor),
    
    SE_mz_srswor = sqrt(VAR_mz_srswor)
  ))
}



PROP_ALLOC <- function(population_data, vec_strata, n, seed){
  
  #population_data = BCI
  #vec_strata =  crfFe_ID
  
  population_data$stratum <- vec_strata
  population_size <- nrow(population_data)
  Nh <- tapply(population_data$stratum, INDEX = population_data$stratum, FUN = length)
  wh <- Nh / sum(Nh)
  sigma2z_h <- tapply(population_data$Ca, INDEX = population_data$stratum, FUN = var)
  
  ###### Proportional allocation
  nh_prop <- round_preserve_sum(n * wh, n)
  
  # control sum(nh_prop) = n
  if (sum(nh_prop) != n) {
    stop("The rounded prop allocation do not sum up to n!")
  }
  
  # we need to ensure that all allocated sample sizes are greater than 2 in order to estimate variances
  if (!all(nh_prop >= 2)){
    stop("Some of the allocated sample sizes contain fewer than 2 elements!")
  }
  
  ###### Sampling of grid cells
  ord <- unique(population_data$stratum)
  #set.seed(271199)
  set.seed(seed)
  
  units_prop <- sampling::strata(
    population_data, stratanames = "stratum", size = nh_prop[ord], method = "srswor") # we select units by simple random sampling WITHOUT replacement
  mysample_prop <- getdata(population_data, units_prop) %>%
    select(-ID_unit, -Stratum) 
  
  ###### Estimation of population mean by stratified simple random sampling (without replacement) with proportional allocation
  
  mz_h_prop <- tapply(mysample_prop$Ca, INDEX = mysample_prop$stratum, FUN = mean)
  mz_prop <- sum(wh * mz_h_prop)
  
  est_sigma2z_h_prop <- tapply(mysample_prop$Ca, INDEX = mysample_prop$stratum, FUN = var)
  
  fpc_prop <- (Nh - nh_prop)/Nh
  est_VAR_mz_prop <- sum(wh^2 * fpc_prop * est_sigma2z_h_prop / nh_prop)
  VAR_mz_prop <- sum(wh^2 * fpc_prop * sigma2z_h / nh_prop)
  
  # Alternative computation of mz_prop and sqrt(est_VAR_mz_prop), confirming that the manual computation above is correct
  
  labels_prop <- sort(unique(mysample_prop$stratum))
  lut_prop <- data.frame(stratum = labels_prop, weight = Nh/nh_prop, fpc = Nh)
  design_stsi_prop <- svydesign(id = ~ 1, strata = ~ stratum, weight = ~ weight, fpc = ~fpc, data = merge(x = mysample_prop, y = lut_prop))
  res_prop <- svymean(~ Ca, design_stsi_prop, deff = TRUE) # is the same to print(c(mz_prop, sqrt(est_VAR_mz_prop)))
  
  # control manual computation with automatic computation
  if (check_not_equal(as.numeric(mean(res_prop)),  mz_prop) || check_not_equal(as.numeric(SE(res_prop)), sqrt(est_VAR_mz_prop))){
    stop("In stratified simple random sampling with proportional allocation, the manual and automatic estimates of the mean/std. error are not equal!")
  }
  
  ### Confidence interval for the estimate of the population mean using stratified simple random sampling (without replacement) with proportional allocation
  conf_int_prop <- confint(svymean(~ Ca, design_stsi_prop), df = degf(design_stsi_prop), level = 0.95)
  
  ### Design effect of stratified random sampling (without replacement) with proportional allocation 
  
  fpc_srswor <- (population_size - n)/population_size
  sigma2z <- var(population_data$Ca)
  VAR_mz_srswor <- fpc_srswor*sigma2z/n
  
  # We can compute the TRUE value of the design effect as
  true_desef_prop <- VAR_mz_prop/VAR_mz_srswor
  
  # "svymean" uses the following estimator of the design effect, where the population variance is estimated using some weighting method
  est_desef_prop <- unname(est_VAR_mz_prop/(svyvar(~ Ca, design_stsi_prop)[1]*fpc_srswor/n))
  
  # control manual computation with automatic computation
  if (check_not_equal(attr(res_prop, "deff")[1], est_desef_prop)){
    stop("In stratified simple random sampling (without replacement) with proportional allocation, the manual and automatic estimates of design effect are not equal!")
  }
  
  # weighted estimation of the population variance, svyvar(~ Ca, design_stsi_prop)[1], can be also computed as 
  s2(mysample_prop$Ca, w = 1/mysample_prop$Prob)
  svyvar(~ Ca, design_stsi_prop)[1]
  # which is different from 
  var(mysample_prop$Ca)
  
  
  return(list(
    Nh = Nh,
    wh = wh,
    
    #sample_prop = mysample_prop,
    
    conf_int_prop = conf_int_prop,
    
    nh_prop = nh_prop,
    
    mz_h_prop = mz_h_prop,
    
    true_desef_prop = true_desef_prop,
    est_desef_prop = est_desef_prop,
    
    #est_sigma2z_h_prop = est_sigma2z_h_prop,
    
    fpc_prop = fpc_prop,
    
    mz_prop = mz_prop,
    
    est_VAR_mz_prop = est_VAR_mz_prop,
    
    VAR_mz_prop = VAR_mz_prop,
    
    est_SE_mz_prop = sqrt(est_VAR_mz_prop),
    
    SE_mz_prop = sqrt(VAR_mz_prop)
  ))
}



RUN_PARALLEL_SIMULATIONS <- function(population_data, SEQ_n, H, SIM_SEQ){
  
  crfFe_strata <- strata.cumrootf(
    x = population_data$Fe, n = 100, Ls = H, nclass = 500)
  crfFe_ID <- as.factor(crfFe_strata$stratumID)
  
  population_data_geo <- population_data
  gridded(population_data_geo) <- ~ x + y # here we are defining the grid structure based on the coordinates x and y
  set.seed(271199) # We should set a seed because each time the k-means algorithm runs, it starts with new random starting points
  geo_strata <- stratify(
    object = population_data_geo, nStrata = H, nTry = 100, equalArea = FALSE)
  geo_ID <- as.factor(geo_strata@stratumId)
  
  library(foreach)
  library(doParallel)
  
  # Set up parallel backend
  num_cores <- parallel::detectCores() - 1
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  n_srswor_list <- list()
  n_geo_list <- list()
  n_crfFe_list <- list()
  
  # Outer loop in parallel
  results_list <- foreach(i = seq_along(SEQ_n), .packages = c(
    "sswr", "sampling", "dplyr", "stratification", "ggplot2", "viridis",
    "spcosa", "sp", "survey", "surveyplanning", "patchwork"
  )) %dopar% {
    # Each worker must re-source necessary functions
    source("final_project_functions.R")
    
    n <- SEQ_n[i]
    
    DICT_srswor_prop <- list(mz = rep(NA, length(SIM_SEQ)), SE = c(NA), est_SE = rep(NA, length(SIM_SEQ)))
    DICT_geo_prop <- list(mz = rep(NA, length(SIM_SEQ)), SE = c(NA), est_SE = rep(NA, length(SIM_SEQ)), desef = c(NA), est_desef = rep(NA, length(SIM_SEQ)))
    DICT_crfFe_prop <- list(mz = rep(NA, length(SIM_SEQ)), SE = c(NA), est_SE = rep(NA, length(SIM_SEQ)), desef = c(NA), est_desef = rep(NA, length(SIM_SEQ)))
    
    for (count in SIM_SEQ){
      
      srswor_res <- SRSWOR(population_data, n, count)
      
      DICT_srswor_prop$mz[count] <- srswor_res$mz_srswor
      DICT_srswor_prop$est_SE[count] <- srswor_res$est_SE_mz_srswor
      
      geo_PROP_res <- PROP_ALLOC(population_data, geo_ID, n, count)
      
      DICT_geo_prop$mz[count] <- geo_PROP_res$mz_prop
      DICT_geo_prop$est_SE[count] <- geo_PROP_res$est_SE_mz_prop
      DICT_geo_prop$est_desef[count] <- geo_PROP_res$est_desef_prop
      
      crfFe_PROP_res <- PROP_ALLOC(population_data, crfFe_ID, n, count)
      
      DICT_crfFe_prop$mz[count] <- crfFe_PROP_res$mz_prop
      DICT_crfFe_prop$est_SE[count] <- crfFe_PROP_res$est_SE_mz_prop
      DICT_crfFe_prop$est_desef[count] <- crfFe_PROP_res$est_desef_prop
      
      
      if (count == 1){
        DICT_srswor_prop$SE[count] <- srswor_res$SE_mz_srswor
        
        DICT_geo_prop$SE[count] <- geo_PROP_res$SE_mz_prop
        DICT_geo_prop$desef[count] <- geo_PROP_res$true_desef_prop
        
        DICT_crfFe_prop$SE[count] <- crfFe_PROP_res$SE_mz_prop
        DICT_crfFe_prop$desef[count] <- crfFe_PROP_res$true_desef_prop
      }
    }
    
    list(
      srswor = DICT_srswor_prop,
      geo = DICT_geo_prop,
      crfFe = DICT_crfFe_prop
    )
  }  
  
  stopCluster(cl)
  
  # Unpack results
  for (i in seq_along(results_list)) {
    n_srswor_list[[i]] <- results_list[[i]]$srswor
    n_geo_list[[i]] <- results_list[[i]]$geo
    n_crfFe_list[[i]] <- results_list[[i]]$crfFe
  }
  
  return(list(
    srswor = n_srswor_list,
    geo = n_geo_list,
    crfFe = n_crfFe_list
  ))
}



round_preserve_sum <- function(x, total) {
  rx <- round(x)
  rdiff <- total - sum(rx)
  if (rdiff > 0) {
    # Add 1 to the smallest elements
    indices <- order(rx)[1:rdiff]
    rx[indices] <- rx[indices] + 1
  } else if (rdiff < 0) {
    # Subtract 1 from the largest elements
    indices <- order(rx, decreasing = TRUE)[1:abs(rdiff)]
    rx[indices] <- rx[indices] - 1
  }
  return(rx)
}



check_not_equal <- function(a, b) {
  return(!isTRUE(all.equal(a, b)))
}





