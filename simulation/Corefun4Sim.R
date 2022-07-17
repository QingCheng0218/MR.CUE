
findindex = function(block_inf, ind){
  indxAL = NULL;
  if(length(ind)==1){
    indxAL = block_inf[ind, 1]:block_inf[ind, 2];
  }else{
    for(i in 1:length(ind)){
      tmp = block_inf[ind[i], 1]:block_inf[ind[i], 2]
      indxAL = append(indxAL, tmp);
    }
  }
  
  return(indxAL)
}


ld_prune_cormat <- function(R, snp_names, p_vals,  p_val_thresh, r2_thresh){
  stopifnot(nrow(R) == length(snp_names))
  stopifnot(length(snp_names) == length(p_vals))
  
  ix <- which(p_vals < p_val_thresh)
  
  if(length(ix) == 0) return(c())
  
  snp_names <- snp_names[ix]
  R <- R[ix, ix, drop=FALSE]
  p_vals <- p_vals[ix]
  o_c <- order(p_vals, decreasing=FALSE)
  keep <- c()
  while(length(o_c) > 0){
    keep <- c(keep, snp_names[o_c[1]])
    myld <- R[o_c, o_c, drop=FALSE]
    remove_ix <- which(myld[,1]^2 > r2_thresh)
    o_c <- o_c[-remove_ix]
  }
  return(keep)
}

new_cause_data1 <- function(x = data.frame()){
  stopifnot(inherits(x, "data.frame"))
  # x <- validate_cause_data(x)
  #stopifnot(all(c("snp", "beta_hat_1", "seb1", "beta_hat_2", "seb2") %in% names(x)))
  structure(x, class = c("cause_data", "data.frame"))
}
cause_sims <- function(dat, param_ests, sigma_g, qalpha=1, qbeta=10,
                       no_ld = FALSE,
                       thresh = 1e-3){
  
  if(no_ld) dat <- process_dat_nold(dat)
  
  vars <- filter(dat, ld_prune == TRUE & p_value < thresh) %>% with(., snp)
  X <- new_cause_data(dat)
  
  #get sigma_g from data
  if(missing(sigma_g)) sigma_g <- cause:::eta_gamma_prior(X, vars)
  if(is.na(sigma_g)) sigma_g <- cause:::eta_gamma_prior(X, vars)
  
  res <- cause::cause(X=X, variants = vars, param_ests = param_ests, sigma_g = sigma_g,
                      qalpha = qalpha, qbeta = qbeta, force=TRUE)
  return(res)
}

#this is a helper function for switching from with ld to no ld data
process_dat_nold <- function(dat){
  dat <- dat %>%  select(-beta_hat_1, -beta_hat_2, -p_value, -ld_prune) %>%
    rename(beta_hat_1 = beta_hat_1_nold,
           beta_hat_2 = beta_hat_2_nold,
           p_value = p_value_nold) %>%
    mutate(ld_prune = TRUE)
  return(dat)
}


# ======= #
# GRAPPLE
# ======= #
my_GRAPPLE_fit <- function(gammah, Gammah, se1, se2, p.thres = NULL) {
  gammah = as.numeric(gammah)
  Gammah = as.numeric(Gammah)
  se1 = as.numeric(se1)
  se2 = as.numeric(se2)
  p_value = 2*pnorm(-abs(gammah/se1))
  
  data <- data.frame(
    SNP = 1:length(gammah),
    gamma_exp1 = gammah,
    gamma_out1 = Gammah,
    se_exp1 = se1,
    se_out1 = se2,
    selection_pvals = p_value
  )
  
  # p.thres = 1e-4
  # diagnosis <- GRAPPLE::findModes(data, p.thres = p.thres)
  # this method require to run findModes first, when only one mode is found, then we can run grappleRobustEst
  # elsewise we should find the confounding trait
  
  fit = GRAPPLE::grappleRobustEst(data, p.thres = p.thres, plot.it=F)
  
  beta.hat = fit$beta.hat
  beta.se = sqrt(fit$beta.var)
  beta.p.value = 2*pnorm(abs(beta.hat / beta.se), lower.tail=F)
  
  output = list(
    beta.hat = beta.hat,
    beta.se = beta.se,
    beta.p.value = beta.p.value
  )
  return(output)
}

# ======= #
# mrclust
# ======= #

my_mrclust_fit <- function(gammah, Gammah, se1, se2) {
  bx = as.numeric(gammah)
  by = as.numeric(Gammah)
  bxse = as.numeric(se1)
  byse = as.numeric(se2)
  ratio_est = by/bx
  ratio_est_se = byse/abs(bx)
  obs_names <- paste0("snp_", 1:length(bx))
  
  res_em = mrclust::mr_clust_em(
    theta = ratio_est, theta_se = ratio_est_se, bx = bx,
    by = by, bxse = bxse, byse = byse
  )
  
  bfit = res_em$results$best
  opt = as.numeric(names(which.max(table(bfit$cluster))))
  id = match(bfit$observation[bfit$cluster == opt], obs_names)
  mr.obj = mr_input(bx = bx[id], bxse = bxse[id], by = by[id], byse = byse[id])
  fit1 = mr_ivw(mr.obj)
  
  beta.hat = fit1$Estimate
  beta.se = fit1$StdError
  beta.p.value = fit1$Pvalue
  
  output = list(
    beta.hat = beta.hat,
    beta.se = beta.se,
    beta.p.value = beta.p.value
  )
  return(output)
}

# ======= #
# MRMix
# ======= #

my_MRMix_fit <- function(gammah, Gammah, se1, se2) {
  gammah = as.numeric(gammah)
  Gammah = as.numeric(Gammah)
  se1 = as.numeric(se1)
  se2 = as.numeric(se2)
  
  fit = MRMix::MRMix(gammah, Gammah, se1, se2, profile = TRUE)
  
  beta.hat = fit$theta
  beta.se = fit$SE_theta
  beta.p.value = fit$pvalue_theta
  
  output = list(
    beta.hat = beta.hat,
    beta.se = beta.se,
    beta.p.value = beta.p.value
  )
  return(output)
}