rm(list=ls());
library("mvtnorm");
library("dplyr")
library(cause);
library(dplyr);
library(MR.LDP);
library(mr.raps);
library(gsmr);
library(MR.Corr2);
library(MR.CUE);
library(MendelianRandomization);
library(MRcML);


# ------------------------------------------------------------------ #
Garate <- 0.05;
Thrate <- 1;
bratio <- 4;
rho <- 0.4;
L <- 200;
h2x <- 0.3;
h2y <- 0.25;

filename <- paste0("Garate", Garate,  "Thrate", Thrate, "bratio", bratio,
                   "rho", rho, "L", L,  "h2x", h2x,
                   "h2y", h2y, "LD4cMLMAJuly20.Rdata");
M = 10; p = M*L; id4ld <- seq(1, p, 5);
n1 = 50000; n2 = 50000;  n3 = 4000;

maf = runif(p, 0.05, 0.5);
G = genRawGeno(maf, L, M, rho, n1 + n2 + n3);
G1 = G[1:n1,];G2 = G[(n1+1):(n1+n2),]; G12 = G[1:(n1+n2),];
G3 = G[(n1+n2+1):(n1+n2+n3),];

block_inf <- cbind(seq(1, p, M), seq(M, p, M));
block_inf1 <- block_inf - 1;
lambda = 0.85
R = Cal_block_SimR(block_inf1, G3, lambda)


sigma2g = 1;
sigma2t = 1;
coreNum = 20;


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
# ------------------------------------------------------------------ #
nrep = 100;
bhat0 = pvalue0 = matrix(NA, nrow = nrep, ncol = 10); # for reverse causal.
bhat1 = pvalue1 = matrix(NA, nrow = nrep, ncol = 10); # for causal.

for(irep in 1:nrep){
  # ---------------------------------------------------------------- #
  # for the exposure x.
  # sparse setting for gamma.
  gano = floor(L*Garate);
  indga = sample(1:L,gano);
  if(M==1){
    idga = indga;
  }else{
    idga = findindex(block_inf, indga);
  }
  
  gamma = numeric(p);
  
  gamma.nz = rnorm(length(idga))*sqrt(sigma2g);
  gamma[idga] = gamma.nz;
  
  G12g = G12%*%gamma;
  xall = G12g + (as.vector(var(G12g)*(1-h2x)/h2x))^0.5*rnorm(n1+n2);
  # check the heritability
  var(G12g)/var(xall);
  # ---------------------------------------------------------------- #
  # for the outcome y.
  theta = numeric(p);
  idth = 1:p;
  theta.nz = rnorm(length(idth))*sqrt(sigma2t);
  theta[idth] = theta.nz;
  G12t = G12%*%theta;
  
  # set beta0 to contorl var(beta0*G12g)/var(G12t) = 1:4, bratio = 4;
  b0 = as.numeric(sqrt(var(G12t)/(bratio*var(G12g))));
  b0*b0*var(G12g)/var(G12t)
  y0 = b0*G12g + G12t;
  yall = y0 + (as.vector(var(y0)*(1-h2y)/h2y))^0.5*rnorm(n1+n2);
  # check
  var(y0)/var(yall)
  var(b0*G12g)/var(yall)
  var(G12t)/var(yall)
  # ---------------------------------------------------------------- #
  x = xall[1:n1];
  y = yall[(n1+1):(n1+n2)];
  # ------------------------------------
  # create summary statistics
  gammaSS = fastSigLm(x, G1);
  GammaSS = fastSigLm(y, G2);
  gammah = gammaSS$coef;
  se1 = gammaSS$std;
  Gammah = GammaSS$coef;
  se2 = GammaSS$std;
  
  # rename summary statistics due to reverse causation.
  Gammah_r = gammah;
  gammah_r = Gammah;
  se1_r = se2;
  se2_r = se1;
  
  #*******************************************************************************#
  # For Correlated SNPs
  #*******************************************************************************#
  #####----   MR-CUE  ----#####
  # --------------------------#
  #-------------------------------------------------------#
  # for causal.
  resM31 = MRCUESim(gammah, Gammah, se1, se2, 0,  R, block_inf1,  coreNum);

  M3beta1 = resM31$beta.hat;
  M3pvalue1 = resM31$beta.p.value;
  bhat1[irep, 1] = M3beta1;
  pvalue1[irep, 1] = M3pvalue1;
  rm(resM31,  M3beta1, M3pvalue1);
  #-------------------------------------------------------#
  # for reverse causal.
  resM30 = MRCUESim(gammah_r, Gammah_r, se1_r, se2_r, 0,  R, block_inf1, coreNum);

  M3beta0 = resM30$beta.hat;
  M3pvalue0 = resM30$beta.p.value;
  bhat0[irep, 1] = M3beta0;
  pvalue0[irep, 1] = M3pvalue0;
  rm(resM30, M3beta0, M3pvalue0);
  #-------------------------------------------------------#
  # --------------------------#
  #####----   MR-LDP  ----#####
  # --------------------------#
  # initial value.
  epsStopLogLik <- 1e-7; maxIter <- 10000; beta0 <- 0;
  gamma <- rep(0, p); alpha <- rep(0, p); sgga2 <- 0.01; sgal2 <- 0.01;
  #-------------------------------------------------------#
  # for causal.
  SimMRLDP_Hb = MRLDP_SimPXvb(gammah, Gammah, se1, se2, gamma, alpha,
                              beta0, sgga2, sgal2, R,
                              0, epsStopLogLik, maxIter, model = 2);

  beta_hat2 = SimMRLDP_Hb$beta0;
  SimMRLDP_H0 = MRLDP_SimPXvb(gammah, Gammah, se1, se2, gamma, alpha,
                              beta0, sgga2, sgal2, R,
                              1, epsStopLogLik, maxIter, model = 2);
  tstat = 2*(SimMRLDP_Hb$tstat - SimMRLDP_H0$tstat);
  pval2 = pchisq(tstat, 1, lower.tail = F);

  bhat1[irep, 2] = beta_hat2;
  pvalue1[irep, 2] = pval2;
  rm(SimMRLDP_Hb, SimMRLDP_H0, beta_hat2, pval2);
  #-------------------------------------------------------#
  # for reverse causal.
  SimMRLDP_Hbr = MRLDP_SimPXvb(gammah_r, Gammah_r, se1_r, se2_r, gamma,
                               alpha,  beta0, sgga2, sgal2, R,
                               0, epsStopLogLik, maxIter, model = 2);

  beta_hat2r = SimMRLDP_Hbr$beta0;
  SimMRLDP_H0r = MRLDP_SimPXvb(gammah_r, Gammah_r, se1_r, se2_r, gamma,
                               alpha,  beta0, sgga2, sgal2, R,
                               1, epsStopLogLik, maxIter, model = 2);
  tstatr = 2*(SimMRLDP_Hbr$tstat - SimMRLDP_H0r$tstat);
  pval2r = pchisq(tstatr, 1, lower.tail = F);

  bhat0[irep, 2] = beta_hat2r;
  pvalue0[irep, 2] = pval2r;
  rm(SimMRLDP_Hbr, SimMRLDP_H0r, beta_hat2r, pval2r);
  #------------------------------------------------------#
  # --------------------------#
  ####---- MR-Corr2  ----#####
  # --------------------------#
  opt = list(agm = 0.001, bgm = 0.001, aal = 0.001, bal = 0.001,
             a = 1, b = L, maxIter = 4000, thin = 10, burnin = 1000);

  #-------------------------------------------------------#
  # for causal.
  SimRes = MRCorr2Sim(gammah, Gammah, se1, se2, R, block_inf1, coreNum, opt);
  beta_hat3 = mean(SimRes$Beta0res);
  beta_se3 = sd(SimRes$Beta0res);
  beta_pvalue3 = 2*pnorm(abs(beta_hat3/beta_se3), lower.tail=FALSE);
  bhat1[irep, 3] = beta_hat3;
  pvalue1[irep, 3] = beta_pvalue3;
  rm(SimRes, beta_hat3, beta_pvalue3);
  #-------------------------------------------------------#
  # for reverse causal.
  SimResr = MRCorr2Sim(gammah_r, Gammah_r, se1_r, se2_r, R, block_inf1, coreNum, opt);
  beta_hat3r = mean(SimResr$Beta0res);
  beta_se3r = sd(SimResr$Beta0res);
  beta_pvalue3r = 2*(1 - pnorm(abs(beta_hat3r/beta_se3r)));
  bhat0[irep, 3] = beta_hat3r;
  pvalue0[irep, 3] =beta_pvalue3r;
  rm(SimResr, beta_hat3r, beta_pvalue3r);
  #*******************************************************************************#
  # For Independent SNPs
  #*******************************************************************************#
  # --------------------------#
  #####------ CAUSE  -----#####
  # --------------------------#
  #-------------------------------------------------------#
  # for causal.
  snp = 1:p;
  ld_prune_pval_thresh = 1;
  r2_thresh = 1;
  ld_prune = rep(5==1, p);
  ld_prune[id4ld]= TRUE;

  p_value = 2*pnorm(-abs(gammah/se1));
  df <- data.frame(snp = snp,
                   beta_hat_1 = as.vector(gammah),
                   beta_hat_2 = as.vector(Gammah),
                   seb1 = as.vector(se1),
                   seb2 = as.vector(se2),
                   ld_prune = ld_prune,
                   p_value = p_value);

  X1 = data.frame(snp = snp,
                  beta_hat_1 = as.vector(gammah),
                  beta_hat_2 = as.vector(Gammah),
                  seb1 = as.vector(se1),
                  seb2 = as.vector(se2))

  X1 = new_cause_data1(X1);
  param = est_cause_params(X1, X1$snp, null_wt = 10, max_candidates = Inf);
  cause_res <- try(cause_sims(df, param, no_ld = FALSE));
  if(class(cause_res) != 'try-error'){
    qs = summary(cause_res)$quants;
    z <- -1*summary(cause_res)$z;
    bhat1[irep, 4] = qs[[2]][1, 1];
    pvalue1[irep, 4] = pnorm(-z);
    rm(z);
    rm(qs);
    rm(cause_res)
  }

  #-------------------------------------------------------#
  # for reverse causal.
  p_valuer = 2*pnorm(-abs(gammah_r/se1_r));

  dfr <- data.frame(snp = snp,
                    beta_hat_1 = as.vector(gammah_r),
                    beta_hat_2 = as.vector(Gammah_r),
                    seb1 = as.vector(se1_r),
                    seb2 = as.vector(se2_r),
                    ld_prune = ld_prune,
                    p_value = p_valuer);

  X1 = data.frame(snp = snp,
                  beta_hat_1 = as.vector(gammah_r),
                  beta_hat_2 = as.vector(Gammah_r),
                  seb1 = as.vector(se1_r),
                  seb2 = as.vector(se2_r))

  X1 = new_cause_data1(X1);
  param = est_cause_params(X1, X1$snp, null_wt = 10, max_candidates = Inf);
  cause_res <- try(cause_sims(dfr, param, no_ld = FALSE));
  if(class(cause_res) != 'try-error'){
    qs = summary(cause_res)$quants;
    z <- -1*summary(cause_res)$z;
    bhat0[irep, 4] = qs[[2]][1, 1];
    pvalue0[irep, 4] = pnorm(-z);
    rm(qs, z);

  }

  #-------------------------------------------------------#
  # --------------------------#
  #####------ GSMR  ------#####
  # --------------------------#
  cat("start GSMR procedure:", "\n")
  p_val_thresh  = 1
  # id4gsmr = which(p_value < p_val_thresh & ld_prune=="TRUE");
  id4gsmr = id4ld;
  Rgsmr = R[id4gsmr, id4gsmr];
  Rgsmr <- data.frame(Rgsmr)
  names(Rgsmr) <- seq_along(id4gsmr)
  rownames(Rgsmr) <- seq_along(id4gsmr);

  #-------------------------------------------------------#
  # for causal.
  dat <- df[id4gsmr,] %>%
    mutate(p_value2 = 2*pnorm(-abs(beta_hat_2/seb2)))

  resgsmr1 <-  try(with(dat, gsmr(beta_hat_1, seb1, p_value,
                                  beta_hat_2, seb2, bzy_pval = p_value2,
                                  ldrho=Rgsmr, snpid=seq_along(id4gsmr),
                                  n_ref = 1, nsnps_thresh=1, gwas_thresh=p_val_thresh)))


  if(class(resgsmr1)!="try-error") {
    bhat1[irep, 5] <- resgsmr1$bxy;
    pvalue1[irep, 5] <- resgsmr1$bxy_pval;
  }

  #-------------------------------------------------------#
  # for reverse causal.
  datr <- dfr[id4gsmr,] %>%
    mutate(p_value2 = 2*pnorm(-abs(beta_hat_2/seb2)))

  resgsmr0 <-  try(with(datr, gsmr(beta_hat_1, seb1, p_value,
                                   beta_hat_2, seb2, bzy_pval = p_value2,
                                   ldrho=Rgsmr, snpid=seq_along(id4gsmr),
                                   n_ref = 1, nsnps_thresh=1, gwas_thresh=p_val_thresh)))


  if(class(resgsmr0)!="try-error") {
    bhat0[irep, 5] <- resgsmr0$bxy;
    pvalue0[irep, 5] <- resgsmr0$bxy_pval;
  }
  # --------------------------#
  #####------ Raps  ------#####
  # --------------------------#

  #-------------------------------------------------------#
  # for causal.
  raps = mr.raps(gammah[id4ld], Gammah[id4ld], se1[id4ld], se2[id4ld], over.dispersion = TRUE);
  bhat1[irep, 6] = raps$beta.hat;
  pvalue1[irep, 6] = raps$beta.p.value;
  #-------------------------------------------------------#
  # for reverse causal.
  rapsr = mr.raps(gammah_r[id4ld], Gammah_r[id4ld], se1_r[id4ld], se2_r[id4ld], over.dispersion = TRUE);
  bhat0[irep, 6] = rapsr$beta.hat;
  pvalue0[irep, 6] = rapsr$beta.p.value;

  
  
  #-------------------------------------------------------#
  #####------------ MR-cML(DP) for causal ------------#####
  #-------------------------------------------------------#
  cML_result_DP = mr_cML_DP(gammah[id4ld], Gammah[id4ld], se1[id4ld], se2[id4ld], n = n1);
  bhat1[irep, 7] = cML_result_DP$MA_BIC_theta;
  pvalue1[irep, 7] = cML_result_DP$MA_BIC_p;
  bhat1[irep, 8] = cML_result_DP$MA_BIC_DP_theta;
  pvalue1[irep, 8] = cML_result_DP$MA_BIC_DP_p;
  
  rm(cML_result_DP);
  #-------------------------------------------------------#
  #####-------- MR-cML(DP) for reverse causal --------#####
  #-------------------------------------------------------#
  cML_result_DP_r = mr_cML_DP(gammah_r[id4ld], Gammah_r[id4ld], se1_r[id4ld], se2_r[id4ld], n = n1);
  bhat0[irep, 7] = cML_result_DP_r$MA_BIC_theta;
  pvalue0[irep, 7] = cML_result_DP_r$MA_BIC_p;
  bhat0[irep, 8] = cML_result_DP_r$MA_BIC_DP_theta;
  pvalue0[irep, 8] = cML_result_DP_r$MA_BIC_DP_p;
  rm(cML_result_DP_r);
  # --------------------------#
  #####------- IVW  ------#####
  # --------------------------#
  mr.obj = mr_input(bx = as.vector(gammah[id4ld]), bxse = as.vector(se1[id4ld]),
                    by = as.vector(Gammah[id4ld]), byse = as.vector(se2[id4ld]));

  mr.objr = mr_input(bx = as.vector(gammah_r[id4ld]), bxse = as.vector(se1_r[id4ld]),
                     by = as.vector(Gammah_r[id4ld]), byse = as.vector(se2_r[id4ld]));

  IVW = mr_ivw(mr.obj);
  IVWr = mr_ivw(mr.objr);

  bhat1[irep, 9] = IVW$Estimate;
  pvalue1[irep, 9] = IVW$Pvalue;
  bhat0[irep, 9] = IVWr$Estimate;
  pvalue0[irep, 9] = IVWr$Pvalue;
  # --------------------------#
  #####---- MR-Egger----#####
  # --------------------------#
  Egger = try(mr_egger(mr.obj));
  Eggerr = try(mr_egger(mr.objr));

  if(class(Egger) != 'try-error'){
    bhat1[irep, 10] = Egger$Estimate;
    pvalue1[irep, 10] = Egger$Pvalue.Est;
  }
  if(class(Eggerr) != 'try-error'){
    bhat0[irep, 10] = Eggerr$Estimate;
    pvalue0[irep, 10] = Eggerr$Pvalue.Est;
  }
  
  if(irep%%10==0) {
    print(irep);
    print(Sys.time());
    save(pvalue0, bhat0, pvalue1, bhat1, file = filename);
  }
  
}