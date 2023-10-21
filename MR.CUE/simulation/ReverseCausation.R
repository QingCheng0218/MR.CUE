rm(list=ls());
source("Corefun4Sim.R")
library("mvtnorm");
library("dplyr")
library(cause);
library(dplyr);
library(MR.LDP);
library(mr.raps);
library(MR.CUE);
library(MRcML);
library(GRAPPLE);
library(mrclust);
library(MRMix);
library(MendelianRandomization);
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
  # ========================= #
  # MR-CUE 
  # ========================= #
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
  # ========================= #
  # MR-LDP
  # ========================= #
  # initial value.
  epsStopLogLik <- 1e-7; maxIter <- 10000; beta0 <- 0;
  gamma <- rep(0, p); alpha <- rep(0, p); sgga2 <- 0.01; sgal2 <- 0.01;
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

  #*******************************************************************************#
  # For Independent SNPs
  #*******************************************************************************#
  # ========================= #
  # CAUSE
  # ========================= #
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
    bhat1[irep, 3] = qs[[2]][1, 1];
    pvalue1[irep, 3] = pnorm(-z);
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
    bhat0[irep, 3] = qs[[2]][1, 1];
    pvalue0[irep, 3] = pnorm(-z);
    rm(qs, z);

  }

  # ========================= #
  # Raps
  # ========================= #
  #-------------------------------------------------------#
  # for causal.
  raps = mr.raps(gammah[id4ld], Gammah[id4ld], se1[id4ld], se2[id4ld], over.dispersion = TRUE);
  bhat1[irep, 4] = raps$beta.hat;
  pvalue1[irep, 4] = raps$beta.p.value;
  #-------------------------------------------------------#
  # for reverse causal.
  rapsr = mr.raps(gammah_r[id4ld], Gammah_r[id4ld], se1_r[id4ld], se2_r[id4ld], over.dispersion = TRUE);
  bhat0[irep, 4] = rapsr$beta.hat;
  pvalue0[irep, 4] = rapsr$beta.p.value;

  # ========================= #
  # MR-cML(DP)
  # ========================= #
  #-------------------------------------------------------#
  # for causal.
  cML_result_DP = mr_cML_DP(gammah[id4ld], Gammah[id4ld], se1[id4ld], se2[id4ld], n = n1);
  bhat1[irep, 5] = cML_result_DP$MA_BIC_DP_theta;
  pvalue1[irep, 5] = cML_result_DP$MA_BIC_DP_p;
  
  rm(cML_result_DP);
  #-------------------------------------------------------#
  # for reverse causal.
  cML_result_DP_r = mr_cML_DP(gammah_r[id4ld], Gammah_r[id4ld], se1_r[id4ld], se2_r[id4ld], n = n1);
  bhat0[irep, 5] = cML_result_DP_r$MA_BIC_DP_theta;
  pvalue0[irep, 5] = cML_result_DP_r$MA_BIC_DP_p;
  rm(cML_result_DP_r);
  # ========================= #
  # IVW 
  # ========================= #
  #-------------------------------------------------------#
  # for causal.
  mr.obj = mr_input(bx = as.vector(gammah[id4ld]), bxse = as.vector(se1[id4ld]),
                    by = as.vector(Gammah[id4ld]), byse = as.vector(se2[id4ld]));

  #-------------------------------------------------------#
  # for reverse causal.
  mr.objr = mr_input(bx = as.vector(gammah_r[id4ld]), bxse = as.vector(se1_r[id4ld]),
                     by = as.vector(Gammah_r[id4ld]), byse = as.vector(se2_r[id4ld]));

  IVW = mr_ivw(mr.obj);
  IVWr = mr_ivw(mr.objr);

  bhat1[irep, 6] = IVW$Estimate;
  pvalue1[irep, 6] = IVW$Pvalue;
  bhat0[irep, 6] = IVWr$Estimate;
  pvalue0[irep, 6] = IVWr$Pvalue;
  # ========================= #
  # MR-Egger
  # ========================= #
  #-------------------------------------------------------#
  # for causal.
  Egger = try(mr_egger(mr.obj));
  #-------------------------------------------------------#
  # for reverse causal.
  Eggerr = try(mr_egger(mr.objr));

  if(class(Egger) != 'try-error'){
    bhat1[irep, 7] = Egger$Estimate;
    pvalue1[irep, 7] = Egger$Pvalue.Est;
  }
  if(class(Eggerr) != 'try-error'){
    bhat0[irep, 7] = Eggerr$Estimate;
    pvalue0[irep, 7] = Eggerr$Pvalue.Est;
  }
  
  
  # =========================#
  # GRAPPLE
  # =========================#
  #-------------------------------------------------------#
  # for causal.
  GRAPPLE_fit <- my_GRAPPLE_fit(gammah[id4ld], Gammah[id4ld], se1[id4ld], se2[id4ld])
  bhat1[irep,   8] = GRAPPLE_fit$beta.hat;
  pvalue1[irep, 8] = GRAPPLE_fit$beta.p.value;
  
  rm(GRAPPLE_fit);
  
  
  GRAPPLE_fitr <- my_GRAPPLE_fit(gammah_r[id4ld], Gammah_r[id4ld], se1_r[id4ld], se2_r[id4ld])
  bhat0[irep,   8] = GRAPPLE_fitr$beta.hat;
  pvalue0[irep, 8] = GRAPPLE_fitr$beta.p.value;
  
  # ========================= #
  # MR-Clust
  # ========================= #
  #-------------------------------------------------------#
  # for causal.
  mrclust_fit <- my_mrclust_fit(gammah[id4ld], Gammah[id4ld], se1[id4ld], se2[id4ld])
  bhat1[irep,   9] = mrclust_fit$beta.hat;
  pvalue1[irep, 9] = mrclust_fit$beta.p.value;
  rm(mrclust_fit);
  
  mrclust_fitr <- my_mrclust_fit(gammah_r[id4ld], Gammah_r[id4ld], se1_r[id4ld], se2_r[id4ld])
  bhat0[irep,   9] = mrclust_fitr$beta.hat;
  pvalue0[irep, 9] = mrclust_fitr$beta.p.value;
  
  # ========================= #
  # MRMix
  # ========================= #
  #-------------------------------------------------------#
  # for causal.
  MRMix_fit <- my_MRMix_fit(gammah[id4ld], Gammah[id4ld], se1[id4ld], se2[id4ld])
  bhat1[irep,   10] = MRMix_fit$beta.hat;
  pvalue1[irep, 10] = MRMix_fit$beta.p.value;
  rm(MRMix_fit);
  
  #-------------------------------------------------------#
  # for reverse causal.
  MRMix_fitr <- my_MRMix_fit(gammah_r[id4ld], Gammah_r[id4ld], se1_r[id4ld], se2_r[id4ld])
  bhat0[irep,   10] = MRMix_fitr$beta.hat;
  pvalue0[irep, 10] = MRMix_fitr$beta.p.value;
  
  
  if(irep%%10==0) {
    print(irep);
    print(Sys.time());
    save(pvalue0, bhat0, pvalue1, bhat1, file = filename);
  }
  
}