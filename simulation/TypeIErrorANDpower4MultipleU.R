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
H2a.pa = c(0.05, 0.1);
H2t.pa = c(0.02, 0.05);
B.pa = c(0.1, 0);
Rho.pa = c(0.4, 0.8);
L.pa = 100;
AL.pa = 0.1;

Sim = "TypeIerror";
if(Sim=="TypeIerror"){
  H2g.pa = 0;
}else if(Sim=="Power"){
  aa = seq(1, 3.5, 0.5);
  H2g.pa = c(0, 10^(-aa));
}

para.all <- expand.grid(H2a.pa, H2t.pa, Rho.pa, L.pa, AL.pa, H2g.pa)
for(i in 1:nrow(para.all)){
  assign(paste0("res",i), para.all[i, ])
}

h2a <- as.numeric(res1[1]);
h2t <- as.numeric(res1[2]);
rho <- as.numeric(res1[3]);
L <- as.numeric(res1[4]);
Alrate <- as.numeric(res1[5]);
h2g <- as.numeric(res1[6]);
# Garate <- 1;
if(h2g==0){
  b0 = 0;
}else{
  b0 = 0.1;
}
M = 10; rho_ag = 0.2; n1 = 50000; n2 = 50000; n3 = 4000; lam = 0.85;


p = M*L;
block_inf <- cbind(seq(1, p, M), seq(M, p, M));
block_inf1 <- block_inf - 1;

coreNum = 5;
if(h2g==0){
  h2gIndex = 0;
}else{
  h2gIndex = -log10(h2g)
}
filename <- paste0("r", 10*rho,  "h2a", h2a*100, "h2t", h2t*100,"h2g", h2gIndex,
                   "p", p,  "AG", rho_ag,
                   "Al", Alrate, "lam", lam*100, "MutipleU5.Rdata");
if(Sim=="TypeIerror"){nrep = 1000;}else{nrep = 500}


maf = runif(p,0.05,0.5);
x = genRawGeno(maf, L, M, rho, n1 + n2 + n3);
#-------------------------------------------------------------------------#
x1 = x[1:n1,];
x2 = x[(n1+1):(n1+n2),];
x12 = x[1:(n1+n2),];
x3 = x[(n1+n2+1):(n1+n2+n3),];
R = Cal_block_SimR(block_inf1, x3, lam);
#-------------------------------------------------------------------------#

m <- p <- M*L;
nblocks = L;
block_inf <- cbind(seq(1, p, M), seq(M, p, M));
block_inf1 <- block_inf - 1;
coreNum =24;
#-------------------------------------------------------------------------#
if(rho==0.4){
  id4ld <- seq(1, p, 5);
  m = p;
  index = id4ld;
}else if(rho==0.8){
  id4ld <- seq(1, p, 20);
  m = p;
  index = id4ld;
}else if(rho==0){
  index = 1:p;
}

bhat = matrix(NA, nrow = nrep, ncol = 10);
pvalue = matrix(NA, nrow = nrep, ncol = 10);


Delta <- c(0.02, 0.05, 0.1, 0.5, 0.8)
SigA2 <- c(0.02, 0.01, 0.01, 0.005, 0.005);
tmpM = matrix(1:p, ncol = p / L, nrow = L, byrow = TRUE);

if(file.exists(filename)) {
  da = load(filename)
  ini_irep = max(which(!is.na(bhat[, 1])))+1
} else {
  ini_irep = 1
}


for(irep in ini_irep:nrep){
  q = 50
  u = matrix(rnorm( (n1+n2) * q),ncol=q);
  sigma2g = 1;
  sigma2a = 1;
  
  S_ag = matrix(c(sigma2a, rho_ag, rho_ag, sigma2g), nrow=2);
  AG = rmvnorm(p, mean=rep(0, 2), sigma = S_ag,method="chol");
  alpha = AG[,1]; gamma = AG[,2];
  
  Garate = 1;
  if(Garate!=1){
    gano = floor(p*(1 - Garate));
    indxGA = sample(1:p,gano);
    gamma[indxGA] = 0;
  }
  
  sigma2g = 0.005;
  gamma = rnorm(p)*sqrt(sigma2g);
  
  # ------------------------------------------------------------------------
  # The correlated horizontal pleiotropy (alpha).
  # First:generate alpha titlde and delta to obtain alpha.
  alpha <- rep(NA, p);
  a = sample(rep(1:5, L/5), L,rep = F); # split the blocks into 5 groups.
  for(g in 1:5){
    gL = which(a == g);
    idL = tmpM[gL, ];
    alphaT = rnorm(length(idL))*sqrt(SigA2[g]);
    alpha[idL] = Delta[g]*gamma[idL] + alphaT;
  }
  
  # ------------------------------------------------------------------------
  Su = matrix(c(1,0.8,0.8,1),nrow=2)
  bu = rmvnorm(q,mean=rep(0,2), sigma = Su,method="chol")
  by = bu[,1]; bz = bu[,2];
  uby = u%*%by; ubz = u%*%bz;
  uby = uby/sqrt(as.numeric(var(uby)/0.6));
  ubz = ubz/sqrt(as.numeric(var(ubz)/0.2));
  
  x12g = x12%*%gamma;
  
  var1 = (h2g*b0*b0 + h2g)/(b0*b0*(1 - (h2g + h2a + h2t)));
  var2 = (h2a*b0*b0 + h2a)/(1 - (h2g + h2a + h2t));
  var3 = (h2t*b0*b0 + h2t)/(1 - (h2g + h2a + h2t));
  
  if(b0!=0){
    gamma0 = gamma/sqrt(as.numeric(var(x12g)/var1));
    x12g = x12%*%gamma0;
  }
  
  yall = x12g + uby + rnorm(n1+n2)*as.numeric(sqrt(1-var(uby)));
  
  #-------------------------------------------------------------------------#
  # remark that alpha is sparse
  if(h2a==0){
    alpha = rep(0, p);
    x12a = x12%*%alpha;
  }else{
    if(Alrate!=1){
      if(Alrate!=1){
        alno = floor(L*Alrate);
        ind = sample(1:L,alno);
        if(length(ind)==1){
          indxAL = block_inf[ind, 1]:block_inf[ind, 2]
          alpha[-indxAL] = 0;
        }else{
          indxAL = NULL;
          for(i in 1:length(ind)){
            tmp = block_inf[ind[i], 1]:block_inf[ind[i], 2]
            indxAL = append(indxAL, tmp);
          }
          alpha[-indxAL] = 0;
        }
        x12a = x12%*%alpha;
        alpha0 = alpha/sqrt(as.numeric(var(x12a)/(var2)));
        x12a = x12%*%alpha0;
      }
    }
    
  }
  
  # The independent horizontal pleiotropy (theta).
  sigma2t <- 0.005;
  if(h2t==0){
    theta0 = rep(0, m);
    x12t = x12%*%theta0;
  }else{
    theta = rnorm(p)*sqrt(sigma2t);
    x12t = x12%*%theta;
    # theta0 = theta/sqrt(as.numeric(var(x12t)/(h2al)));
    theta0 = theta/sqrt(as.numeric(var(x12t)/var3));
    x12t = x12%*%theta0;
  }
  sgth2 = var(theta0);
  
  
  #-------------------------------------------------------------------------#
  resz = ubz + rnorm(n1+n2)*as.numeric(sqrt(1-var(ubz)));
  zall = b0*yall  + x12a +  x12t + resz;
  
  y = yall[1:n1];
  z = zall[(n1+1):(n1+n2)];
  
  x1 = x12[1:n1, ];
  x2 = x12[(n1+1):(n1+n2), ]
  var(b0*yall)/var(zall)
  
  var(x12a)/var(zall)
  var(x12t)/var(zall)
  #-------------------------------------------------------------------------#
  # create summary statistics
  gammaSS = fastSigLm(y, x1);
  GammaSS = fastSigLm(z, x2);
  gammah = gammaSS$coef;
  se1 = gammaSS$std;
  Gammah = GammaSS$coef;
  se2 = GammaSS$std;
  
  #*******************************************************************************#
  # For Correlated SNPs
  #*******************************************************************************#
  # ========================= #
  # MR-CUE
  # ========================= #
  resM3 = MRCUESim(gammah, Gammah, se1, se2, 0, R, block_inf1, coreNum);
  
  M3beta = resM3$beta.hat;
  M3pvalue = resM3$beta.p.value;
  bhat[irep, 1] = M3beta;
  pvalue[irep, 1] = M3pvalue;
  
  # ========================= #
  # MR-LDP
  # ========================= #
  # initial value.
  epsStopLogLik <- 1e-7; maxIter <- 10000; beta0 <- 0;
  gamma <- rep(0, p);
  alpha <- rep(0, p);
  sgga2 <- 0.01;
  sgal2 <- 0.01;
  SimMRLDP_Hb = MRLDP_SimPXvb(gammah, Gammah, se1, se2, gamma, alpha,  beta0, sgga2, sgal2, R,
                              0, epsStopLogLik, maxIter, model = 2);
  
  beta_hat2 = SimMRLDP_Hb$beta0;
  
  SimMRLDP_H0 = MRLDP_SimPXvb(gammah, Gammah, se1, se2, gamma, alpha,  beta0, sgga2, sgal2, R,
                              1, epsStopLogLik, maxIter, model = 2);
  
  tstat = 2*(SimMRLDP_Hb$tstat - SimMRLDP_H0$tstat);
  pval2 = pchisq(tstat, 1, lower.tail = F);
  
  bhat[irep, 2] = beta_hat2;
  pvalue[irep, 2] = pval2;
 
  #*******************************************************************************#
  # For Independent SNPs
  #*******************************************************************************#
  # ========================= #
  # CAUSE
  # ========================= #
  ld_prune_pval_thresh = 1e-3;
  r2_thresh = 0.001;
  
  snp = 1:p;
  p_value = 2*pnorm(-abs(gammah/se1));
  ld_prune = rep(5==1, p);
  ld_prune[id4ld]= TRUE;
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
    bhat[irep, 3] = qs[[2]][1, 1];
    pvalue[irep, 3] = pnorm(-z);
  }
  
  # ========================= #
  # Raps
  # ========================= #
  raps = mr.raps(gammah[id4ld], Gammah[id4ld], se1[id4ld], se2[id4ld], over.dispersion = TRUE);
  bhat[irep, 4] = raps$beta.hat;
  pvalue[irep, 4] = raps$beta.p.value;
  
  # ========================= #
  # MR-cML(DP)
  # ========================= #
  cML_result_DP = mr_cML_DP(gammah[id4ld], Gammah[id4ld], se1[id4ld], se2[id4ld], n = n1);
  bhat[irep, 5] = cML_result_DP$MA_BIC_DP_theta;
  pvalue[irep, 5] = cML_result_DP$MA_BIC_DP_p;
  
  # ========================= #
  # IVW 
  # ========================= #
  mr.obj = mr_input(bx = as.vector(gammah[id4ld]), bxse = as.vector(se1[id4ld]),
                    by = as.vector(Gammah[id4ld]), byse = as.vector(se2[id4ld]));
  
  IVW = mr_ivw(mr.obj);
  bhat[irep, 6] = IVW$Estimate;
  pvalue[irep, 6] = IVW$Pvalue;
  
  # ========================= #
  # MR-Egger
  # ========================= #
  Egger = try(mr_egger(mr.obj));
  
  if(class(Egger) != 'try-error'){
    bhat[irep, 7] = Egger$Estimate;
    pvalue[irep, 7] = Egger$Pvalue.Est;
  }
  
  # =========================#
  # GRAPPLE
  # =========================#
  GRAPPLE_fit <- my_GRAPPLE_fit(gammah[id4ld], Gammah[id4ld], se1[id4ld], se2[id4ld])
  bhat[irep,   8] = GRAPPLE_fit$beta.hat;
  pvalue[irep, 8] = GRAPPLE_fit$beta.p.value;
  
  # ========================= #
  # MR-Clust
  # ========================= #
  mrclust_fit <- my_mrclust_fit(gammah[id4ld], Gammah[id4ld], se1[id4ld], se2[id4ld])
  bhat[irep, 9] = mrclust_fit$beta.hat;
  pvalue[irep, 9] = mrclust_fit$beta.p.value;
  
  # ========================= #
  # MRMix
  # ========================= #
  MRMix_fit <- my_MRMix_fit(gammah[id4ld], Gammah[id4ld], se1[id4ld], se2[id4ld])
  bhat[irep, 10] = MRMix_fit$beta.hat;
  pvalue[irep, 10] = MRMix_fit$beta.p.value;
  
  
  if(irep%%50==0){
    print(irep);
    print(Sys.time());
    save(bhat, pvalue, file = filename);
  }
}