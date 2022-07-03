rm(list = ls());
library("mvtnorm");
library(MR.CUE);

H2a.pa = c(0.05, 0.1);
H2t.pa = c(0.02, 0.05);
B.pa = c(0, 0.1);
Rho.pa = c(0.4, 0.8);
L.pa = c(100, 200);

para.all <- expand.grid(H2a.pa, H2t.pa, B.pa, Rho.pa, L.pa)
for(i in 1:nrow(para.all)){
  assign(paste0("res",i), para.all[i, ])
}

h2a <- as.numeric(res1[1]);
h2t <- as.numeric(res1[2]);
b0 <- as.numeric(res1[3]);
rho <- as.numeric(res1[4]);
L <- as.numeric(res1[5]);

Alrate = 0.1;
M = 10;

rho_ag = 0.2;
n1 = 50000;
n2 = 50000;
n3 = 4000;
lam = 0.85;
maxIter = 4000;
burnin = 1000;
thin = 10;

# -------------------------------------------------------
h2g <- 0.1;  Garate <- 1;
# -------------------------------------------------------
p = M*L;
block_inf <- cbind(seq(1, p, M), seq(M, p, M));
block_inf1 <- block_inf - 1;

p = M*L;

filename <- paste0("r", 10*rho,  "h2a", h2a*100, "h2t", h2t*10,
                   "p", p, "b", b0,  "AG",rho_ag,
                   "Ga", Garate, "Al", Alrate,"N",n1/10000, "Weta.Rdata");


maf = runif(p,0.05,0.5);
x = genRawGeno(maf, L, M, rho, n1 + n2 + n3);
#-------------------------------------------------------------------------#

x1 = x[1:n1,];
x2 = x[(n1+1):(n1+n2),];
x12 = x[1:(n1+n2),];
x3 = x[(n1+n2+1):(n1+n2+n3),];
R = Cal_block_SimR(block_inf1, x3, lam)
#-------------------------------------------------------------------------#
coreNum = 20;
nrep = 100;
bres11 = se11 = pva11 = bres2L = se2L = pva2L= rep(NA, nrep);
TrueEta = EstEta11 = EstEta2L = matrix(NA, nrow = L, ncol = nrep);

for(irep in 1:nrep){
  q = 50
  u = matrix(rnorm( (n1+n2) * q),ncol=q);
  sigma2g = 1;
  sigma2a = 1;
  
  S_ag = matrix(c(sigma2a, rho_ag, rho_ag, sigma2g), nrow=2);
  AG = rmvnorm(p, mean=rep(0, 2), sigma = S_ag,method="chol")
  alpha = AG[,1]; gamma = AG[,2];
  
  
  if(Garate!=1){
    gano = floor(p*(1 - Garate));
    indxGA = sample(1:p,gano);
    gamma[indxGA] = 0;
  }
  
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
  
  # ------------------------------------------------------------------------
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
  # ------------------------------------------------------------------------
  # The correlated horizontal pleiotropy (alpha).
  
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
  
  eta = rep(0, L);
  eta[ind] = 1;
  TrueEta[, irep] = eta;
  # ------------------------------------------------------------------------
  
  
  resz = ubz + rnorm(n1+n2)*as.numeric(sqrt(1-var(ubz)));
  zall = b0*yall  + x12a +  x12t + resz;
  
  # var(x12t) /var(zall);
  # var(x12a) /var(zall);
  # var(b0*x12g) /var(zall);
  
  
  y = yall[1:n1];
  z = zall[(n1+1):(n1+n2)];
  
  x1 = x12[1:n1, ];
  x2 = x12[(n1+1):(n1+n2), ]
  # ------------------------------------------------------------------------
  # create summary statistics
  gammaSS = fastSigLm(y, x1);
  GammaSS = fastSigLm(z, x2);
  gammah = gammaSS$coef;
  se1 = gammaSS$std;
  Gammah = GammaSS$coef;
  se2 = GammaSS$std;
  # ------------------------------------------------------------------------
  opt11 = list(agm = 0, bgm = 0, atau1 = 0, btau1 = 0,
             atau2 = 0, btau2 = 0,
             a = 1, b = 1, maxIter = maxIter, thin = 10, burnin = burnin);
  opt2L = list(agm = 0, bgm = 0, atau1 = 0, btau1 = 0,
               atau2 = 0, btau2 = 0,
               a = 2, b = L, maxIter = maxIter, thin = 10, burnin = burnin);
  
  resM11 = MRCUESim(gammah, Gammah, se1, se2, 0, R, block_inf1, coreNum, opt11);
  bres11[irep] = resM11$beta.hat;
  se11[irep] = resM11$beta.se;
  pva11[irep] = resM11$beta.p.value;
  EstEta11[, irep] = resM11$EtaIterRate;
  
  resM2L = MRCUESim(gammah, Gammah, se1, se2, 0, R, block_inf1, coreNum, opt2L);
  bres2L[irep] = resM2L$beta.hat;
  se2L[irep] = resM2L$beta.se;
  pva2L[irep] = resM2L$beta.p.value;
  EstEta2L[, irep] = resM2L$EtaIterRate;
  # ------------------------------------------------------------------------
  
  if(irep%%10==0) {
    print(irep);
    print(Sys.time());
    save(bres11, se11, pva11, bres2L, se2L, pva2L, TrueEta, EstEta11, EstEta2L, file = filename);
  }
}