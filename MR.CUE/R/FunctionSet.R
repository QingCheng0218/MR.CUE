genRawGeno <- function(maf, L, M, rho, n){
  SIGMA = matrix(nrow=M,ncol=M)
  for (i in 1:M){
    for (j in 1:M){
      SIGMA[i,j] = rho^(abs(i-j));
    }
  }
  
  nsnp = L*M;
  X = NULL;
  for ( l in 1:L ){
    
    index = (M*(l-1)+1): (M*l);
    AAprob = maf[index]^2.;
    Aaprob = 2*maf[index]*(1-maf[index]);
    quanti = matrix(c(1-Aaprob-AAprob, 1- AAprob),M,2);
    Xt = rmvnorm(n, mean=rep(0,M), sigma=SIGMA, method="chol")
    Xt2 = matrix(0,n,M);
    for (j in 1:M){
      cutoff = qnorm(quanti[j,]);
      Xt2[Xt[,j] < cutoff[1],j] = 0;
      Xt2[Xt[,j] >= cutoff[1] & Xt[,j] < cutoff[2],j] = 1;  ## attention
      Xt2[Xt[,j] >= cutoff[2],j] = 2;
    }
    X <- cbind(X,Xt2);
  }
  return(X)
}


genSumStat <- function(x12, n1, n2, M, L, b1, rho_ag, Alrate, h2a, h2t, h2g){
  
  p = M*L;
  block_inf <- cbind(seq(1, p, M), seq(M, p, M));
  
  q = 50
  u = matrix(rnorm( (n1+n2) * q),ncol=q);
  
  # ------------------------------------------------------------------------
  sigma2g = 1;
  sigma2a = 1;
  
  S_ag = matrix(c(sigma2a, rho_ag, rho_ag, sigma2g), nrow=2);
  AG = rmvnorm(p, mean=rep(0, 2), sigma = S_ag,method="chol")
  alpha = AG[,1]; gamma = AG[,2];
  
  # Assume that gamma is dense.
  Garate = 1;
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
  
  var1 = (h2g*b1*b1 + h2g)/(b1*b1*(1 - (h2g + h2a + h2t)));
  var2 = (h2a*b1*b1 + h2a)/(1 - (h2g + h2a + h2t));
  var3 = (h2t*b1*b1 + h2t)/(1 - (h2g + h2a + h2t));
  
  if(b1!=0){
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
      if(M!=1){
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
      }else{
        alno = floor(L*Alrate);
        ind = sample(1:L, alno);
        alpha[-ind] = 0;
        x12a = x12%*%alpha;
        alpha0 = alpha/sqrt(as.numeric(var(x12a)/(var2)));
        x12a = x12%*%alpha0;
      }
    }
    
  }
  
  # ------------------------------------------------------------------------
  
  
  resz = ubz + rnorm(n1+n2)*as.numeric(sqrt(1-var(ubz)));
  zall = b1*yall  + x12a +  x12t + resz;
  
  
  y = yall[1:n1];
  z = zall[(n1+1):(n1+n2)];
  
  x1 = x12[1:n1, ];
  x2 = x12[(n1+1):(n1+n2), ]
  # ------------------------#
  # create summary statistics
  gammaSS = fastSigLm(y, x1);
  GammaSS = fastSigLm(z, x2);
  gammah = gammaSS$coef;
  se1 = gammaSS$std;
  Gammah = GammaSS$coef;
  se2 = GammaSS$std;
  # ------------------------#
  return(list(gammah = gammah, se1 = se1, se2 = se2, Gammah = Gammah, CHPindex = ind))
}


mhcstart = 28477797;
mhcend = 33448354;
summaryQC = function(mhcstart, mhcend, bh1, bh2, s12, s22, bp, chr,
                     rsname, avbIndex, idx4panel, xbound, ybound){
  # remove SNPs in MHC region:
  idxchr6 = which(chr==6);
  idxcut = idxchr6[which(bp[idxchr6]>=mhcstart & bp[idxchr6]<=mhcend)];
  pmhc = length(idxcut);
  
  if(pmhc!=0){
    bh1Rmhc = bh1[-idxcut];
    bh2Rmhc = bh2[-idxcut];
    s12Rmhc = s12[-idxcut];
    s22Rmhc = s22[-idxcut];
    bpRmhc = bp[-idxcut];
    chrRmhc = chr[-idxcut];
    rsnameRmhc = rsname[-idxcut];
    avbIndexRmhc = avbIndex[-idxcut];
    
    tmp0 = 1:length(bh1);
    tmp = tmp0[-idxcut];
    if(length(idx4panel)!=0){
      idx4panelRmhc = match(avbIndex[intersect((idx4panel + 1), tmp)], avbIndexRmhc) -1
    }else{
      idx4panelRmhc = idx4panel;
    }
    
  }else{
    bh1Rmhc = bh1;
    bh2Rmhc = bh2;
    s12Rmhc = s12;
    s22Rmhc = s22;
    bpRmhc = bp;
    chrRmhc = chr;
    rsnameRmhc = rsname;
    avbIndexRmhc = avbIndex;
    idx4panelRmhc = idx4panel;
  }
  
  
  # remove SNPs(exposure) with chi-square >80
  idx = which((bh1Rmhc/s12Rmhc)^2>xbound);
  px = length(idx);
  if(px!=0){
    bh1Rmhc_x = bh1Rmhc[-idx];
    bh2Rmhc_x = bh2Rmhc[-idx];
    s12Rmhc_x = s12Rmhc[-idx];
    s22Rmhc_x = s22Rmhc[-idx];
    bpRmhc_x = bpRmhc[-idx];
    chrRmhc_x = chrRmhc[-idx];
    rsnameRmhc_x = rsnameRmhc[-idx];
    avbIndexRmhc_x = avbIndexRmhc[-idx];
    # idx4panelRmhc_x = idx4panelRmhc[-idx];
    
    
    tmp0 = 1:length(bh1Rmhc);
    tmp = tmp0[-idx];
    if(length(idx4panel)!=0){
      idx4panelRmhc_x = match(avbIndexRmhc[intersect((idx4panelRmhc + 1), tmp)], avbIndexRmhc_x) -1;
    }else{
      idx4panelRmhc_x = idx4panel;
    }
    
    
  }else{
    bh1Rmhc_x = bh1Rmhc;
    bh2Rmhc_x = bh2Rmhc;
    s12Rmhc_x = s12Rmhc;
    s22Rmhc_x = s22Rmhc;
    bpRmhc_x = bpRmhc;
    chrRmhc_x = chrRmhc;
    rsnameRmhc_x = rsnameRmhc;
    avbIndexRmhc_x = avbIndexRmhc;
    idx4panelRmhc_x = idx4panelRmhc;
  }
  
  # remove SNPs(outcome) with chi-square >80
  idy = which((bh2Rmhc_x/s22Rmhc_x)^2>ybound);
  py = length(idy);
  if(py!=0){
    bh1Rmhc_xy = bh1Rmhc_x[-idy];
    bh2Rmhc_xy = bh2Rmhc_x[-idy];
    s12Rmhc_xy = s12Rmhc_x[-idy];
    s22Rmhc_xy = s22Rmhc_x[-idy];
    bpRmhc_xy = bpRmhc_x[-idy];
    chrRmhc_xy = chrRmhc_x[-idy];
    rsnameRmhc_xy = rsnameRmhc_x[-idy];
    avbIndexRmhc_xy = avbIndexRmhc_x[-idy];
    # idx4panelRmhc_xy = idx4panelRmhc_x[-idy];
    
    tmp0 = 1:length(bh1Rmhc_x);
    tmp = tmp0[-idx];
    
    if(length(idx4panel)!=0){
      idx4panelRmhc_xy = match(avbIndexRmhc_x[intersect((idx4panelRmhc_x + 1), tmp)], avbIndexRmhc_xy) -1;
    }else{
      idx4panelRmhc_xy = idx4panel;
    }
    
  }else{
    bh1Rmhc_xy = bh1Rmhc_x;
    bh2Rmhc_xy = bh2Rmhc_x;
    s12Rmhc_xy = s12Rmhc_x;
    s22Rmhc_xy = s22Rmhc_x;
    bpRmhc_xy = bpRmhc_x;
    chrRmhc_xy = chrRmhc_x;
    rsnameRmhc_xy = rsnameRmhc_x;
    avbIndexRmhc_xy = avbIndexRmhc_x;
    idx4panelRmhc_xy = idx4panelRmhc_x;
  }
  return(list(bh1new = bh1Rmhc_xy, bh2new = bh2Rmhc_xy, s12new = s12Rmhc_xy, s22new = s22Rmhc_xy,
              bpnew = bpRmhc_xy, chrnew = chrRmhc_xy, rsnamenew = rsnameRmhc_xy,
              avbIndexnew = avbIndexRmhc_xy, idx4panelnew = idx4panelRmhc_xy, pmhc = pmhc, px = px, py = py))
}

traceplot <- function(bhatpoint){
  y <- bhatpoint
  x <- 1:length(bhatpoint);
  da <- cbind(x, y);
  dat <- data.frame(da);
  p1 <- ggplot(data = dat, aes(x= x, y = y))+  geom_line()  +
    labs( x = paste0("GibbsSampleIndex"), y =  expression(hat(beta[0])));
  p1 = p1 + theme(axis.title.x = element_text(size=10,face = "bold"),
                  axis.text.x = element_text(size=12,face = "bold"),
                  axis.title.y = element_text(size=10,face = "bold"),
                  axis.text.y = element_text(size=12,face = "bold"));
  return(p1);
}



# ---------------------------------------------------------
EstRhofun <- function(fileexposure, fileoutcome, stringname3,
                      ld_r2_thresh, lam, pth){
  
  # Estimate the rho
  res = matchsnp(fileexposure, fileoutcome, stringname3, FALSE);
  bh1 = as.numeric(res$bh1);
  bh2 = as.numeric(res$bh2);
  s12 = as.numeric(res$s12);
  s22 = as.numeric(res$s22);
  chr = as.numeric(res$chr);
  bp = res$bp;
  rsname = res$rsname
  avbIndex = res$idxin;
  idx4panel = res$idx4panel;
  QCresult = summaryQC(mhcstart, mhcend, bh1, bh2, s12, s22, bp,
                       chr, rsname, avbIndex, idx4panel, Inf, Inf);
  
  
  bh1new = QCresult$bh1new;
  bh2new = QCresult$bh2new;
  s12new = QCresult$s12new;
  s22new = QCresult$s22new;
  bpnew = QCresult$bpnew;
  chrnew = QCresult$chrnew;
  avbIndexnew = QCresult$avbIndexnew;
  rsnamenew = QCresult$rsnamenew;
  pmhc = QCresult$pmhc;
  px = QCresult$px;
  py = QCresult$py;
  
  coreNum = 24;
  
  IndSumRes = IndepSummary(bpnew, chrnew, avbIndexnew - 1, block_file, stringname3, 
                           bh1new, bh2new, s12new, s22new, coreNum, 
                           lam, ld_r2_thresh);
  bh1_ind = IndSumRes$bh1_ind;
  bh2_ind = IndSumRes$bh2_ind;
  se1_ind = IndSumRes$se1_ind;
  se2_ind = IndSumRes$se2_ind;
  
  z1_ind = bh1_ind / se1_ind;
  z2_ind = bh2_ind / se2_ind;
  
  # a = rep(-pth, 2);
  # b = rep(pth, 2);
  # z1_new = z1_ind [which(abs(z1_ind) < pth&abs(z2_ind) < pth)];
  # z2_new = z2_ind[which(abs(z1_ind) < pth&abs(z2_ind) < pth)];
  # rhores = truncEstfun(a, b, z1_new, z2_new, 4000, 1000, 10)
  # rhohat = mean(rhores);
  # p1 = length(z1_new);
  # pvalue = testR(rhohat, p1);
  # 
  maxIter = 4000;
  thin = 10;
  burnin = 1000;
  
  nsave = maxIter / thin;
  
  if(length(pth)==1){
    Rhores = rep(NA, nsave);
    a = rep(-pth, 2);
    b = rep(pth, 2);
    z1_new = z1_ind [which(abs(z1_ind) < pth&abs(z2_ind) < pth)];
    z2_new = z2_ind[which(abs(z1_ind) < pth&abs(z2_ind) < pth)];
    rhores = truncEstfun(a, b, z1_new, z2_new, maxIter, burnin, thin)
    rhohat = mean(rhores);
    p1 = length(z1_new);
    pvalue = testR(rhohat, p1);
    Rhores = rhores;
    pres = p1;
  }else{
    rhohat = pvalue = rep(NA, length(pth));
    Rhores = matrix(NA, nrow = nsave, ncol = length(pth));
    pres = rep(NA, length(pth));
    for(i in 1:length(pth)){
      pth1 = pth[i];
      a = rep(-pth1, 2);
      b = rep(pth1, 2);
      z1_new = z1_ind [which(abs(z1_ind) < pth1&abs(z2_ind) < pth1)];
      z2_new = z2_ind[which(abs(z1_ind) < pth1&abs(z2_ind) < pth1)];
      rhores = truncEstfun(a, b, z1_new, z2_new, 4000, 1000, 10)
      rhohat[i] = mean(rhores);
      p1 = length(z1_new);
      pvalue[i] = testR(rhohat[i], p1);
      Rhores[, i] = rhores;
      pres[i] = p1;
    }
  }
  # ---------------------------------------------------------
  return(list(rhohat = rhohat, pvalue = pvalue, pres = pres, Rhores = Rhores));
  
}
