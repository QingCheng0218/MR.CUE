# ------------------------------------------------------------ #
# Remark
# The sigle SNP should be in one block (June 18, 2022, CQ)
# Remove the MHC region. (June 26, 2022. CQ)
# ------------------------------------------------------------ #

mhcstart = 28477797;
mhcend = 33448354;
RemoveMHCCfun = function(mhcstart, mhcend, bh1, bh2, s12, s22, bp, chr,
                         rsname,  xbound, ybound){
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
    
    
  }else{
    bh1Rmhc = bh1;
    bh2Rmhc = bh2;
    s12Rmhc = s12;
    s22Rmhc = s22;
    bpRmhc = bp;
    chrRmhc = chr;
    rsnameRmhc = rsname;
    
  }
  
  # If xbound = Inf & ybound = Inf, skip this procedure.
  # remove SNPs(exposure) with chi-square > xbound
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
    
  }else{
    bh1Rmhc_x = bh1Rmhc;
    bh2Rmhc_x = bh2Rmhc;
    s12Rmhc_x = s12Rmhc;
    s22Rmhc_x = s22Rmhc;
    bpRmhc_x = bpRmhc;
    chrRmhc_x = chrRmhc;
    rsnameRmhc_x = rsnameRmhc;
  }
  
  # remove SNPs(outcome) with chi-square >ybound
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
  }
  return(list(bh1new = bh1Rmhc_xy, bh2new = bh2Rmhc_xy, s12new = s12Rmhc_xy, s22new = s22Rmhc_xy,
              bpnew = bpRmhc_xy, chrnew = chrRmhc_xy, rsnamenew = rsnameRmhc_xy,
              pmhc = pmhc, px = px, py = py))
}

ReadSummaryStat <- function(fileexp, fileout, filepan, snpinfo, pva_cutoff, lambad) {
  daexp = read.table(fileexp, header = T)
  daout = read.table(fileout, header = T)
  
  daexp = daexp[daexp$pvalue < pva_cutoff, ]
  daexp$a1 = 0
  daexp$a2 = 0
  daexp$a1[daexp$A1=="C"] = 1
  daexp$a1[daexp$A1=="G"] = 1
  daexp$a2[daexp$A2=="C"] = 1
  daexp$a2[daexp$A2=="G"] = 1
  daexp = daexp[daexp$a1+daexp$a2==1,]
  daexp[,c("A1", "A2", "a2", "pvalue", "BP")] = NULL
  daexp$chr = as.numeric(daexp$chr)
  daexp = unique(daexp)
  rownames(daexp) = daexp$SNP
  
  daout = daout[daout$SNP %in% daexp$SNP, ]
  rownames(daout) = daout$SNP
  daout$a1 = 0
  daout$a2 = 0
  daout$a1[daout$A1=="C"] = 1
  daout$a1[daout$A1=="G"] = 1
  daout$a2[daout$A2=="C"] = 1
  daout$a2[daout$A2=="G"] = 1
  daout = daout[daout$a1+daout$a2==1,]
  bp = daout$BP;
  daout[,c("A1", "A2", "a2", "pvalue", "BP")] = NULL
  daout$chr = as.numeric(daout$chr)
  daout = unique(daout)
  
  daexp = daexp[daout$SNP,]
  # CQ: daout and daexp do not need match the A1 !!!
  # ===================================================== #
  # Remove MHC region.
  QCres = RemoveMHCCfun(mhcstart, mhcend, daexp$beta, daout$beta, daout$beta, daout$se, bp, daout$chr,
                        daout$SNP,  Inf, Inf)
  rsnameRMHC = QCres$rsnamenew
  if(length(rsnameRMHC)!=dim(daout)[1]){
    daexp = daexp[daexp$SNP%in%rsnameRMHC, ];
    daout = daout[daout$SNP%in%rsnameRMHC, ];
  }
  
  # ===================================================== #
  daRefsnpinfo = readRDS(snpinfo);
  daRefsnpinfosame = daRefsnpinfo[daRefsnpinfo$SNP%in%daexp$SNP, ];
  daRefsnpinfosame$a1 = 0
  daRefsnpinfosame$a1[daRefsnpinfosame$A1=="C"] = 1
  daRefsnpinfosame$a1[daRefsnpinfosame$A1=="G"] = 1
  daRefsnpinfosame[, "A1"] = NULL
  # ===================================================== #
  
  ResF4gamma   = vector("list", 0)
  ResF4gammase = vector("list", 0)
  ResF4Gamma   = vector("list", 0)
  ResF4Gammase = vector("list", 0)
  ResF4SNPchr  = vector("list", 0)
  ResF4R       = vector("list", 0)
  Res_id  = c()
  # CQ: In case that not all 22 chromosomes are included.
  comchrs = sort(unique(daexp$chr));
  
  for (chrk in comchrs) {
    if (is.null(filepan[[chrk]])) next
    if (!file.exists(filepan[[chrk]])) next
    dapan = readRDS(filepan[[chrk]])
    dapan[] <- lapply(dapan, function(x) if (is.factor(x)) as.character(x) else {x})
    
    daexp_k = daexp[daexp$chr==chrk,]
    daout_k = daout[daexp_k$SNP,]
    # ===================================================== #
    RefSNP_k = daRefsnpinfosame[daRefsnpinfosame$CHR==chrk, ]
    snpint = intersect(union(dapan$SNP1, dapan$SNP2), daout_k$SNP);
    stopifnot(sum(RefSNP_k$SNP%in%snpint)==length(snpint))
    # ===================================================== #
    
    dapan = dapan[dapan$SNP1 %in% daout_k$SNP,]
    dapan = dapan[dapan$SNP2 %in% daout_k$SNP,]
    snppan = union(dapan$SNP1, dapan$SNP2)
    
    if(length(snppan)!=length(snpint)){
      # CQ: In case that each block just has one snp.
      sigleSNP = snpint[!snpint%in%snppan];
      if(length(snppan)!=0){
        Sigpan = data.frame(rep(chrk, length(sigleSNP)),paste0(chrk*10, max(dapan$BlockID)+seq(1:length(sigleSNP))), sigleSNP, sigleSNP, rep(1, length(sigleSNP)))
        colnames(Sigpan) = colnames(dapan)
        dapan = rbind(dapan, Sigpan)
      }else{
        Sigpan = data.frame(rep(chrk, length(sigleSNP)),paste0(chrk*10, seq(1:length(sigleSNP))), sigleSNP, sigleSNP, rep(1, length(sigleSNP)))
        colnames(Sigpan) = colnames(dapan);
        dapan = rbind(dapan, Sigpan);
      }
      
    }
    
    if(nrow(dapan)==0) next
    
    # move row 345-347 to here, as dapan may expand in line 349-355, ZX
    snp = union(dapan$SNP1, dapan$SNP2)
    daexp_k = daexp_k[daexp_k$SNP %in% snp, ]
    daout_k = daout_k[daexp_k$SNP,]
    
    
    arma_blc_id = unique(dapan$BlockID)
    
    nid = length(arma_blc_id)
    F4gamma   = vector("list", nid)
    F4gammase = vector("list", nid)
    F4Gamma   = vector("list", nid)
    F4Gammase = vector("list", nid)
    F4SNPchr  = vector("list", nid)
    F4R       = vector("list", nid)
    snpset = rep(FALSE, nid)
    for (k in 1:nid) {
      uid = arma_blc_id[k]
      da2 = dapan[dapan$BlockID==uid,]
      comsnp = union(da2$SNP1, da2$SNP2)
      
      daRef = RefSNP_k[RefSNP_k$SNP%in%comsnp, ];
      da0 = daexp_k[comsnp, ]; 
      da1 = daout_k[comsnp, ];
      # ===================================================== #
      # match Ref panel data.
      da0$beta[da0$a1!=daRef$a1] = -1*da0$beta[da0$a1!=daRef$a1];
      da1$beta[da0$a1!=daRef$a1] = -1*da1$beta[da0$a1!=daRef$a1];
      # ===================================================== #
      
      F4gamma[[k]] = da0$beta
      F4gammase[[k]] = da0$se
      F4Gamma[[k]] = da1$beta
      F4Gammase[[k]] = da1$se
      F4SNPchr[[k]] = comsnp
      
      M = diag(rep(1,length(comsnp)))
      rownames(M) = comsnp
      colnames(M) = comsnp
      if(length(comsnp)!=1){
        for (j in 1:nrow(da2)) {
          M[da2$SNP1[j], da2$SNP2[j]] = da2$r[j]*lambad
          M[da2$SNP2[j], da2$SNP1[j]] = da2$r[j]*lambad
        }
      }
      rownames(M) = NULL
      colnames(M) = NULL
      F4R[[k]] = M
    }
    
    ResF4gamma    = append(ResF4gamma, F4gamma)
    ResF4gammase  = append(ResF4gammase, F4gammase)
    ResF4Gamma    = append(ResF4Gamma, F4Gamma)
    ResF4Gammase = append(ResF4Gammase, F4Gammase)
    ResF4SNPchr   = append(ResF4SNPchr, F4SNPchr)
    ResF4R        = append(ResF4R, F4R)
    Res_id        = c(Res_id, arma_blc_id)
  }
  
  return(
    list(
      ResF4gammah = ResF4gamma,
      ResF4se1 = ResF4gammase,
      ResF4Gammah = ResF4Gamma,
      ResF4se2 = ResF4Gammase,
      ResF4SNPchr = ResF4SNPchr,
      ResF4Rblock = ResF4R,
      arma_blc_id = Res_id
    )
  )
}


EstRho <- function(fileexp, fileout, filepan, snpinfo, ld_r2_thresh, lambad, pth) {
  
  daexp = read.table(fileexp, header = T)
  daout = read.table(fileout, header = T)
  
  daexp$a1 = 0
  daexp$a2 = 0
  daexp$a1[daexp$A1=="C"] = 1
  daexp$a1[daexp$A1=="G"] = 1
  daexp$a2[daexp$A2=="C"] = 1
  daexp$a2[daexp$A2=="G"] = 1
  daexp = daexp[daexp$a1+daexp$a2==1,]
  daexp[,c("A1", "A2", "a2", "pvalue", "BP")] = NULL
  daexp$chr = as.numeric(daexp$chr)
  daexp = daexp[!is.na(daexp$chr),]
  daexp = unique(daexp)
  rownames(daexp) = daexp$SNP
  
  daout = daout[daout$SNP %in% daexp$SNP, ]
  daout = unique(daout)
  rownames(daout) = daout$SNP
  daout$a1 = 0
  daout$a2 = 0
  daout$a1[daout$A1=="C"] = 1
  daout$a1[daout$A1=="G"] = 1
  daout$a2[daout$A2=="C"] = 1
  daout$a2[daout$A2=="G"] = 1
  daout = daout[daout$a1+daout$a2==1,]
  bp = daout$BP;
  daout[,c("A1", "A2", "a2", "pvalue", "BP")] = NULL
  daout$chr = as.numeric(daout$chr)
  daout = daout[!is.na(daout$chr),]
  
  daexp = daexp[daout$SNP,]
  # ===================================================== #
  # Remove MHC region.
  QCres = RemoveMHCCfun(mhcstart, mhcend, daexp$beta, daout$beta, daexp$se, daout$se, bp, daout$chr,
                        daout$SNP,  Inf, Inf)
  rsnameRMHC = QCres$rsnamenew;
  
  if(length(rsnameRMHC)!=dim(daout)[1]){
    daexp = daexp[daexp$SNP%in%rsnameRMHC, ];
    daout = daout[daout$SNP%in%rsnameRMHC, ];
  }
  
  
  # ===================================================== #
  daRefsnpinfo = readRDS(snpinfo);
  daRefsnpinfosame = daRefsnpinfo[daRefsnpinfo$SNP%in%daexp$SNP, ];
  daRefsnpinfosame$a1 = 0
  daRefsnpinfosame$a1[daRefsnpinfosame$A1=="C"] = 1
  daRefsnpinfosame$a1[daRefsnpinfosame$A1=="G"] = 1
  daRefsnpinfosame[, "A1"] = NULL
  
  # ===================================================== #
  # Random select SNPs in case that the number of SNPs are too large.
  if(dim(daRefsnpinfosame)[1] > 1e6){
    idxrand = sort(sample(1:dim(daRefsnpinfosame)[1], 1000000, replace = FALSE));
    daRefsnpinfosame = daRefsnpinfosame[idxrand, ];
    daexp = daexp[daexp$SNP%in%daRefsnpinfosame$SNP, ];
    daout = daexp[daout$SNP%in%daRefsnpinfosame$SNP, ];
  }
  
  # ===================================================== #
  ResF4gamma   = vector("list", 0)
  ResF4gammase = vector("list", 0)
  ResF4Gamma   = vector("list", 0)
  ResF4Gammase = vector("list", 0)
  # ResF4SNPchr  = vector("list", 0)
  # ResF4R       = vector("list", 0)
  # Res_id  = c()
  comchrs = sort(unique(daexp$chr));
  for (chrk in comchrs) {
    if (is.null(filepan[[chrk]])) next
    if (!file.exists(filepan[[chrk]])) next
    dapan = readRDS(filepan[[chrk]])
    dapan[] <- lapply(dapan, function(x) if (is.factor(x)) as.character(x) else {x})
    
    daexp_k = daexp[daexp$chr==chrk,]
    daout_k = daout[daexp_k$SNP,]
    
    # ===================================================== #
    RefSNP_k = daRefsnpinfosame[daRefsnpinfosame$CHR==chrk, ]
    snpint = intersect(union(dapan$SNP1, dapan$SNP2), daout_k$SNP);
    stopifnot(sum(RefSNP_k$SNP%in%snpint)==length(snpint))
    # ===================================================== #
    
    dapan = dapan[dapan$SNP1 %in% daout_k$SNP,]
    dapan = dapan[dapan$SNP2 %in% daout_k$SNP,]
    
    # snp = union(dapan$SNP1, dapan$SNP2)
    # daexp_k = daexp_k[daexp_k$SNP %in% snp, ]
    # daout_k = daout_k[daexp_k$SNP,]
    
    snppan = union(dapan$SNP1, dapan$SNP2)
    
    if(length(snppan)!=length(snpint)){
      # In case that each block just has one snp.
      sigleSNP = snpint[!snpint%in%snppan];
      if(length(snppan)!=0){
        Sigpan = data.frame(rep(chrk, length(sigleSNP)),paste0(chrk*10, max(dapan$BlockID)+seq(1:length(sigleSNP))), sigleSNP, sigleSNP, rep(1, length(sigleSNP)))
        colnames(Sigpan) = colnames(dapan)
        dapan = rbind(dapan, Sigpan)
      }else{
        Sigpan = data.frame(rep(chrk, length(sigleSNP)),paste0(chrk*10, seq(1:length(sigleSNP))), sigleSNP, sigleSNP, rep(1, length(sigleSNP)))
        colnames(Sigpan) = colnames(dapan);
        dapan = rbind(dapan, Sigpan);
      }
      
    }
    
    # move row 345-347 to here, as dapan may expand in line 349-355, ZX
    snp = union(dapan$SNP1, dapan$SNP2)
    daexp_k = daexp_k[daexp_k$SNP %in% snp, ]
    daout_k = daout_k[daexp_k$SNP,]
    
    
    arma_blc_id = unique(dapan$BlockID)
    
    nid = length(arma_blc_id)
    F4gamma   = vector("list", nid)
    F4gammase = vector("list", nid)
    F4Gamma   = vector("list", nid)
    F4Gammase = vector("list", nid)
    # F4SNPchr  = vector("list", nid)
    # F4R       = vector("list", nid)
    snpset = rep(FALSE, nid)
    for (k in 1:nid) {
      uid = arma_blc_id[k]
      da2 = dapan[dapan$BlockID==uid,]
      comsnp = union(da2$SNP1, da2$SNP2)
      
      daRef = RefSNP_k[RefSNP_k$SNP%in%comsnp, ];
      da0 = daexp_k[comsnp, ]; 
      da1 = daout_k[comsnp, ];
      # ===================================================== #
      # match Ref panel data.
      da0$beta[da0$a1!=daRef$a1] = -1*da0$beta[da0$a1!=daRef$a1];
      da1$beta[da0$a1!=daRef$a1] = -1*da1$beta[da0$a1!=daRef$a1];
      # ===================================================== #
      
      M = diag(rep(1,length(comsnp)))
      rownames(M) = comsnp
      colnames(M) = comsnp
      for (j in 1:nrow(da2)) {
        M[da2$SNP1[j], da2$SNP2[j]] = da2$r[j]*lambad
        M[da2$SNP2[j], da2$SNP1[j]] = da2$r[j]*lambad
      }
      rownames(M) = NULL
      colnames(M) = NULL
      
      id = LDclump(M*1.0, ld_r2_thresh) + 1
      
      F4gamma[[k]] = da0$beta[id]
      F4gammase[[k]] = da0$se[id]
      F4Gamma[[k]] = da1$beta[id]
      F4Gammase[[k]] = da1$se[id]
      # F4SNPchr[[k]] = comsnp[id]
      # F4R[[k]] = M[id,id]
    }
    
    ResF4gamma    = append(ResF4gamma, F4gamma)
    ResF4gammase  = append(ResF4gammase, F4gammase)
    ResF4Gamma    = append(ResF4Gamma, F4Gamma)
    ResF4Gammase = append(ResF4Gammase, F4Gammase)
    # ResF4SNPchr   = append(ResF4SNPchr, F4SNPchr)
    # ResF4R        = append(ResF4R, F4R)
    # Res_id        = c(Res_id, arma_blc_id)
  }
  
  try(
    rm(daexp, daout, dapan, daexp_k, daout_k, da0, da1, da2)
  )
  
  bh1_ind = unlist(ResF4gamma)
  bh2_ind = unlist(ResF4Gamma)
  se1_ind = unlist(ResF4gammase)
  se2_ind = unlist(ResF4Gammase)
  
  z1_ind = bh1_ind / se1_ind;
  z2_ind = bh2_ind / se2_ind;
  
  maxIter = 4000;
  thin = 10;
  burnin = 1000;
  
  nsave = maxIter / thin;
  
  if(length(pth)==1){
    Rhores = rep(NA, nsave);
    a = rep(-pth, 2);
    b = rep(pth, 2);
    # z1_new may has length 0, ZX
    z1_new = z1_ind [which(abs(z1_ind) < pth&abs(z2_ind) < pth)];
    z2_new = z2_ind[which(abs(z1_ind) < pth&abs(z2_ind) < pth)];
    # ZX
    if (length(z1_new) == 0) {
      stop("No significant SNPs were found and the program terminated.")
    }
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
      # ZX
      if (length(z1_new) == 0) {
        warning(paste0("No significant SNPs were found with pth = ", pth, " and the program skip it."))
        next
      }
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
