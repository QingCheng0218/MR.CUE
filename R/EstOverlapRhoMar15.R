


EstRhoNewfun <- function(Namesave, fileexposure, fileoutcome, stringname3,
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
  # In case that the number of SNP is too large.
  if(length(bh1new) > 1000000){
    p = length(bh1new);
    idxrand = sort(sample(1:p, 1000000, replace = FALSE));
    
    bh1new1 = bh1new[idxrand];
    bh2new1 = bh2new[idxrand];
    s12new1 = s12new[idxrand];
    s22new1 = s22new[idxrand];
    bpnew1 = bpnew[idxrand];
    chrnew1 = chrnew[idxrand];
    avbIndexnew1 = avbIndexnew[idxrand];
    rsnamenew1 = rsnamenew[idxrand];
    
    # generate new panel data
    write.table(rsnamenew1,  row.names = FALSE,
                col.names = FALSE,
                quote = FALSE,
                file = paste0(Namesave,"RhoUK10K_match.snplist"));
    
  }else{
    bh1new1 = bh1new;
    bh2new1 = bh2new;
    s12new1 = s12new;
    s22new1 = s22new;
    bpnew1 = bpnew;
    chrnew1 = chrnew;
    avbIndexnew1 = avbIndexnew;
    rsnamenew1 = rsnamenew;
    
    # generate new panel data
    write.table(rsnamenew1,  row.names = FALSE,
                col.names = FALSE,
                quote = FALSE,
                file = paste0(Namesave,"RhoUK10K_match.snplist"));
  }
  
  cmd = paste("plink --bfile", stringname3,
              "--extract ", paste0(Namesave,"RhoUK10K_match.snplist"),
              "--make-bed",
              "--out", paste0(Namesave,"RhoUK10K_match"))
  
  system(cmd);
  
  newpanel = paste0(Namesave,"RhoUK10K_match");
  Index = 1:length(bpnew1) - 1;
  IndSumRes = IndepSummary(bpnew1, chrnew1, Index, block_file, newpanel, 
                           bh1new1, bh2new1, s12new1, s22new1, coreNum, 
                           lam, ld_r2_thresh);
  
  system(paste0("rm -rf ", Namesave, "RhoUK10K_match.*"));
  #-------------------------------------------#
  
  
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

