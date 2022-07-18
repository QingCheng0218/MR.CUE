rm(list = ls());
source("Coreplotfun.R")
library(ggplot2);
# ************************************************************ #
# All the results can be download from: https://drive.google.com/drive/folders/1eeXxxp2H4Nv9EJJNUh9_l22-zOYgEvl3
# ************************************************************ #
# ************************************************************ #
# Figure 2 
# ************************************************************ #
da = load("Figure2results.Rdata");
MethOrder <- MethOrderlabel <- c("MR-CUE", "CAUSE","GRAPPLE" ,"cML-MA-BIC-DP", 
                                 "RAPS" , "IVW", "MR-Egger", "MRMix","MR-Clust",  "MR-LDP");
color = c("#d81f2a","#232121","#0033a1","#ff9900","#9ea900",
          "#6ec9e0","#007ea3","#631d76","#9e4770","#8d6cab");

# ====================================== #
# Figure 2a and 2c  for single confounder
# ====================================== #
Fig2a_SingleUTypeIE = TypeIfun4main(TypeIE4SingleU, color, alvalue = 0.7, MethOrder, MethOrderlabel)
Fig2c_SingleUpower = powerplotfun(Power4SingleU, color, MethOrder, MethOrderlabel);

# ====================================== #
# Figure 2b and 2d for multiple confounders
# ====================================== #
Fig2b_MultiUTypeIE = TypeIfun4main(TypeIE4MultipleU, color, alvalue = 0.7, MethOrder, MethOrderlabel);
Fig2d_MultiUpower = powerplotfun(Power4MultipleU, color, MethOrder, MethOrderlabel);

# ====================================== #
# Figure 2e and 2f for Negative Controls and Positive Controls (with existing methods).
# ====================================== #
Fig2e_NegCon4Meth10 = Negfigurefun(NegMeth105e4, color);
Fig2f_PosCon4Meth10 = Posfigurefun(PosMeth105e4, color);

# ====================================== #
# Figure 2g for Positive Controls (with MR-CUE-Indep).
# ====================================== #
Fig2g_PosCon4CorrIndep = PosfigureCorrIndepfun(Neg4CorrIndep);

# ====================================== #
# Figure 2h for Reverse Causation.
# ====================================== #
library(plotROC);
Fig2h_RC = RCplotfun(RCres, color);


# ************************************************************ #
# Figure 3 
# ************************************************************ #
# ====================================== #
# Figure 3a (left panel) for Genetic correlations.
# ====================================== #
da = load("Figure3results4IL6.Rdata");
library(reshape2)
Fig3aL_IL6GenCorr = IL6GenCorrfun(lower_tri_IL6GenCorr, lower_pva_IL6GenCorr);

# ====================================== #
# Figure 3a (right panel) for CHP.
# ====================================== #
library(ComplexHeatmap);
library(circlize);
Fig3aR_IL6hm = IL6heatmapplotfun(bhat4IL6, pval4IL6, PvaMatfix4IL6, chridx4IL6);

# ====================================== #
# Figure 3b for cis-genes.
# ====================================== #
Fig3b_IL6cisgenes = ILcisgenesplotfun(Chrcisgenes, Pvacisgenes, cisgenes, xcolorcisgenes, 5, 6)

# ************************************************************ #
# Figure 4 
# ************************************************************ #
da = load("Figure4results4T2D.Rdata");
# ====================================== #
# Figure 4a: Heatmap for T2D in European.
# ====================================== #
Fig4a_hmEur = T2DEurhmplotfun(bhat4T2DEur, log10pval4T2DEur, PvaMatfix4T2DEur, chridx4T2DEur);

# ====================================== #
# Figure 4b: Heatmap for T2D in East Asian.
# ====================================== #
Fig4b_hmAsian = T2DAisanhmplotfun(bhat4T2DAsian, log10pval4T2DAsian, PvaMatfix4T2DAsian, chridx4T2DAsian)

# ====================================== #
# Figure 4c: Pathway for T2D in European.
# ====================================== #
Fig4c_pathwayEur = T2DEurpathwayfun(daEur, daAsian, T2DEurdata, T2DEurPVname, T2DEurPvaMat, UnionPath)

# ====================================== #
# Figure 4d: Pathway for T2D in East Asian.
# ====================================== #
Fig4d_pathwayAsian = T2DAsianpathwayfun(daEur, daAsian, T2DAsiandata, T2DAsianPVName)
