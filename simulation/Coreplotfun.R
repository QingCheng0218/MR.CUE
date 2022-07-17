# ========================================= #
# Plot functions used in Simulation Studies.
# ========================================= #
TypeIfun4main <- function(dat1, color, alvalue, MethOrder, MethOrderlabel){
  datplot <- data.frame(
    Method =  ordered(dat1[, "Method"], 
                      levels = MethOrder), 
    bhat = as.numeric(dat1[, "bhat"]),
    power = as.numeric(dat1[, "power"]),
    h2a = factor(dat1[, "H2a"], labels = c(expression(h[alpha]^2==0.05),  expression(h[alpha]^2==0.1))),
    h2t = factor(dat1[, "H2t"], labels = c(expression(h[theta]^2==0.02), expression(h[theta]^2==0.05))),
    # p = factor(dat1[, "p"]),
    Alrate = factor(dat1[, "Alrate"]),
    Rho = factor(dat1[, "Rho"]),
    # hhz = factor(dat1[, "H2z"], labels = c(expression(h[alpha]^2==0.05), expression(h[alpha]^2==0.1))),
    rho_ag = factor(dat1[, "rho_ag"]));
  
  p1 <- ggplot(data=datplot, aes(x=Rho, y=power, fill = Method)) +
    geom_bar(stat="identity", position=position_dodge())+
    facet_grid(h2a~h2t, labeller = label_parsed, scales = "free_y")+
    # geom_hline(yintercept=0.035,colour="red",linetype="dashed")+
    geom_hline(yintercept=0.05,colour="red",linetype="dashed")
  # p1 = p1 + scale_y_continuous(breaks=seq(0,1,0.05))
  
  p1 <- p1 + theme_bw() + theme(panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                plot.title = element_text(size=23,face = "bold",hjust=0.5),
                                axis.title.x = element_text(size=23,face = "bold"),
                                strip.text.x = element_text(size = 23,face = "bold"),
                                strip.text.y = element_text(size = 23,face = "bold"),
                                axis.text.x = element_text(size=23,face = "bold"),
                                axis.title.y = element_text(size=25.5,face = "bold"),
                                axis.text.y = element_text(size=20,face = "bold"),
                                legend.position="right",
                                legend.title=element_text(size=20,face = "bold"),
                                legend.text=element_text(size=20,face = "bold"));
  p1 = p1 + labs(x = expression(r), y = "Type I error rate");
  
  
  p1 <- p1 + scale_fill_manual( values= alpha(color, alvalue), labels = MethOrderlabel)
  
  
  
  p1 = p1 + geom_vline(xintercept=1+.5,color="grey", alpha = 0.5, size = 1);
  
  return(p1);
}

powerplotfun <- function(dat, color, MethOrder, MethOrderlabel){
  dat[which(dat[, "h2gIndex"]=="0"), "h2gIndex"] = -4;
  
  dapower = data.frame(
    Method = ordered(dat[, "Method"], levels = MethOrder),
    Power = as.numeric(dat[,"power"]),
    h2g = as.numeric(dat[, "h2gIndex"])
  );
  
  
  p1 <- ggplot(dapower, aes(x = h2g, y = Power, color = Method, shape = Method))+ 
    geom_line(aes(linetype= Method),size = 1)+geom_point(size = 4)
  
  p1 = p1 + scale_shape_manual(values=c(17,16,15, 18, 1,2,3,7,8, 10));
  
  
  # color = c("#FF0000","#0000FF","darkmagenta", "darkslateblue", "#33638DFF", "#228B22", "#599119", "#909811", "#C79E08" ,"#FFA500");
  p1 = p1 + scale_color_manual(values=color,labels = MethOrderlabel )+
    scale_shape_manual(values=c(17,16,15, 18, 1,2,3,7,8, 10),labels = MethOrderlabel)
  
  p1 = p1 + scale_x_continuous(breaks = seq(-4, -1, length.out = 7),
                               labels = c(0,expression(bold(10^{-3.5})),expression(bold(10^{-3}))
                                          ,expression(bold(10^{-2.5})),expression(bold(10^{-2}))
                                          ,expression(bold(10^{-1.5})),expression(bold(10^{-1}))));
  p1 = p1 + scale_y_continuous(breaks=c(0.0, 0.05,0.25, 0.5, 0.75, 1), limits = c(-0.02, 1));
  
  lt = c("solid", "F1", "dashed", "12345678", "4C88C488",  "dotdash", "longdash", "twodash", "4C88C488",   "12345678",   "dotted")
  p1 = p1 +  scale_linetype_manual( labels = MethOrderlabel,
                                    values = lt)
  
  ycol = c("black", "red", rep("black", 4))
  
  
  
  p1 <- p1+theme_bw() + theme(legend.key.width= unit(1.5, 'cm'),
                              panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              plot.title = element_text(size=20,face = "bold",hjust=0.5),
                              axis.title.x = element_text(size=20,face = "bold"),
                              strip.text.x = element_text(size = 20,face = "bold"),
                              strip.text.y = element_text(size = 20,face = "bold"),
                              axis.text.x = element_text(size=20,face = "bold"),
                              axis.title.y = element_text(size=25,face = "bold"),
                              axis.text.y = element_text(size=20,face = "bold", color = ycol),
                              legend.position="right",
                              legend.title=element_text(size=19,face = "bold"),
                              legend.text=element_text(size=19,face = "bold"));
  
  p1 = p1 + annotate("rect", xmin = -4.1, xmax = -3.9, ymin = -0.02, ymax = 1,
                     alpha = .1,fill = "mediumblue");
  
  p1 = p1 + labs(x = expression(bold(h[gamma]^2)), y = "Power/Type I error rate");
  
  p1 = p1 + guides(fill=guide_legend(nrow=2,byrow=TRUE))
  
  return(p1);
}

Posfigurefun <- function(dat, color){
  plotdata = data.frame(
    Method = ordered(dat[, "Method"], levels = c("MR-CUE", "CAUSE","GRAPPLE" ,"cML-MA-BIC-DP", 
                                                 "RAPS" , "IVW", "MR-Egger", "MRMix","MR-Clust",  "MR-LDP")), 
    observed = as.numeric(dat[,"observed"]),
    expected = as.numeric(dat[,"expected"]),
    clower = as.numeric(dat[, "clower"]),
    cupper = as.numeric(dat[, "cupper"])
  )
  
  p1 = ggplot(plotdata, aes(x = expected, y = observed, shape = Method))+ 
    geom_ribbon(aes(ymax = cupper, ymin = clower), fill = "grey70", alpha = 0.1)+ geom_point(aes(color = Method), size = 3, alpha = 0.9)
  p1 = p1 + geom_abline(slope=1, intercept=0);
  # color = c("#FF0000","#0000FF", "#33638DFF", "#228B22", "#599119", "#909811", "#C79E08" ,"#FFA500");
  # color = c("#FF0000","#0000FF","darkmagenta", "darkslateblue", "#33638DFF", "#228B22", "#599119", "#909811", "#C79E08" ,"#FFA500");
  
  # color = c("#FF0000", "deeppink1","darkmagenta","#0000FF", 
  #           "#33638DFF", "#537AC2",  "steelblue2",
  #           "#B6D5B9","#FFD200","darkorange3");
  
  p1 = p1 + scale_color_manual(values=color,labels = c("MR-CUE", "CAUSE","GRAPPLE" ,"cML-MA-BIC-DP", 
                                                       "RAPS" , "IVW", "MR-Egger", "MRMix","MR-Clust",  "MR-LDP"))+
    scale_shape_manual(values=c(17,16,15, 18, 11,2,13,7,8, 10),labels =c("MR-CUE", "CAUSE","GRAPPLE" ,"cML-MA-BIC-DP", 
                                                                         "RAPS" , "IVW", "MR-Egger", "MRMix","MR-Clust",  "MR-LDP"));
  
  # p1 = p1 + scale_color_manual(values=color)+ scale_shape_manual(values=c(17,0,1,2,3,7,8, 10));
  p1 = p1  + xlab(expression(bold(paste("Expected (",-log[10], " p-value)")))) + ylab(expression(bold(paste("Observed (",-log[10], " p-value)"))))
  
  p1 <- p1+theme_bw() + theme(panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              plot.title = element_text(size=24,face = "bold",hjust=0.5),
                              axis.title.x = element_text(size=24,face = "bold"),
                              strip.text.x = element_text(size = 24,face = "bold"),
                              strip.text.y = element_text(size = 24,face = "bold"),
                              axis.text.x = element_text(size=24,face = "bold"),
                              axis.title.y = element_text(size=24,face = "bold"),
                              axis.text.y = element_text(size=24,face = "bold"),
                              legend.position="right",
                              legend.title=element_text(size=19,face = "bold"),
                              legend.text=element_text(size=19,face = "bold"));
  
  return(p1); 
  
  
  
}

Negfigurefun <- function(dat, color){
  
  plotdata = data.frame(
    Method = ordered(dat[, "Method"], levels = c("MR-CUE", "CAUSE","GRAPPLE" ,"cML-MA-BIC-DP", 
                                                 "RAPS" , "IVW", "MR-Egger", "MRMix","MR-Clust",  "MR-LDP")), 
    observed = as.numeric(dat[,"observed"]),
    expected = as.numeric(dat[,"expected"]),
    clower = as.numeric(dat[, "clower"]),
    cupper = as.numeric(dat[, "cupper"])
  )
  
  p1 = ggplot(plotdata, aes(x = expected, y = observed, shape = Method))+ 
    geom_ribbon(aes(ymax = cupper, ymin = clower), fill = "grey70", alpha = 0.1)+ geom_point(aes(color = Method), size = 3, alpha = 0.9)
  p1 = p1 + geom_abline(slope=1, intercept=0);
  # color = c("#FF0000","#000000", "#6495ED","#228B22", "#599119", "#909811", "#C79E08" ,"#FFA500");
  
  
  
  p1 = p1 + scale_color_manual(values=color)+ scale_shape_manual(values=c(17,16,15, 18, 11,2,13,7,8, 10));
  p1 = p1 + ylim(c(0, 4)) + xlab(expression(bold(paste("Expected (",-log[10], " p-value)")))) + ylab(expression(bold(paste("Observed (",-log[10], " p-value)"))))
  
  
  p1 <- p1+theme_bw() + theme(panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              plot.title = element_text(size=26,face = "bold",hjust=0.5),
                              axis.title.x = element_text(size=26,face = "bold"),
                              strip.text.x = element_text(size = 26,face = "bold"),
                              strip.text.y = element_text(size = 26,face = "bold"),
                              axis.text.x = element_text(size=26,face = "bold"),
                              axis.title.y = element_text(size=26,face = "bold"),
                              axis.text.y = element_text(size=26,face = "bold"),
                              legend.position="none",
                              legend.title=element_text(size=19,face = "bold"),
                              legend.text=element_text(size=19,face = "bold"));
  
  return(p1) 
  
}

PosfigureCorrIndepfun <- function(dat){
  
  fillshape = c(17, 2, 19, 13, 15, 7,18,5, 20, 10);
  color2 = c( "#d81f2a","#1B2431", "#FE46A5","#048243", "#FDB0C0",
              "#0E87CC","#D46A7E", "#13BBAF","#FDB915", "#214761");
  label = c(bquote(bold(MR~"-"~CUE)~"("~5*x*10^{-4}~")"),bquote(bold(MR~"-"~CUE~"-"~Indep)~"("~5*x*10^{-4}~")"),
            bquote(bold(MR~"-"~CUE)~"("~10^{-4}~")"), bquote(bold(MR~"-"~CUE~"-"~Indep)~"("~10^{-4}~")"),
            bquote(bold(MR~"-"~CUE)~"("~5*x*10^{-5}~")") ,bquote(bold(MR~"-"~CUE~"-"~Indep)~"("~5*x*10^{-5}~")"),
            bquote(bold(MR~"-"~CUE)~"("~10^{-5}~")"), bquote(bold(MR~"-"~CUE~"-"~Indep)~"("~10^{-5}~")"),
            bquote(bold(MR~"-"~CUE)~"("~5*x*10^{-6}~")") ,bquote(bold(MR~"-"~CUE~"-"~Indep)~"("~5*x*10^{-6}~")"));
  
  
  plotdata = data.frame(
    Method = ordered(dat[, "Method"], levels = c("MR-CUE (5e-4)", "MR-CUE-Indep (5e-4)",
                                                 "MR-CUE (1e-4)", "MR-CUE-Indep (1e-4)",
                                                 "MR-CUE (5e-5)" ,"MR-CUE-Indep (5e-5)",
                                                 "MR-CUE (1e-5)", "MR-CUE-Indep (1e-5)",
                                                 "MR-CUE (5e-6)" , "MR-CUE-Indep (5e-6)")), 
    observed = as.numeric(dat[,"observed"]),
    expected = as.numeric(dat[,"expected"]),
    clower = as.numeric(dat[, "clower"]),
    cupper = as.numeric(dat[, "cupper"])
  )
  
  p1 = ggplot(plotdata, aes(x = expected, y = observed, shape = Method))+ 
    geom_ribbon(aes(ymax = cupper, ymin = clower), fill = "grey70", alpha = 0.1)+ geom_point(aes(color = Method), size = 3, alpha = 0.9)
  p1 = p1 + geom_abline(slope=1, intercept=0);
  # color = c("#FF0000","#0000FF", "#33638DFF", "#228B22", "#599119", "#909811", "#C79E08" ,"#FFA500");
  # color = c("#FF0000","#0000FF","darkmagenta", "darkslateblue", "#33638DFF", "#228B22", "#599119", "#909811", "#C79E08" ,"#FFA500");
  
  # color = c("#FF0000", "deeppink1","darkmagenta","#0000FF", 
  #           "#33638DFF", "#537AC2",  "steelblue2",
  #           "#B6D5B9","#FFD200","darkorange3");
  
  p1 = p1 + scale_color_manual(values=color,labels = label)+
    scale_shape_manual(values=fillshape,labels = label);
  
  # p1 = p1 + scale_color_manual(values=color)+ scale_shape_manual(values=c(17,0,1,2,3,7,8, 10));
  p1 = p1  + xlab(expression(bold(paste("Expected (",-log[10], " p-value)")))) + ylab(expression(bold(paste("Observed (",-log[10], " p-value)"))))
  
  p1 <- p1+theme_bw() + theme(panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              plot.title = element_text(size=23,face = "bold",hjust=0.5),
                              axis.title.x = element_text(size=23,face = "bold"),
                              strip.text.x = element_text(size = 23,face = "bold"),
                              strip.text.y = element_text(size = 23,face = "bold"),
                              axis.text.x = element_text(size=23,face = "bold"),
                              axis.title.y = element_text(size=23,face = "bold"),
                              axis.text.y = element_text(size=23,face = "bold"),
                              legend.position="right",
                              legend.title=element_text(size=18,face = "bold"),
                              legend.text=element_text(size=18,face = "bold"));
  
  p1 = p1+scale_x_continuous(limits = c(0, 2))
  p1 = p1 +guides(fill=guide_legend(nrow=2,byrow=TRUE))
  return(p1); 
  
  
  
}

RCplotfun <- function(RCres, color){
  modellevel = c("MR-CUE", "CAUSE","GRAPPLE" ,"cML-MA-BIC-DP", 
                 "RAPS" , "IVW", "MR-Egger", "MRMix","MR-Clust",  "MR-LDP");
  
  
  
  
  colnames(RCres) = c("Method", "pva", "btrue")
  # colnames(res1) = c("Method", "pva", "btrue")
  pva.df <- data.frame(Method = ordered(RCres[, "Method"], 
                                        levels = modellevel),
                       pva = as.numeric(RCres[, "pva"]),
                       Btrue = as.numeric(RCres[, "btrue"]));
  
  p1 = ggplot(pva.df, aes(m = pva, d = Btrue, col = Method, linetype =Method, shape = Method)) +
    geom_roc(n.cuts = 10, linealpha = 0.9, size = 1, pointsize =1, labels = FALSE);
  
  
  
  color = c(
    "#d81f2a","#232121","#0033a1","#ff9900","#9ea900",
    "#6ec9e0","#007ea3","#631d76","#9e4770","#8d6cab"
  )
  p1 = p1 + scale_color_manual(values=color,labels = modellevel)+
    scale_shape_manual(values=c(17,16,15, 18, 11,2,13,7,8, 10),labels = modellevel)
  
  lt = c("solid", "F1", "dashed", "12345678", "4C88C488",  "dotdash", "longdash", "twodash", "4C88C488",   "12345678",   "dotted")
  p1 = p1 +  scale_linetype_manual( labels = modellevel,
                                    values = lt)
  
  p1 <- p1+theme_bw() + theme(legend.key.width= unit(1.5, 'cm'),
                              panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              plot.title = element_text(size=23,face = "bold",hjust=0.5),
                              axis.title.x = element_text(size=24,face = "bold"),
                              strip.text.x = element_text(size = 23,face = "bold"),
                              strip.text.y = element_text(size = 23,face = "bold"),
                              axis.text.x = element_text(size=23,face = "bold"),
                              axis.title.y = element_text(size=24,face = "bold"),
                              axis.text.y = element_text(size=23,face = "bold"),
                              legend.position="right",
                              legend.title=element_text(size=19,face = "bold"),
                              legend.text=element_text(size=19,face = "bold"));
  p1 = p1 + labs(x = "False positive rate", y = "True positive rate");
  
  
  p1 = p1+ guides(colour = guide_legend(override.aes = list(size=2)))
  
  return(p1);
  
}
# ========================================= #
# Plot functions used in IL-6 Analysis.
# ========================================= #
IL6GenCorrfun <- function(lower_tri, lower_pva){
  melted_cormat <- melt(lower_tri, na.rm = TRUE);
  melted_pva <- melt(lower_pva, na.rm = TRUE);
  dim(melted_cormat);
  dim(melted_pva)
  melted_cormat[which(melted_cormat$value<0), "value"] = NA;
  
  
  # check the order for genetic correlation and pvalue
  sum(melted_cormat$Var1==melted_pva$Var1)
  sum(melted_cormat$Var2==melted_pva$Var2)
  dim(melted_cormat)
  which(is.na(melted_cormat$value))
  # construct the plot data.
  plot.data <- melted_cormat
  plot.data$pvalue = melted_pva$value;
  melted_pva[which(melted_pva$value<0.05), ];
  melted_cormat[which(melted_pva$value<0.05), ];
  
  plot.data$stars <- cut(plot.data$pvalue, breaks=c(-Inf,0.05, Inf), 
                         label=c("*", ""))  # Create column of significance labels
  
  
  
  plow = ggplot(data = plot.data, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "white", high = "steelblue",, limit = c(0, 1), space = "Lab", na.value = "grey",
                         name=expression(paste("|", rho[g], "|"))) +  
    geom_text(aes(label=stars), color="black", size=5)+
    scale_y_discrete(expand = c(0, 0), position = "left")+
    scale_x_discrete(expand = c(0, 0), position = "top")+
    theme_bw()+ 
    theme(  axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            panel.grid.major = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.ticks = element_blank(),
            legend.position = "left",
            legend.title = element_text(angle = -90),
            legend.text = element_text(angle = -90),
            legend.key.height = unit(1, 'cm'), 
            axis.text.y.left = element_text(vjust=0.5),
            axis.text.x.top = element_text(angle = 270, vjust = 0.5,hjust=1)) +coord_fixed()
  
  
  
  return(plow)
}

IL6heatmapplotfun <- function(bhat, pva1, PvaMatfix, chridx){
  pch = rep(NA, length(bhat));
  # is_sig = which(pva<0.05);
  pch[which(sign(bhat)==-1)] = 6;
  pch[which(sign(bhat)==1)] = 2;
  
  # col_fun1 = colorRamp2(c(0, 5, 10, 15, 20),brewer.pal(5, "Greys"));
  col_fun1 = colorRamp2(c(1, 2, 3),c("grey", "#CDCD00", "#CD6600"));
  
  
  col_fun2 = colorRamp2(c(0, 0.3, 0.6, 0.9), c( "#F7FCFD","#3F3FD9", "#7B7EE5", "#0000CD"))
  
  # 
  # col_fun2 = colorRamp2(seq(0, .9, length.out = length(bhat)), 
  #                       colorRampPalette(brewer.pal(9,"BuPu"))(length(bhat)))
  row_ha = rowAnnotation(foo2 = anno_simple(pva1, col = col_fun1, 
                                            pt_gp = gpar(col = "#0000CD"),pch = pch,
                                            pt_size = unit(2, "mm")),
                         annotation_label = expression(italic('p') ["-value"]))
  
  #-----------------------------------------------------------------------------#
  # legend for pvalue
  lgd_pvalue = Legend(title = expression(italic('p') ["-value"]),
                      labels = c(">0.05", "<0.05", "<0.001"),
                      at = c(1, 2, 3),
                      legend_gp = gpar(fill =  c("grey", "#CDCD00", "#CD6600")),
                      labels_gp = gpar(fontsize = 15), title_gp = gpar(footsize = 25))
  # lgd_pvalue = Legend(at = 1:3, title = "foo", legend_gp = gpar(fill = 1:3, col = c("grey", "blue", "green")))
  # lgd_sig = Legend(pch = 8, legend_gp = gpar(col = "red"), type = "points",
  #                  labels =  expression(italic('p') ["-value"]<0.05));
  
  
  lgd_bhatsign = Legend(labels = c( "Positive", "Negative"), title = expression(hat(beta)), 
                        type = "points", 
                        pch = c(2, 6), legend_gp = gpar(col = c("#0000CD", "#0000CD")), background = "white",
                        labels_gp = gpar(fontsize = 15), title_gp = gpar(footsize = 25))
  
  
  
  ht = Heatmap(PvaMatfix,column_split = factor(chridx, levels = paste0("chr",1:22)), column_title_rot = 90,
               column_title_gp = gpar( fontsize = 5),
               heatmap_legend_param = 
                 list(title = expression(bold(paste(-log[10],"(",1-"Pr("*eta==1*")",")"))),
                      at = c(0,  0.3, 0.6, 0.9), 
                      labels = c("0", "0.3", "0.6", "0.9"),
                      direction = "horizontal",
                      legend_width = unit(6, "cm"),
                      labels_gp = gpar(fontsize = 15),
                      title_gp = gpar(footsize = 25),
                      title_position = "leftcenter"),
               cluster_rows = FALSE, 
               cluster_column_slices = FALSE,
               na_col = "white",row_names_side = "left",
               col = col_fun2,right_annotation = row_ha);
  
  
  ht1 = draw(ht, annotation_legend_list = list(lgd_pvalue, lgd_bhatsign), 
             heatmap_legend_side = "bottom", annotation_legend_side = "right")
  
  return(ht1);
  
}

ILcisgenesplotfun <- function(ChrRes, PvaminRes, genesymbol, xcolor, xsize,ysize){
  melted_Num <- melt(PvaminRes);
  chrinfo3 = rep(NA, dim(melted_Num)[1]);
  chrinfo3 =  ChrRes[match(melted_Num$Var1, genesymbol)]
  melted_Num$Chr = chrinfo3;
  
  
  genchr3= paste0("chr", ChrRes,":", genesymbol)
  melted_Num$genNloc = paste0("chr",chrinfo3, ":",  melted_Num$Var1);
  melted_Num$genNloc <- factor(melted_Num$genNloc, levels = genchr3)
  
  
  p2 = ggplot(data = melted_Num, aes(genNloc, Var2, fill = value))+
    geom_tile(color = "white")+labs(y = "Traits", x = "Genes")+
    scale_fill_gradient2(low = "white", high = "red", mid = "#FFE2E2", 
                         space = "Lab", na.value = "Gray",
                         name="-log10(p-value)") +
    theme_bw()+ 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = xsize, hjust = 1, face="bold",
                                     color = xcolor),
          axis.text.y=element_text(size = ysize,face="bold"))+ coord_fixed();
  return(p2);
}

# ========================================= #
# Plot functions used in T2D Analysis.
# ========================================= #
T2DEurhmplotfun <- function(bhat, log10pva, PvaMatfix, chridx){
  pch = rep(NA, length(bhat));
  # is_sig = which(pva<0.05);
  pch[which(sign(bhat)==-1)] = 6;
  pch[which(sign(bhat)==1)] = 2;
  
  # col_fun1 = colorRamp2(c(0, 5, 10, 15, 20),brewer.pal(5, "Greys"));
  col_fun1 = colorRamp2(c(0, 5, 10, 15, 20),c("#FFFFFF", "#F2F2BF", "#E6E67F", "#D9D93F", "#CDCD00"));
  col_fun2 = colorRamp2(c(0, 0.3, 0.6, 0.9), c( "#F7FCFD", "#7F7FE6", "#3F3FD9", "#0000CD"))
  
  # 
  # col_fun2 = colorRamp2(seq(0, .9, length.out = length(bhat)), 
  #                       colorRampPalette(brewer.pal(9,"BuPu"))(length(bhat)))
  row_ha = rowAnnotation(foo2 = anno_simple(log10pva, col = col_fun1, 
                                            pt_gp = gpar(col = "red"),pch = pch,
                                            pt_size = unit(2, "mm")),
                         annotation_label = expression(italic('p') ["-value"]))
  
  ########################################
  # legend
  lgd_pvalue = Legend(title = expression(italic('p') ["-value"]), col_fun = col_fun1, at = c(0, 5, 10, 15, 20), 
                      labels = c("1", "1e-5", "1e-10", "1e-15", "1e-20"),
                      labels_gp = gpar(fontsize = 15), title_gp = gpar(footsize = 25))
  
  # lgd_sig = Legend(pch = 8, legend_gp = gpar(col = "red"), type = "points",
  #                  labels =  expression(italic('p') ["-value"]<0.05));
  
  
  lgd_bhatsign = Legend(labels = c( "Positive", "Negative"), title = expression(hat(beta)), 
                        type = "points", 
                        pch = c(2, 6), legend_gp = gpar(col = c("red", "red")), background = "white",
                        labels_gp = gpar(fontsize = 15), title_gp = gpar(footsize = 25))
  
  
  
  ht = Heatmap(PvaMatfix,column_split = factor(chridx, levels = paste0("chr",1:22)), column_title_rot = 90,
               column_title_gp = gpar( fontsize = 5),
               heatmap_legend_param = 
                 list(title = expression(bold(paste(-log[10],"(",1-"Pr("*eta==1*")",")"))),
                      at = c(0,  0.3, 0.6, 0.9), 
                      labels = c("0", "0.3", "0.6", "0.9"),
                      direction = "horizontal",
                      legend_width = unit(6, "cm"),
                      labels_gp = gpar(fontsize = 15),
                      title_gp = gpar(footsize = 25),
                      title_position = "leftcenter"),
               cluster_rows = FALSE, 
               cluster_column_slices = FALSE,
               na_col = "white",row_names_side = "left",
               col = col_fun2,right_annotation = row_ha);
  
  
  ht1 = draw(ht, annotation_legend_list = list(lgd_pvalue, lgd_bhatsign), 
             heatmap_legend_side = "bottom", annotation_legend_side = "right")
  
  return(ht1);
  
}

T2DAisanhmplotfun <- function(bhat, log10pva, PvaMatfix, chridx){
  pch = rep(NA, length(bhat));
  # is_sig = which(pva<0.05);
  pch[which(sign(bhat)==-1)] = 6;
  pch[which(sign(bhat)==1)] = 2;
  
  # col_fun1 = colorRamp2(c(0, 5, 10, 15, 20),brewer.pal(5, "Greys"));
  col_fun1 = colorRamp2(c(-1, 0, 5, 10, 15, 20),c("grey", "#FFFFFF", "#F2F2BF", "#E6E67F", "#D9D93F", "#CDCD00"));
  col_fun2 = colorRamp2(c(0, 0.3, 0.6, 0.9), c("#F7FCFD", "#7F7FE6", "#3F3FD9", "#0000CD"))
  
  # 
  # col_fun2 = colorRamp2(seq(0, .9, length.out = length(bhat)), 
  #                       colorRampPalette(brewer.pal(9,"BuPu"))(length(bhat)))
  row_ha = rowAnnotation(foo2 = anno_simple(log10pva, col = col_fun1, 
                                            pt_gp = gpar(col = "red"),pch = pch,
                                            pt_size = unit(2, "mm")),
                         annotation_label = expression(italic('p') ["-value"]))
  
  ########################################
  # legend
  lgd_pvalue = Legend(title = expression(italic('p') ["-value"]), col_fun = col_fun1, at = c( 0, 5, 10, 15, 20), 
                      labels = c("1", "1e-5", "1e-10", "1e-15", "1e-20"),
                      labels_gp = gpar(fontsize = 15), title_gp = gpar(footsize = 25))
  
  # lgd_sig = Legend(pch = 8, legend_gp = gpar(col = "red"), type = "points",
  #                  labels =  expression(italic('p') ["-value"]<0.05));
  
  
  lgd_bhatsign = Legend(labels = c( "Positive", "Negative"), title = expression(hat(beta)), 
                        type = "points", 
                        pch = c(2, 6), legend_gp = gpar(col = c("red", "red")), background = "white",
                        labels_gp = gpar(fontsize = 15), title_gp = gpar(footsize = 25))
  
  
  
  ht = Heatmap(PvaMatfix,column_split = factor(chridx, levels = paste0("chr",1:22)), column_title_rot = 90,
               column_title_gp = gpar( fontsize = 5),
               heatmap_legend_param = 
                 list(title = expression(bold(paste(-log[10],"(",1-"Pr("*eta==1*")",")"))),
                      at = c(0,  0.3, 0.6, 0.9), 
                      labels = c("0", "0.3", "0.6", "0.9"),
                      direction = "horizontal",
                      legend_width = unit(6, "cm"),
                      labels_gp = gpar(fontsize = 15),
                      title_gp = gpar(footsize = 25),
                      title_position = "leftcenter"),
               cluster_rows = FALSE, 
               cluster_column_slices = FALSE,
               na_col = "white",row_names_side = "left",
               col = col_fun2,right_annotation = row_ha);
  
  
  ht1 = draw(ht, annotation_legend_list = list(lgd_pvalue, lgd_bhatsign), 
             heatmap_legend_side = "bottom", annotation_legend_side = "right")
  
  return(ht1);
  
}

T2DEurpathwayfun <- function(daEur, daAsian, data.plot, PVName, PvaMat, UnionPath){

  pathcom = intersect(as.vector(unique(daAsian$Pathway)), as.vector(unique(daEur$Pathway)));
  Eurpath = unique(daEur$Pathway)
  idxRedEur = which(!is.na(match(Eurpath[order(Eurpath)], pathcom)))
  # length(as.vector(unique(daEur$Pathway)))
  ycol = rep("black", length(as.vector(unique(daEur$Pathway))))
  ycol[1:sum(!is.na(match(Eurpath, pathcom)))] = "blue"
  
  Eurpathorder = c(sort(pathcom), as.vector(sort(Eurpath[which(is.na(match(Eurpath, pathcom)))])))
  
  NewName1 = c("HIP",  "WC", "WHR", "BMI", "BFP",
               "BL", "BW", "HR", "HF","HDL-C",  "TG", 
               "SBP", "DBP", "PP", "AIS", 
               "MDD", "UA",
               "RBC", "WBC", "PLT","Glu", "Gly","DHA","OA",
               "ISI","FG","FI","IR","HbA1c");
  
  dataorderEur = data.frame(
    Pathway = ordered(daEur$Pathway, levels = Eurpathorder),
    variable = ordered(daEur$newVar, levels = NewName1),
    value = daEur$value
  )
  p1 <-  ggplot(dataorderEur, aes(x =variable,y = Pathway)) +
    geom_tile(aes(fill = value),color="black",size=0.1) +
    scale_fill_gradient2(low = "white",high = "red",mid =  "white", midpoint = 1.5, na.value = "Gray",limit = c(1, 6),breaks = c(2, 4,6),
                         guide = guide_colorbar(frame.colour = "black", ticks.colour = "white"))+
    scale_y_discrete(expand = c(0,0))+
    scale_x_discrete(expand = c(0,0),position = "top")+
    theme_bw()+
    theme(axis.text.x.top = element_text(angle = 90, vjust = 0.5, hjust=0),
          axis.text.x=element_text(size=20,face="bold"), 
          axis.text.y=element_text(size=20,face="bold", color = ycol),  
          axis.title.x = element_text(size=25,face = "bold"),
          axis.title.y = element_text(size=25,face = "bold"),
          # axis.text.y = element_text(size=25,face = "bold", color = ycol),
          # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
          legend.title=element_text(size=25),
          legend.position = "none",
          text = element_text(size=20))+
    # rotate_x_text(180) +
    # rotate_y_text(20)+
    xlab("Trait")+
    labs(fill = "-log10(p-value)");
  
  p1 <- p1 + coord_fixed(ratio=1)+guides(fill = guide_colourbar(title = expression(bold(paste(-log[10], "(p-value)")))));
  
  
  
}

T2DAsianpathwayfun <- function(daEur, daAsian, data.plot, PVName){
  pathcom = intersect(as.vector(unique(daAsian$Pathway)), as.vector(unique(daEur$Pathway)));
  
  Aisanpath = unique(daAsian$Pathway)
  idxRedAsian = which(!is.na(match(Aisanpath[order(Aisanpath)], pathcom)))
  
  ycol = rep("black", length(as.vector(unique(daAsian$Pathway))))
  ycol[idxRedAsian] = "blue"
  
  
  p2 <- ggplot(daAsian, aes(x =varAsian,y = Pathway)) +
    geom_tile(aes(fill = value),color="black",size=0.1) +
    scale_fill_gradient2(low = "white",high = "red",mid =  "white", midpoint = 1.5,  na.value = "Gray",limit = c(0, 6),breaks = c(2, 4,6),
                         guide = guide_colorbar(frame.colour = "black", ticks.colour = "white"))+
    scale_y_discrete(expand = c(0,0))+
    scale_x_discrete(expand = c(0,0),position = "top")+
    theme_bw()+
    theme(
      axis.text.x.top = element_text(angle = 90, vjust = 0.5, hjust=0),
      axis.text.x=element_text(size=21,face="bold"), 
      axis.text.y=element_text(size=21,face="bold", color = ycol),  
      axis.title.x = element_text(size=26,face = "bold"),
      axis.title.y = element_text(size=26,face = "bold"),
      # axis.text.y = element_text(size=25,face = "bold", color = ycol),
      # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
      legend.title=element_text(size=26),
      text = element_text(size=21))+
    ylab("")+
    xlab("Trait")+
    labs(fill = "-log10(p-value)");
  
  
  p2 = p2 + coord_fixed(ratio=1)
  p2 = p2 +guides(fill = guide_colourbar(title = expression(bold(paste(-log[10], "(p-value)")))))
  return(p2);
}




