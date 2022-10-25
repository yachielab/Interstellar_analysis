library(Seurat)
library(tidyverse)
library(ggpubr)
#library(VennDiagram)
library(effsize)
library(ggforce)

fnc_draw_scatter <- function(df_rand,
                             colset,
                             outprefix,
                             both_rand_and_keep=TRUE,
                             df_keep=NULL,
                             show_legend=TRUE,
                             w=3.8,
                             h=3.8,
                             show_label=TRUE,
                             mgn=c(.3,.3,.3,.3),
                             legpos="right",
                             ptsize=.1,
                             bin=40,
                             dpi=600,
                             carib_x=FALSE,
                             carib_y=FALSE){
  write.table(paste("RandomUMI:",cor(df_rand$original,df_rand$target,method = "pearson")),paste0(outprefix,"_stat.txt"),quote=F,append = F,row.names = F,col.names = F)
  write.table(paste("RandomUMI P-value:",cor.test(df_rand$original,df_rand$target,method = "pearson")$p.value),paste0(outprefix,"_stat.test.txt"),quote=F,append = F,row.names = F,col.names = F)
  
  df_rand$method <- "random_UMI"
  if(both_rand_and_keep&&!is.null(df_keep)){
    df_keep$method <- "keep_UMI"
    df_merge <- rbind(df_rand,df_keep)
    write.table(paste("keepUMI:",cor(df_keep$original,df_keep$target,method = "pearson")),paste0(outprefix,"_stat.txt"),quote=F,append = T,row.names = F,col.names = F)
    write.table(paste("keepUMI P-value:",cor.test(df_rand$original,df_rand$target,method = "pearson")$p.value),paste0(outprefix,"_stat.test.txt"),quote=F,append = T,row.names = F,col.names = F)
    
  }else{
    df_merge <- df_rand
  }
  
  g <- ggplot(df_merge,aes(x=original,y=target,col=method))+
    theme_classic()+
    geom_point(size=ptsize,shape=16)+
    # geom_point(size=ptsize,shape=".")+
    # coord_fixed(ratio = 1)+
    theme(axis.title = element_blank(),
          axis.text = element_text(size=7),
          plot.margin= unit(mgn, "lines"),
          axis.line = element_line(size=0.3481),
          axis.ticks = element_line(size=0.3481))
          # axis.line = element_line(size=.3),
          # axis.ticks = element_line(size=.3))
  
  if(carib_x!=F){
    g <- g+annotate("text",x=carib_x,y=0,label="")+scale_x_continuous(limits = c(0,carib_x),breaks = seq(0,carib_x,bin))
  }
  
  if(carib_y!=F){
    g <- g+annotate("text",x=0,y=carib_y,label="")+scale_y_continuous(limits = c(0,carib_y),breaks = seq(0,carib_y,bin))
  }

  if(!show_legend){
    g <- g + NoLegend()
  }else{
    g <- g+
      guides(color=guide_legend(override.aes = list(alpha=1,size=4)))+
      theme(legend.position = legpos,
            legend.title = element_blank(),
            legend.text = element_text(size=10))
  }
  
  if(!show_label){
    g <- g+theme(axis.title = element_blank())
  }
  
  if(both_rand_and_keep){
    g <- g+
      scale_color_manual(values = colset,labels=c(random_UMI="Random UMI",keep_UMI="Conserved UMI"))
  }else{
    g <- g+scale_color_manual(values = colset[2])
  }
  
  ggsave(paste0(outprefix,"_CellCellDist_scatter.jpeg"),g,width = w,height = h,units="cm",dpi=600)
  return(g)
}

diff_rank_sina <- function(rand_ccdist_path,
                           keep_ccdist_path="",
                           colset,
                           keep=FALSE,
                           h=1.4,
                           w=1.4,
                           ymax=50000,
                           outdir,
                           outname){
  df_rand <- read.table(rand_ccdist_path,sep="\t",row.names = 1)
  df_rand$orig_rank <- rank(df_rand$original)
  df_rand$target_rank <- rank(df_rand$target)
  df_rand$rank_diff <- abs(df_rand$orig_rank - df_rand$target_rank)
  df_rand$type <- "random_UMI"
  
  df_shuf <- read.table(rand_ccdist_path,sep="\t",row.names = 1)
  set.seed(1)
  df_shuf$orig_rank <- rank(df_shuf$original) %>% sample
  set.seed(0)
  df_shuf$target_rank <- rank(df_shuf$target) %>% sample
  df_shuf$rank_diff <- abs(df_shuf$orig_rank - df_shuf$target_rank)
  df_shuf$type <- "shuffle"
  
  df_rand <- rbind(df_rand,df_shuf)
  
  if(!keep){
    merge.df <- df_rand
    merge.df$type <- factor(merge.df$type,levels = c("random_UMI","shuffle"))
    write.table(median(filter(merge.df,type=="random_UMI")$rank_diff),paste0(outdir,"/MEDIAN_",outname,".txt"),col.names = F,row.names = F,quote = F)
  }else{
    df_keep <- read.table(keep_ccdist_path,sep="\t",row.names = 1)
    df_keep$orig_rank <- rank(df_keep$original)
    df_keep$target_rank <- rank(df_keep$target)
    df_keep$rank_diff <- abs(df_keep$orig_rank - df_keep$target_rank)
    df_keep$type <- "keep_UMI"
    
    merge.df <- rbind(df_rand,df_keep)
    merge.df$type <- factor(merge.df$type,levels = c("random_UMI","keep_UMI","shuffle"))
  }
  
  col_trans <- adjustcolor(colset,alpha.f = 0.1)
  names(col_trans) <- names(colset)
  
  if(keep){
    v_rand <- (merge.df %>% filter(type=="random_UMI"))$rank_diff
    v_keep <- (merge.df %>% filter(type=="keep_UMI"))$rank_diff
    v_shuf <- (merge.df %>% filter(type=="shuffle"))$rank_diff
    
    tes_r <- wilcox.test(v_rand,v_shuf)
    tes_k <- wilcox.test(v_keep,v_shuf)
    write.table(paste0("random UMI:",tes_r$p.value,"\n","keep UMI:",tes_k$p.value),paste0(outdir,"/",outname,".Utest.txt"),col.names = F,row.names = F,quote = F)
    
    tes <- wilcox.test(v_rand,v_keep)
    
    eff <- cliff.delta(v_rand,v_keep)
    # print(eff)
    
    # print(median(v_keep))
    # print(median(v_rand))
    # return(list(rand=v_rand,keep=v_keep))
    
    write.table(tes$p.value,paste0(outdir,"/",outname,".rand_vs_keep.Utest.txt"),col.names = F,row.names = F,quote = F)
    write.table(eff$estimate,paste0(outdir,"/",outname,".rand_vs_keep.EFF.txt"),col.names = F,row.names = F,quote = F)
    write.table(median(v_keep) / median(v_rand),paste0(outdir,"/",outname,".rand_vs_keep.FC.txt"),col.names = F,row.names = F,quote = F)
    df_temp <- data.frame(cat=c("rand","keep"),median=c(median(v_rand), median(v_keep)))
    write.table(df_temp,paste0(outdir,"/MEDIAN_",outname,".txt"),col.names = F,row.names = F,quote = F)
  }else{
    v_rand <- (merge.df %>% filter(type=="random_UMI"))$rank_diff
    v_shuf <- (merge.df %>% filter(type=="shuffle"))$rank_diff
    
    tes_r <- wilcox.test(v_rand,v_shuf)
    write.table(paste0("random UMI:",tes_r$p.value),paste0(outdir,"/",outname,".Utest.txt"),col.names = F,row.names = F,quote = F)
  }
  
  g <- ggplot(merge.df,aes(x=type,y=rank_diff,color=type,fill=type))+
    theme_classic()+
    geom_sina(size=0.001,shape=16,scale="width")+
    # geom_hline(yintercept = 1,linetype="dotted",size=0.2*1.3)+
    stat_summary(fun = "median",geom = "crossbar",size=.1,color="black")+
    theme(axis.title = element_blank(),
          axis.text.x = element_blank(),
          axis.text = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.y = element_line(size=0.2),
          axis.ticks.y = element_line(size=0.2),
          # axis.text = element_text(size=7),
          plot.background = element_rect(fill = "transparent",color = NA),
          plot.margin = unit(c(1,0,0,0),"mm"))+
    NoLegend()+
    scale_color_manual(values = col_trans)+
    scale_y_continuous(breaks = c(0,ymax))
  # }
  # scale_color_manual(values = c(random="black",keep="black"))+
  
  # if(ymax==0){
  #   g <- g+scale_y_continuous(breaks = c(0,1))
  # }else{
  #   g <- g+scale_y_continuous(breaks = c(0,1),limits = c(-0.01,ymax))
  # }
  ggsave(paste0(outdir,"/",outname,".jpeg"),g,width = w,height = h,units = "cm",dpi=600)
  return(g)
}


vln_plot <- function(rand_rds_path,
                     rand_rds_scr_path,
                     keep_rds_path="",
                     keep_rds_scr_path="",
                     colset,
                     keep=FALSE,
                     h=1.4,
                     w=1.4,
                     ymax=0,
                     outdir,
                     outname){
  diff.df <- readRDS(rand_rds_path)
  diff.df.scr <- readRDS(rand_rds_scr_path)
  merge.rand.df <- merge(diff.df,diff.df.scr,by.x=1,by.y=1)
  merge.rand.df.norm <- data.frame(norm=merge.rand.df$case_diff/merge.rand.df$scramble_diff,
                                   type="random_UMI",
                                   stringsAsFactors = F)
  if(!keep){
    merge.df <- merge.rand.df.norm
    write.table(median(merge.df$norm),paste0(outdir,"/MEDIAN_",outname,".txt"),col.names = F,row.names = F,quote = F)
  }else{
    diff.df.keep <- readRDS(keep_rds_path)
    diff.df.scr.keep <- readRDS(keep_rds_scr_path)
    merge.keep.df <- merge(diff.df.keep,diff.df.scr.keep,by.x=1,by.y=1)
    merge.keep.df.norm <- data.frame(norm=merge.keep.df$case_diff/merge.keep.df$scramble_diff,
                                     type="keep_UMI",
                                     stringsAsFactors = F)
    
    merge.df <- rbind(merge.rand.df.norm,merge.keep.df.norm)
    merge.df$type <- factor(merge.df$type,levels = c("random_UMI","keep_UMI"))
  }
  
  col_trans <- adjustcolor(colset,alpha.f = 0.1)
  names(col_trans) <- names(colset)
  
  if(keep){
    v_rand <- (merge.df %>% filter(type=="random_UMI"))$norm
    v_keep <- (merge.df %>% filter(type=="keep_UMI"))$norm
    
    tes <- wilcox.test(v_rand,v_keep,alternative = "greater")
    print(tes)
    
    eff <- cliff.delta(v_rand,v_keep)
    print(eff)
    
    print(median(v_keep))
    print(median(v_rand))
    
    write.table(tes$p.value,paste0(outdir,"/",outname,".Utest.txt"),col.names = F,row.names = F,quote = F)
    write.table(eff$estimate,paste0(outdir,"/",outname,".EFF.txt"),col.names = F,row.names = F,quote = F)
    write.table(median(v_keep) / median(v_rand),paste0(outdir,"/",outname,".FC.txt"),col.names = F,row.names = F,quote = F)
    df_temp <- data.frame(cat=c("rand","keep"),median=c(median(v_rand), median(v_keep)))
    write.table(df_temp,paste0(outdir,"/MEDIAN_",outname,".txt"),col.names = F,row.names = F,quote = F)
  }
  
  g <- ggplot(merge.df,aes(x=type,y=norm,color=type))+
    theme_classic()+
    geom_jitter(size=0.001,shape=16)+
    geom_hline(yintercept = 1,linetype="dotted",size=0.2*1.3)+
    stat_summary(fun = "median",geom = "crossbar",size=.1,color="black")+
    theme(axis.title = element_blank(),
          axis.text.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.y = element_line(size=0.2),
          axis.ticks.y = element_line(size=0.2),
          axis.text = element_text(size=7),
          plot.background = element_rect(fill = "transparent",color = NA),
          plot.margin = unit(c(1,0,0,0),"mm"))+
    NoLegend()+
    scale_color_manual(values = col_trans)
    # scale_color_manual(values = c(random="black",keep="black"))+
  
  if(ymax==0){
    g <- g+scale_y_continuous(breaks = c(0,1))
  }else{
    g <- g+scale_y_continuous(breaks = c(0,1),limits = c(-0.01,ymax))
  }
  ggsave(paste0(outdir,"/",outname,".jpeg"),g,width = w,height = h,units = "cm",dpi=600)
  return(g)
}


fnc_vln <- function(df_diff,
                    outprefix,
                    factor_level=c("to_10x","to_10x_keepUMI","to_drop","to_drop_keepUMI","dwn90","dwn80","dwn70","scramble"),
                    w=8,
                    h=3,
                    colvec,
                    fig_label=c("Barista\nCR","Barista\nCR-UMI","Barista\nDST","Barista\nDST-UMI","dwn90","dwn80","dwn70","scramble"),
                    dpi=600){
  df_diff_long <- pivot_longer(df_diff,
                               colnames(df_diff),
                               names_to = "destination",
                               values_to = "diff")
  df_diff_long$diff <- log(df_diff_long$diff+1)
  df_diff_long$destination <- factor(df_diff_long$destination,levels = factor_level)
  
  g <- ggplot(df_diff_long,aes(x=destination,y=diff,fill=destination))+
    geom_violin(scale = "width",size=.1)+
    theme_classic()+
    geom_boxplot(fill="white",width=0.2,outlier.shape = NA,size=.1)+
    scale_fill_manual(values = colvec)+
    theme(axis.title.x=element_blank(),
          axis.title.y = element_text(size=10),
          axis.ticks.length = unit(.1,"cm"),
          axis.text = element_text(size=5))+
    scale_y_continuous(breaks=seq(0,10,1))+
    ylab("Log(1+Diff)")+
    NoLegend()+
    scale_x_discrete(labels=fig_label)
  
  ggsave(paste0(outprefix,"_CellCellDist_violin.jpeg"),g,width = w,height = h,units="cm",dpi=dpi)
  return(g)
}



fnc_compare_markers <- function(dir,target_prefix,topn=10){
  
  files <- list.files(dir,pattern = paste0("^",target_prefix),full.names = T)
  for(f in files){
    print(f)
    if(str_detect(f,".orig.tsv")){
      orig   <- read.table(f,header = T,sep="\t",stringsAsFactors = F)
    }else if(str_detect(f,".to10x.tsv")){
      to10x  <- read.table(f,header = T,sep="\t",stringsAsFactors = F)
    }else if(str_detect(f,".todrop.tsv")){
      todrop <- read.table(f,header = T,sep="\t",stringsAsFactors = F)
    }
  }
  
  clusters <- unique(orig$cluster)
  totalmarkers.orig <- c()
  totalmarkers.to10x <- c()
  totalmarkers.todrop <- c()
  
  inter_orig_vs_10x <- c()
  inter_orig_vs_drop <- c()
  
  for(cls in clusters){
    topn.genes.orig.now   <- subset(orig,  cluster==cls)$gene[1:topn] %>% na.omit %>% as.vector
    topn.genes.to10x.now  <- subset(to10x, cluster==cls)$gene[1:topn] %>% na.omit %>% as.vector
    topn.genes.todrop.now <- subset(todrop,cluster==cls)$gene[1:topn] %>% na.omit %>% as.vector
    
    inter_orig_vs_10x <- c(inter_orig_vs_10x,
                           length(intersect(topn.genes.to10x.now,topn.genes.orig.now))/length(topn.genes.orig.now))
    inter_orig_vs_drop <- c(inter_orig_vs_drop,
                           length(intersect(topn.genes.todrop.now,topn.genes.orig.now))/length(topn.genes.orig.now))
    
    totalmarkers.orig   <- c(totalmarkers.orig,  topn.genes.orig.now)
    totalmarkers.to10x  <- c(totalmarkers.to10x, topn.genes.to10x.now)
    totalmarkers.todrop <- c(totalmarkers.todrop,topn.genes.todrop.now)
  }  
  
  L <- list(total_orig = totalmarkers.orig,
            total_10x=totalmarkers.to10x,
            total_drop=totalmarkers.todrop,
            inter_orig_vs_10x=inter_orig_vs_10x,
            inter_orig_vs_drop=inter_orig_vs_drop)
  
  return(L)
}

fnc_draw_venn <- function(L,
                          outdir,
                          outprefix,
                          h=10,
                          w=10,
                          cex=c(1.5,1.5,1.5),
                          cat.pos=c(0,0),
                          cat.dist=c(0.03,0.03),
                          cat.cex=c(1.5,1.5),
                          res=400){
  
  outname_10x  <- paste0(outdir,"/",outprefix,"_to10x_venn.png")
  outname_drop <- paste0(outdir,"/",outprefix,"_todrop_venn.png")
  
  L_to10x_venn=list(orig=unique(L[["total_orig"]]),to10x=unique(L[["total_10x"]]))
  rot=0
  if(length(L_to10x_venn$to10x)<=length(L_to10x_venn$orig)){rot=180}
  venn.diagram(L_to10x_venn, 
               filename=outname_10x,
               category.names = c("Original","Barista-CR"),
               imagetype="png", 
               fontfamily = "sans",
               cat.fontfamily = "sans",
               rotation.degree=rot,
               height=h, 
               width=w, 
               units = "cm",
               fill=c("grey","lightgreen"), 
               lty=0,
               scaled=F, 
               cex=cex, 
               cat.pos=cat.pos, 
               cat.dist=cat.dist, 
               cat.cex=cat.cex,
               resolution = res)
  
  L_todrop_venn=list(orig=unique(L[["total_orig"]]),to10x=unique(L[["total_drop"]]))
  rot=0
  if(length(L_todrop_venn$to10x)<=length(L_todrop_venn$orig)){rot=180}
  venn.diagram(L_todrop_venn, 
               filename=outname_drop,
               category.names = c("Original","Barista-DST"),
               imagetype="png", 
               fontfamily = "sans",
               cat.fontfamily = "sans",
               rotation.degree=rot,
               height=h, 
               width=w, 
               units = "cm",
               fill=c("grey","gold"), 
               lty=0,
               scaled=F, 
               cex=cex, 
               cat.pos=cat.pos, 
               cat.dist=cat.dist, 
               cat.cex=cat.cex,
               resolution = res)
  
  logfiles <- list.files(outdir,pattern = paste0(".log$"),full.names = T)
  file.remove(logfiles)
  return()
}


GeomSplitViolin <- ggproto(
  "GeomSplitViolin",
  GeomViolin,
  draw_group = function(self, data, ..., draw_quantiles = NULL) {
    data <- transform(data,
                      xminv = x - violinwidth * (x - xmin),
                      xmaxv = x + violinwidth * (xmax - x)
    )
    grp <- data[1, "group"]
    newdata <- plyr::arrange(
      transform(data, x = if (grp %% 2 == 1) xminv else xmaxv),
      if (grp %% 2 == 1) y else -y
    )
    newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
    newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
    if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
      stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
      quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
      aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
      aesthetics$alpha <- rep(1, nrow(quantiles))
      both <- cbind(quantiles, aesthetics)
      quantile_grob <- GeomPath$draw_panel(both, ...)
      ggplot2:::ggname(
        "geom_split_violin",
        grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob)
      )
    } else {
      ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
    }
  }
)

geom_split_violin <- function(mapping = NULL,
                              data = NULL,
                              stat = "ydensity",
                              position = "identity", ...,
                              draw_quantiles = NULL,
                              trim = TRUE,
                              scale = "area",
                              na.rm = FALSE,
                              show.legend = NA,
                              inherit.aes = TRUE) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomSplitViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      trim = trim,
      scale = scale,
      draw_quantiles = draw_quantiles,
      na.rm = na.rm, ...
    )
  )
}



