library(tidyverse)
library(RColorBrewer)
library(Seurat)

L <- list()
for(i in c("10x","split")){
  in_file_todrop <- paste0("~/work/yachie/droParser/tool_test/paper_data/scRNA/efficiency/",i,"_todrop_UMIcount.tsv")
  in_file_todrop_noalloc <- paste0("~/work/yachie/droParser/tool_test/paper_data/scRNA/efficiency_220307//",i,"_todrop_UMIcount.tsv")
  
  df_to10x <- read.table(in_file_todrop,sep="\t",header=F)
  df_to10x <- filter(df_to10x,V3>200)
  
  conv_eff_to10x <- (df_to10x[,3]-df_to10x[,4])/df_to10x[,3]
  df_eff_to10x <- data.frame(dest=rep("With value space opt.",length(conv_eff_to10x)),loss=conv_eff_to10x)

  df_todrop <- read.table(in_file_todrop_noalloc,sep="\t",header=F)
  df_todrop <- filter(df_todrop,V3>200)
  
  conv_eff_todrop <- (df_todrop[,3]-df_todrop[,4])/df_todrop[,3]
  df_eff_todrop <- data.frame(dest=rep("Without value space opt.",length(conv_eff_todrop)),loss=conv_eff_todrop)
  
  tes <- wilcox.test(conv_eff_to10x,conv_eff_todrop)
  
  df_eff <- rbind(df_eff_to10x,df_eff_todrop)
  df_eff$dest <- factor(df_eff$dest,levels = c("With value space opt.","Without value space opt."))
  L[[i]] <- df_eff
  
  g <- ggplot(df_eff,aes(x=dest,y=loss,fill=dest))+
    theme_classic()+
    geom_violin(size=0.15,scale = "width",fill="forestgreen",col="forestgreen")+
    theme(axis.title = element_blank(),
          axis.text.y = element_text(size=7),
          plot.margin = margin(unit(c(0,0,0,1),unit="mm")),
          axis.line.y = element_line(size = 0.3481),
          axis.line.x=element_blank(),
          axis.ticks.x = element_blank(),
          plot.background = element_rect(fill = "transparent",color = NA))+
    NoLegend()
  ggsave(paste0(outdir,"/",i,"drop_loss_vln.pdf"),g,width =4,height=4,units = "cm")
  
}

