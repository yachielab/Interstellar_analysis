library(RColorBrewer)
library(viridis)
library(ggthemes)

###Validation of coordinate expansion with N=5
library(RColorBrewer)
library(viridis)
library(ggthemes)
library(tidyverse)

paths <- list(cereb2=c("data/zoom-in/cereb2_exp1_aligned.tsv","data/zoom-in/cereb2_exp5_aligned.tsv","data/zoom-in/cereb2_exp10_aligned.tsv"),
              cereb1=c("data/zoom-in/cereb1_exp1_aligned.tsv","data/zoom-in/cereb1_exp5_aligned.tsv","data/zoom-in/cereb1_exp10_aligned.tsv"),
              kidney=c("data/zoom-in/kidney_exp1_aligned.tsv","data/zoom-in/kidney_exp5_aligned.tsv","data/zoom-in/kidney_exp10_aligned.tsv"),
              liver=c("data/zoom-in/liver_exp1_aligned.tsv","data/zoom-in/liver_exp5_aligned.tsv","data/zoom-in/liver_exp10_aligned.tsv"),
              hippo=c("data/zoom-in/hippo_exp1_aligned.tsv","data/zoom-in/hippo_exp5_aligned.tsv","data/zoom-in/hippo_exp10_aligned.tsv"))

table_pileup <- function(tbl1,tbl5,tbl10){
  out_df <- data.frame(Expansion_rate=c(),N_plet=c(),Freq=c())
  for(rate in c("1","2","3","4<")){
    if(rate != "4<"){
      out_df <- rbind(out_df,data.frame(Expansion_rate="1x", N_plet=rate, Freq=ifelse(is.na(tbl1[rate]/sum(tbl1)),0,tbl1[rate]/sum(tbl1)) %>% unname))
      out_df <- rbind(out_df,data.frame(Expansion_rate="5x", N_plet=rate, Freq=ifelse(is.na(tbl5[rate]/sum(tbl5)),0,tbl5[rate]/sum(tbl5)) %>% unname))
      out_df <- rbind(out_df,data.frame(Expansion_rate="10x",N_plet=rate, Freq=ifelse(is.na(tbl10[rate]/sum(tbl10)),0,tbl10[rate]/sum(tbl10)) %>% unname))
    }else{
      out_df <- rbind(out_df,data.frame(Expansion_rate="1x", N_plet=rate, Freq=ifelse(is.na((sum(tbl1)-sum(tbl1[1:3]))/sum(tbl1)),0,(sum(tbl1)-sum(tbl1[1:3]))/sum(tbl1)) %>% unname))
      out_df <- rbind(out_df,data.frame(Expansion_rate="5x", N_plet=rate, Freq=ifelse(is.na((sum(tbl5)-sum(tbl5[1:3]))/sum(tbl5)),0,(sum(tbl5)-sum(tbl5[1:3]))/sum(tbl5)) %>% unname))
      out_df <- rbind(out_df,data.frame(Expansion_rate="10x",N_plet=rate, Freq=ifelse(is.na((sum(tbl10)-sum(tbl10[1:3]))/sum(tbl10)),0,(sum(tbl10)-sum(tbl10[1:3]))/sum(tbl10)) %>% unname))
    }
  }
  out_df$Freq[is.na(out_df$Freq)] <- 0
  return(out_df)
}

summary_df <- data.frame()
for(path in paths){
  samplename <- strsplit(basename(path),"_")[[1]][1]
  coord_aligned_1  <- read.table(path[1],sep="\t",row.names = 1,header = T)
  coord_aligned_5  <- read.table(path[2],sep="\t",row.names = 1,header = T)
  coord_aligned_10 <- read.table(path[3],sep="\t",row.names = 1,header = T)
  
  tbl.1  <- table(table(coord_aligned_1$d_ind_bc))
  tbl.5  <- table(table(coord_aligned_5$d_ind_bc))
  tbl.10 <- table(table(coord_aligned_10$d_ind_bc))
  
  tmp_df <- table_pileup(tbl.1,tbl.5,tbl.10)
  tmp_df$sample <- samplename
  summary_df <- rbind(summary_df,tmp_df)
}

summary_df$N_plet <- as.character(summary_df$N_plet) %>% factor(levels = c("1","2","3","4<"), labels = c("1 puck","2 pucks","3 pucks","4< pucks"))
summary_df$Expansion_rate <- factor(summary_df$Expansion_rate,levels = c("1x","5x","10x"),labels = c("1x1","5x5","10x10"))
summary_df$xcoord <- rep(c(c(1,6,11),c(1,6,11)+1,c(1,6,11)+2,c(1,6,11)+3),5)

summary_df_tmp <- summary_df
summary_df_tmp$cat <- unite(summary_df_tmp[,1:2],col=unite)$unite
summary_df_mean<- aggregate(summary_df_tmp$Freq,list(summary_df_tmp$cat),mean)
summary_df_mean <- summary_df_mean %>% separate(col = Group.1,sep="_",into=c("Expansion_rate","N_plet"))
summary_df_mean$N_plet <- as.character(summary_df_mean$N_plet) %>% factor(levels = c("1 puck","2 pucks","3 pucks","4< pucks"), labels = c("1 puck","2 pucks","3 pucks","4< pucks"))
summary_df_mean$Expansion_rate <- factor(summary_df_mean$Expansion_rate,levels = c("1x1","5x5","10x10"),labels = c("1x1","5x5","10x10"))
summary_df_mean$xcoord <- c(11:14,1:4,6:9)

library(plotrix)
summary_df_se <- aggregate(summary_df_tmp$Freq,list(summary_df_tmp$cat),std.error)
summary_df_se <- summary_df_se %>% separate(col = Group.1,sep="_",into=c("Expansion_rate","N_plet"))
summary_df_se$N_plet <- as.character(summary_df_se$N_plet) %>% factor(levels = c("1 puck","2 pucks","3 pucks","4< pucks"), labels = c("1 puck","2 pucks","3 pucks","4< pucks"))
summary_df_se$Expansion_rate <- factor(summary_df_se$Expansion_rate,levels = c("1x1","5x5","10x10"),labels = c("1x1","5x5","10x10"))
summary_df_se$xcoord <- c(11:14,1:4,6:9)
colnames(summary_df_se) <- c("Expansion_rate","N_plet","SE","xcoord")
summary_df_se <- merge(summary_df_se,summary_df_mean,by = c(1,2,4))


library(Seurat)

g <- ggplot(summary_df,aes(x=xcoord,y=Freq))+
  geom_bar(data=summary_df_se,aes(x=xcoord,y=x,fill=Expansion_rate),
           stat="identity",
           position="dodge",
           size=1)+
  geom_jitter(size=.05,width = .2)+
  geom_errorbar(stat = "summary", fun.data = "mean_se", width=.2, size=.25)+
  theme_classic()+
  scale_fill_manual(values = brewer.pal(3,"Paired"))+
  NoLegend()+
  theme(axis.title = element_blank(),
        axis.text.y = element_text(size=7),
        axis.text.x = element_text(size=7,angle=45,hjust =.7,vjust = 1),
        axis.line.x = element_blank(),
        axis.line.y = element_line(size=.34),
        axis.ticks.y = element_line(size=0.34),
        axis.ticks.x = element_blank())+
  scale_x_continuous(breaks = c(1:4,6:9,11:14),labels = rep(c("1 puck","2 pucks","3 pucks","4< pucks"),3))
g
ggsave("frequency_bar.pdf",g,width = 10,height = 4,units = "cm")  

