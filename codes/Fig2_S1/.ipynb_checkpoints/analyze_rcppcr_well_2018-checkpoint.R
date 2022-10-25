# Output of analyze_rcppcr.py
rcppcr_AD_result <- read.table("AD_thresh05_summary.tsv",sep="\t",header = T,stringsAsFactors = T)
rcppcr_DB_result <- read.table("DB_thresh05_summary.tsv",sep="\t",header = T,stringsAsFactors = T)

outdir <- "./"
library(Seurat)
library(tidyverse)

# Quality check function
qc <- function(df,dtype){
  if(dtype=="all"){
    for(i in 1:nrow(df)){
      if(df$Quality_filter_bc1[i]=="-" || df$Quality_filter_bc2[i]=="-" || df$Quality_filter_LoxP[i]=="-" || df$Quality_filter_Lox2272[i]=="-"){
        df$Quality_check[i] <- "Not_detected"
      }else if(df$Quality_filter_bc1[i]=="Drop" || df$Quality_filter_bc2[i]=="Drop"){
        df$Quality_check[i] <- "Low_quality_barcode"
      }else if(df$Quality_filter_LoxP[i]=="Drop" || df$Quality_filter_Lox2272[i]=="Drop"){
        df$Quality_check[i] <- "Low_quality_Lox"
      }else if(df$Correspondence_BC_Lox[i]=="Drop"){
        df$Quality_check[i] <- "Conflict_Lox_barcode"
      }
    }
    df$Quality_check <- factor(df$Quality_check,levels = c("Pass","Low_quality_barcode","Low_quality_Lox","Conflict_Lox_barcode","Not_detected"))
    
  }else if(dtype=="bc"){
    for(i in 1:nrow(df)){
      if(df$Quality_filter_bc1[i]=="-" || df$Quality_filter_bc2[i]=="-"){
        df$Quality_check[i] <- "Not_detected"
      }else if(df$Quality_filter_bc1[i]=="Drop" || df$Quality_filter_bc2[i]=="Drop"){
        df$Quality_check[i] <- "Low_quality_barcode"
      }
    }
    df$Quality_check <- factor(df$Quality_check,levels = c("Pass","Low_quality_barcode","Not_detected"))
    
  }else if(dtype=="lox"){
    for(i in 1:nrow(df)){
      if(df$Quality_filter_LoxP[i]=="-" || df$Quality_filter_Lox2272[i]=="-"){
        df$Quality_check[i] <- "Not_detected"
      }else if(df$Quality_filter_LoxP[i]=="Drop" || df$Quality_filter_Lox2272[i]=="Drop"){
        df$Quality_check[i] <- "Low_quality_Lox"
      }
    }
    df$Quality_check <- factor(df$Quality_check,levels = c("Pass","Low_quality_Lox","Not_detected"))
    
  }else{
    stop("dtype should be all, bc or lox.")
  }
  
  return(df)
}

    
# Table generating function
make_df <- function(raw_df,dtype,pname){
  df <- subset(raw_df,Plate_name==pname)
  if(dtype=="all"){
    df <- df[,c("Row",
                "Column",
                "Quality_filter_bc1",
                "Quality_filter_bc2",
                "Correspondence_BC_Lox",
                "Quality_filter_LoxP",
                "Quality_filter_Lox2272")]
  }else if(dtype=="bc"){
    df <- df[,c("Row",
                "Column",
                "Quality_filter_bc1",
                "Quality_filter_bc2")]
  }else if(dtype=="lox"){
    df <- df[,c("Row",
                "Column",
                "Quality_filter_LoxP",
                "Quality_filter_Lox2272")]
  }else{
    stop("dtype should be all, bc or lox.")
  }
  
  df$Row <- factor(df$Row,levels = rev(levels(df$Row)))
  df$Quality_check <- "Pass"
  df <- qc(df,dtype)
  
  return(df)
}

# Sub function
add_non_detect <- function(df){
  for(r in LETTERS[1:16]){
    for(col in 1:24){
      if(nrow(filter(df,Row==r&Column==col))==0){
        df <- rbind(df,c(r,as.integer(col),"Drop","Drop","Not_detected"))
      }
    }
  }
  return(df)
}

# Sub function
add_non_detect_merge <- function(df){
  for(r in LETTERS[1:16]){
    for(col in 1:24){
      if(nrow(filter(df,Row==r&Column==col))==0){
        df <- rbind(df,c(r,as.integer(col),rep("Drop",5),"Not_detected"))
      }
    }
  }
  return(df)
}

group.colors <- c("Pass" = "blue" ,
                  "Low_quality_barcode" ="red",
                  "Low_quality_Lox" =  "red",
                  "Conflict_Lox_barcode"="gold",
                  "Not_detected"="black")
group.names <- c("Pass" = "Pass" ,
                 "Low_quality_barcode" ="Failed",
                 "Low_quality_Lox" =  "Failed",
                 "Conflict_Lox_barcode"="Conflict tag calling",
                 "Not_detected"="n.d.")


library(ggthemes)
    
#DHFR12, bc
for(i in 10:12){
  pname <- paste0("DHFR[12]0",i)
  df <- make_df(rcppcr_AD_result,"bc",pname)
  df <- add_non_detect(df)
  df$Column <- as.integer(df$Column)
  g <- ggplot(df,aes(x=Column,y=Row,col=Quality_check))+
    theme_few()+
    geom_point(size=.175)+
    scale_x_continuous(position = 'top',breaks = seq(1,24))+
    coord_fixed(ratio=0.9)+
    scale_color_manual(values = group.colors)+
    theme(axis.title = element_blank(),
          plot.title = element_text(hjust = 0.5,size=10),
          axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
  g <- g+NoLegend()
  ggsave(paste0(outdir,"/",pname,"_barcode.pdf"),g,width = 2.5,height = 2,units = "cm")
}

#DHFR12, lox
for(i in 10:12){
  pname <- paste0("DHFR[12]0",i)
  df <- make_df(rcppcr_AD_result,"lox",pname)
  df <- add_non_detect(df)
  df$Column <- as.integer(df$Column)
  g <- ggplot(df,aes(x=Column,y=Row,col=Quality_check))+
    theme_few()+
    geom_point(size=.175)+
    scale_x_continuous(position = 'top',breaks = seq(1,24))+
    coord_fixed(ratio=0.9)+
    scale_color_manual(values = group.colors)+
    theme(axis.title = element_blank(),
          plot.title = element_text(hjust = 0.5,size=10),
          axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
  g <- g+NoLegend()
  ggsave(paste0(outdir,"/",pname,"_lox.pdf"),g,width = 2.5,height = 2,units = "cm")
}

#DHFR12, merge
L_AD <- list()
for(i in 10:12){
  pname <- paste0("DHFR[12]0",i)
  df <- make_df(rcppcr_AD_result,"all",pname)
  df <- add_non_detect_merge(df)
  df$Column <- as.integer(df$Column)
  g <- ggplot(df,aes(x=Column,y=Row,col=Quality_check))+
    theme_few()+
    geom_point(size=.175)+
    scale_x_continuous(position = 'top',breaks = seq(1,24))+
    coord_fixed(ratio=0.9)+
    scale_color_manual(values = group.colors)+
    theme(axis.title = element_blank(),
          plot.title = element_text(hjust = 0.5,size=10),
          axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
  g <- g+NoLegend()
  L_AD[[i]] <- df
}

#DHFR3, bc
for(i in 10:13){
  pname <- paste0("DHFR[3]0",i)
  df <- make_df(rcppcr_DB_result,"bc",pname)
  df <- add_non_detect(df)
  df$Column <- as.integer(df$Column)
  g <- ggplot(df,aes(x=Column,y=Row,col=Quality_check))+
    theme_few()+
    geom_point(size=.175)+
    scale_x_continuous(position = 'top',breaks = seq(1,24))+
    coord_fixed(ratio=0.9)+
    scale_color_manual(values = group.colors)+
    theme(axis.title = element_blank(),
          plot.title = element_text(hjust = 0.5,size=10),
          axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
  g <- g+NoLegend()
  ggsave(paste0(outdir,"/",pname,"_barcode.pdf"),g,width = 2.5,height = 2,units = "cm")
}

#DHFR, lox
for(i in 10:13){
  pname <- paste0("DHFR[3]0",i)
  df <- make_df(rcppcr_DB_result,"lox",pname)
  df <- add_non_detect(df)
  df$Column <- as.integer(df$Column)
  g <- ggplot(df,aes(x=Column,y=Row,col=Quality_check))+
    theme_few()+
    geom_point(size=.175)+
    scale_x_continuous(position = 'top',breaks = seq(1,24))+
    coord_fixed(ratio=0.9)+
    scale_color_manual(values = group.colors)+
    theme(axis.title = element_blank(),
          plot.title = element_text(hjust = 0.5,size=10),
          axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
  g <- g+NoLegend()
  ggsave(paste0(outdir,"/",pname,"_lox.pdf"),g,width = 2.5,height = 2,units = "cm",dpi=600)
}

#DHFR, merge
L_DB <- list()
for(i in 10:13){
  pname <- paste0("DHFR[3]0",i)
  df <- make_df(rcppcr_DB_result,"all",pname)
  df <- add_non_detect_merge(df)
  df$Column <- as.integer(df$Column)
  g <- ggplot(df,aes(x=Column,y=Row,col=Quality_check))+
    theme_few()+
    geom_point(size=.175)+
    scale_x_continuous(position = 'top',breaks = seq(1,24))+
    coord_fixed(ratio=0.9)+
    scale_color_manual(values = group.colors)+
    theme(axis.title = element_blank(),
          plot.title = element_text(hjust = 0.5,size=10),
          axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
  g <- g+NoLegend()
  L_DB[[i]] <- df
}
    

#Legend
i=13
pname <- paste0("DHFR[3]0",i)
df <- make_df(rcppcr_DB_result,"all",pname)
g <- ggplot(df,aes(x=Column,y=Row,col=Quality_check))+
  theme_bw(base_rect_size = .5)+
  geom_point(size=.9)+
  scale_x_continuous(position = 'top',breaks = seq(1,24))+
  coord_fixed(ratio=0.9)+
  scale_color_manual(values = group.colors,
                     labels = group.names)+
  # NoLegend()+
  theme(axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5,size=10),
        axis.line = element_blank(),
        axis.text = element_text(size=3.4),
        axis.ticks = element_line(size=.1),
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        legend.position = "bottom")+
  guides(colour=guide_legend(override.aes = list(size=3)))
ggsave(paste0(outdir,"/legend_",pname,"_merge.jpeg"),g,width = 20,height = 8,units = "cm")

    
    
#Well select for sanger sequencing
seed <- 3
rcppcr_AD_result <- read.table("AD_thresh05_summary.tsv",sep="\t",header = T)
rcppcr_DB_result <- read.table("DB_thresh05_summary.tsv",sep="\t",header = T)
outdir <- "./"

library(tidyverse)
set.seed(seed)
AD_selected <- filter(rcppcr_AD_result,Quality_filter_total=="Pass") %>%
  group_by(Plate_name) %>%
  sample_n(4)
set.seed(seed)
DB_selected <- filter(rcppcr_DB_result,Quality_filter_total=="Pass") %>%
  group_by(Plate_name) %>%
  sample_n(3)

write.table(AD_selected,paste0(outdir,"DHFR12.selected.seed1.tsv"),quote=F,row.names = F,col.names = T,sep="\t")
write.table(DB_selected,paste0(outdir,"DHFR3.selected.seed3.tsv"),quote=F,row.names = F,col.names = T,sep="\t")
