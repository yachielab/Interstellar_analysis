library(tidyverse)
library(RColorBrewer)
library(ggforce)

setwd("pacrat/MSH2/")


# Version 221022
onlyABP.cntdist <- read.table("BCcount_distribution.onlyABP.txt",header=F) %>% mutate(norm_cnt = V2/sum(V2),name="ABP",src="ABP")
onlyIS.cntdist <- read.table("BCcount_distribution.onlyIS.txt",header=F) %>% mutate(norm_cnt = V2/sum(V2),name="IS",src="IS")
commonIS.cntdist <- read.table("BCcount_distribution.common.txt",header=F) %>% mutate(norm_cnt = V2/sum(V2),name="common",src="IS")
commonABP.cntdist <- read.table("BCcount_distribution.common.fromPacRat.txt",header=F) %>% mutate(norm_cnt = V2/sum(V2),name="common",src="ABP")

all.dist <- rbind(onlyABP.cntdist,onlyIS.cntdist,commonIS.cntdist,commonABP.cntdist)
all.dist <- all.dist %>% mutate(src=factor(src,levels = c("IS","ABP")),
                                name=factor(name,levels = c("common","IS","ABP")))
col_vec <- c("common"=brewer.pal(3,"Set2")[2],
             "IS"=brewer.pal(3,"Set2")[3],
             "ABP"=brewer.pal(3,"Set2")[1])

g1 <- ggplot(all.dist %>% filter(src=="IS"),aes(V1,norm_cnt,fill=name))+
  theme_bw()+
  geom_bar(stat = "identity", position = "stack")+
  scale_fill_manual(values = col_vec)+
  theme(axis.title = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border=element_rect(colour = "black", fill=NA, size=.68),
        axis.line = element_blank(),
        axis.ticks = element_line(size=.34),
        axis.text = element_text(size=7),
        strip.text = element_blank(),
        strip.background = element_rect(fill = "white",size = .68,color="black"),
        legend.position = "none")
ggsave("bc_cnt_distribution.IS.pdf",g1,width = 6,height = 4.25,units = "cm")

g1 <- ggplot(all.dist %>% filter(src=="ABP"),aes(V1,norm_cnt,fill=name))+
  theme_bw()+
  geom_bar(stat = "identity", position = "stack")+
  scale_fill_manual(values = col_vec)+
  theme(axis.title = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border=element_rect(colour = "black", fill=NA, size=.68),
        axis.line = element_blank(),
        axis.ticks = element_line(size=.34),
        axis.text = element_text(size=7),
        strip.text = element_blank(),
        strip.background = element_rect(fill = "white",size = .68,color="black"),
        legend.position = "none")
ggsave("bc_cnt_distribution.ABP.pdf",g1,width = 6,height = 4.25,units = "cm")

g2 <- ggplot(all.dist %>% filter(src=="IS"),aes(V1,norm_cnt,fill=name))+
  theme_bw()+
  geom_bar(stat = "identity", position = "stack")+
  scale_fill_manual(values = col_vec)+
  scale_x_continuous(limits = c(0,100))+
  scale_y_continuous(limits = c(0,0.025),oob = squish)+
  theme(axis.title = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border=element_rect(colour = "black", fill=NA, size=.68),
        axis.line = element_blank(),
        axis.ticks = element_line(size=.34),
        axis.text = element_text(size=7),
        strip.background = element_rect(fill = "white",size = .68,color="black"),
        legend.position = "none")
ggsave("bc_cnt_distribution.crop.IS.pdf",g2,width = 4.25,height = 3,units = "cm")

g2 <- ggplot(all.dist %>% filter(src=="ABP"),aes(V1,norm_cnt,fill=name))+
  theme_bw()+
  geom_bar(stat = "identity", position = "stack")+
  scale_fill_manual(values = col_vec)+
  scale_x_continuous(limits = c(0,100))+
  scale_y_continuous(limits = c(0,0.025),oob = squish)+
  theme(axis.title = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border=element_rect(colour = "black", fill=NA, size=.68),
        axis.line = element_blank(),
        axis.ticks = element_line(size=.34),
        axis.text = element_text(size=7),
        strip.background = element_rect(fill = "white",size = .68,color="black"),
        legend.position = "none")
ggsave("bc_cnt_distribution.crop.ABP.pdf",g2,width = 4.25,height = 3,units = "cm")




# length distribution 221022
len_dist <- read.table("analysis_221013/length_distribution.withFromPacRat.tsv",sep="\t",header=F)
len_dist <- len_dist %>% group_by(V3) %>% 
  mutate(norm_cnt = V2/sum(V2),V3=factor(V3,levels = c("thresh0_onlyIS",
                                                       "thresh0_onlyABP",
                                                       "thresh0_common",
                                                       "thresh0_common_fromPacRat",
                                                       "thresh1_common_fromPacRat",
                                                       "thresh1_common",
                                                       "thresh1_onlyABP",
                                                       "thresh1_onlyIS")))

g <- ggplot(len_dist %>% filter(V1>2910 & V1<2976),aes(V1,norm_cnt,fill=V3))+
  theme_bw()+
  geom_bar(stat = "identity", position = "identity",alpha=1)+
  facet_grid(rows = vars(V3))+
  theme(axis.title = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border=element_rect(colour = "black", fill=NA, size=.68),
        axis.line = element_blank(),
        axis.ticks = element_line(size=.34),
        axis.text = element_text(size=7),
        strip.text = element_blank(),
        strip.background = element_rect(fill = "white",size = .34,color="white"),
        legend.position = "none")+
  scale_fill_manual(values = brewer.pal(6,"Set2")[c(3,1,2,2,6,6,4,5)])+
  scale_x_continuous(breaks = c(2928,2950,2975))
ggsave("analysis_221013/len_dist.withFromPacRat.pdf",g,width = 4.6,height = 10,units = "cm")





# CDS hits 221022
t0_ABP <- read.table("analysis_221013/summary/thresh0_onlyABP_summary.tsv",header=T)
t0_IS <- read.table("analysis_221013/summary/thresh0_onlyIS_summary.tsv",header=T)
t0_common <- read.table("analysis_221013/summary/thresh0_common_summary.tsv",header=T)
t0_common_ABP <- read.table("analysis_221013/summary/thresh0_common_fromPacRat_summary.tsv",header=T)
t1_ABP <- read.table("analysis_221013/summary/thresh1_onlyABP_summary.tsv",header=T)
t1_IS <- read.table("analysis_221013/summary/thresh1_onlyIS_summary.tsv",header=T)
t1_common <- read.table("analysis_221013/summary/thresh1_common_summary.tsv",header=T)
t1_common_ABP <- read.table("analysis_221013/summary/thresh1_common_fromPacRat_summary.tsv",header=T)

t0_ABP <- mutate(t0_ABP, name="t0_onlyABP",src="t0_onlyABP")
t0_IS <- mutate(t0_IS, name="t0_onlyIS",src="t0_onlyIS")
t0_common <- mutate(t0_common, name="t0_common",src="t0_common")
t0_common_ABP <- mutate(t0_common, name="t0_common",src="t0_common_ABP")
t1_ABP <- mutate(t1_ABP, name="t1_onlyABP",src="t1_onlyABP")
t1_IS <- mutate(t1_IS, name="t1_onlyIS",src="t1_onlyIS")
t1_common <- mutate(t1_common, name="t1_common",src="t1_common")
t1_common_ABP <- mutate(t1_common, name="t1_common",src="t1_common_ABP")

summary_table <- rbind(t0_ABP,t0_IS,t0_common,t0_common_ABP,t1_ABP,t1_IS,t1_common,t1_common_ABP)

cds_table <- read.csv("MSH2_Expected_alleles.csv",header=T)
aa_conv <- read.csv("aa_conv.csv",header=F)
aa_conv_vec <- aa_conv$V2; names(aa_conv_vec) <- aa_conv$V1
cds_table <- mutate(cds_table, Alt1char=aa_conv_vec[cds_table$Mutant_AA])
cds_table <- mutate(cds_table, pat = paste0(Ref_AA,AA_number,Alt1char))

summary_table <- mutate(summary_table, 
                        pat_hit = ifelse(mut_pat %in% cds_table$pat, "yes", "no"))
cds_hit_aggr <- aggregate(summary_table$pat_hit,as.list(summary_table[,c("name","src")]),table) %>% as.matrix() %>% as.data.frame
colnames(cds_hit_aggr) <- c("Group","src","No","Yes")
cds_hit_aggr.long <- pivot_longer(cds_hit_aggr,cols = colnames(cds_hit_aggr)[c(-1,-2)])
cds_hit_aggr.long$value <- as.integer(cds_hit_aggr.long$value)

cds_hit_aggr.long$Group <- factor(cds_hit_aggr.long$Group, levels = c("t0_onlyIS",
                                                                      "t0_onlyABP",
                                                                      "t0_common",
                                                                      "t1_common",
                                                                      "t1_onlyABP",
                                                                      "t1_onlyIS"))
cds_hit_aggr.long$src <- factor(cds_hit_aggr.long$src, levels = c("t0_onlyIS",
                                                                  "t0_onlyABP",
                                                                  "t0_common",
                                                                  "t0_common_ABP",
                                                                  "t1_common_ABP",
                                                                  "t1_common",
                                                                  "t1_onlyABP",
                                                                  "t1_onlyIS"))

g <- ggplot(cds_hit_aggr.long, aes(Group,value,fill=name))+
  theme_classic()+
  facet_grid(rows = vars(src),space = "free",scales = "free")+
  geom_bar(stat = "identity",position = "fill",alpha=.5)+
  scale_fill_manual(values = c(No="grey",Yes="red"))+
  theme(axis.title = element_blank(),
        axis.line = element_line(size=.34),
        axis.ticks = element_line(size=.34),
        axis.text = element_text(size=7),
        strip.text = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "none")+
  coord_flip()

ggsave("analysis_221013/CDS_hit_stat.withFromABP.pdf",g,width = 4.2,height = 10,units = "cm")


