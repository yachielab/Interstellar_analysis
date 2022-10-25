library(ggplot2)
library(RColorBrewer)
library(tidyverse)
my_col <- brewer.pal(4,"Paired")[1:2]; names(my_col) <- c("keep_UMI","random_UMI")
my_col_drop <- brewer.pal(4,"Paired")[3:4]; names(my_col_drop) <- c("keep_UMI","random_UMI")
col_vln_base <- c("paleturquoise3","blue3","gold2","thistle2","grey30","grey50","grey70","black")
my_col["shuffle"] <- "grey"
my_col_drop["shuffle"] <- "grey"

source("~/work/yachie/droParser/not_use/barista_fig2_source.R")
setwd("/Users/yusukekijima/work/yachie/droParser/tool_test/paper_data/scRNA")

#10x_to_all
rand_10xto10x <- read.table("CCdist_220306/Mar2022_10x_to_10x_randReassign_CCDist.tsv",header = T)
rand_10xtodrop <- read.table("CCdist_220306/Mar2022_10x_to_drop_randReassign_CCDist.tsv",header = T)
keep_10xto10x <- read.table("CCdist_220306/Mar2022_10x_to_10x_keep_CCDist.tsv",header = T)
noAlloc_10xtodrop <- read.table("CCdist_220306/Mar2022_10x_to_drop_NoAlloc_CCDist.tsv",header = T)

g <- fnc_draw_scatter(df_rand = rand_10xto10x,
                      colset = my_col, 
                      outprefix = "fig_CCdist_220306/10x_to_10x",
                      both_rand_and_keep = T,
                      df_keep = keep_10xto10x,
                      show_legend = F,
                      show_label=F,
                      carib_x=210,
                      carib_y=210)

g <- fnc_draw_scatter(df_rand = rand_10xtodrop,
                      colset=my_col_drop,
                      outprefix = "fig_CCdist_220306/10x_to_drop",
                      both_rand_and_keep = F,
                      show_legend = F,
                      show_label=F,
                      carib_x=210,
                      carib_y=210)

g <- fnc_draw_scatter(df_rand = noAlloc_10xtodrop,
                      colset=my_col_drop,
                      outprefix = "fig_CCdist_220306/10x_to_drop_NoAlloc",
                      both_rand_and_keep = F,
                      show_legend = F,
                      show_label=F,
                      carib_x=210,
                      carib_y=210)
# 
# g <- vln_plot(rand_rds_path = "CCdist_0429/10x_to10x.diff.rds",
#               rand_rds_scr_path = "CCdist_0429/10x_to10x.scr.diff.rds",
#               keep = T,
#               keep_rds_path = "CCdist_0429/10x_to10x_keep.diff.rds",
#               keep_rds_scr_path = "CCdist_0429/10x_to10x_keep.scr.diff.rds",
#               colset = my_col,
#               outdir = "Utest_jitter//",
#               outname="10x_to_10x",
#               w = 1.06,
#               h = 1.06)
# 
# g <- vln_plot(rand_rds_path = "CCdist_0429/10x_todrop.diff.rds",
#               rand_rds_scr_path = "CCdist_0429/10x_todrop.scr.diff.rds",
#               keep = F,
#               colset = my_col_drop[2],
#               outdir = "Utest_jitter//",
#               outname="10x_to_drop",
#               w=0.76,
#               h = 1.06)

g <- diff_rank_sina(rand_ccdist_path = "CCdist_220306/Mar2022_10x_to_10x_randReassign_CCDist.tsv",
              keep_ccdist_path = "CCdist_220306//Mar2022_10x_to_10x_keep_CCDist.tsv",
              keep = T,
              colset = my_col,
              outdir = "Utest_jitter_220307/",
              outname="10x_to_10x",
              w=0.76,
              h = 1.06)

g <- diff_rank_sina(rand_ccdist_path = "CCdist_220306/Mar2022_10x_to_drop_randReassign_CCDist.tsv",
                    keep = F,
                    colset = my_col_drop,
                    outdir = "Utest_jitter_220307/",
                    outname="10x_to_drop",
                    w=0.76,
                    h = 1.06)

#drop_to_all
rand_dropto10x <- read.table("CCdist_0429//drop_to10x_CCDist.tsv",header = T)
rand_droptodrop<- read.table("CCdist_0429/drop_todrop_CCDist.tsv",header = T)
keep_dropto10x <- read.table("CCdist_0429/drop_to10x_keep_CCDist.tsv",header = T)
keep_droptodrop<- read.table("CCdist_0429/drop_todrop_keep_CCDist.tsv",header = T)
seq2seq_droptodrop<- read.table("CCdist_220306/Mar2022_processed_drop_to_drop_SEQ2SEQ_CCDist.tsv",header = T)
source("~/work/yachie/droParser/not_use/barista_fig2_source.R")
g <- fnc_draw_scatter(df_rand = rand_dropto10x,
                      colset = my_col, 
                      outprefix = "fig_CCdist_220306//drop_to_10x",
                      both_rand_and_keep = T,
                      df_keep = keep_dropto10x,
                      show_legend = F,
                      show_label=F,
                      carib_x=210,
                      carib_y=210)

g <- fnc_draw_scatter(df_rand = rand_droptodrop,
                      colset=my_col_drop,
                      outprefix = "fig_CCdist_220306///drop_to_drop",
                      both_rand_and_keep = T,
                      df_keep = keep_droptodrop,
                      show_legend = F,
                      show_label=F,
                      # legpos = "bottom",
                      carib_x=210,
                      carib_y=210)

g <- fnc_draw_scatter(df_rand = seq2seq_droptodrop,
                      colset="gold",
                      outprefix = "fig_CCdist_220306/drop_to_drop_SEQ2SEQ",
                      show_legend = F,
                      show_label=F,
                      # legpos = "bottom",
                      carib_x=210,
                      carib_y=210)

# g <- vln_plot(rand_rds_path = "CCdist_0429/drop_to10x.diff.rds",
#               rand_rds_scr_path = "CCdist_0429/drop_to10x.scr.diff.rds",
#               keep = T,
#               keep_rds_path = "CCdist_0429/drop_to10x_keep.diff.rds",
#               keep_rds_scr_path = "CCdist_0429/drop_to10x_keep.scr.diff.rds",
#               colset = my_col,
#               outdir = "Utest_jitter///",
#               outname="drop_to_10x",
#               w = 1.06,
#               h = 1.06)
# 
# g <- vln_plot(rand_rds_path = "CCdist_0429/drop_todrop.diff.rds",
#               rand_rds_scr_path = "CCdist_0429/drop_todrop.scr.diff.rds",
#               keep = T,
#               keep_rds_path = "CCdist_0429/drop_todrop_keep.diff.rds",
#               keep_rds_scr_path = "CCdist_0429/drop_todrop_keep.scr.diff.rds",
#               colset = my_col_drop,
#               outdir = "Utest_jitter//",
#               outname="drop_to_drop",
#               w = 1.06,
#               h = 1.06)

g <- diff_rank_sina(rand_ccdist_path = "CCdist_0429/drop_to10x_CCDist.tsv",
                    keep_ccdist_path = "CCdist_0429/drop_to10x_keep_CCDist.tsv",
                    keep = T,
                    colset = my_col,
                    outdir = "Utest_jitter_220307/",
                    outname="drop_to_10x",
                    w=0.76,
                    h = 1.06)

g <- diff_rank_sina(rand_ccdist_path = "CCdist_0429/drop_todrop_CCDist.tsv",
                    keep_ccdist_path = "CCdist_0429/drop_todrop_keep_CCDist.tsv",
                    keep = T,
                    colset = my_col_drop,
                    outdir = "Utest_jitter_220307/",
                    outname="drop_to_drop",
                    w=0.76,
                    h = 1.06)




#Quartz_to_all
rand_qto10x <- read.table("CCdist_0429/qtz_to10x_CCDist.tsv",header = T)
rand_qtodrop<- read.table("CCdist_0429/qtz_todrop_CCDist.tsv",header = T)
keep_qto10x <- read.table("CCdist_0429/qtz_to10x_keep_CCDist.tsv",header = T)
keep_qtodrop<- read.table("CCdist_0429/qtz_todrop_keep_CCDist.tsv",header = T)

g <- fnc_draw_scatter(df_rand = rand_qto10x,
                      colset = my_col, 
                      outprefix = "fig_CCdist_220306//quartz_to_10x",
                      both_rand_and_keep = T,
                      df_keep = keep_qto10x,
                      show_legend = F,
                      show_label=F,
                      carib_x=210,
                      carib_y=210)

g <- fnc_draw_scatter(df_rand = rand_qtodrop,
                      colset=my_col_drop,
                      outprefix = "fig_CCdist_220306///quartz_to_drop",
                      both_rand_and_keep = T,
                      df_keep = keep_qtodrop,
                      show_legend = F,
                      show_label=F,
                      carib_x=210,
                      carib_y=210)

# g <- vln_plot(rand_rds_path = "CCdist_0429/qtz_to10x.diff.rds",
#               rand_rds_scr_path = "CCdist_0429/qtz_to10x.scr.diff.rds",
#               keep = T,
#               keep_rds_path = "CCdist_0429/qtz_to10x_keep.diff.rds",
#               keep_rds_scr_path = "CCdist_0429/qtz_to10x_keep.scr.diff.rds",
#               colset = my_col,
#               outdir = "Utest_jitter/",
#               outname="quartz_to_10x",
#               w = 1.06,
#               h = 1.06)
# 
# g <- vln_plot(rand_rds_path = "CCdist_0429/qtz_todrop.diff.rds",
#               rand_rds_scr_path = "CCdist_0429/qtz_todrop.scr.diff.rds",
#               keep = T,
#               keep_rds_path = "CCdist_0429/qtz_todrop_keep.diff.rds",
#               keep_rds_scr_path = "CCdist_0429/qtz_todrop_keep.scr.diff.rds",
#               colset = my_col_drop,
#               outdir = "Utest_jitter/",
#               outname="quartz_to_drop",
#               w = 1.06,
#               h = 1.06)

g <- diff_rank_sina(rand_ccdist_path = "CCdist_0429/qtz_to10x_CCDist.tsv",
                    keep_ccdist_path = "CCdist_0429/qtz_to10x_keep_CCDist.tsv",
                    keep = T,
                    colset = my_col,
                    outdir = "Utest_jitter_220307/",
                    outname="qtz_to_10x",
                    w=0.76,
                    h = 1.06)

g <- diff_rank_sina(rand_ccdist_path = "CCdist_0429/qtz_todrop_CCDist.tsv",
                    keep_ccdist_path = "CCdist_0429/qtz_todrop_keep_CCDist.tsv",
                    keep = T,
                    colset = my_col_drop,
                    outdir = "Utest_jitter_220307/",
                    outname="qtz_to_drop",
                    w=0.76,
                    h = 1.06)



#split_to_all
rand_splitto10x <- read.table("CCdist_0429/split_to10x_CCDist.tsv",header = T)
rand_splittodrop <- read.table("CCdist_0429/split_todrop_CCDist.tsv",header = T)
keep_splitto10x <- read.table("CCdist_0429/split_to10x_keep_CCDist.tsv",header = T)
noaaloc_splittodrop <- read.table("CCdist_220306/Mar2022_SPLiT_to_drop_randUMI_NoAlloc_CCDist.tsv",header = T)


g <- fnc_draw_scatter(df_rand = rand_splitto10x,
                      colset = my_col, 
                      outprefix = "fig_CCdist_220306///split_to_10x",
                      both_rand_and_keep = T,
                      df_keep = keep_splitto10x,
                      show_legend = F,
                      show_label=F,
                      carib_x=210,
                      carib_y=210)

g <- fnc_draw_scatter(df_rand = rand_splittodrop,
                      colset=my_col_drop,
                      outprefix = "fig_CCdist_220306//split_to_drop",
                      both_rand_and_keep = F,
                      show_legend = F,
                      show_label=F,
                      carib_x=210,
                      carib_y=210)

g <- fnc_draw_scatter(df_rand = noaaloc_splittodrop,
                      colset=my_col_drop,
                      outprefix = "fig_CCdist_220306///split_to_drop_noAlloc",
                      both_rand_and_keep = F,
                      show_legend = F,
                      show_label=F,
                      carib_x=210,
                      carib_y=210)

# g <- vln_plot(rand_rds_path = "CCdist_0429/split_to10x.diff.rds",
#               rand_rds_scr_path = "CCdist_0429/split_to10x.scr.diff.rds",
#               keep = T,
#               keep_rds_path = "CCdist_0429/split_to10x_keep.diff.rds",
#               keep_rds_scr_path = "CCdist_0429/split_to10x_keep.scr.diff.rds",
#               colset = my_col,
#               outdir = "Utest_jitter/",
#               outname="split_to_10x",
#               w=1.06,
#               h = 1.06)
# 
# g <- vln_plot(rand_rds_path = "CCdist_0429/split_todrop.diff.rds",
#               rand_rds_scr_path = "CCdist_0429/split_todrop.scr.diff.rds",
#               keep = F,
#               colset = my_col_drop[2],
#               outdir = "Utest_jitter/",
#               outname="split_to_drop",
#               w=0.76,
#               h = 1.06)

g <- diff_rank_sina(rand_ccdist_path = "CCdist_0429/split_to10x_CCDist.tsv",
                    keep_ccdist_path = "CCdist_0429/split_to10x_keep_CCDist.tsv",
                    keep = T,
                    colset = my_col,
                    outdir = "Utest_jitter_220307/",
                    outname="split_to_10x",
                    w=0.76,
                    h = 1.06)

g <- diff_rank_sina(rand_ccdist_path = "CCdist_0429//split_todrop_CCDist.tsv",
                    keep = F,
                    colset = my_col_drop,
                    outdir = "Utest_jitter_220307/",
                    outname="split_to_drop",
                    w=0.76,
                    h = 1.06)



#crazy 10x
rand_10xto10x1 <- read.table("CCdist_0429/crazy_10x_heart1_10x_randUMI_samegene_0428_CCDist.tsv",header = T)
rand_10xto10x2 <- read.table("CCdist_0429/crazy_10x_heart2_10x_randUMI_samegene_0428_CCDist.tsv",header = T)
rand_10xto10x3 <- read.table("CCdist_0429/crazy_10x_neuron1_10x_randUMI_samegene_0428_CCDist.tsv",header = T)
rand_10xto10x4 <- read.table("CCdist_0429/crazy_10x_neuron2_10x_randUMI_samegene_0428_CCDist.tsv",header = T)

L <- list(crazy_heart1=rand_10xto10x1,
          crazy_heart2=rand_10xto10x2,
          crazy_neuron1=rand_10xto10x3,
          crazy_neuron2=rand_10xto10x4)

for(i in names(L)){
  g <- fnc_draw_scatter(df_rand = L[[i]],
                        colset = my_col, 
                        outprefix = paste0("fig_CCdist_220306//",i),
                        both_rand_and_keep = F,
                        show_legend = F,
                        show_label=F,
                        carib_x=210,
                        carib_y=210)
}

L2 <- list(crazy_heart1=c("CCdist_0429/crazy_10x_heart1.diff.rds","CCdist_0429/crazy_10x_heart1.scr.diff.rds"),
           crazy_heart2=c("CCdist_0429/crazy_10x_heart2.diff.rds","CCdist_0429/crazy_10x_heart2.scr.diff.rds"),
           crazy_neuron1=c("CCdist_0429/crazy_10x_neuron1.diff.rds","CCdist_0429/crazy_10x_neuron1.scr.diff.rds"),
           crazy_neuron2=c("CCdist_0429/crazy_10x_neuron2.diff.rds","CCdist_0429/crazy_10x_neuron2.scr.diff.rds"))

# for(i in names(L2)){
#   g <- vln_plot(rand_rds_path = L2[[i]][1],
#                 rand_rds_scr_path = L2[[i]][2],
#                 keep = F,
#                 colset = my_col["random_UMI"],
#                 outdir = "Utest_jitter///",
#                 outname=i,
#                 ymax=1.1,
#                 w=0.78,
#                 h = 0.91)
#   
# }

g <- diff_rank_sina(rand_ccdist_path = "CCdist_0429//crazy_10x_heart1_10x_randUMI_samegene_0428_CCDist.tsv",
                    keep = F,
                    colset = my_col,
                    outdir = "Utest_jitter_220307/",
                    outname="crazy_10x_heart1",
                    w=0.76,
                    h = 1.06)

g <- diff_rank_sina(rand_ccdist_path = "CCdist_0429//crazy_10x_heart2_10x_randUMI_samegene_0428_CCDist.tsv",
                    keep = F,
                    colset = my_col,
                    outdir = "Utest_jitter_220307/",
                    outname="crazy_10x_heart2",
                    w=0.76,
                    h = 1.06)

g <- diff_rank_sina(rand_ccdist_path = "CCdist_0429//crazy_10x_neuron1_10x_randUMI_samegene_0428_CCDist.tsv",
                    keep = F,
                    colset = my_col,
                    outdir = "Utest_jitter_220307/",
                    outname="crazy_10x_neuron1",
                    w=0.76,
                    h = 1.06)

g <- diff_rank_sina(rand_ccdist_path = "CCdist_0429//crazy_10x_neuron2_10x_randUMI_samegene_0428_CCDist.tsv",
                    keep = F,
                    colset = my_col,
                    outdir = "Utest_jitter_220307/",
                    outname="crazy_10x_neuron2",
                    w=0.76,
                    h = 1.06)


# drop to drop, SEQ2SEQ
g <- diff_rank_sina(rand_ccdist_path = "CCdist_220306/Mar2022_processed_drop_to_drop_SEQ2SEQ_CCDist.tsv",
                    keep = F,
                    colset = c("random_UMI"="gold","shuffle"="grey"),
                    outdir = "Utest_jitter_220307/",
                    outname="drop_to_drop_SEQ2SEQ",
                    w=0.76,
                    h = 1.06)

# No allocation
g <- diff_rank_sina(rand_ccdist_path = "CCdist_220306/Mar2022_10x_to_drop_NoAlloc_CCDist.tsv",
                    keep = F,
                    colset = my_col_drop,
                    outdir = "Utest_jitter_220307/",
                    outname="10x_to_drop_noalloc",
                    w=0.76,
                    h = 1.06)

g <- diff_rank_sina(rand_ccdist_path = "CCdist_220306/Mar2022_SPLiT_to_drop_randUMI_NoAlloc_CCDist.tsv",
                    keep = F,
                    colset = my_col_drop,
                    outdir = "Utest_jitter_220307/",
                    outname="split_to_drop_noalloc",
                    w=0.76,
                    h = 1.06)


# Round trip conversion, revision round 1
roundtrip_10x <- read.table("~/work/yachie/droParser/tool_test/paper_data/REVISE/roundtrip/CCdist/Aug2022_10x_roundtrip_CCDist.tsv",header = T)
roundtrip_drop <- read.table("~/work/yachie/droParser/tool_test/paper_data/REVISE/roundtrip/CCdist/Aug2022_Drop_roundtrip_CCDist.tsv",header = T)
roundtrip_quartz <- read.table("~/work/yachie/droParser/tool_test/paper_data/REVISE/roundtrip/CCdist/Aug2022_Quartz_roundtrip_CCDist.tsv",header = T)
roundtrip_split <- read.table("~/work/yachie/droParser/tool_test/paper_data/REVISE/roundtrip/CCdist/Aug2022_SPLiT_to_roundtrip_CCDist.tsv",header = T)

L <- list(roundtrip_10x=roundtrip_10x,
          roundtrip_drop=roundtrip_drop,
          roundtrip_quartz=roundtrip_quartz,
          roundtrip_split=roundtrip_split)

for(i in names(L)){
  g <- fnc_draw_scatter(df_rand = L[[i]],
                        colset = my_col, 
                        outprefix = paste0("~/work/yachie/droParser/tool_test/paper_data/REVISE/roundtrip/CCdist/",i),
                        both_rand_and_keep = F,
                        show_legend = F,
                        show_label=F,
                        carib_x=210,
                        carib_y=210)
}

g <- diff_rank_sina(rand_ccdist_path = "~/work/yachie/droParser/tool_test/paper_data/REVISE/roundtrip/CCdist/Aug2022_10x_roundtrip_CCDist.tsv",
                    keep = F,
                    colset = my_col,
                    outdir = "~/work/yachie/droParser/tool_test/paper_data/REVISE/roundtrip/sina/",
                    outname="10x",
                    w=0.76,
                    h = 1.06)

g <- diff_rank_sina(rand_ccdist_path = "~/work/yachie/droParser/tool_test/paper_data/REVISE/roundtrip/CCdist/Aug2022_Drop_roundtrip_CCDist.tsv",
                    keep = F,
                    colset = my_col,
                    outdir = "~/work/yachie/droParser/tool_test/paper_data/REVISE/roundtrip/sina/",
                    outname="drop",
                    w=0.76,
                    h = 1.06)

g <- diff_rank_sina(rand_ccdist_path = "~/work/yachie/droParser/tool_test/paper_data/REVISE/roundtrip/CCdist/Aug2022_Quartz_roundtrip_CCDist.tsv",
                    keep = F,
                    colset = my_col,
                    outdir = "~/work/yachie/droParser/tool_test/paper_data/REVISE/roundtrip/sina//",
                    outname="quartz",
                    w=0.76,
                    h = 1.06)

g <- diff_rank_sina(rand_ccdist_path = "~/work/yachie/droParser/tool_test/paper_data/REVISE/roundtrip/CCdist/Aug2022_SPLiT_to_roundtrip_CCDist.tsv",
                    keep = F,
                    colset = my_col,
                    outdir = "~/work/yachie/droParser/tool_test/paper_data/REVISE/roundtrip/sina//",
                    outname="split",
                    w=0.76,
                    h = 1.06)














#supple markers
# source("~/work/yachie/droParser/not_use/barista_fig2_source.R")
# 
# to10x_scatter_df  <- data.frame(Platform=c(),Ratio=c())
# todrop_scatter_df <- data.frame(Platform=c(),Ratio=c())
# 
# for(i in c("to10x","Drop","Quartz","SPLiT")){
#   L <- fnc_compare_markers("~/work/yachie/droParser/tool_test/markers/",i,topn=50)
#   fnc_draw_venn(L,"~/work/yachie/droParser/tool_test/markers/venn/",outprefix = i)
#   print(names(L))
#   
#   to10x_scatter_df <- rbind(to10x_scatter_df,data.frame(Platform=i,Ratio=L$inter_orig_vs_10x))
#   todrop_scatter_df <- rbind(todrop_scatter_df,data.frame(Platform=i,Ratio=L$inter_orig_vs_drop))
# }
# 
# g <- ggplot(to10x_scatter_df,aes(x=Platform,y=Ratio,col=Platform))+
#   theme_classic()+
#   geom_jitter(size=2)+
#   scale_color_manual(values = c("lightpink3","turquoise1","goldenrod2","lightgreen"))+
#   scale_y_continuous(breaks = seq(0,1,0.1),limits = c(-0.01,1.01))+
#   scale_x_discrete(labels=c("10X Chromium v3","Drop-seq","Quartz-seq v3.1","SPLiT-seq"))+
#   theme(text = element_text(size=10),
#         axis.text.y = element_text(size=6),
#         axis.text.x = element_text(size=10))+
#   NoLegend()
# ggsave("~/work/yachie/droParser/tool_test/markers/venn/jitter_percluster_to10x.jpeg",
#        g,
#        height = 5,
#        width = 12,
#        units = "cm",
#        dpi=500)
# 
# g <- ggplot(todrop_scatter_df,aes(x=Platform,y=Ratio,col=Platform))+
#   theme_classic()+
#   geom_jitter(size=2)+
#   scale_color_manual(values = c("lightpink3","turquoise1","goldenrod2","lightgreen"))+
#   scale_y_continuous(breaks = seq(0,1,0.1),limits = c(-0.01,1.01))+
#   scale_x_discrete(labels=c("10X Chromium v3","Drop-seq","Quartz-seq v3.1","SPLiT-seq"))+
#   theme(text = element_text(size=10),
#         axis.text.y = element_text(size=6),
#         axis.text.x = element_text(size=10))+
#   NoLegend()
# ggsave("~/work/yachie/droParser/tool_test/markers/venn/jitter_percluster_todrop.jpeg",
#        g,
#        height = 5,
#        width = 12,
#        units = "cm",
#        dpi=500)

# L_to10x_venn=list(orig=L[["total_orig"]],to10x=L[["total_10x"]])
# venn.diagram(L_to10x_venn, 
#              filename="~/work/yachie/droParser/tool_test/markers/venn/Qtz_to_10x_venn.png",
#              category.names = c("Original","Barista-CR"),
#              imagetype="png", 
#              fontfamily = "sans",
#              cat.fontfamily = "sans",
#              height=10, 
#              width=10, s
#              units = "cm",
#              fill=c("grey","lightgreen"), 
#              lty=0,
#              scaled=F, 
#              cex=c(1.5,1.5,1.5), 
#              cat.pos=c(0,0), 
#              cat.dist=c(0.03,0.03), 
#              cat.cex=c(1.5,1.5))
# 
# L_todrop_venn=list(orig=L[["total_orig"]],to10x=L[["total_drop"]])
# venn.diagram(L_todrop_venn, 
#              filename="~/work/yachie/droParser/tool_test/markers/venn/Qtz_to_drop_venn.png",
#              category.names = c("Original","Barista-DST"),
#              imagetype="png", 
#              fontfamily = "sans",
#              cat.fontfamily = "sans",
#              height=10, 
#              width=10, 
#              units = "cm",
#              fill=c("grey","gold2"), 
#              lty=0,
#              scaled=F, 
#              cex=c(1.5,1.5,1.5), 
#              cat.pos=c(0,0), 
#              cat.dist=c(0.03,0.03), 
#              cat.cex=c(1.5,1.5))
# 

# seu2 <- seu.4dpf
# v <- rep(100,20)
# v <- c(v,rep(11111,ncol(seu2)-20))
# names(v) <- colnames(seu2)
# 
# Idents(seu2) <- v
# FindAllMarkers(seu2)
