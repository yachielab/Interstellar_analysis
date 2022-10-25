#0828
source("~/work/yachie/droParser/not_use/barista_fig2_source.R")
my_col <- adjustcolor(c("blue3","blue3"),alpha.f = 0.1)
my_col2 <- adjustcolor(c("blue3","blue3"),alpha.f = 1)

orig <- readRDS("~/work/yachie/droParser/tool_test/sciATAC/data_Aug2020/sciATAC_fly_original_umap.rds")
to10x <- readRDS("~/work/yachie/droParser/tool_test/sciATAC/data_Aug2020/sciATAC_fly_to_10x_umap.rds")

ggsave("~/work/yachie/droParser/tool_test/sciATAC/data_Aug2020/sciATAC_fly_original_umap.local.jpeg",
       orig,
       height = 8,
       width = 8,
       units = "cm",
       dpi=500)
ggsave("~/work/yachie/droParser/tool_test/sciATAC/data_Aug2020/sciATAC_fly_to10x_umap.local.jpeg",
       to10x,
       height = 8,
       width = 8,
       units = "cm",
       dpi=500)
ccdist <- read.table("~/work/yachie/droParser/tool_test/sciATAC/data_Aug2020/sciATAC_fly_CellCellDist.tsv",header = T)

g <- fnc_draw_scatter(df_rand = ccdist,
                      colset = my_col, 
                      outprefix = "~/work/yachie/droParser/tool_test/sciATAC//data_Aug2020/sciATAC_to_10x",
                      both_rand_and_keep = F,
                      show_legend = F,
                      show_label=F,
                      ptsize=.2,
                      w = 5.4,
                      h=5.4)


ccdist2 <- read.table("~/work/yachie/droParser/tool_test/sciATAC/data_Aug2020/sciATAC_fly_rawsignal_scaletest_CellCellDist.tsv",header = T)

g <- fnc_draw_scatter(df_rand = ccdist2,
                      colset = my_col2, 
                      outprefix = "~/work/yachie/droParser/tool_test/sciATAC//data_Aug2020/sciATAC_to_10x_rawsignal_scaletest",
                      both_rand_and_keep = F,
                      show_legend = F,
                      show_label=T,
                      ptsize=.2,
                      w = 5.4,
                      h=5.4)


mef2_ftp <- readRDS("~/work/yachie/droParser/tool_test/sciATAC/data_Aug2020/sciATAC_fly_Mef2_ftp.rds")
mef2_ftp <- mef2_ftp + theme(title = element_text(size=10),axis.title = element_text(size=10))
ggsave("~/work/yachie/droParser/tool_test/sciATAC//data_Aug2020/sciATAC_to_10x_Mef2.ftp.jpeg",width = 5.96,height = 5.96,units="cm",dpi=500)

mef2_cvp <- readRDS("~/work/yachie/droParser/tool_test/sciATAC/data_Aug2020/sciATAC_fly_Mef2_cvp.rds")
ggsave("~/work/yachie/droParser/tool_test/sciATAC//data_Aug2020/sciATAC_to_10x_Mef2.cvp.jpeg",mef2_cvp,width = 8,height = 16,units="cm",dpi=500)

ftz_cvp <- readRDS("~/work/yachie/droParser/tool_test/sciATAC/data_Aug2020/sciATAC_fly_ftz_cvp.rds")
ftz_cvp <- ftz_cvp + scale_x_continuous(breaks = c(6854000,6860000,6866000,6871000))
ggsave("~/work/yachie/droParser/tool_test/sciATAC//data_Aug2020/sciATAC_to_10x_ftz.cvp.jpeg",ftz_cvp,width = 8,height = 16,units="cm",dpi=500)


ftp <- readRDS("~/work/yachie/droParser/tool_test/sciATAC/data_Aug2020/sciATAC_fly_ftz_ftz_ftp.rds")
ftp <- ftp + theme(plot.title = element_blank(),
                   legend.text = element_text(size=10),
                   legend.key.size = grid::unit(0.6, "lines"),
                   axis.title = element_blank(),
                   axis.text = element_blank(),
                   axis.line = element_blank(),
                   axis.ticks = element_blank())+
        NoLegend()
ggsave("~/work/yachie/droParser/tool_test/sciATAC//data_Aug2020/sciATAC_to_10x_ftz_noleg_0916.ftp.jpeg",ftp,width = 5,height = 5,units="cm",dpi=500)


biological_dimp <- readRDS("~/work/yachie/droParser/tool_test/sciATAC//data_Aug2020/umap_6to8_biological.rds")
biological_dimp <- biological_dimp+
        scale_color_manual(values = brewer.pal(10,"Set3"))+
        theme(axis.title = element_blank(),
              axis.text = element_blank(),
              axis.line = element_blank(),
              axis.ticks = element_blank(),
              legend.text = element_text(size=10))+
        NoLegend()
ggsave("~/work/yachie/droParser/tool_test/sciATAC//data_Aug2020/umap_6to8_biological.jpeg",
       biological_dimp,
       width = 10,
       height = 10,
       units = "cm",
       dpi=500)
#Drawing 1005
for(g in c("CR44226","CG11249","CR44677","GATAe","btl","scrt","CG6415","CG4393","CG5225","wor")){
        ftp <- readRDS(paste0("~/work/yachie/droParser/tool_test/sciATAC/data_Aug2020/ftp/sciATAC_fly_label_",g,"_ftp.rds"))
        ftp <- ftp + theme(title = element_text(size=10),
             axis.title = element_blank(),
             axis.text = element_blank(),
             axis.line = element_blank(),
             axis.ticks = element_blank()) 
        ggsave(paste0("~/work/yachie/droParser/tool_test/sciATAC//data_Aug2020/sciATAC_to_10x_",g,".label.ftp.jpeg"),width = 9,height = 9,units="cm",dpi=500)
}
for(g in c("Mef2","ftz")){
        ftp <- readRDS(paste0("~/work/yachie/droParser/tool_test/sciATAC/data_Aug2020/ftp/sciATAC_fly_label_",g,"_ftp.rds"))
        ftp1 <- ftp + theme(plot.title = element_blank(),
                           legend.text = element_text(size=10),
                           legend.key.size = grid::unit(0.6, "lines"),
                           axis.title = element_blank(),
                           axis.text = element_blank(),
                           axis.line = element_blank(),
                           axis.ticks = element_blank()) + NoLegend()
        ggsave(paste0("~/work/yachie/droParser/tool_test/sciATAC//data_Aug2020/sciATAC_to_10x_",g,".label.ftp.jpeg"),ftp1,width = 8,height = 8,units="cm",dpi=500)
        leg_ftp <- ftp + theme(title = element_text(size=10),axis.title = element_text(size=10),legend.key.size = grid::unit(0.6, "lines"))
        ggsave(paste0("~/work/yachie/droParser/tool_test/sciATAC//data_Aug2020/sciATAC_to_10x_",g,".legend.ftp.jpeg"),leg_ftp,width = 4,height = 4,units="cm",dpi=500)
        
}


#Feb2021
setwd("~/work/yachie/droParser/tool_test/paper_data/scRNA/")
rand_10xto10x1 <- read.table("CCdist_0429/fly_atac_to10x_CCDist.tsv",header = T)

L <- list(fly_atac_to10x=rand_10xto10x1)

for(i in names(L)){
        g <- fnc_draw_scatter(df_rand = L[[i]],
                              colset = my_col, 
                              outprefix = paste0("../scATAC/fig_CCdist_May2021/",i),
                              both_rand_and_keep = F,
                              show_legend = F,
                              show_label=F,
                              w=6.21,
                              h=6.21)
}

L2 <- list(fly_atac_to10x=c("CCdist_0429/fly_atac_to10x.diff.rds","CCdist_0429/fly_atac_to10x.scr.diff.rds"))

for(i in names(L2)){
        g <- vln_plot(rand_rds_path = L2[[i]][1],
                      rand_rds_scr_path = L2[[i]][2],
                      keep = F,
                      colset = my_col[1],
                      outdir = "../scATAC/Utest_jitter/",
                      outname=i,
                      w=1.23,
                      h=2.5)
        
}

x <- readRDS("CCdist_0429/fly_atac_to10x.diff.rds")
x$case_diff %>% median

y <- readRDS("CCdist_0429/fly_atac_to10x.scr.diff.rds")
y$scramble_diff %>% median


#2022 March
g <- diff_rank_sina(rand_ccdist_path = "CCdist_0429/fly_atac_to10x_CCDist.tsv",
                    keep = F,
                    colset = my_col,
                    outdir = "Utest_jitter_220307/",
                    outname="sciATAC_to_10x",
                    w=0.76,
                    h = 1.06)











