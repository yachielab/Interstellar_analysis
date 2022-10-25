#0828
source("../Fig4_S3_S4/Fig4_source.R")
my_col <- adjustcolor(c("blue3","blue3"),alpha.f = 0.1)
my_col2 <- adjustcolor(c("blue3","blue3"),alpha.f = 1)



# scatter plot
rand_10xto10x1 <- read.table("./fly_atac_to10x_CCDist.tsv",header = T)

L <- list(fly_atac_to10x=rand_10xto10x1)

for(i in names(L)){
        g <- fnc_draw_scatter(df_rand = L[[i]],
                              colset = my_col, 
                              outprefix = paste0("./",i),
                              both_rand_and_keep = F,
                              show_legend = F,
                              show_label=F,
                              w=6.21,
                              h=6.21)
}

# sina plot
g <- diff_rank_sina(rand_ccdist_path = "./fly_atac_to10x_CCDist.tsv",
                    keep = F,
                    colset = my_col,
                    outdir = "./",
                    outname="sciATAC_to_10x",
                    w=0.76,
                    h = 1.06)











