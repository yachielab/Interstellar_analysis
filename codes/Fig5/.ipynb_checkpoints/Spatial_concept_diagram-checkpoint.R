setwd("~/work/yachie/YUSUKE_KIJIMA.LAB/Experiments/Projects/droParser/tool_test/slide_to_visium/")
library(RColorBrewer)
library(viridis)
library(ggthemes)

visium.orig <- read.table("visium_spot.tsv",sep="\t",header=F)
g.orig <- ggplot(visium.orig,aes(V2,V3))+
  theme_classic()+
  geom_point(size=.01,col="black",shape=16)+
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank())


aln <- read.table("20200507_bc_expansion/slide_to_vis_exp10_aligned.tsv",sep="\t",header=T)
cvec <- rep("0",length(aln$d_ind_bc))
names(cvec) <- aln$d_ind_bc
cvec[filter(aln,x>40000&x<42500&y>20000&y<22500)$d_ind_bc %>% as.character] <- "1"
cvec <- as.character(cvec)
aln$col <- cvec
colset <- c("0"="grey30","1"="blue")
ggplot(aln,aes(x,y,col=col))+
  geom_point()+
  theme_classic()+
  NoLegend()+
  scale_color_manual(values = colset)

g_aligned <- ggplot(aln %>% filter(col=="1"),aes(x,y,col=col))+
  geom_point(size=1,shape=16)+
  theme_classic()+
  NoLegend()+
  scale_color_manual(values = colset)+
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank())+
  NoLegend()
ggsave("Fig5_concept/slide_region5_2.pdf",g,width = 6,height = 6,units = "cm")


bc_use <- filter(aln,col=="1")$barcode %>% as.character()
slide.hippo.orig <- read.table("20200507_bc_expansion/slide_to_vis_exp10_expanded.tsv",header=T)
slide.hippo.orig$grp <- "not_use"
slide.hippo.orig$grp[slide.hippo.orig$barcode %in% bc_use] <- "use"
g.hippo <- ggplot(slide.hippo.orig,aes(x,y,col=grp))+
  theme_classic()+
  geom_point(size=.01,shape=16)+
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank())+
  scale_color_manual(values = c("navy","deeppink2"))+
  NoLegend()


# ggsave("Fig5_concept/visium_orig_coord.pdf",g.orig,width = 6,height = 6,units = "cm")
ggsave("Fig5_concept/hippo_orig_coord.jpeg",g.hippo,width = 8,height = 8,units = "cm",dpi=600)

# g.hippo.filt <- ggplot(slide.hippo.orig %>% filter(grp=="use"),aes(x,y))+
#   theme_classic()+
#   geom_point(size=1,shape=16,col="salmon")+
#   theme(axis.line = element_blank(),
#         axis.ticks = element_blank(),
#         axis.title = element_blank(),
#         axis.text = element_blank())+
#   # scale_color_manual(values = c("grey","salmon"))+
#   NoLegend()

get_nearest_coord

df.filt <- slide.hippo.orig %>% filter(grp=="use")
df.filt <- filter(aln,col=="1")


tiled <- read.table("test/test_hippo_exp10_tiled_visium.tsv")
visium.orig.moved <- filter(tiled,V1<=max(df.filt$x)&V1>=min(df.filt$x)&V2<=max(df.filt$y)&V2>=min(df.filt$y))

g.orig <- ggplot(visium.orig.moved,aes(V1,V2))+
  theme_classic()+
  geom_point(size=2,col="black",shape=16)+
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank())

df.merged <- data.frame(x=c(visium.orig.moved$V1,df.filt$x),
                        y=c(visium.orig.moved$V2,df.filt$y),
                        category=c(rep("Visium",nrow(visium.orig.moved)),
                                   rep("Slide-seq",nrow(df.filt))))
df.merged$category <- factor(df.merged$category,levels = c("Slide-seq","Visium"))
g <- ggplot(df.merged,aes(x,y,col=category,size=category,shape=category))+
  theme_classic()+
  geom_point(data=filter(df.merged,category=="Slide-seq"))+
  geom_point(data=filter(df.merged,category=="Visium"))+
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank())+
  scale_color_manual(values = c("orange","black"))+
  scale_size_manual(values = c(1,2))+
  scale_shape_manual(values = c(16,21))+
  NoLegend()
# ggsave("Fig5_concept/visium_slide_region5_2.pdf",g,width = 6,height = 6,units = "cm")
ggsave("Fig5_concept/visium_slide_region5_2.adjusted.pdf",g,width = 6,height = 6,units = "cm")



ggplot(tiled %>% filter(V1>7500&V1<10000&V2<11000),aes(x=V1,y=V2,col=V3))+
  theme_classic()+
  geom_point(size=1)+NoLegend()+coord_fixed()+
  scale_color_manual(values = adjustcolor(c("red","blue","yellow","pink"),alpha.f = 0.5))
filter(tiled,V1>40240) %>% View
