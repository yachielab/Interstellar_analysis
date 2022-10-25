library(RColorBrewer)
library(viridis)
library(ggthemes)

visium.orig <- read.table("data/visium_spot.tsv",sep="\t",header=F)
g.orig <- ggplot(visium.orig,aes(V2,V3))+
  theme_classic()+
  geom_point(size=.01,col="black",shape=16)+
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank())


aln <- read.table("data/slide_to_vis_exp10_aligned.tsv",sep="\t",header=T)
cvec <- rep("0",length(aln$d_ind_bc))
names(cvec) <- aln$d_ind_bc
cvec[filter(aln,x>40000&x<42500&y>20000&y<22500)$d_ind_bc %>% as.character] <- "1"
cvec <- as.character(cvec)
aln$col <- cvec
colset <- c("0"="grey30","1"="blue")

bc_use <- filter(aln,col=="1")$barcode %>% as.character()
slide.hippo.orig <- read.table("data/slide_to_vis_exp10_expanded.tsv",header=T)
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


ggsave("hippo_orig_coord.jpeg",g.hippo,width = 8,height = 8,units = "cm",dpi=600)


df.filt <- slide.hippo.orig %>% filter(grp=="use")
df.filt <- filter(aln,col=="1")


tiled <- read.table("data/hippo_exp10_tiled_visium.tsv")
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
ggsave("visium_slide_region5_2.adjusted.pdf",g,width = 6,height = 6,units = "cm")
