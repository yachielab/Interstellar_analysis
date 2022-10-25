library(Seurat)
library(Signac)
library(Matrix)
library(tidyverse)
library(rtracklayer)
library(patchwork)
library(RColorBrewer)

#Set unified cell barcode to dge
refresh_dge <- function(raw_dge,vec_source="",vec_dest="",direct=FALSE,header=FALSE,src=1,dest=2,path="",add=""){
    if(direct){
        df <- read.table(path,header=header)
        vec_source=df[,src]
        vec_dest <- df[,dest]
        if(add!=""){
            vec_dest <- paste0(vec_dest,add)
        }
        
        v <- vec_source
        names(v) <- vec_dest
        new_colname <- v[colnames(raw_dge)] %>% unname %>% as.character
        colnames(raw_dge) <- new_colname
    }else{
        if(add!=""){
            vec_source <- paste0(vec_source,add)
            vec_dest   <- paste0(vec_dest,  add)
        }
        v <- vec_source
        names(v) <- vec_dest
        new_colname <- v[colnames(raw_dge)] %>% unname %>% as.character
        colnames(raw_dge) <- new_colname
    }
    return(raw_dge)
}

#Read dge for 10x output
read_dge_10x <- function(target_dir){
    dge <- Read10X(target_dir)
    dge <- dge[,colSums(dge)>0]
    return(dge)
}

#Read dge for DStools output
read_dge_drop <- function(target_dge_txt){
    dge <- read.table(target_dge_txt,header = T)
    dge <- column_to_rownames(dge,var="GENE")
    dge <- Matrix(as.matrix(dge))
    dge <- dge[,colSums(dge)>0]
    return(dge)
}

#get intersect cells
get_cells_used <- function(dge_base,
                            target_dge_list,
                             min_features=200,
                             num_feature_rna_min=200,
                             min_cells = 2,
                             num_feature_rna_max = 5000,
                             mito_pattern="^mt-",
                             percent_mito=5,
                             seurat_processing_script="/home/ha5656/work/script_general/R/sc_analysis.R"){
    #reading processing script
    source(seurat_processing_script)
    
    #generate base seurat object
    base.seu <- seurat_preprocessing(dge_base,
                                   min_features=min_features,
                                   num_feature_rna_min=num_feature_rna_min,
                                   min_cells = min_cells,
                                   num_feature_rna_max = num_feature_rna_max,
                                   mito_pattern=mito_pattern,
                                   percent_mito=percent_mito)
    #check mitochondrial genes
    print(rownames(base.seu)[grep(mito_pattern,rownames(base.seu))])

    base.cells <- colnames(base.seu)
    
    #get intersection
    out.cells <- base.cells
    for(dge in target_dge_list){
        out.cells <- intersect(colnames(dge),out.cells)
    }
    
	message("Parameters:")
	message(paste("min_features:",min_features))
	message(paste("num_feature_rna_min:",num_feature_rna_min))
	message(paste("num_feature_rna_max:",num_feature_rna_max))
	message(paste("num_feature_rna_min",num_feature_rna_min))
	message(paste("mito_pattern",mito_pattern))
	message(paste("percent_mito",percent_mito))
	message("")
	message(paste("Number of cells after QCing original:",length(base.cells)))
	message(paste("Number of cells after taking intersection:",length(out.cells)))
	
    return(out.cells)
}



getVarGenes_untilPCA <- function(seu,var_genes=NULL,num_varFeatures=2500){
    seu <- NormalizeData(seu,normalization.method = "LogNormalize", scale.factor = 10000,verbose = F)
    if(is.null(var_genes)){
        seu <- FindVariableFeatures(seu,selection.method="vst",nfeatures = num_varFeatures,verbose = F)
        return(VariableFeatures(seu))
    }else{
        VariableFeatures(seu) <- var_genes
        seu <- ScaleData(seu,verbose=F)
        seu <- RunPCA(seu,verbose = F)
        return(seu)
    }
}


#Draw and store UMAP
draw_umap <- function(seu_base,seu_to10x,seu_todrop,col_seed,outdir,outname,dot_pt=.05,w=3.5,outrds=TRUE,h=3.5,dpi=400,is_pdf=FALSE,units="cm",is_bmp=FALSE,celltype_v=""){
    library(RColorBrewer)
    
    if(all(!colnames(seu_base)==colnames(seu_to10x))||!all(colnames(seu_base)==colnames(seu_todrop))){
        stop("Colnames of base, to10x and todrop are inconsistent!")
    }

    if(length(celltype_v)==0){
        seu_to10x$cluster_original <- seu_base$`RNA_snn_res.0.6`
        seu_todrop$cluster_original <- seu_base$`RNA_snn_res.0.6`
    }else{
        seu_to10x$cluster_original <- celltype_v
        seu_todrop$cluster_original <- celltype_v
    }
    
    #Color code
    n <- length(unique(seu_to10x$cluster_original))
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    set.seed(col_seed)
    my_col <- sample(col_vector, n)
    
    #Draw
    g_base <- DimPlot(seu_base,group.by='RNA_snn_res.0.6',cols=my_col,pt.size=dot_pt)+
    theme_void()+
    theme(axis.title = element_blank(),axis.ticks = element_blank(),axis.text=element_blank())+
    NoLegend()+
	ggtitle("")

    g_to10x <- DimPlot(seu_to10x,group.by='cluster_original',cols=my_col,pt.size=dot_pt)+
    theme_void()+
    theme(axis.title = element_blank(),axis.ticks = element_blank(),axis.text=element_blank())+
    NoLegend()+
	ggtitle("")


    g_todrop <- DimPlot(seu_todrop,group.by='cluster_original',cols=my_col,pt.size=dot_pt)+
    theme_void()+
    theme(axis.title = element_blank(),axis.ticks = element_blank(),axis.text=element_blank())+
    NoLegend()+
	ggtitle("")

    
    #Export
    outprefix <- paste0(outdir,"/",outname)
    
    if(!is_pdf && !is_bmp){
        ggsave(paste0(outprefix,"_original_umap.jpeg"),g_base,width=w,height=h,units="cm",dpi=dpi)
        ggsave(paste0(outprefix,"_to_10x_umap.jpeg"),g_to10x,width=w,height=h,units="cm",dpi=dpi)
        ggsave(paste0(outprefix,"_to_drop_umap.jpeg"),g_todrop,width=w,height=h,units="cm",dpi=dpi)
    }else if(is_bmp){
        ggsave(paste0(outprefix,"_original_umap.bmp"),g_base,width=w,height=h,units="cm",dpi=400)
        ggsave(paste0(outprefix,"_to_10x_umap.bmp"),g_to10x,width=w,height=h,units="cm",dpi=400)
        ggsave(paste0(outprefix,"_to_drop_umap.bmp"),g_todrop,width=w,height=h,units="cm",dpi=400)
    }else{
        ggsave(paste0(outprefix,"_original_umap.pdf"),g_base,width=w,height=h,units="cm")
        ggsave(paste0(outprefix,"_to_10x_umap.pdf"),g_to10x,width=w,height=h,units="cm")
        ggsave(paste0(outprefix,"_to_drop_umap.pdf"),g_todrop,width=w,height=h,units="cm")
    }
    
    if(outrds){
        saveRDS(g_base,paste0(outprefix,"_original_umap.rds"))
        saveRDS(g_to10x,paste0(outprefix,"_to10x_umap.rds"))
        saveRDS(g_todrop,paste0(outprefix,"_todrop_umap.rds"))
    }
    
    
}

draw_umap_crazy10x <- function(seu_base,seu_to10x,col_seed,outdir,outname,dot_pt=.05,w=3.5,h=3.5,dpi=600,is_pdf=FALSE){
    library(RColorBrewer)
    
    if(all(!colnames(seu_base)==colnames(seu_to10x))){
        stop("Colnames of base and to10x are inconsistent!")
    }
    seu_to10x$cluster_original <- seu_base$`RNA_snn_res.0.6`
    
    #Color code
    n <- length(unique(seu_to10x$cluster_original))
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    set.seed(col_seed)
    my_col <- sample(col_vector, n)
    
    #Draw
    g_base <- DimPlot(seu_base,group.by='RNA_snn_res.0.6',cols=my_col,pt.size=dot_pt)+
    theme_void()+
    theme(axis.title = element_blank(),axis.ticks = element_blank(),axis.text=element_blank())+
    NoLegend()+
	ggtitle("")


    g_to10x <- DimPlot(seu_to10x,group.by='cluster_original',cols=my_col,pt.size=dot_pt)+
    theme_void()+
    theme(axis.title = element_blank(),axis.ticks = element_blank(),axis.text=element_blank())+
    NoLegend()+
	ggtitle("")

    
    #Export
    outprefix <- paste0(outdir,"/",outname)
    
    if(!is_pdf){
        ggsave(paste0(outprefix,"_original_umap.jpeg"),g_base,width=w,height=h,units="cm",dpi=dpi)
        ggsave(paste0(outprefix,"_to_10x_umap.jpeg"),g_to10x,width=w,height=h,units="cm",dpi=dpi)
    }else{
        ggsave(paste0(outprefix,"_original_umap.pdf"),g_base,width=w,height=h,units="cm")
        ggsave(paste0(outprefix,"_to_10x_umap.pdf"),g_to10x,width=w,height=h,units="cm")
    }
    
}

#Distace calculation 
caldist <- function(pos,posvec1,posvec2,df_input){
    res <- tryCatch({
        v1 <- df_input[,posvec1[pos]]
#     print(posvec2[pos])
        v2 <- df_input[,posvec2[pos]]
    },error=function(e){
        print(e)
        print(pos)
        print(posvec1[pos])
        print(posvec2[pos])
        stop()
    })
    
    return(sqrt(sum((v1-v2)^2)))
}

#sampling
sampling_combination <- function(seedval,cells,n_cell_sampling){
    set.seed(seedval)
    x <- cells
    N <- 2
    out <- list()
    while(length(out) < n_cell_sampling) {
      out <- c(out,
               unique(replicate(n_cell_sampling - length(out), sort(sample(x, N)), simplify = FALSE)))
      out <- unique(out)
    }
    sampled_df <- data.frame(out) %>% t 
    colnames(sampled_df) <- c("c1","c2")
    return(sampled_df)
}

scr_sampling.subfun <- function(x,n,pool){
  seed <- strsplit(x,":")[[1]][2] %>% as.integer()
  x <- strsplit(x,":")[[1]][1]
  
  out <- c()
  n_sampling_now <- n
  while(length(out)<n){
    pool <- setdiff(pool,c(out,x))
    n_sampling_now <- n_sampling_now-length(out)
    set.seed(seed)
    out <- c(out,sample(size=n_sampling_now,x=pool,replace = F))
  }
  out <- paste(out,collapse = ":")
  out <- paste(x,out,sep=";")
  return(out)
}

#scrable sampling
scr_sampling <- function(random_50000_pair_1,random_50000_pair_2,cells,n_scramble_sampling,anotherset=FALSE){
#     redundant_v <- c(rep(random_50000_pair_1,each=scramble_sampling),
#                     rep(random_50000_pair_2,each=scramble_sampling))
    df <- data.frame(V1=c(random_50000_pair_1,random_50000_pair_2),stringsAsFactors=F)
#     print(head(df))
#     print(tail(df))
    if(anotherset){
        df$V1 <- paste(df$V1,(nrow(df)+1):(2*nrow(df)),sep=":")
    }else{
        df$V1 <- paste(df$V1,1:nrow(df),sep=":")
    }
   
    df <- apply(df,c(1,2),scr_sampling.subfun,n=n_scramble_sampling,cells)
    df <- str_split_fixed(df,";",2)
    df <- cbind(df[,1],
                rep(paste0(random_50000_pair_1,":",random_50000_pair_2),2),
                c(rep("left",length(random_50000_pair_1)),rep("right",length(random_50000_pair_1))),
                str_split_fixed(df[,2],":",n_scramble_sampling)) %>% as.data.frame()
    df <- pivot_longer(df,col=colnames(df)[4:ncol(df)])
    return(df[,c(1,2,3,5)])
}

#Quantitative analysis
get_dist_pcaspace_scatter <- function(seu_base,
                                      seu_target,
                                      outdir,
                                      outname,
                                      n_cell_sampling=50000,
                                      n_scramble_sampling=0,
                                      scramble=FALSE,
                                      rowcount=TRUE,
                                      pca=FALSE,
                                      lsi=FALSE,
                                      expression=FALSE,
                                      rawsignal=FALSE,
                                      use_dimension_base=c(),
                                     use_dimension_target=c()){
    print("start...")
    
    if(rowcount){
        orig.pcaspace   <- seu_base@assays$RNA@scale.data[VariableFeatures(seu_base),]
        target.pcaspace <- seu_target@assays$RNA@scale.data[VariableFeatures(seu_target),]
    }else if(lsi){
        orig.pcaspace   <- seu_base@reductions$lsi@cell.embeddings        
        target.pcaspace <- seu_target@reductions$lsi@cell.embeddings 
        orig.pcaspace   <- orig.pcaspace[,use_dimension_base]
        target.pcaspace <- target.pcaspace[,use_dimension_target]
    }else if(pca){
        orig.pcaspace   <- seu_base@reductions$pca@cell.embeddings
        target.pcaspace <- seu_target@reductions$pca@cell.embeddings 
    }else if(rawsignal){
        orig.pcaspace   <- seu_base@assays$peaks@counts
        target.pcaspace <- seu_target@assays$peaks@counts
    }
    
    if(!all(colnames(seu_base)==colnames(seu_target))){
        stop("Colnames of base and target are inconsistent!")
    }

    #Random sampling of cells
    if(rowcount || rawsignal){
        print("start sampling...")
        sampled_df <- sampling_combination(1111,colnames(orig.pcaspace),n_cell_sampling) %>% as_tibble
        random_50000_pair_1 <- sampled_df$c1 %>% as.character
        random_50000_pair_2 <- sampled_df$c2 %>% as.character
        
#         sampled_df <- sampling_combination(2222,colnames(orig.pcaspace),n_cell_sampling) %>% as_tibble
#         random_50000_pair_3 <- sampled_df$c1
#         random_50000_pair_4 <- sampled_df$c2
#         set.seed(1212)
#         random_50000_pair_1 <- sample(1:ncol(orig.pcaspace),size = n_cell_sampling,replace = T)
#         set.seed(2121)
#         random_50000_pair_2 <- sample(1:ncol(orig.pcaspace),size = n_cell_sampling,replace = T)
#         set.seed(1111)
#         random_50000_pair_3 <- sample(1:ncol(orig.pcaspace),size = n_cell_sampling,replace = T)
#         set.seed(2222)
#         random_50000_pair_4 <- sample(1:ncol(orig.pcaspace),size = n_cell_sampling,replace = T)
    }else{
        sampled_df <- sampling_combination(1111,rownames(orig.pcaspace),n_cell_sampling) %>% as_tibble
        random_50000_pair_1 <- sampled_df$c1 %>% as.character
        random_50000_pair_2 <- sampled_df$c2 %>% as.character
        
#         sampled_df <- sampling_combination(2222,rownames(orig.pcaspace),n_cell_sampling) %>% as_tibble
#         random_50000_pair_3 <- sampled_df$c1 %>% as.character
#         random_50000_pair_4 <- sampled_df$c2 %>% as.character
        
#         set.seed(1212)
#         random_50000_pair_1 <- sample(1:nrow(orig.pcaspace),size = n_cell_sampling,replace = T)
#         set.seed(2121)
#         random_50000_pair_2 <- sample(1:nrow(orig.pcaspace),size = n_cell_sampling,replace = T)
#         set.seed(1111)
#         random_50000_pair_3 <- sample(1:nrow(orig.pcaspace),size = n_cell_sampling,replace = T)
#         set.seed(2222)
#         random_50000_pair_4 <- sample(1:nrow(orig.pcaspace),size = n_cell_sampling,replace = T)
    }
    
#     print(length(unique(random_50000_pair_1)))
#     print(length(unique(random_50000_pair_2)))
#     print(length(unique(random_50000_pair_3)))
#     print(length(unique(random_50000_pair_4)))
#     print("##")
    #Distance dataframe
    if(((rowcount)&&(!scramble))||((rawsignal)&&(!scramble))){
        df.orig.pcaspace <- data.frame(dist_orig=seq(1,n_cell_sampling))
        df.orig.pcaspace <- apply(df.orig.pcaspace,
                                  c(1,2),
                                  caldist,
                                  random_50000_pair_1,
                                  random_50000_pair_2,
                                  orig.pcaspace)
        df.orig.pcaspace <- as.data.frame(df.orig.pcaspace)
        df.orig.pcaspace$combi <- paste0(random_50000_pair_1,":",random_50000_pair_2)
        
        df.target.pcaspace <- data.frame(dist_target=seq(1,n_cell_sampling))
        df.target.pcaspace <- apply(df.target.pcaspace,
                                   c(1,2),
                                   caldist,
                                   random_50000_pair_1,
                                   random_50000_pair_2,
                                   target.pcaspace)
        df.target.pcaspace <- as.data.frame(df.target.pcaspace)
        df.target.pcaspace$combi <- paste0(random_50000_pair_1,":",random_50000_pair_2)
        
        df.merge <- merge(df.orig.pcaspace,df.target.pcaspace,by.x=2,by.y=2)
        
    }else if(!scramble){
        df.orig.pcaspace <- data.frame(dist_orig=seq(1,n_cell_sampling))
        df.orig.pcaspace <- apply(df.orig.pcaspace,
                                  c(1,2),
                                  caldist,
                                  random_50000_pair_1,
                                  random_50000_pair_2,
                                  t(orig.pcaspace))
        df.orig.pcaspace <- as.data.frame(df.orig.pcaspace)
        df.orig.pcaspace$combi <- paste0(random_50000_pair_1,":",random_50000_pair_2)

        df.target.pcaspace <- data.frame(dist_target=seq(1,n_cell_sampling))
        df.target.pcaspace <- apply(df.target.pcaspace,
                                   c(1,2),
                                   caldist,
                                   random_50000_pair_1,
                                   random_50000_pair_2,
                                   t(target.pcaspace))
        df.target.pcaspace <- as.data.frame(df.target.pcaspace)
        df.target.pcaspace$combi <- paste0(random_50000_pair_1,":",random_50000_pair_2)
        
        df.merge <- merge(df.orig.pcaspace,df.target.pcaspace,by.x=2,by.y=2)

    }else if((scramble && rowcount)||(scramble && rawsignal)||(scramble && lsi)){
        pairs.orig <- scr_sampling(random_50000_pair_1 %>% as.character,
                                   random_50000_pair_2 %>% as.character,
                                   colnames(seu_base),
                                   n_scramble_sampling)
        pairs.target <- scr_sampling(random_50000_pair_1 %>% as.character,
                                     random_50000_pair_2 %>% as.character,
                                     colnames(seu_target),
                                     n_scramble_sampling,
                                     anotherset=T)

        df.orig.pcaspace <- data.frame(dist_orig=seq(1,nrow(pairs.orig)))
        if(lsi){
            df.orig.pcaspace <- apply(df.orig.pcaspace,
                                  c(1,2),
                                  caldist,
                                  pairs.orig$V1 %>% as.character,
                                  pairs.orig$value %>% as.character,
                                  t(orig.pcaspace))
        }else{
            df.orig.pcaspace <- apply(df.orig.pcaspace,
                                  c(1,2),
                                  caldist,
                                  pairs.orig$V1 %>% as.character,
                                  pairs.orig$value %>% as.character,
                                  orig.pcaspace)
        }
        
        df.orig.pcaspace <- as.data.frame(df.orig.pcaspace)
        df.orig.pcaspace$pairs <- pairs.orig$V2 %>% as.character
        df.orig.pcaspace$side <- pairs.orig$V3 %>% as.character
        
        df.orig.left <- filter(df.orig.pcaspace,side=="left")[,c("pairs","dist_orig")]
        colnames(df.orig.left) <- c("pairs","dist_orig_left")
        df.orig.right <- filter(df.orig.pcaspace,side=="right")[,c("pairs","dist_orig")]
        colnames(df.orig.right) <- c("pairs","dist_orig_right")

        df.target.pcaspace <- data.frame(dist_target=seq(1,nrow(pairs.target)))
        if(lsi){
            df.target.pcaspace <- apply(df.target.pcaspace,
                                   c(1,2),
                                   caldist,
                                   pairs.target$V1 %>% as.character,
                                  pairs.target$value %>% as.character,
                                   t(target.pcaspace))
        }else{
            df.target.pcaspace <- apply(df.target.pcaspace,
                                   c(1,2),
                                   caldist,
                                   pairs.target$V1 %>% as.character,
                                  pairs.target$value %>% as.character,
                                   target.pcaspace)
        }
        
        df.target.pcaspace <- as.data.frame(df.target.pcaspace)
        df.target.pcaspace$pairs <- pairs.target$V2 %>% as.character
        df.target.pcaspace$side <- pairs.target$V3 %>% as.character
        
        df.target.left <- filter(df.target.pcaspace,side=="left")[,c("pairs","dist_target")]
        colnames(df.target.left) <- c("pairs","dist_target_left")
        df.target.right <- filter(df.target.pcaspace,side=="right")[,c("pairs","dist_target")]
        colnames(df.target.right) <- c("pairs","dist_target_right")
        
        if(!all((df.orig.left$pairs %>% as.character)==(df.target.left$pairs %>% as.character))){
            stop("inconsistent!")
        }
        if(!all((df.orig.right$pairs %>% as.character)==(df.target.right$pairs %>% as.character))){
            stop("inconsistent!")
        }
        diff_left <- data.frame(pairs=df.orig.left$pairs %>% as.character,
                               diff_left=abs(df.orig.left$dist_orig_left - df.target.left$dist_target_left),
                                stringsAsFactors=F)
        diff_right <- data.frame(pairs=df.orig.right$pairs %>% as.character,
                               diff_right=abs(df.orig.right$dist_orig_right - df.target.right$dist_target_right),
                                stringsAsFactors=F)
        diff_left <- aggregate(diff_left$diff_left,list(diff_left$pairs),mean)
        diff_right<- aggregate(diff_right$diff_right,list(diff_right$pairs),mean)
        
        
        
        df.merge <- merge(diff_left,diff_right,by.x=1,by.y=1)
        df.out <- data.frame(pairs=df.merge$Group.1,scramble_diff=(df.merge$x.x+df.merge$x.y)/2,stringsAsFactors=F)
        
        return(df.out)

    }else{
        df.orig.pcaspace <- data.frame(v=seq(1,n_cell_sampling))
        df.orig.pcaspace <- apply(df.orig.pcaspace,
                                  c(1,2),
                                  caldist,
                                  random_50000_pair_1,
                                  random_50000_pair_2,
                                  t(orig.pcaspace))
        df.target.pcaspace <- data.frame(v=seq(1,n_cell_sampling))
        df.target.pcaspace <- apply(df.target.pcaspace,
                                   c(1,2),
                                   caldist,
                                   random_50000_pair_3,
                                   random_50000_pair_4,
                                   t(target.pcaspace))
    
    }
    
    
    #Set dataframe fow drawing
    df.dist.orig_vs_target  <- data.frame(original=df.merge$dist_orig,
                                          target=df.merge$dist_target)
    rownames(df.dist.orig_vs_target) <- df.merge$combi %>% as.character
        
    #Output
    outprefix <- paste0(outdir,"/",outname)
    write.table(df.dist.orig_vs_target,paste0(outprefix,"_CCDist.tsv"),sep="\t",col.names=T,row.names=T,quote=F)
    
    #For violin
    v.vln <- data.frame(pairs=rownames(df.dist.orig_vs_target),
                        case_diff=abs(df.dist.orig_vs_target$original-df.dist.orig_vs_target$target))
#     rownames(v.vln) <- rownames(df.dist.orig_vs_target)
    
    return(v.vln)
}

#ATAC preprocess
atac_data_build <- function(raw_cnt="",h5_path="",sccsv="",frag="",sci=FALSE){
    if(!sci){
        counts <- Read10X_h5(h5_path)
        metadata <- read.csv(
          file = sccsv,
          header = TRUE,
          row.names = 1
        )

        dm_assay <- CreateChromatinAssay(
          counts = counts,
          sep = c(":", "-"),
          fragments = frag,
          min.cells = 1
        )

        dm <- CreateSeuratObject(
          counts = dm_assay,
          assay = 'peaks',
          project = 'ATAC',
          meta.data = metadata
        )
    }else{
        dm <- CreateSeuratObject(
          counts = raw_cnt,
          assay = 'peaks',
          project = 'ATAC',
          min.cells = 1
        )
    }
        
    return(dm)
}

atac_data_qc <- function(dm,gtf,min_thresh_nucleosome_signal=4,region='2L-1-23513712',min_TSS.enrichment=2){
    print("Start extracting signal...")
    dm <- NucleosomeSignal(object = dm)
    dm$nucleosome_group <- ifelse(dm$nucleosome_signal > min_thresh_nucleosome_signal,paste0('NS>',min_thresh_nucleosome_signal),paste0('NS<',min_thresh_nucleosome_signal))
#     dm$nucleosome_group <- ifelse(dm$nucleosome_signal > min_thresh_nucleosome_signal,"NS>4","NS<4")
    g1 <- FragmentHistogram(object = dm, group.by = 'nucleosome_group', region = region)
    show(g1)
    
    print("Generating reference...")
    gtf <- rtracklayer::import(gtf)
    gene.coords <- gtf[gtf$type == 'gene']
    seqlevelsStyle(gene.coords) <- 'Ensembl'
    gene.ranges <- keepStandardChromosomes(gene.coords, pruning.mode = 'coarse')
    tss.ranges <- GRanges(
      seqnames = seqnames(gene.ranges),
      ranges = IRanges(start = start(gene.ranges), width = 2),
      strand = strand(gene.ranges)
    )
    seqlevelsStyle(tss.ranges) <- 'Ensembl'
    tss.ranges <- keepStandardChromosomes(tss.ranges, pruning.mode = 'coarse')

    print("Start extracting TSS...")
    dm <- TSSEnrichment(dm,tss.positions = tss.ranges, fast = FALSE)
    dm$high.tss <- ifelse(dm$TSS.enrichment > min_TSS.enrichment, 'High', 'Low')
    g2 <- TSSPlot(dm, group.by = 'high.tss') + NoLegend()
    show(g2)
    
    dm$pct_reads_in_peaks <- dm$peak_region_fragments / dm$passed_filters * 100
    dm$blacklist_ratio <- dm$blacklist_region_fragments / dm$peak_region_fragments
    g3 <- VlnPlot(
      object = dm,
      features = c('pct_reads_in_peaks', 'peak_region_fragments',
                   'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
      pt.size = 0.1,
      ncol = 5
    )
    show(g3)
    
    return(dm)
}

atac_data_filtering <- function(dm,
                                min_peak_region_fragments=200,
                                max_peak_region_fragments=100000,
                                min_pct_reads_in_peaks=40,
                                max_blacklist_ratio=0.025,
                                max_nucleosome_signal=4,
                                min_TSS.enrichment=1,
                                max_TSS.enrichment=2){
    dm <- subset(
      x = dm,
      subset = peak_region_fragments > min_peak_region_fragments &
        peak_region_fragments < max_peak_region_fragments &
#         pct_reads_in_peaks > min_pct_reads_in_peaks &
#         blacklist_ratio < max_blacklist_ratio &
        nucleosome_signal < max_nucleosome_signal &
        TSS.enrichment > min_TSS.enrichment 
#         TSS.enrichment < max_TSS.enrichment
    )
    return(dm)
}

atac_dimred_lsi <- function(fly_atac){
    fly_atac <- RunTFIDF(fly_atac)
    fly_atac <- FindTopFeatures(fly_atac, min.cutoff = 'q0')
    fly_atac <- RunSVD(
      object = fly_atac,
      assay = 'peaks',
      reduction.key = 'LSI_',
      reduction.name = 'lsi'
    )
    show(DepthCor(fly_atac))
    return(fly_atac)
}

atac_dimred_umap <- function(fly_atac,use_dims=c(2:30)){
     fly_atac <- RunUMAP(
      object = fly_atac,
      reduction = 'lsi',
      dims = use_dims
    )
    fly_atac <- FindNeighbors(
      object = fly_atac,
      reduction = 'lsi',
      dims = use_dims
    )
    fly_atac <- FindClusters(
      object = fly_atac,
      algorithm = 3,
      resolution = 1,
      verbose = FALSE
    )
    return(fly_atac)
}

filter_thick_col <- function(x){
  val <- col2rgb(x)
  if(sum(val)<600){return(x)}else{return(NA)}
}

draw_atac_umap <- function(seu_base,seu_to10x,col_seed,outdir,outname,dot_pt=.05,w=3.5,h=3.5,dpi=600){
    library(RColorBrewer)
    
    if(!all(colnames(seu_base)==colnames(seu_to10x))){
        stop("Colnames of base and to10x are inconsistent!")
    }
    seu_to10x$cluster_original <- seu_base$seurat_clusters
    
    #Color code
    n <- length(unique(seu_to10x$cluster_original))
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    col_vector <- sapply(col_vector,filter_thick_col)
    col_vector <- na.omit(col_vector)
    col_vector <- unname(col_vector)
    set.seed(col_seed)
    my_col <- sample(col_vector, n)
    
    #Draw
    g_base <- DimPlot(seu_base,group.by='seurat_clusters',cols=my_col,pt.size=dot_pt,label=F)+
    theme_void()+
    theme(axis.title = element_blank(),axis.ticks = element_blank(),axis.text=element_blank())+
    NoLegend()

    g_to10x <- DimPlot(seu_to10x,group.by='cluster_original',cols=my_col,pt.size=dot_pt,label=F)+
    theme_void()+
    theme(axis.title = element_blank(),axis.ticks = element_blank(),axis.text=element_blank())+
    NoLegend()
    
#     plot(g_base)
#     plot(g_to10x)
    #Export
    outprefix <- paste0(outdir,"/",outname)
    
    ggsave(paste0(outprefix,"_original_umap.jpeg"),g_base,width=w,height=h,units="cm",dpi=dpi)
    ggsave(paste0(outprefix,"_to_10x_umap.jpeg"),g_to10x,width=w,height=h,units="cm",dpi=dpi)
    
    
    #Draw
    g_base <- DimPlot(seu_base,group.by='seurat_clusters',cols=my_col,pt.size=dot_pt,label=T)+
    theme_void()+
    theme(axis.title = element_blank(),axis.ticks = element_blank(),axis.text=element_blank())+
    NoLegend()

    g_to10x <- DimPlot(seu_to10x,group.by='cluster_original',cols=my_col,pt.size=dot_pt,label=T)+
    theme_void()+
    theme(axis.title = element_blank(),axis.ticks = element_blank(),axis.text=element_blank())+
    NoLegend()
    
    ggsave(paste0(outprefix,"_original_umap.pdf"),g_base,width=w,height=h,units="cm")
    ggsave(paste0(outprefix,"_to_10x_umap.pdf"),g_to10x,width=w,height=h,units="cm")
    
#     saveRDS(g_base,paste0(outprefix,"_original_umap.rds"))
#     saveRDS(g_to10x,paste0(outprefix,"_to_10x_umap.rds"))
}

atac_genecellmatrix <- function(fly_atac,gtf,frag){
    print("Generating reference...")
    gtf <- rtracklayer::import(gtf)
    gene.coords <- gtf[gtf$type == 'gene']
    seqlevelsStyle(gene.coords) <- 'Ensembl'
    gene.ranges <- keepStandardChromosomes(gene.coords, pruning.mode = 'coarse')
    
    # Extend coordinates upstream to include the promoter
    genebodyandpromoter.coords <- Extend(x = gene.ranges, upstream = 2000, downstream = 0)
    
    print("build a gene by cell matrix with extended promoter region")
    frag <- CreateFragmentObject(frag,cells=colnames(fly_atac))
    # build a gene by cell matrix
    gene.activities <- FeatureMatrix(
      fragments = frag,
      features = genebodyandpromoter.coords,
      cells = colnames(fly_atac)
    )

    # convert rownames from chromsomal coordinates into gene names
    gene.key <- genebodyandpromoter.coords$gene_name
    names(gene.key) <- GRangesToString(grange = genebodyandpromoter.coords)
    rownames(gene.activities) <- make.unique(gene.key[rownames(gene.activities)])
    gene.activities <- gene.activities[rownames(gene.activities)!="",]

    #Add the gene activity matrix to the Seurat object as a new assay, and normalize it
    fly_atac[['RNA']] <- CreateAssayObject(counts = gene.activities)
    fly_atac <- NormalizeData(
      object = fly_atac,
      assay = 'RNA',
      normalization.method = 'LogNormalize',
      scale.factor = median(fly_atac$nCount_RNA)
    )
    
    return(fly_atac)
}

draw_atac_ftp <- function(fly_atac,seu_base,gene,outdir,outname,label.size=4,saveRDS=TRUE){
    if(!all(colnames(seu_base)==colnames(fly_atac))){
        stop("Colnames of base and to10x are inconsistent!")
    }
    DefaultAssay(fly_atac) <- "RNA"
    Idents(fly_atac) <- seu_base$seurat_clusters
    
    ftp <- FeaturePlot(
      object = fly_atac,
      label=F,
#         group.by='cluster_original',
      features = gene,
      pt.size = .1,
#         label.size=label.size,
        cols=c("grey90","deeppink2"),
      max.cutoff = 'q95')+
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          legend.text = element_text(size=7),
          text = element_blank(),
          legend.key.width = unit(0.2,"cm"),
         legend.key.height = unit(0.5,"cm"),
         plot.margin = unit(c(0,0,0,0),"cm"),
         legend.box.margin = unit(c(0,0,0,0),"cm"),
         legend.spacing.x = unit(c(0.25,0,0,0),"cm"))
    ggsave(paste0(outdir,"/",outname,"_",gene,"_ftp.pdf"),ftp,width=8,height=8,units="cm")
    ftp <- ftp+NoLegend()+theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          text = element_blank(),
          plot.margin = unit(c(0,0,0,0),"cm"))
    ggsave(paste0(outdir,"/",outname,"_",gene,"_ftp.jpeg"),ftp,width=6.7,height=8.46,units="cm",dpi=600)
    show(ftp)
    
    if(saveRDS){
        saveRDS(ftp,paste0(outdir,"/",outname,"_",gene,"_ftp.rds"))
    }
} 

draw_atac_cvp <- function(seu_base,seu_to10x,gene,gtf,col_seed,ncluster,outdir,outname,saveRDS=TRUE){
    gtf <- rtracklayer::import(gtf)
    gene.coords <- gtf[gtf$type == 'gene']
    seqlevelsStyle(gene.coords) <- 'Ensembl'
    gene.ranges <- keepStandardChromosomes(gene.coords, pruning.mode = 'coarse')
    
    DefaultAssay(seu_to10x) <- "peaks"
    region1 <- GRangesToString(subset(gene.ranges, gene_name==gene))
    
    seu_to10x$cluster_original <- seu_base$seurat_clusters
    
    cvp <- CoveragePlot(
      object = seu_to10x,
      region = region1,
      group.by="cluster_original",
        annotation=T,
#       sep = c(":", "-"),
#       idents=c("2","6"),
#       annotation = gene.ranges,
#       peaks = StringToGRanges(regions = rownames(seu_to10x), sep = c(":", "-")),
      extend.upstream = 11000,
      extend.downstream = 5000,
#       ncol = 1
    )
    plot(cvp)
    
    #Color code
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    col_vector <- sapply(col_vector,filter_thick_col)
    col_vector <- na.omit(col_vector)
    col_vector <- unname(col_vector)
    set.seed(col_seed)
    my_col <- sample(col_vector, ncluster)
    
    cvp1 <- cvp[[1]]+scale_fill_manual(values=my_col)
    cvp2 <- cvp[[2]]
    cvp3 <- cvp[[3]]
    
    cvp <- cvp1 + cvp2 + cvp3 + plot_layout(ncol = 1, heights = c(10,0.75, 1))
    
    ggsave(paste0(outdir,"/",outname,"_",gene,"_cvp.pdf"),cvp,width=9,height=17,units="cm")
    
    if(saveRDS){
        saveRDS(cvp,paste0(outdir,"/",outname,"_",gene,"_cvp.rds"))
    }
    plot(cvp)
    return(cvp)
}

ReadSPLiT <- function(mtx,meta,genes){
    meta <- read.csv(meta,header=T,stringsAsFactors=F)
    genes <- read.csv(genes,header=T,stringsAsFactors=F)
    mtx <- readMM(mtx)
    mtx <- t(mtx)
    
    cells <- meta$cell_barcode
    genes <- genes$gene_name
    
#     print(nrow(mtx))
#     print(ncol(mtx))
    
#     print(length(cells))
#     print(length(genes))
    rownames(mtx) <- genes
    colnames(mtx) <- cells
    
    return(mtx)
}
