library(Matrix)
library(dplyr)
library(ggplot2)
library(Seurat)
library(RColorBrewer)

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

#scramble sampling
scr_sampling <- function(random_50000_pair_1,random_50000_pair_2,cells,n_scramble_sampling,anotherset=FALSE){
  df <- data.frame(V1=c(random_50000_pair_1,random_50000_pair_2),stringsAsFactors=F)
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
                                      use_dimension_target=c(),
                                      conv_vec=conv_vec){
  print("start...")
  
  if(rowcount){
    orig.pcaspace   <- seu_base@assays$Spatial@scale.data[VariableFeatures(seu_base),]
    target.pcaspace <- seu_target@assays$Spatial@scale.data[VariableFeatures(seu_target),]
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
  
  orig.pcaspace <- orig.pcaspace[,unname(conv_vec)]
  target.pcaspace <- target.pcaspace[,names(conv_vec)]
  colnames(target.pcaspace) <- conv_vec[colnames(target.pcaspace)]
  
  #Random sampling of cells
  if(rowcount || rawsignal){
    print("start sampling...")
    sampled_df <- sampling_combination(1111,colnames(orig.pcaspace),n_cell_sampling) %>% as_tibble
    random_50000_pair_1 <- sampled_df$c1 %>% as.character
    random_50000_pair_2 <- sampled_df$c2 %>% as.character
  }else{
    sampled_df <- sampling_combination(1111,rownames(orig.pcaspace),n_cell_sampling) %>% as_tibble
    random_50000_pair_1 <- sampled_df$c1 %>% as.character
    random_50000_pair_2 <- sampled_df$c2 %>% as.character
  }  
    
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
                               colnames(orig.pcaspace),
                               n_scramble_sampling)
    pairs.target <- scr_sampling(random_50000_pair_1 %>% as.character,
                                 random_50000_pair_2 %>% as.character,
                                 colnames(target.pcaspace),
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


caldist <- function(pos,posvec1,posvec2,df_input){
  res <- tryCatch({
    v1 <- df_input[,posvec1[pos]]
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
filter_thick_col <- function(x){
  val <- col2rgb(x)
  if(sum(val)<600){return(x)}else{return(NA)}
}

combi_sample <- function(n_cell,seed,n_combi=50000){
  set.seed(seed)
  x <- 1:n_cell
  N <- 2
  M <- n_combi
  out <- list()
  maxcombi <- n_cell*(n_cell-1)/N
  while(length(out) < M) {
    out <- c(out,
             unique(replicate(M - length(out), sort(sample(x, N)), simplify = FALSE)))
    out <- unique(out)
    if(length(out)==maxcombi){
      break
    }
  }
  df <- do.call(rbind.data.frame,out)
  colnames(df) <- c("combi1","combi2")
  return(df)
}



orig.hippo <- "original_rds//slide_hippo_orig.rds"
orig.cereb1 <- "original_rds//slide_cereb1_orig.rds"
orig.cereb2 <- "original_rds//slide_cereb2_orig.rds"
orig.kidney <- "original_rds//slide_kidney_orig.rds"
orig.liver <- "original_rds//slide_liver_orig.rds"

hippo.ed0 <- "rds/hippo_ed0_slide.rds"
cereb1.ed0 <- "rds/cereb1_ed0_slide.rds"
cereb2.ed0 <- "rds/cereb2_ed0_slide.rds"
kidney.ed0 <- "rds/kidney_ed0_slide.rds"
liver.ed0 <- "rds/liver_ed0_slide.rds"

correspo.hippo <- "data/zoom_in/hippo_exp10_aligned.tsv"
correspo.cereb1 <- "data/zoom_in/cereb1_exp10_aligned.tsv"
correspo.cereb2 <- "data/zoom_in/cereb2_exp10_aligned.tsv"
correspo.kidney <- "data/zoom_in/kidney_exp10_aligned.tsv"
correspo.liver <- "data/zoom_in/liver_exp10_aligned.tsv"

L <- list(hippo_ed0=c(orig.hippo,hippo.ed0,correspo.hippo,"benchmarking/hippo_ed0.pdf"),
          cereb1_ed0=c(orig.cereb1,cereb1.ed0,correspo.cereb1,"benchmarking/cereb1_ed0.pdf"),
          cereb2_ed0=c(orig.cereb2,cereb2.ed0,correspo.cereb2,"benchmarking/cereb2_ed0.pdf"),
          kidney_ed0=c(orig.kidney,kidney.ed0,correspo.kidney,"benchmarking/kidney_ed0.pdf"),
          liver_ed0=c(orig.liver,liver.ed0,correspo.liver,"benchmarking/liver_ed0.pdf"))


# Get spot-to-spot Euclidean distances
for(l in names(L)){
  print(l)
  orig.seu <- readRDS(L[[l]][1])
  barista.seu <- readRDS(L[[l]][2])
  correspo.table <- read.table(L[[l]][3],header = T,stringsAsFactors = F)
  outfile <- L[[l]][4]
  
  dup_vec <- correspo.table$d_ind_bc %>% table
  d_ind_bc.dedup <- correspo.table$d_ind_bc[!duplicated(correspo.table$d_ind_bc)] %>% as.character()
  d_ind_bc.dedup <- paste0("export_",d_ind_bc.dedup,"-1")
  d_ind_bc.used <- intersect(d_ind_bc.dedup,colnames(barista.seu))
  correspo.table.dedup <- correspo.table[!duplicated(correspo.table$d_ind_bc),]
  rownames(correspo.table.dedup) <- paste0("export_",correspo.table.dedup$d_ind_bc,"-1")
  
  spot_used.df <- correspo.table.dedup[d_ind_bc.used,]
  conv_vec <- spot_used.df$barcode %>% as.character()
  names(conv_vec) <- rownames(spot_used.df)
  
  v.vln.to10x <- get_dist_pcaspace_scatter(orig.seu,
                                           n_cell_sampling = 50000,
                                           barista.seu,
                                           outdir = "./",
                                           outname = l,
                                           conv_vec=conv_vec)
  
  print("scramble")
  df.scramble.to10x <- get_dist_pcaspace_scatter(orig.seu,
                                                 barista.seu,
                                                 n_cell_sampling = 50000,
                                                 n_scramble_sampling=100,
                                                 scramble = T,
                                                 outdir = "./",
                                                 outname = l,
                                                 conv_vec=conv_vec)

}

#Plot
my_col <- brewer.pal(4,"Paired")[1:2]; names(my_col) <- c("keep_UMI","random_UMI")
L <- list(cereb1="cereb1_to10x_CCDist.tsv",
          cereb2="cereb2_to10x_CCDist.tsv",
          hippo="hippo_to10x_CCDist.tsv",
          kidney="kidney_to10x_CCDist.tsv",
          liver="liver_to10x_CCDist.tsv")

for(i in names(L)){
  df_now <- read.table(L[[i]],header=T)
  g <- fnc_draw_scatter(df_rand = df_now,
                        colset = my_col, 
                        outprefix = paste0("./",i),
                        both_rand_and_keep = F,
                        show_legend = F,
                        show_label=F,
                        w=4.48,
                        h=4.48,
                        carib_x=140,
                        carib_y=140)
  
}

my_col["shuffle"] <- "grey"
for(i in names(L)){
  g <- diff_rank_sina(rand_ccdist_path = L[[i]],
                      keep = F,
                      colset = my_col,
                      outdir = "./",
                      outname=paste0(i,"_to_visium"),
                      w=0.76,
                      h = 1.06)
  
}
