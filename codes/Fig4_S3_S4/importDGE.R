library(Seurat)
library(Matrix)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

dge_path    <- args[1]
input_category <- args[2]
outdir  <-  args[3]
outname <- args[4]

source("/home/ha5656/work/barista/paper_data/R_fig_script/src/sc_preprocess.R")

#This script export the maximum size of DGE
#Only for SPLiT dataset, mtx, meta and genes data should be provided in comma-separated manner

if (! input_category %in% c("10x","drop","split")){
        stop("input_category should be chosen from 10x, drop or split.")
    }

outfullpath <- paste0(outdir,"/",outname,".",input_category,".rds")

if(input_category=="10x"){
        dge_10x <- Read10X(dge_path)
        dge_10x <- dge_10x[,colSums(dge_10x)>0]
        saveRDS(dge_10x,outfullpath)
    }else if(input_category=="drop"){
        dge_drop <- read.table(dge_path,header=T)
        dge_drop <- column_to_rownames(dge_drop,var="GENE")
        dge.mat_drop <- Matrix(as.matrix(dge_drop))
        dge_drop <- dge_drop[,colSums(dge_drop)>0]
        saveRDS(dge_drop,outfullpath)
    }else{
        mtx <- strsplit(dge_path,",")[[1]][1]
        meta <-strsplit(dge_path,",")[[1]][2]
        genes<-strsplit(dge_path,",")[[1]][3]

        dge_split <- ReadSPLiT(mtx,meta,genes)
        dge_split <- dge_split[,colSums(dge_split)>0]
        saveRDS(dge_split,outfullpath)
    }

