setwd("path")
# load packages used in this vignette 


suppressMessages(library(dsb))
suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
suppressMessages(library(SeuratDisk))


# STEP1. Load RNA and protein alignment (Cell Ranger) data and define metadata

## list all input folders (c ould adapt with args PATH and tissue, cond could be a txt file)
#PATH=args[1]


## do not remove the bar at the end!
tissue = args[2]

OUTPUT=paste0(tissue,"_","combined_results/")



### include here your samples, for testing select first two

cond = c("C","ED","LD","PL","IC")

res_names=paste(cond,tissue, sep="-")
res_names

#md.all=readRDS(paste0(OUTPUT,tissue,"combined_results/WA.md.all.rds"))
s=readRDS(paste0(OUTPUT,tissue,"_seuratObj_step1.rds"))

### check if cells with hastags x cells with hastags > 2.147e+9, choose other script
if(ncol(s)^2> 2.147e+9){
  print("The dataset is too large, use instead 2_combined_hdist_long_vectors.R script")
}else{
  print("The dataset size is suitable for this script")
}


out=readRDS(paste0(OUTPUT,tissue,".scaled_prot_DSB.rds"))
out=out[,which(colnames(out)%in%colnames(s))]

head(out)
colnames(out)[1:5]

### add condition to metadata
condition = do.call(rbind, strsplit(colnames(s),"-1."))
rownames(condition)=colnames(s)
head(condition)
s=AddMetaData(s, condition[,2], col.name="conditionName")

## calculate max hash signal
hashMax=apply(out[1:3,],2,max)

s=AddMetaData(s, hashMax, col.name="hashMax")

dim(s@assays$HTO)


## scaled val at least 2 standard dev from neg distr
hashDet=out[1:3,]>2
isDet=apply(hashDet,2,any)
nDet=apply(hashDet,2,sum)
print("cell detected positive for Nhash: ")
print(summary(nDet))
print(length(which(isDet)))
print("from :")
print(ncol(hashDet))


s=AddMetaData(s, nDet, col.name="nHashPos")


### make the table ready
out=out[,which(colnames(out)%in%colnames(s))]
s=s[,colnames(out)]
dim(out)
s=s[,isDet]

options(future.globals.maxSize= 14194304000)
# calculate eucledian distances
h_dist = dist(t(s@assays$HTO@data))

h_dist = as.matrix(h_dist)

s <- RunTSNE(s, distance.matrix = h_dist, perplexity =50)

## hash signal separation

pdf(paste0(OUTPUT,tissue,"_all_tsne_hash_plots_celltype.pdf"))
print(DimPlot(s,reduction="tsne"))
#dev.off()

# pdf(paste0(res_name,"_",sub,".tsne_hash_signals.pdf"))
# FeaturePlot(s, features = paste0("adt_",rownames(GetAssayData(s,assay="ADT"))), combine=F, order=T, cols=c("midnightblue","seagreen","chartreuse1", "yellow"))
print(FeaturePlot(s, reduction="tsne",features = paste0("hto_",rownames(GetAssayData(s,assay="HTO"))), combine=F, order=T, cols=c("midnightblue","seagreen","chartreuse1", "yellow")))

# dev.off()

# pdf(paste0(res_name,".tsne_hash_QC_clus10.pdf"))
print(FeaturePlot(s, features = c("ngene", "rna_size","propmt","prot_size","nHashPos")))
dev.off()

getwd()


#-------------#


shash=s

# s <- readRDS("CA_WA_tervaseurat_hashClu_clus5.rds")
## can continue with this object to de-multiplex data by hashtag signals
#saveRDS(shash, file="seurat_hashClu_WA_C-LD.rds")
#write.table(t(out), file="WA.all.signal_hash.txt", sep="\t", quote=F)

#saveRDS(s,file="seurat_hashClu_WA_C_LD.rds")
## -------------
## for specific ab can process all cells, resume from the object s1 stored above
#s=s1
##  signal separation
pdf(paste0(OUTPUT,tissue,"all.tsne_scaledOnly_terva_celltype.pdf"))
print(DimPlot(s, reduction = "tsne"))
dev.off()


saveRDS(s, file=paste0(OUTPUT,tissue,"_seurat_scaled_step2.rds"))



