setwd("/research/groups/sysgen/PROJECTS/terva/HTO_classification/hist_distribution/")
### set tissue
#tissue=args[1]
tissue="MO"

## do not remove the bar at the end!
OUTPUT=paste0(tissue,"_","combined_results/")



### check for each tissue if the droplet identification is better or similar to Seurat
getwd()
s=readRDS(paste0(OUTPUT,tissue,".sObj.with.dsb.classif_step3.rds"))
table(s$HTO_classification.global)
#### print signals with 

##a) DSB classsification

pdf(paste0(OUTPUT,tissue,".all.tsne_hash_signals_bycell_dsb.pdf"))
print(FeaturePlot(s, reduction="tsne", features = paste0("hto_",rownames(GetAssayData(s,assay="HTO"))), combine=F, order=T, cols=c("midnightblue","seagreen","chartreuse1", "yellow")))
dev.off()

##b) SEURAT

## create a vector with seurat HTO files per tissue

cond= c("C","ED","LD","PL","IC")
res_names=paste(cond,tissue, sep="-")


path.hto="/research/groups/sysgen/PROJECTS/terva/cellranger/"
path.hto="/research/groups/sysgen/PROJECTS/terva/HTO_classification/HTO_classif_seurat/"
seurat_hto_samples = paste0(path.hto,res_names,"_primirs_HTO_info.txt")
seurat_hto_samples = paste0(path.hto,res_names,"_HTO_info.txt")

list.DFs <- lapply(seurat_hto_samples,read.table)
names(list.DFs)=res_names

makeObj = function(x,par,N){
  
  if(par=="rows"){
    ## combined data frames
    all.df = do.call(rbind.data.frame, x)
    namescells=do.call(rbind,strsplit(rownames(all.df),"\\."))
    rownames(all.df)=paste(namescells[,2],namescells[,1],sep=".")
  }else {
    all.df=do.call(cbind.data.frame,x)
    namescells=do.call(rbind,strsplit(colnames(all.df),"\\."))
    colnames(all.df)=paste(namescells[,2],namescells[,1],sep=".")
  }
  ## rewrite rownames to match with seurat object
  
  
  ## save metatata for later
  saveRDS(all.df,file=paste0(OUTPUT,tissue,N))
  
  return(all.df)
}


seurat_hto = makeObj(list.DFs, "rows",paste0("seurat_hto_",tissue,".rds"))
table(seurat_hto$HTO_classification.global)
rownames(seurat_hto)[1:5]

s <- AddMetaData(object = s, metadata = seurat_hto)


Idents(s) <- s@meta.data$HTO_classification.global
table(Idents(s))
pdf(paste0(tissue,"_all.tsne_hash_clustering_classified_SEURAT_celltype.pdf"),width = 10, height = 7)
p = DimPlot(s, reduction = "tsne")
print(p + labs(title = "tsne_hash_clustering_classified_SEURAT"
))
dev.off()


#DSB by cluster HTO
s = saves
newclss = do.call(rbind.data.frame, listseurat)
names = do.call(rbind,strsplit(rownames(newclss),"\\."))
head(names)
rownames(newclss)=paste(names[,2],names[,3],sep=".")
head(newclss)
s <- AddMetaData(object = s, metadata = newclss)
Idents(s) <- s@meta.data$HTO_classification

pdf("wa_tsne_hash_clustering_classified_DSB_byclus_celltype.pdf",width = 10, height = 7)
p = DimPlot(s)
print(p + labs(title = "WA all tsne_hash_clustering_classified_DSB by clus v3"
))
dev.off()

save.image("dsb_by_cluster_c-wa_ld.rda")

#### plot things with ditto seq

table(Idents(s))
#BiocManager::install("dittoSeq")
library(Seurat)
s =FindVariableFeatures(s)
library("dittoSeq")
dittoBarPlot(s, "celltype", group.by = "ConditionName") 
## subset for example clusters 

FeaturePlot(s, features = "Cd3e")

Idents(s) = s$leiden_r08

## t-cell cluster 11,10 and 7
RidgePlot(s, features = c("Cd3e","Cd4"), ncol = 3)

## NK cluster
RidgePlot(s, features = c("Cd27"), ncol = 2)

## myeloid markers
RidgePlot(s,features=c("Ccr2","Cd14","Cd36"),ncol=3)

### for C-MO select lympyoid clusters 11,10 and 7 and myeloid 9
table(Idents(s))
sL =s[, WhichCells(s, idents = c("11","7","10"))]
sM = s[,WhichCells(s,idents = "9")]


table(sL$HTO_classification.global)
library(dittoSeq)
dev.off()
dittoScatterPlot(
  object = s,
  x.var = "n_counts", y.var = "n_genes",
  color.var = "HTO_classification.global")

### write down the cell names single and cell names for doublet, plot scaled matrix

s$HTO_classification.global=="Doublet"
doubletL = sL[["HTO"]]@data[s$HTO_classification.global=="Doublet"]
singletL = sL[["HTO"]]@data[s$HTO_classification.global=="Singlet"]

plot(doubletL,singletL)

dsb = readRDS("classification_global_cellranger.rds")
s = readRDS(input)
s = AddMetaData(s, metadata=dsb)
s= AddMetaData(s,metadata =obs )
Idents(s) = s
VlnPlot(s, features = "n_counts") 

VlnPlot(s,features="nFeature_RNA")
saveRDS(s,file="seuratObj_dsb_byclus.rds")
saveC = s

### repeat with seurat HTO
pdf(file="plots_tocheck_Dsb_wa_5.pdf")
#s= AddMetaData(s, metadata= classauto)
s = AddMetaData(s, metadata=class5)
Idents(s) = s$HTO_classification.global
VlnPlot(s, features="nCount_RNA", pt.size=0)
VlnPlot(s, features="nFeature_RNA", pt.size=0)
FeatureScatter(s, "n_counts","n_genes")
dev.off()

saveRDS(s, file="seurat_object_celltype_FigxS1.rds")

dittoBarPlot(s, "HTO_classification.global", group.by = "celltype.pred") 

