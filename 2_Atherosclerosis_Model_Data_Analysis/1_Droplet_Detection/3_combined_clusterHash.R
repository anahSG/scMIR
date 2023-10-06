commandArgs()
# load packages used in this vignette 
AUTOINSTALL=function(){
  # check libraries, autoinstall:
  
  # bioconductor packages
  list.of.packages <- c("tidyverse" ,"Seurat", "dsb", "Matrix","cluster","fitdistrplus")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  
  if(length(new.packages)){
    BiocManager::install(new.packages)
  }
  
  list.of.packages <- c("tidyverse", "Seurat", "dsb","Matrix","cluster","fitdistrplus")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages, dependencies = T, repos='http://cran.us.r-project.org')
  
}

# Autoinstall packages
cat("Checking R packages, installing if missing...", sep="\n\n")
AUTOINSTALL()
cat("Checking R packages, installing if missing... Done", sep="\n\n")


suppressMessages(library(dsb))
suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
suppressMessages(library(Matrix))
suppressMessages(library(cluster))
suppressMessages(library(fitdistrplus))


MaxN <- function(x, N = 2){
  len <- length(x)
  if (N > len) {
    warning('N greater than length(x).  Setting N=length(x)')
    N <- length(x)
  }
  sort(x, partial = len - N + 1)[len - N + 1]
}


### set wd

setwd("/research/groups/sysgen/PROJECTS/terva/HTO_classification/hist_distribution/")
#setwd(getwd)


### set tissue
#tissue=args[1]
tissue="MO"

## do not remove the bar at the end!
OUTPUT=paste0(tissue,"_","combined_results/")

## read seurat objects and  OBS (metadata from scanpy)

#obs.path=args[2]
#obs.path="/research/groups/bioinformaticians/internship/celtundi/outs/"

# obs=read.delim(paste0(obs.path,tissue,"/",tissue,"_metadata_human_01.csv"),
#               sep=",", header=T, row.names=1)
#obs=read.delim(paste0(obs.path,tissue,"/",tissue,"_metadata.csv"),
#               sep=",", header=T, row.names=1)


##AO
 obs.path="/research/groups/bioinformaticians/internship/celtundi/outs/integrated_datasets/"

obs=read.delim(paste0(obs.path,tissue,"_metadata.csv"),
               sep=",", header=T, row.names=1)




s=readRDS(paste0(OUTPUT,tissue,"_seurat_scaled_step2.rds"))
table(obs$predicted.celltype.l1)
table(obs$celltype.pred)
rownames(obs)[1:5]
rownames(obs)[13385:13392]

cleanNames=function(df){
  df$names=gsub("-1.*","-1",rownames(df)) 
  rownames(df)=paste0(df$names,".",df$ConditionName,"-", df$Tissue)
  return(df)
}

obs = cleanNames(obs)
rownames(obs)[1:5]


# add cell type predictions
s = AddMetaData(s, metadata=obs)
#s=saves
 # s$celltype.pred=s$predicted.celltype.l1
 # obs$celltype.pred=obs$predicted.celltype.l1
 # 
table(s$celltype.pred)
table(s$ConditionName)
#s$celltype.pred[(s$leiden_r08%in%"14"|s$leiden_r08%in%"13")]="13_14group"

table(s$celltype.pred)


assay = "HTO"
data <- GetAssayData(object = s, assay = assay)

positive.quantile = 0.99
init = nrow(x = data) + 1 + choose(nrow(x = data),2)

nstarts = 100
kfunc = "clara"
nsamples = 100
seed = 42
set.seed(seed)

## originally HTOdemux takes the counts but we try to adapt this to DSB output
ncenters=init


### clean list to start the loop

listseurat=NULL
dev.off()
data=NULL

pdf(paste0(OUTPUT,tissue,"_all.tse_hash_clustering_classification_sep_celltype.pdf"))

s$celltype.pred[is.na(s$celltype.pred)] = "Neg"
table(s$celltype.pred)

### add a piece of code that in case the are cells under 30 events to call them neg

library("dplyr")

### function to substitute those cell types that have less than 100 cells

too_few= function(x){
  test=x %>%
    count(celltype.pred)
 return( test$celltype.pred[test$n<50])
}

too_few_cell=too_few(obs)
s$celltype.pred[s$celltype.pred%in%too_few_cell]="Neg"


table(s$celltype.pred)
saves=s
table(saves$celltype.pred)
cell.pred=unique(s$celltype.pred)

listseurat=NULL
for (j in cell.pred) {

 s=saves

 Idents(s) <- s$celltype.pred
# Idents(s) = s$leiden_r08
 table(Idents(s))
 #j="mesenchymal stem cell of adipose"
 #rm(j)
  g <- WhichCells(s, idents = j)
 # g = WhichCells(s,idents = cell.pred[1])
  s <- s[,g]
  data <- GetAssayData(object = s, assay = assay)
  
  ### k -means clusters
  init.clusters.k <- kmeans(x = t(data),centers = ncenters,nstart = nstarts)
  
  Idents(object = s, cells = names(x = init.clusters.k$cluster)) <-init.clusters.k$cluster 
#pdf(file="cluster14_WA_test.pdf")
  p= DimPlot(s, reduction="tsne")
  print(p + labs(title = "tsne_hash_clustering_k-means",
                 subtitle = paste("cluster",j)
  ))

  #### clara clusters
  
  init.clusters.c <- clara(x = t(data),k = ncenters,samples = nsamples)
  Idents(object = s, cells = names(x = init.clusters.c$cluster)) <-init.clusters.c$cluster 

  p = DimPlot(s,reduction="tsne")
  
  print(p + labs(title = "tsne_hash_clustering_clara",
                 subtitle = paste("cluster",j)))
  

##average hto signals per cluster, default AverageExpression will assume data is log normalized, averages in non-log scale and returns that
#average.expression <- AverageExpression(object = s,assays = assay,verbose = FALSE)[[assay]]



 # Idents(object = s, cells = names(x = init.clusters.k$cluster)) <-init.clusters.k$cluster 

average.expression<-matrix(data=0,nrow=nrow(data), ncol=length(unique(init.clusters.c$clustering)))
rownames(average.expression)=rownames(data)

for(i in 1:ncol(average.expression)){
	average.expression[,i]=apply(data[,init.clusters.c$clustering==unique(init.clusters.c$clustering)[i]],1,mean)
}

print(average.expression)
rm(discrete)
#create a matrix to store classification result
discrete <- matrix(data=0, nrow=nrow(data), ncol=ncol(data))
rownames(discrete)=rownames(data)
colnames(discrete)=colnames(data)
## for each HTO, we will use the cluster with min value for fitting
## originally this fits neg binom (nbinom) to counts (not normalized?)
## but here we use the DSB output that should give us values resembling normal distribution
for (iter in rownames(x = data)) {
    values <- data[iter, ]

	print(which.min(x = average.expression[iter, ]))
    values.use <- values[WhichCells(object = s,idents = levels(x=Idents(object = s))[[which.min(x = average.expression[iter, ])]])]
	print(summary(values.use))
	print(quantile(values.use,positive.quantile))
	## if for some reason we would have a cluster of non-hash labeled cells
	## DSB still would give some low value per cell, so we don't need to skip fitting
    fit <- suppressWarnings(expr = fitdist(data = values.use, distr = "norm"))
    cutoff <- as.numeric(x = quantile(x = fit, probs = positive.quantile)$quantiles[1])
    if (cutoff<2.5){
     cutoff=2.5
     print("Cutoff under 2.5, will be replace with 2.5")
    }else{

      }

    
    discrete[iter, values > cutoff] <- 1
    ## same as : discrete[iter, names(x = which(x = values > cutoff))] <- 1
    print(cutoff)
   }

#### plots ####


# now assign cells to HTO based on discretized values
npositive <- colSums(x = discrete)
print("Number of positive hash per cell: ")
summary(npositive)

classification.global <- npositive
classification.global[npositive == 0] <- "Negative"
classification.global[npositive == 1] <- "Singlet"
classification.global[npositive > 1] <- "Doublet"

table(classification.global)
#classification.global
# Doublet Negative  Singlet
#     572      561     2127

donor.id = rownames(x = data)
hash.max <- apply(X = data, MARGIN = 2, FUN = max)
hash.maxID <- apply(X = data, MARGIN = 2, FUN = which.max)
hash.second <- apply(X = data, MARGIN = 2, FUN = MaxN, N = 2)
hash.maxID <- as.character(x = donor.id[sapply(X = 1:ncol(x = data),
    FUN = function(x) {return(which(x = data[, x] == hash.max[x])[1])}
  )])
hash.secondID <- as.character( x = donor.id[sapply(X = 1:ncol(x = data),
    FUN = function(x) {return(which(x = data[, x] == hash.second[x])[1])}
  )])
  
hash.margin <- hash.max - hash.second

print("Hash separation summary:")
summary(hash.margin)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.00182  2.24355  5.11126  5.50259  8.09166 60.99850

doublet_id <- sapply(
    X = 1:length(x = hash.maxID),
    FUN = function(x) {
      return(paste(sort(x = c(hash.maxID[x], hash.secondID[x])), collapse = "_"))
    }
)


classification <- classification.global
classification[classification.global == "Negative"] <- "Negative"
classification[classification.global == "Singlet"] <- hash.maxID[which(x = classification.global == "Singlet")]
classification[classification.global == "Doublet"] <- doublet_id[which(x = classification.global == "Doublet")]

## prepare metadata
classification.metadata <- data.frame(
    hash.maxID,
    hash.secondID,
    hash.margin,
    classification,
    classification.global
  )

colnames(x = classification.metadata) <- paste(assay,
    c('maxID', 'secondID', 'margin', 'classification', 'classification.global'),sep = '_')
 
## add to Seurat object
#class=classification.metadata

s <- AddMetaData(object = s, metadata = classification.metadata)
Idents(s) <- paste0(assay, '_classification')
table(Idents(s))

FeatureScatter(s, 'anti-mouse-Hashtag-4', 'anti-mouse-Hashtag-5')
FeatureScatter(s, 'anti-mouse-Hashtag-4', 'anti-mouse-Hashtag-6')
FeatureScatter(s, 'anti-mouse-Hashtag-5', 'anti-mouse-Hashtag-6')

doublets <- rownames(x = s[[]])[which(s[[paste0(assay, "_classification.global")]] == "Doublet")]
singlets <- rownames(x = s[[]])[which(s[[paste0(assay, "_classification.global")]] == "Singlet")]

#length(doublets)>0

if (length(doublets)>0)  {
  Idents(object = s, cells = doublets) <- 'Doublet'
} else {
  
}

p = DimPlot(s,reduction="tsne")

print(p + labs(title = "tsne_hash_clustering_classified",
                   subtitle = paste("cluster",j)))
          

s$hash.ID <- Idents(object = s)


## viz also what cells got paired in doublets

Idents(s) <- paste0(assay, '_classification')
Idents(object = s, cells = singlets) <- 'Singlet'


p = DimPlot(s, reduction="tsne")
print(p + labs(title = "tsne_hash_clustering_classified.vizDoubletPairs",
                   subtitle = paste("cluster",j)))

    
clasif = cbind(HTO_classification=s$HTO_classification,HTO_classification.global=s$HTO_classification.global)    

FeaturePlot(s, features = paste0("hto_",rownames(GetAssayData(s,assay="HTO"))), 
            combine=F, order=T, cols=c("midnightblue","seagreen","chartreuse1", "yellow")
           # , max.cutoff = 9
)

FeaturePlot(s,features=paste0("hto_",rownames(GetAssayData(s, assay="HTO"))),
            combine=F, order=T, cols=c("midnightblue","seagreen","chartreuse1","yellow"),
            max.cutoff = 9)

#Idents(s) = s$leiden_r08
#DimPlot(s, cols=DiscretePalette(n=19,palette="alphabet"))

listseurat[[j]]  = as.data.frame(clasif)


}
dev.off()

names(listseurat)

s=saves
#saveRDS(s, file=paste0(OUTPUT,tissue,"_by_celltype_demux.rds"))
saveRDS(listseurat,paste0(OUTPUT,"list_classification_",tissue,"_by_cell.rds"))

### build the table and print it

buildTable =function(x){
  newclss=do.call(rbind.data.frame,listseurat)
  names = do.call(rbind,strsplit(rownames(newclss),"\\."))
  rownames(newclss)=paste(names[,2],names[,3],sep=".")
  #newclss$names=names
  saveRDS(newclss,file=paste0(OUTPUT,tissue,".hto.classific.dsb.rds"))
  return(newclss  )
}

newclss=buildTable(newclss)

### DSB ###
s <- AddMetaData(object = s, metadata = newclss)
#s <- AddMetaData(object = s, metadata = seurat_hto)


#s <- AddMetaData(object = s, metadata = newclss)
Idents(s) <- s@meta.data$HTO_classification
doublets <- rownames(x = s[[]])[which(s[[paste0(assay, "_classification.global")]] == "Doublet")]
if (length(doublets)>0)  {
  Idents(object = s, cells = doublets) <- 'Doublet'
} else {
  
}
table(Idents(s))




saveRDS(s, file=paste0(OUTPUT,tissue,".sObj.with.dsb.classif_step3.rds"))
table(Idents(s))

### save table in the same format as obs

obs[,c("HTO_classification","HTO_classification.global")] = newclss[match(rownames(obs),rownames(newclss)),c("HTO_classification","HTO_classification.global")]
obs[is.na(obs)]= "Negative"

out.path="/research/groups/bioinformaticians/internship/celtundi/outs/hto_classifications/"
write.table(obs,file=paste0(out.path,tissue,".dsb_HTO_all_celltype_integrated.txt"), sep="\t", quote = F)

### print classification
pdf(paste0(OUTPUT,tissue,"_tsne_hash_clustering_classified_DSB_byclus_celltype.pdf"),width = 10, height = 7)
p = DimPlot(s,reduction="tsne")
print(p + labs(title = paste(tissue,"integrated tsne_hash_clustering_classified_DSB by celltype")
))
dev.off()


### print dimplot 
