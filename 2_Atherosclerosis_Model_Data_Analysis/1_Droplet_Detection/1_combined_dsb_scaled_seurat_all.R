commandArgs()

#### To run DSB droplet detection you need to run:

# 1_combined_dsb_scaled_seurat_all.R --> this script
# 2_combined_hdist.R 
# 3_combined_clusterHash.R

#### This script are meant to have individual tissues where several conditions are pooled together
# in order to perform better in the droplet detection by cell type. 
## usage: nohup Rscript 1_combined_dsb_scaled_seurat_all.R path2files tissue samplestable.txt > log_droplet_detection_date.log
## example: nohup Rscript 1_combined_dsb_scaled_seurat_all.R /user/folder/ WA samplestable.txt > log_droplet_detection_date.log

### there are several options that will allow you to run this in the terminal
setwd(getwd())

# load packages used in this vignette 
AUTOINSTALL=function(){
  # check libraries, autoinstall:
  
  # bioconductor packages
  list.of.packages <- c("tidyverse" ,"Seurat", "dsb", "Matrix")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  
  if(length(new.packages)){
    BiocManager::install(new.packages)
  }
  
  list.of.packages <- c("tidyverse", "Seurat", "dsb","Matrix")
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


# STEP1. Load RNA and protein alignment (Cell Ranger) data and define metadata

## list all input folders (could adapt with args PATH and tissue, cond could be a txt file)

#........................................................................................#
PATH=args[1]

#### samples are combined by tissue
tissue = args[2]

## do not remove the bar at the end!
OUTPUT=paste0(tissue,"_","combined_results/")
dir.create(paste0(getwd(),"/",OUTPUT))

### include here your samples, if run in bash write your samples in a column in a .txt file
cond=read.table(args[3])[1,]
# cond = c("C","ED","LD","PL","IC")

### if different set of samples and tissues, make a vector with the samples to be analyzed
res_names=paste(cond,tissue, sep="-")
res_names

## Let's first read the RNA data and see if we can plot n features and n counts per dataset
## here we could used feature matrix instead of raw 

## this script will produce two lists with scaled datasets and will combine those tables at the end
## preliminary plots are done sample wise

## md.all is a list that stores the metadata for each sample
md.all=NULL

## dsb_norm_list is a list that stores the DSB normalized counts for HTO assay
dsb_norm_list=NULL

for ( i in res_names){
  # read raw data using the Seurat function "Read10X" 
  raw = Seurat::Read10X(file.path(PATH,i,"outs/raw_feature_bc_matrix/"))
  cells = Seurat::Read10X(file.path(PATH,i,"outs/filtered_feature_bc_matrix/"))
  
  # define a vector of cell-containing barcodes and remove them from unfiltered data 
  stained_cells = colnames(cells$`Gene Expression`)
  background = setdiff(colnames(raw$`Gene Expression`), stained_cells)
  
  # split the data into separate matrices per assay 
  prot = raw$'Antibody Capture'
  rna = raw$'Gene Expression'
  
  # create metadata of droplet QC stats used in standard scRNAseq processing
  rna_size = log10(Matrix::colSums(rna))
  prot_size = log10(Matrix::colSums(prot))
  ngene = Matrix::colSums(rna > 0)
  mtgene = grep(pattern = "^mt-", rownames(rna), value = TRUE)
  propmt = Matrix::colSums(rna[mtgene, ]) / Matrix::colSums(rna)
  md = as.data.frame(cbind(propmt, rna_size, ngene, prot_size))
  md$bc = rownames(md)
  md$droplet_class = ifelse(test = md$bc %in% stained_cells, yes = 'cell', no = 'background')
  
  # filter barcodes to only include those with data for both assays 
  md = md %>% dplyr::filter(rna_size > 0 & prot_size > 0 )
  md.all[[i]]=md
  
  # Quality control on cell-containing and background droplets
  #The plot below shows the number of detected genes vs the protein library size for cells vs 
  #background drops. One can also define the cells vs background drops directly by thresholding 
  #on a plot like this if using a different count aligner that does not provide filtered cell output.
  
  pdf(paste0(i,"_rawData_dsb_plots.pdf"))
  p = ggplot(md, aes(x = log10(ngene), y = prot_size )) +
    theme_bw() + 
    geom_bin2d(bins = 300) + 
    scale_fill_viridis_c(option = "C") + 
    facet_wrap(~droplet_class) 
  print(p + labs(title = paste(i,"rawData DSB"
  )))
  
  dev.off()
  
  # We next further filter cells based on thresholds calculated from quality control metrics as in 
  # any standard scRNAseq analysis, e.g. see Luecken et. al. 2019 Mol Syst Biol.
  pdf(paste0(i,"_rawData_dsb_plots_hist.pdf"))
  
  cellmd = md %>% filter(droplet_class == 'cell')
  plot_aes = list(theme_bw(), geom_point(shape = 21 , stroke = 0, size = 0.7), scale_fill_viridis_c(option = "C"))
  p1 = ggplot(cellmd, aes(x = rna_size )) + geom_histogram(bins = 50) + theme_bw() + xlab("log10 RNA library size")
  p2 = ggplot(cellmd, aes(x = propmt)) + geom_histogram(bins = 50) + theme_bw() + xlab("mitochondrial read proportion")
  p3 = ggplot(cellmd, aes(x = log10(ngene), y = rna_size, fill = propmt )) + plot_aes
  p4 = ggplot(cellmd, aes(x = ngene, y = prot_size, fill = propmt )) + plot_aes
  print(p1+p2+p3+p4)
  dev.off()
  #Note cells with the smaller library size are mostly naive CD4 T cells which are small in size and naturally have less mRNA content.
  
  # calculate statistical thresholds for droplet filtering.
  positive_cells = colnames(cells$`Gene Expression`)
  cells_mtx_rawprot = as.matrix(prot[ , positive_cells])
 
  #Sanity check: are the number of cells passing QC in line with the expected recovery from the experiment?
  
  length(positive_cells)
  
  ## background
  brna_size_min = median(md$rna_size[md$droplet_class=="background"]) - (2*mad(md$rna_size[md$droplet_class=="background"]))
  brna_size_max = median(md$rna_size[md$droplet_class=="background"]) + (2*mad(md$rna_size[md$droplet_class=="background"]))
   print(paste0(" rna size lower cutoff in background: ", brna_size_min))
  print(paste0(" rna size upper cutoff in background: ", brna_size_max))
  
  pdf(paste0(res_names,"_raw_plots_cells_param.pdf"))
  ggplot(md, aes(x = log10(ngene), y = rna_size, fill = propmt, group = propmt )) +
    theme_bw() +
    geom_bin2d(bins = 300) +
    scale_fill_viridis_c(option = "C") +
    facet_wrap(~droplet_class) +
     geom_hline(yintercept = brna_size_max, color = 'red', size = 1) +
    geom_hline(yintercept = brna_size_min, color = 'red', size = 1)
  dev.off()
  
  isNeg=md$rna_size < brna_size_max &md$rna_size>brna_size_min &md$propmt < 0.25
  
  # define a vector of background / empty droplet barcodes based on protein library size and mRNA content  
  background_drops = md[isNeg, ]$bc
  negative_mtx_rawprot = prot[ , background_drops] %>%  as.matrix()
  
  #normalize protein data for the cell containing droplets with the dsb method. 
  
  dsb_norm_prot = DSBNormalizeProtein(
    cell_protein_matrix = cells_mtx_rawprot, 
    empty_drop_matrix = negative_mtx_rawprot, 
    denoise.counts = FALSE 
  
  )
  
  dsb_norm_list[[i]]=as.data.frame(dsb_norm_prot)
  
  ## diagnostic plots 
  pdf(paste0(i, ".histScaled.pdf"))
  ## check some hashtags
  hist(cells_mtx_rawprot[1,], breaks=100, main="hash1 raw signal across cells")
  hist(dsb_norm_prot[1,], breaks=100, main="hash1 scaled signal across cells")
  hist(cells_mtx_rawprot[2,], breaks=100, main="hash2 raw signal across cells")
  hist(dsb_norm_prot[2,], breaks=100, main="hash2 scaled signal across cells")
  hist(cells_mtx_rawprot[3,], breaks=100, main="hash3 raw signal across cells")
  hist(dsb_norm_prot[3,], breaks=100, main="hash3 scaled signal across cells")

  dev.off()
  
   
}

names(dsb_norm_list) = res_names

# STEP2. Refine matrices to be compatible with each other

#........................................................................................#
makeObj = function(x,par,N){

#' @x = metadata dataframe
#' @par = "rows" or "columns"
#' @N = outputfilename.rds
  
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

##### A) METADATA ######
### create a dataframe with the metadata for each sample
md.all.df=makeObj(md.all,"rows","md.all.df.rds")

##### B) scaled DSB  ######
### create a dataframe with the scaled DSB data for each sample
out=makeObj(dsb_norm_list,"columns",".scaled_prot_DSB.rds")

### load Seurat object gene expression
d10x.data.rna <- sapply(res_names, function(i){
  d10x <- Read10X(file.path(PATH,i,"outs/filtered_feature_bc_matrix/"))
  d10x = d10x$"Gene Expression"
  colnames(d10x) <- paste(colnames(d10x),i,sep=".")
  d10x
})

## rename accordingly each list
names(d10x.data.rna)=res_names

### bind the rna matrices
rna <- do.call("cbind", d10x.data.rna)

### check naming consistency
colnames(rna)[1:5]
colnames(out)[1:5]
rownames(md.all.df)[1:5]

######### -------------- START WITH THE SEURAT PACKAGE -------------------###########
dim(md.all.df[colnames(rna),])
dim(rna)
table(is.na(colnames(rna)))
table(is.na(md.all.df[colnames(rna),]))


#### fix rna if needed, rna = rna[,which(colnames(rna)%in%rownames(md.all.df))]
s = CreateSeuratObject(counts = rna, meta.data = md.all.df[colnames(rna),], assay = "RNA")

## scaled hashtag data
s[["HTO"]] <- CreateAssayObject(data = as.matrix(out[,which(colnames(out)%in%colnames(s))]))
saveRDS(s, file=paste0(OUTPUT,tissue,"_seuratObj_step1.rds"))


## Session info
print(paste("Last compiled on", Sys.time()))
sessionInfo()
