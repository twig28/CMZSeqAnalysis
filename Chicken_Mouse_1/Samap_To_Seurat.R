#File to convert samap outputs into different formats in R for further analysis
install = TRUE
annDataPickleFileName = "Samap_adata.pkl"
annDataRawPickleFileName = "samap_raw_adata.pkl"
samapObjectPickleFilename = "samap_object.pkl"
h5adFileName = "samap_adata"

if(install){
  remove.packages('reticulate')
  install.packages('reticulate')
  install.packages("Seurat")
  install.packages('devtools')
  devtools::install_github(repo = "mojaveazure/loomR", ref = "develop")
  devtools::install_github("cellgeni/sceasy")
  library("Seurat")
  require("reticulate")
  reticulate::py_install("samap")
  reticulate::py_install("sam-algorithm")
  reticulate::py_install("scanpy==1.9.3")
  reticulate::py_install("colorlover")
  reticulate::py_install("plotly")
  reticulate::py_install("ipywidgets")
  reticulate::py_install("ipyevents")
  reticulate::py_install("numpy==1.23.5")
}

#This Will Use a Lot of memory as it loads the entirety of the samap implementation/object
samap_pickle_data <- reticulate::py_load_object(paste0("./", samapObjectPickleFilename))

anndata_pickle_data <- reticulate::py_load_object(paste0("./", annDataPickleFileName))

anndata_raw_pickle_data <- reticulate::py_load_object(paste0("./", annDataRawPickleFileName))

#Load samap anndata into seurat using h5ad
sceasy::convertFormat(paste0(h5adFileName, ".h5ad"), from="anndata", to="seurat",
                      outFile=paste0(h5adFileName, ".rds"))
seuratObjecth5ad = readRDS(paste0(h5adFileName, ".rds"))

#Print Gene Names
gene_List <- rownames(seuratObjecth5ad@assays$RNA@counts)
print(gene_List)




#---For Sparse Data---






#Convert Sparse R Object into h5ad for Samap Sparse Run
h5adFileName <- "samap_adata_sparse"
sparse_RDS <- readRDS("chMG_NMsubset.rds")

#add cell type label and subset out FI and NFI
print(sparse_RDS@active.ident)
sparse_RDS <- AddMetaData(
  object = sparse_RDS,
  metadata = Idents(sparse_RDS),
  col.name = "active_ident"
)
print(unique(sparse_RDS@meta.data$treatment))

sparse_RDS <- subset(
  sparse_RDS,
  subset = treatment != "FI" & treatment != "NFI"
)

table(sparse_RDS@meta.data$treatment)

#Save for Samap Run
sceasy::convertFormat(sparse_RDS, from="seurat", to="anndata",
                      outFile=paste0(h5adFileName, ".h5ad"))

#---SAMAP RUN----

#Load Samap Run
h5adFileNameNew = "samap_results_h5ad_sparse"
sceasy::convertFormat(paste0(h5adFileNameNew, ".h5ad"), from="anndata", to="seurat",
                      outFile=paste0(h5adFileNameNew, ".rds"))
seuratObjectSparse = readRDS(paste0(h5adFileNameNew, ".rds"))
