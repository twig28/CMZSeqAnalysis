install = TRUE
if(install){
  #remove.packages('reticulate')
  #install.packages('reticulate')
  install.packages("Seurat")
  install.packages('devtools')
  devtools::install_github(repo = "mojaveazure/loomR", ref = "develop")
  devtools::install_github("cellgeni/sceasy")
  library("Seurat")
  require("reticulate")
  #reticulate::py_install("samap")
  #reticulate::py_install("sam-algorithm")
  #reticulate::py_install("scanpy==1.9.3")
  #reticulate::py_install("colorlover")
  #reticulate::py_install("plotly")
  #reticulate::py_install("ipywidgets")
  #reticulate::py_install("ipyevents")
  #reticulate::py_install("numpy==1.23.5")
}

#Convert Sparse R Object into h5ad for Samap Sparse Run
h5adName <- "samap_adata_sparse_mouse"
seuratObjectName <- "mouse.rds"


sparse_RDS <- readRDS(seuratObjectName)
sparse_RDS = UpdateSeuratObject(object = sparse_RDS)

sparse_RDS <- AddMetaData(
  object = sparse_RDS,
  metadata = Idents(sparse_RDS),
  col.name = "active_ident"
)

DimPlot(sparse_RDS)

#sparse_RDS <- subset(sparse_RDS, subset = treatment != "FI" & treatment != "NFI")

table(sparse_RDS$active_ident)
table(sparse_RDS$orig.ident)

#Save for Samap Run
sceasy::convertFormat(sparse_RDS, from="seurat", to="anndata",
                      outFile=paste0(h5adName, ".h5ad"))

#RUN SAMAP HERE

#Load Seurat Object After Samap Run
h5adFileNameNew = "samap_results_h5ad_sparse"
sceasy::convertFormat(paste0(h5adFileNameNew, ".h5ad"), from="anndata", to="seurat",
                      outFile=paste0(h5adFileNameNew, ".rds"))
seuratObjectSparse = readRDS(paste0(h5adFileNameNew, ".rds"))

gene_List <- rownames(seuratObjecth5ad@assays$RNA@counts)
print(gene_List)