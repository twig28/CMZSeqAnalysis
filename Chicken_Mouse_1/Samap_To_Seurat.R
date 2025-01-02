#File to convert samap outputs into different formats in R for further analysis
install = TRUE
annDataPickleFileName = "samap_adata.pkl"
samapObjectPickleFilename = "samap_object.pkl"
h5adFileName = "Samap_adata"

if(install){
  remove.packages('reticulate')
  install.packages('reticulate')
  install.packages("Seurat")
  install.packages('devtools')
  devtools::install_github(repo = "mojaveazure/loomR", ref = "develop")
  devtools::install_github("cellgeni/sceasy")
  library("Seurat")
  library(loomR)
  require("reticulate")
  reticulate::py_install("samap")
  reticulate::py_install("sam-algorithm")
  reticulate::py_install("scanpy==1.9.3")
  reticulate::py_install("colorlover")
  reticulate::py_install("plotly")
  reticulate::py_install("ipywidgets")
  reticulate::py_install("ipyevents")
  reticulate::py_install("numpy==1.23.5")
  reticulate::py_install("loompy")
}

#This Will Use a Lot of memory as it loads the entirety of the samap implementation/object
samap_pickle_data <- reticulate::py_load_object(paste0("./", samapObjectPickleFilename))

anndata_pickle_data <- reticulate::py_load_object(paste0("./", annDataPickleFileName))

#Load into seurat using h5ad
sceasy::convertFormat(paste0(h5adFileName, ".h5ad"), from="anndata", to="seurat",
                      outFile=paste0(h5adFileName, ".rds"))
seuratObjecth5ad = readRDS(paste0(h5adFileName, ".rds"))

gene_List <- rownames(seuratObjecth5ad@assays$RNA@counts)
print(gene_List)

print(Idents(seuratObjecth5ad))

#Load into seurat using h5ad
#sceasy::convertFormat("chMG.h5ad", from="anndata", to="seurat",
#                      outFile=paste0("testchMG", ".rds"))
#seuratObjecth5adTEST = readRDS(paste0(h5adFileName, ".rds"))
