#Here we process ICA BM dataset with Seurat, compute diffusion maps, generate 3D FLE with scSVAtools and save the results as hdf5 (loom) file
#Use with docker (Dockerfile in the same directory) 
#Run the script with scSVA on the Google Cloud Platform
require(Matrix)
require(Seurat)
require(rhdf5)
require(data.table)
require(future.apply)

print("Reading ICA (bone marrow) dataset")
data        <-h5read("ica_bone_marrow_h5.h5","GRCh38/data")
indices     <-h5read("ica_bone_marrow_h5.h5","GRCh38/indices")
indptr      <-h5read("ica_bone_marrow_h5.h5","GRCh38/indptr")
shape       <-h5read("ica_bone_marrow_h5.h5","GRCh38/shape")
gene_names  <-h5read("ica_bone_marrow_h5.h5","GRCh38/gene_names")
barcodes    <-h5read("ica_bone_marrow_h5.h5","GRCh38/barcodes")

print("Converting to sparse matrix")
sparse_exp_mat <- sparseMatrix(i           = indices + 1, 
                               p           = indptr,
                               x           = as.numeric(data),
                               dims        = shape,
                               giveCsparse = F)
colnames(sparse_exp_mat)<-barcodes
rownames(sparse_exp_mat)<-make.unique(gene_names)

#Run Seurat
print("Creating Seurat object")
ica <- CreateSeuratObject(raw.data  = sparse_exp_mat, 
                          min.cells = 50)
print("Filtering cells")
ica <- FilterCells(ica, 
                   subset.names    = "nGene", 
                   low.thresholds  = 500, 
                   high.thresholds = Inf)
rm(sparse_exp_mat)

print("Normalizing data")
ica <- NormalizeData(ica)

print("Finding variable genes")
ica <- FindVariableGenes(ica, 
                         mean.function       = ExpMean, 
                         dispersion.function = LogVMR,  
                         do.plot             = T)
hv_genes <- head(rownames(ica@hvg.info), 1000)

print("Scale expression matrix of variable genes")
ica <- ScaleData(ica, 
                 display.progress = T,
                 genes.use        = hv_genes,
                 num.cores        = availableCores())

print("Running PCA")
ica <- RunPCA(object      = ica, 
              pc.genes    = hv_genes, 
              pcs.compute = 50, 
              do.print    = FALSE)
gc()
print("Finding Clusters")
ica <- FindClusters(object         = ica, 
                    reduction.type = "pca", 
                    dims.use       = 1:50, 
                    resolution     = 0.5, 
                    print.output   = T,
                    n.start        = 10, 
                    nn.eps         = 0.5)

print("Computing UMAP embedding")
ica <- RunUMAP(object        = ica, 
               reduction.use = "pca", 
               dims.use      = 1:50)

print("Computing tSNE embedding")
ica <- RunTSNE(object         = ica, 
               reduction.use  = "pca", 
               dims.use       = 1:50, 
               tsne.method    = "FIt-SNE", 
               nthreads       = availableCores(), 
               reduction.name = "FItSNE", 
               reduction.key  = "FItSNE_", 
               fast_tsne_path = "/FIt-SNE/bin/fast_tsne", 
               max_iter       = 2000)

source("/scSVAtools/R/FUNCTIONS.R")

print("Computing diffusion maps")
dmap<-ComputeDMap(t(as.matrix(ica@data[hv_genes,])),
                  nNN             = 10,
                  k               = 50,
                  nLocalsigma     = 5L,
                  nThreads        = availableCores(),
                  M               = 10,
                  efC             = 100,
                  efS             = 100,
                  AnnMethod       = "Nmslib",
                  EigDecompMethod = "Irlba")

print("Getting nearest neighbors")
NNs<-GetANNs(dmap$vectors,
             nNN            = 10,
             nThreads       = availableCores(),
             returnDistance = FALSE,
             M              = 10,
             efC            = 100,
             efS            = 100,
             AnnMethod      = "Nmslib")
NNs<-as.data.frame(cbind(1:nrow(NNs$NN.ind),NNs$NN.ind))
fwrite(x         = NNs,
       file      = "ica_NNs.csv",
       quote     = F,
       row.names = F,
       col.names = F,
       nThread   = availableCores())

print("Computing approximate nearest neighbor graph in the diffusion component space")
mem<-paste0(as.integer(system(command="awk '{if(NR==2) print $2/1e6}' /proc/meminfo",intern=T)),"g")
GetFLE("ica_NNs.csv",
       "/gephi/",
       memory              = mem,
       nsteps              = 5000,
       nThreads            = availableCores(),
       scalingRatio        = 1,
       barnesHutTheta      = 2,
       barnesHutUpdateIter = 10,
       updateCenter        = F,
       barnesHutSplits     = 2,
       restart             = F
)
FLE<-fread("ica_NNs_FLE.txt")
FLE<-FLE[order(FLE$id),]
FLE<-FLE[,-1]
colnames(FLE)<-c("FLE_1","FLE_2","FLE_3")
dmap<-dmap$vectors
colnames(dmap)<-paste0("DC",1:ncol(dmap))

#Attach diffusion components and FLE coordinates to Seurat's metadata
ica@meta.data<-cbind(ica@meta.data,dmap,FLE)

print("Saving results as h5 (loom) file")
source("/scSVAtools/R/Seurat.R")
GenerateInputFiles_fromSeurat(ica,
                              "ICA_BM.h5",
                              ChunkSize        = 1000,
                              CompressionLevel = 7)
print("Done")
