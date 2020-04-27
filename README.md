# MCVR
Implementation of [molecular cross validation](https://www.biorxiv.org/content/10.1101/786269v1) technique to determine ideal number of PCs for scRNA-seq in R.

The intent is for this to be used as part of a Seurat workflow to identify the ideal number of principal components to use in a rigorous way.  Refer to comments in the code for more detailed documentation.  

A simple example of this in a Seurat (V2) workflow would be:


```
source('https://raw.githubusercontent.com/constantAmateur/MCVR/master/code.R')
#Load data and do basic Seurat 
srat = Read10X(srcs)
srat = CreateSeuratObject(srat)
srat = NormalizeData(srat)
srat = FindVariableGenes(srat,do.plot=FALSE)
#Work out the number of PCs
mcv = molecularCrossValidation(srat,normalisation=minimalSeurat)
plotMCV(mcv)
nPCs = mcv$minimas[1]
#Downstream analysis
srat = ScaleData(srat)
srat = RunPCA(srat,pcs.compute=nPCs)
...

```
