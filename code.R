library(Matrix)
library(Seurat)

#' Create test and training data splits.
#'
#' Splits sample of full molecular pool into two samples from larger molecular pool (for test and training), accounting for the overlap between samples (p).
#'
#' @param src The input source matrix.
#' @param alpha Ratio of training to test fraction.
#' @param p Fraction of total molecules sampled.
#' @param tFrac Titration fraction.
#' @return Two sparse matrices, one for the test, one for the training.
partitionData = function(src,alpha,p,tFrac){
  #Do we need to titrate?  If so do that first.
  src = as(src,'dgTMatrix')
  if(tFrac<1)
    src@x = rbinom(length(src@x),src@x,tFrac)
  #Adjust p by titration fraction
  p = p *tFrac
  #Now work out p' and p''
  p1 = (1+alpha-sqrt((1+alpha)**2-4*alpha*p))/(2*alpha)
  p2 = p1*alpha
  #Construct the training sample
  train = src
  train@x = rbinom(length(train@x),train@x,p1/p)
  #Construct the test sample, including an overlap adjustment.
  tst = src
  tst@x = tst@x-train@x + rbinom(length(train@x),train@x,p2)
  #Reformat the training data into the matrix
  w = which(train@x>0)
  train = sparseMatrix(x=train@x[w],i=train@i[w]+1,j=train@j[w]+1,dims=dim(train),dimnames=dimnames(train))
  #Reformat the test data into the matrix
  w = which(tst@x>0)
  tst = sparseMatrix(x=tst@x[w],i=tst@i[w]+1,j=tst@j[w]+1,dims=dim(tst),dimnames=dimnames(tst))
  return(list(tst,train))
}

#' Minimal Seurat processing
#'
#' Does the steps up to scaling the matrix in Seurat.
#'
#' @param dat Sparse matrix to process.
#' @return The scaled matrix.
minimalSeurat = function(dat){
  ss = CreateSeuratObject(dat)
  ss = NormalizeData(ss,display.progress=FALSE)
  ss = ScaleData(ss,display.progress=FALSE)
  return(ss@scale.data)
}

#' Finds the optimal number of PCs to use
#'
#' Uses molecular cross validation method (https://www.biorxiv.org/content/10.1101/786269v1), which must be applied to raw count data, to determine the optimal number of PCs to use for a given data-set.  This is intended to be run as part of a Seurat v2 workflow, but is written so that it can be used in a general context.  If supplying a seurat object, FindVariableGenes must have been run, otherwise a set of genes on which to perform the principal component analysis must be supplied.
#'
#' Arbitrary normalisation requires equal splits of data.  To check that 50/50 split does not unde-estimate the number of PCs, it is useful to perform a series of "titrations" of tha data, to check the number of PCs is not sensative to the sampling depth.  There is no relationship between the optimum number of PCs for un-normalised data and normalised data, so it is best to live with the limitations of normalising data (50/50 split, need to titrate) than to do find the optimum for a type of normalisation you will not use in practice.
#'
#' @param srat Seurat object which contains the raw data and has been processed up to "FindVariableGenes" or a sparse matrix with raw UMI counts, with cells as columns and genes as rows.  If the latter, varGenes must also be specified.
#' @param varGenes Variable genes to use to perform the principal component analysis.  If missing, \code{srat} is assumed to be a Seurat object and the variable genes are extracted from it.
#' @param trainFrac Fraction of data to be used for training.
#' @param p Fraction of total molecules sampled in experiment.
#' @param tFracs Titration fractions.  A vector indicating how to sub-sample the data to test if results are sensative to the depth of sampling.  If NULL, just run one version with all the data.
#' @param nSplits The number of random splits of the data to performe.
#' @param maxPCs Check up to this many PCs.
#' @param approxPCA Use irlba instead of prcomp to speed up PCA at the expense of some accuracy.  Other things take so much longer in the analysis, please don't use approximate PCA.
#' @param errorMetric Type of error to calculate.  Options are mse for mean squared erorr, or poisson, which calculates something proportional to the mean poisson log likelihood.
#' @param nSEs To quantify the sampling error, this function also returns those PCs within this many standard errors of the global minimum mean error. 1.96 corresponds to 95% confidence interval.
#' @param doPlot Make the standard plots?
#' @return A list.  Contains the average error across all cells for each titration and split of the data and various summaries thereof. 
molecularCrossValidation = function(srat,varGenes,trainFrac=0.5,p=0.01,tFracs=c(1,0.9,0.8,0.5,0.1),nSplits=5,normalisation=minimalSeurat,maxPCs=100,approxPCA=FALSE,errorMetric=c('mse','poisson'),nSEs=1.96,doPlot=TRUE){
  errorMetric = match.arg(errorMetric)
  if(missing(varGenes))
    varGenes = srat@var.genes
  #If it's a seurat object, extract data
  if(is.null(dim(srat))){
    dat = srat@raw.data
  }else{
    dat = srat
  }
  if(is.null(tFracs))
    tFracs=1
  titrate = length(tFracs)>1 || tFracs<1
  #If normalisation is on, have to do equal splits so we don't need to calculate complicated scale factor
  if(!identical(normalisation,identity)){
    #Arbitrary normalisation, have to have 50/50 split.
    if(trainFrac!=0.5){
      warning("Arbitrary normalisation requires equal splits of data. Setting trainFrac=0.5")
      trainFrac=0.5
    }
    if(!titrate){
      warning("When performing arbitrary normalisation it is useful to perform several titrations (by setting tFracs) to ensure results are not sensative to depth.  Consider re-running with tFracs set.")
    }
  }
  #Work out alpha
  alpha = (1-trainFrac)/trainFrac
  titrates = list()
  for(tFrac in tFracs){
    if(titrate)
      message(sprintf("Running with %d%% of data",tFrac*100))
    #Calculate correction factor alpha
    mse = matrix(NA,nrow=nSplits,ncol=maxPCs-1)
    colnames(mse) = seq(2,maxPCs)
    rownames(mse) = paste0("Split",seq(nSplits))
    for(i in seq(nSplits)){
      message(sprintf("Performing split %d of %d",i,nSplits)) 
      tst = partitionData(dat,alpha,p,tFrac)
      #Normalise data
      train = normalisation(tst[[2]])[varGenes,]
      tst = normalisation(tst[[1]])[varGenes,]
      #Normalise data
      #Run PCA on training data.
      message("Running PCA and cross-validating")
      if(approxPCA){
        pca = irlba(train,nv=maxPCs)
        rot = t(pca$v)
        embed = pca$u
      }else{
        pca = prcomp(train,center=FALSE,rank.=maxPCs)
        rot = t(pca$rotation)
        embed = pca$x
      }
      pb = txtProgressBar(min = 2,max=maxPCs,style=1)
      for(k in seq(maxPCs,2)){
        setTxtProgressBar(pb,maxPCs-k+2)
        #message(sprintf("Cross-validating for %d components",k))
        tmp = embed[,seq(k)] %*% rot[seq(k),]
        #Mean squared error
        if(errorMetric=='mse'){
          mse[i,k-1] = mean((tmp*alpha-tst)**2)
        }else if(errorMetric=='poisson'){
          a = as.vector(tmp*alpha)
          b = as.vector(tst)
          #Ignore the entries that would give infinite likelihood (where mean<=0 and tst>0)
          w = which(a>0 | b==0)
          a=a[w]
          b=b[w]
          #The ones where both are 0 need separate processing
          ww = which(a<=0 & b==0)
          nww = which(!(a<=0 & b==0))
          mse[i,k-1] = (sum(a[ww]) + sum(a[nww] - b[nww]*log(a[nww])))/length(w)
        }else{
          stop('Oh No!')
        }
      }
      close(pb)
    }
    titrates[[length(titrates)+1]] = mse
  }
  #Now work out which PC did the best
  #Work out mean and limits
  lower = do.call(cbind,lapply(titrates,apply,2,min))
  means = do.call(cbind,lapply(titrates,apply,2,mean))
  upper = do.call(cbind,lapply(titrates,apply,2,max))
  #Get the minima for each
  mins = seq(2,maxPCs)[apply(means,2,which.min)]
  #Get the standard error of each estimate
  ses = lapply(titrates,function(e) apply(e,2,sd)/sqrt(nSplits))
  ses = do.call(cbind,ses)
  #Find out the PCs that are within x standard errors of the minimum mean error
  seMins = t(t(means-nSEs*ses) < apply(means,2,min))
  seMins = apply(seMins,2,which)
  if(doPlot){
    #Work out normalisation factors
    normFacs = apply(lower,2,min)
    #Define plot area
    yRange = c(1,max(t(upper)/normFacs))
    layout(matrix(c(1,1,2),nrow=3))
    plot(0,0,type='n',xlab='Number of PCs',ylab='avg error / min avg error',ylim=yRange,xlim=c(2,maxPCs),frame.plot=FALSE,main='Error profiles')
    #Make lower,middle and upper plots
    for(i in seq_along(titrates)){
      lines(seq(2,maxPCs),lower[,i]/normFacs[i],col=i,lty=2)
      lines(seq(2,maxPCs),means[,i]/normFacs[i],col=i)
      lines(seq(2,maxPCs),upper[,i]/normFacs[i],col=i,lty=2)
    }
    #abline(v=mins,col=seq(ncol(x)))
    legend(mean(mins),yRange[2],legend=paste0(round(tFracs*100),'% Data, ',round(trainFrac*100),'% train (Min = ',mins,')'),col=seq_along(titrates),lty=1)
    plot(tFracs*100,mins,pch=19,xlab='Data percentage',ylab='Optimal #PCs',main='Convergance',ylim=range(seq(2,maxPCs)[unlist(seMins)]))
    lines(tFracs*100,mins)
    #Add error bars
    arrows(tFracs*100, seq(2,maxPCs)[sapply(seMins,min)], tFracs*100,seq(2,maxPCs)[sapply(seMins,max)], length=0.05, angle=90, code=3)
  }
  if(length(titrates)==1)
    titrates=titrates[[1]]
  out = list(errors = titrates,
             means = means,
             lower = lower,
             upper = upper,
             minimas = mins,
             mins1sd = seMins,
             titrations = tFracs
             )
  return(out)
}


