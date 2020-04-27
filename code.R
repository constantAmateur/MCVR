library(Matrix)
library(Seurat)
if('parallel' %in% rownames(installed.packages()))
  library(parallel)

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
  x = src@x
  i = src@i
  j = src@j
  if(tFrac<1)
    x = rbinom(length(x),x,tFrac)
  #Adjust p by titration fraction
  p = p *tFrac
  #Now work out p' and p''
  p1 = (1+alpha-sqrt((1+alpha)**2-4*alpha*p))/(2*alpha)
  p2 = p1*alpha
  #Construct the training sample
  trx = rbinom(length(x),x,p1/p)
  #Construct the test sample, including an overlap adjustment.
  tsx = x-trx + rbinom(length(trx),trx,p2)
  #Reformat the training data into the matrix
  w = which(trx>0)
  train = sparseMatrix(x=trx[w],i=i[w]+1,j=j[w]+1,dims=dim(src),dimnames=dimnames(src))
  #Reformat the test data into the matrix
  w = which(tsx>0)
  tst = sparseMatrix(x=tsx[w],i=i[w]+1,j=j[w]+1,dims=dim(src),dimnames=dimnames(src))
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

#' Same as above, but for Seurat V3.
minimalSeuratV3 = function(dat){
  ss = CreateSeuratObject(dat)
  ss = NormalizeData(ss,verbose=FALSE)
  ss = ScaleData(ss,verbose=FALSE)
  return(ss@assays$RNA@scale.data)
}

#' Round stochastically to integers
#' 
#' Round each count up to an integer in a sparse matrix with probability equal to fractional part.
#'
#' @param dat An input matrix of things to round.
#' @return A sparse martix.
roundToInt = function(dat){
  dat = as(dat,'dgTMatrix')
  dat = sparseMatrix(x = floor(dat@x) + rbinom(length(dat@x),1,dat@x%%1),
                     i = dat@i+1,
                     j = dat@j+1,
                     dims = dim(dat),
                     dimnames = dimnames(dat))
  return(dat)
}




#' Finds the optimal number of PCs to use
#'
#' Uses molecular cross validation method (https://www.biorxiv.org/content/10.1101/786269v1), which must be applied to raw count data, to determine the optimal number of PCs to use for a given data-set.  This is intended to be run as part of a Seurat workflow, but is written so that it can be used in a general context.  If supplying a seurat object, FindVariableGenes must have been run, otherwise a set of genes on which to perform the principal component analysis must be supplied.
#'
#' Arbitrary normalisation requires equal splits of data.  To check that 50/50 split does not unde-estimate the number of PCs, it is useful to perform a series of "titrations" of tha data, to check the number of PCs is not sensative to the sampling depth.  There is no relationship between the optimum number of PCs for un-normalised data and normalised data, so it is best to live with the limitations of normalising data (50/50 split, need to titrate) than to do find the optimum for a type of normalisation you will not use in practice.
#'
#' @param dat Data matrix with rows being genes, columns being cells, and counts being integers.
#' @param normalisation The normalisation that will be applied to the counts to produced a matrix of counts that will be fed directly into PCA.  For example, the included function minimalSeurat (or minimalSeuratV3 for Seurat V3) does the standard Seurat pre-processing steps run before PCA.
#' @param varGenes Variable genes to use to perform the principal component analysis.  If NULL all genes are used.
#' @param trainFrac Fraction of data to be used for training.
#' @param p Fraction of total molecules sampled in experiment.
#' @param tFracs Titration fractions.  A vector indicating how to sub-sample the data to test if results are sensative to the depth of sampling.  If NULL, just run one version with all the data.
#' @param nSplits The number of random splits of the data to perform.
#' @param maxPCs Check up to this many PCs.
#' @param approxPCA Use irlba instead of prcomp to speed up PCA at the expense of some accuracy.  Other things take so much longer in the analysis, please don't use approximate PCA.
#' @param errorMetric Type of error to calculate.  Options are mse for mean squared erorr, or poisson, which calculates something proportional to the mean poisson log likelihood.
#' @param poissonMin If the PCs predict a number of counts below this value, use this instead to prevent infinite log-likelihood and penalise predicting zeros/negative values. 
#' @param confInt Used in quantifying the sampling error.  The function returns any PC that is within this confidence interval of the  minimum mean error.
#' @param nCores Number of parallel processes to use.
#' @param ... Extra parameters passed to normalisation function.
#' @return A list.  Contains the average error across all cells for each titration and split of the data and various summaries thereof. 
#'
#' @examples
#' #Assuming srat is a Seurat v2.x object that has had FindVariableGenes run
#' mcv = molecularCrossValidation(srat@raw.data,minimalSeurat,srat@var.genes)
#' #If it is a v3.x object
#' mcv = molecularCrossValidation(srat@assays$RNA@count,minimalSeuratV3,srat@assays$RNA@var.features)
molecularCrossValidation = function(dat,normalisation,varGenes=NULL,trainFrac=0.5,p=0.01,tFracs=c(1,0.9,0.8,0.5),nSplits=5,maxPCs=100,approxPCA=FALSE,errorMetric=c('mse','poisson'),poissonMin=1e-6,confInt=0.95,nCores=1,...){
  if(nCores>1 && (!'parallel' %in% rownames(installed.packages()))){
    warning("Package 'parallel' is required for multi-core execution.  Reverting to single threaded mode.")
    nCores=1
    mclapply=lapply
  }
  #Set option so mclapply acts like lapply
  options(mc.cores=nCores)
  errorMetric = match.arg(errorMetric)
  if(is.null(varGenes))
    varGenes = seq(nrow(dat))
  if(is.null(tFracs))
    tFracs=1
  titrate = length(tFracs)>1 || tFracs<1
  #Convert to number std. errors
  nSEs = qnorm(confInt/2 +0.5)
  #Convert to sparse Matrix of a useful type
  dat = as(dat,'dgTMatrix')
  #Check that input data are integer counts
  if(any(dat@x%%1!=0))
    stop("Input matrix must contain only integers.")
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
    if(errorMetric=='poisson')
      warning("Poisson error is not appropriate for non-count data.  Ensure your normalisation produces counts if you wish to proceed with this error profile.")
  }
  #Work out alpha
  alpha = (1-trainFrac)/trainFrac
  titrates = list()
  itrs = data.frame(tFrac = rep(tFracs,nSplits),
                    split = rep(seq(nSplits),each=length(tFracs))
                    )
  titrates = mclapply(seq(nrow(itrs)),
                      FUN = function(e,...){
    tFrac = itrs$tFrac[e]
    i = itrs$split[e]
    mse = rep(NA,maxPCs-1)
  #for(tFrac in tFracs){
    if(titrate)
      message(sprintf("Running with %d%% of data",tFrac*100))
    #Calculate correction factor alpha
    #mse = matrix(NA,nrow=nSplits,ncol=maxPCs-1)
    #colnames(mse) = seq(2,maxPCs)
    #rownames(mse) = paste0("Split",seq(nSplits))
    #tt = mclapply(seq(nSplits),function(i) {
    #for(i in seq(nSplits)){
      message(sprintf("Performing split %d of %d",i,nSplits)) 
      tst = partitionData(dat,alpha,p,tFrac)
      #Normalise data
      train = normalisation(tst[[2]],...)[varGenes,]
      tst = normalisation(tst[[1]],...)[varGenes,]
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
          #mse[i,k-1] = mean((tmp*alpha-tst)**2)
          mse[k-1] = mean((tmp*alpha-tst)**2)
        }else if(errorMetric=='poisson'){
          a = as.vector(tmp*alpha)
          b = as.vector(tst)
          #Fix those that are below minimum
          a[a<poissonMin]=poissonMin
          #Calculate poisson negative log likelihood
          #mse[i,k-1] = -1*mean(dpois(b,a,log=TRUE))
          mse[k-1] = -1*mean(dpois(b,a,log=TRUE))
        }
      }
      close(pb)
      #return(mse)
    #})
    return(mse)
    #titrates[[length(titrates)+1]] = mse
  },...)
  #Merge into titrates
  tmp = lapply(tFracs,function(e) seq(nrow(itrs))[itrs$tFrac==e])
  names(tmp) = tFracs
  titrates = lapply(tmp,function(e) do.call(rbind,titrates[e]))
  names(titrates) = names(tmp)
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
  seMins = t(t(means-nSEs*ses) <= apply(means,2,min))
  seMins = apply(seMins,2,which)
  #Ensure formating the same
  if(length(titrates)==1){
    seMins = list(as.numeric(rownames(seMins)))
  }else{
    seMins = lapply(seMins,function(e) as.numeric(names(e)))
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

#' Plots summary of MCV run
#' 
#' Plots normalised loss curves for each titration of data as returned by molecularCrossValidation.
#' 
#' @param mcv Output of molecularCrossValidation.
#' @param cols Colours to apply to each titration line.
#' @return Nothing.  Produces a plot.
plotMCV = function(mcv,cols=seq(ncol(mcv$lower))){
  pcs = as.numeric(rownames(mcv$means))
  titrates = mcv$errors
  if(!is.list(titrates))
    titrates = list(titrates)
  #Work out normalisation factors
  normFacs = apply(mcv$lower,2,min)
  #Define plot area
  yRange = c(1,max(t(mcv$upper)/normFacs))
  layout(matrix(c(1,1,2),nrow=3))
  plot(0,0,type='n',xlab='Number of PCs',ylab='avg error / min avg error',ylim=yRange,xlim=range(pcs),frame.plot=FALSE,main='Error profiles')
  #Make lower,middle and upper plots
  for(i in seq_along(titrates)){
    lines(pcs,mcv$lower[,i]/normFacs[i],col=cols[i],lty=2)
    lines(pcs,mcv$means[,i]/normFacs[i],col=cols[i])
    lines(pcs,mcv$upper[,i]/normFacs[i],col=cols[i],lty=2)
  }
  #abline(v=mcv$minimas,col=seq(ncol(x)))
  legend(mean(mcv$minimas),yRange[2],legend=paste0(round(mcv$titrations*100),'% Data (Min = ',mcv$minimas,')'),col=cols,lty=1)
  plot(mcv$titrations*100,mcv$minimas,pch=19,xlab='Data percentage',ylab='Optimal #PCs',main='Convergance',ylim=range(pcs[unlist(mcv$mins1sd)]))
  lines(mcv$titrations*100,mcv$minimas)
  #Add error bars
  arrows(mcv$titrations*100, pcs[sapply(mcv$mins1sd,min)], mcv$titrations*100,pcs[sapply(mcv$mins1sd,max)], length=0.05, angle=90, code=3)
}
