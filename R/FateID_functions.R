#' @title Feature selection based on differentially expressed genes
#'
#' @description This function performs a feature selection based on the inference of differentially expressed genes between each target cluster and all remaining cells.
#' @param x expression data frame with genes as rows and cells as columns. Gene IDs should be given as row names and cell IDs should be given as column names. This can be a reduced expression table only including the features (genes) to be used in the analysis.
#' @param y clustering partition. A vector with an integer cluster number for each cell. The order of the cells has to be the same as for the columns of x.
#' @param tar vector of integers representing target cluster numbers. Each element of \code{tar} corresponds to a cluster of cells committed towards a particular mature state. One cluster per different cell lineage has to be given and is used as a starting point for learning the differentiation trajectory.
#' @param fpv p-value cutoff for calling differentially expressed genes. This is a cutoff for the Benjamini-Hochberg corrected false discovery rate. Default value is 0.05.
#' @param  ... additional arguments to be passed to the low level function \code{diffexpnb}.
#' @details The function determines differentially expressed between the cells in each of the target clusters in comparison to the remaining cells by using \code{diffexpnb} function.
#' @return A filtered expression table with features extracted based on differentially expressed genes.
#' @examples
#' x <- intestine$x
#' y <- intestine$y
#' tar <- c(6,9,13)
#' xf <- getFeat(x,y,tar,fpv=.05)
#' @export
getFeat <- function(x,y,tar,fpv=.05,...){
  g    <- c()
  for ( i in tar ){
    n <- names(y)
    A <- names(y)[y != i]
    B <- names(y)[y == i]
    d <- diffexpnb(x,A,B,...)
    ga <- rownames(d$res)[d$res$foldChange > 1 & d$res$padj < fpv]
    if (length(ga)>0){
      g <- c(g,ga[!ga %in% g])
    }
  }
  xf <- x[g,]
  xf
}

#' @title Inference of a cell type partition
#'
#' @description This function performs an inference of a cell type partition based on the expression of marker genes.
#' @param x expression data frame with genes as rows and cells as columns. Gene IDs should be given as row names and cell IDs should be given as column names. This can be a reduced expression table only including the features (genes) to be used in the analysis.
#' @param FMarker list of vectors of gene IDs corresponding to valid rownames of \code{x}. The gene IDs within each component of \code{FMarker} are considered as marker genes of one of the cell types in the dataset. The aggregated expression of the genes for each component is compared to a threshold defined by the input argument \code{fthr} or \code{n}. All cells exceeding this threshold are assigned to a cluster representing cells with expression of the respective marker genes.
#' @param fthr vector of real positive numbers. This vector has to have the same length as the list \code{FMarker} and contains a threshold for the aggregated expression of all genes in the corresponding component of \code{FMarker}. If NULL then a threshold is inferred from the \code{n} top-expressing cells for the genes in the respective component of \code{FMarker}.
#' @param n positive integer number. For each component of \code{FMarker} the expression of all genes is aggregated in every cell and the \code{n} top-expressing cells are extracted. The average expression across these cell defines the expression threshold applied to infer the partitioning. Default value is 25.
#' @return A list with the following three components:
#'   \item{part}{A vector with a partitioning, i. e. cluster assignment for each cell.}
#'   \item{tar}{A vector with the numbers of target clusters. Cluster 1 comprises all cells with no enrichment of marker genes. The remaining clusters correspond to cell types up-regulating the markers in the list \code{FMarker} in the same order as in this list.}
#'   \item{thr}{A vector with expression threshold values applied for each component in the list \code{FMarker} in the same order as in this list.}
#' @examples
#' x <- intestine$x
#' y <- intestine$y
#' 
#' FMarker <- list(c("Defa20__chr8","Defa24__chr8"), "Clca3__chr3", "Alpi__chr1")
#' xf <- getPart(x,FMarker,fthr=NULL,n=5)
#'
#' @importFrom utils head
#' @export
getPart <- function(x,FMarker,fthr=NULL,n=25){
  y <- rep(1,ncol(x))
  names(y) <- names(x)
  if ( is.null(fthr) ){ flag <- TRUE}
  for ( i in 1:length(FMarker)){
    if ( flag ){
      if ( length(FMarker[[i]]) == 1 ){
        u <- t(x[FMarker[[i]],])
      }else{
        u <- apply(x[FMarker[[i]],],2,sum)
      }
      u <- u[order(u,decreasing=TRUE)]
      fthr[i] <- mean(head(u,n))
    }
    if ( length(FMarker[[i]]) == 1 ){
      y[t(x[FMarker[[i]],]) > fthr[i]] <- i+1
    }else{
      y[apply(x[FMarker[[i]],],2,sum) > fthr[i]] <- i+1
    }
  }
  tar <- 2:(length(FMarker) + 1 )
  return(list(part=y,tar=tar,fthr=fthr))
}

#' @title Reclassification of cells
#'
#' @description This function attempts to reassign additional cells in the dataset to one of the target clusters.
#' @param x expression data frame with genes as rows and cells as columns. Gene IDs should be given as row names and cell IDs should be given as column names. This can be a reduced expression table only including the features (genes) to be used in the analysis.
#' @param y clustering partition. A vector with an integer cluster number for each cell. The order of the cells has to be the same as for the columns of x.
#' @param tar vector of integers representing target cluster numbers. Each element of \code{tar} corresponds to a cluster of cells committed towards a particular mature state. One cluster per different cell lineage has to be given and is used as a starting point for learning the differentiation trajectory.
#' @param z Matrix containing cell-to-cell distances to be used in the fate bias computation. Default is \code{NULL}. In this case, a correlation-based distance is computed from \code{x} by \code{1 - cor(x)}
#' @param clthr real number between zero and one. This is the threshold for the fraction of random forest votes required to assign a cell not contained within the target clusters to one of these clusters. The value of this parameter should be sufficiently high to only reclassify cells with a high-confidence assignment. Default value is 0.9.
#' @param nbfactor positive integer number. Determines the number of trees grown for each random forest. The number of trees is given by the number of columns of th training set multiplied by \code{nbfactor}. Default value is 5.
#' @param use.dist logical value. If \code{TRUE} then the distance matrix is used as feature matrix (i. e. \code{z} if not equal to \code{NULL} and \code{1-cor(x)} otherwise). If \code{FALSE}, gene expression values in \code{x} are used. Default is \code{FALSE}.
#' @param seed integer seed for initialization. If equal to \code{NULL} then each run will yield slightly different results due to the radomness of the random forest algorithm. Default is \code{NULL}
#' @param nbtree integer value. If given, it specifies the number of trees for each random forest explicitely. Default is \code{NULL}.
#' @param q real value between zero and one. This number specifies a threshold used for feature selection based on importance sampling. A reduced expression table is generated containing only features with an importance larger than the q-quantile for at least one of the classes (i. e. target clusters). Default value is 0.75.
#' @param  ... additional arguments to be passed to the low level function \code{randomForest}.
#' @details The function uses random forest based supervised learning to assign cells not contained in the target clusters to one of these clusters. All cells not within any of the target clusters which receive a fraction of votes larger than \code{clthr} for one of the target clusters will be reassigned to this cluster. Since this function is developed to reclassify cells only if they can be assigned with high confidence, a high value of \code{clthr} (e. g. > 0.75) should be applied.
#' @return A list with the following three components:
#'   \item{part}{A vector with the revised cluster assignment for each cell in the same order as in the input argument \code{y}.}
#'   \item{rf}{The random forest object generated for the reclassification, with enabled importance sampling (see \pkg{randomForest}).}
#'   \item{xf}{A filtered expression table with features extracted based on the important samples, only features with an importance larger than the q-quantile are for at least one of the classes are retained.}
#' @examples
#' x <- intestine$x
#' y <- intestine$y
#' tar <- c(6,9,13)
#' rc <- reclassify(x,y,tar,z=NULL,nbfactor=5,use.dist=FALSE,seed=NULL,nbtree=NULL,q=.9)
#' @importFrom stats quantile cor
#' @importFrom randomForest randomForest
#' @export
reclassify <- function(x,y,tar,z=NULL,clthr=.75,nbfactor=5,use.dist=FALSE,seed=NULL,nbtree=NULL,q=.9,...){
  if (!is.null(seed) ) set.seed(seed)
  
  # create cell-to-cell distances if not supplied
  di <- if ( is.null(z) ) as.data.frame( 1 - cor(x) ) else as.data.frame(as.matrix(z))

  # order of columns/rows in z has to be the same as order of columns in x
  names(di) <- rownames(di) <- colnames(x)
  names(y)  <- colnames(x)
  # assign whether distances or expression values should be used for random forest classification
  if ( use.dist ){
    d <- di
  }else{
    d <- x
  }
  trt <- names(y)[y %in% tar]
  n   <- names(y)[!y%in% tar]
  if ( use.dist ){
      # reference set
    xr <- di[trt,trt]
      # test set
    xt <- di[n,trt]
  }else{
    xr <- t(d[,trt])
    xt <- t(d[,n])
  }
  pr <- y[trt]

  if ( is.null(nbtree) ) nbtree = nrow(xr)*nbfactor
  rf <- randomForest(xr,as.factor(pr),xt,nbtree=nbtree,norm.votes=TRUE,importance=TRUE,...)
  
  tv <- as.data.frame(rf$test$votes)
  names(tv) <- paste("t",colnames(tv),sep="")
  for ( k in names(tv) ) y[rownames(tv)[tv[,k] > clthr]] <- as.numeric(sub("t","",k))

  g <- c()
  for ( i in tar ){
    k <- rf$importance[,as.character(i)]
    k <- k[order(k,decreasing=TRUE)]
    g <- append(g,names(k)[k > quantile(k,q)])
  }
  g <- unique(g)
  return(list(part=y,rf=rf,xf=x[g,]))
}

#' @title Computation of fate bias
#'
#' @description This function computes fate biases for single cells based on expression data from a single cell sequencing experiment. It requires a clustering partition and a target cluster representing a commited state for each trajectory. 
#' @param x expression data frame with genes as rows and cells as columns. Gene IDs should be given as row names and cell IDs should be given as column names. This can be a reduced expression table only including the features (genes) to be used in the analysis.
#' @param y clustering partition. A vector with an integer cluster number for each cell. The order of the cells has to be the same as for the columns of x.
#' @param tar vector of integers representing target cluster numbers. Each element of \code{tar} corresponds to a cluster of cells committed towards a particular mature state. One cluster per different cell lineage has to be given and is used as a starting point for learning the differentiation trajectory.
#' @param z Matrix containing cell-to-cell distances to be used in the fate bias computation. Default is \code{NULL}. In this case, a correlation-based distance is computed from \code{x} by \code{1 - cor(x)}.
#' @param minnr integer number of cells per target cluster to be selected for classification (test set) in each iteration. For each target cluster, the \code{minnr} cells with the highest similarity to a cell in the training set are selected for classification. If \code{z} is not \code{NULL} it is used as the similarity matrix for this step. Otherwise, \code{1-cor(x)} is used. Default value is \code{NULL} and \code{minnr} is estimated as the minimum of and 20 and half the median of target cluster sizes. 
#' @param minnrh integer number of cells from the training set used for classification. From each training set, the \code{minnrh} cells with the highest similarity to the training set are selected. If \code{z} is not \code{NULL} it is used as the similarity matrix for this step. Default value is \code{NULL} and \code{minnrh} is estimated as the maximum of and 20 and half the median of target cluster sizes.
#' @param adapt logical. If \code{TRUE} then the size of the test set for each target cluster is adapted based on the classification success in the previous iteration. For each target cluster, the number of successfully classified cells is determined, i.e. the number of cells with a minimum fraction of votes given by the \code{confidence} parameter for the target cluster, which gave rise to the inclusion of the cell in the test set (see \code{minnr}). Weights are then derived by dividing this number by the maximum across all clusters after adding a pseudocount of 1. The test set size \code{minnr} is rescaled for each cluster by the respective weight in the next iteration. Default is \code{TRUE}.
#' @param confidence real number between 0 and 1. See \code{adapt} parameter. Default is 0.75.
#' @param nbfactor positive integer number. Determines the number of trees grown for each random forest. The number of trees is given by the number of columns of th training set multiplied by \code{nbfactor}. Default value is 5.
#' @param use.dist logical value. If \code{TRUE} then the distance matrix is used as feature matrix (i. e. \code{z} if not equal to \code{NULL} and \code{1-cor(x)} otherwise). If \code{FALSE}, gene expression values in \code{x} are used. Default is \code{FALSE}.
#' @param seed integer seed for initialization. If equal to \code{NULL} then each run will yield slightly different results due to the radomness of the random forest algorithm. Default is \code{NULL}
#' @param nbtree integer value. If given, it specifies the number of trees for each random forest explicitely. Default is \code{NULL}.
#' @param  ... additional arguments to be passed to the low level function \code{randomForest}.
#' @details The bias is computed as the ratio of the number of random forest votes for a trajectory and the number of votes for the trajectory with the second largest number of votes. By this means only the trajectory with the largest number of votes will receive a bias >1. The siginifcance is computed based on counting statistics on the difference in the number of votes. A significant bias requires a p-value < 0.05. Cells are assigned to a trajectory if they exhibit a significant bias >1 for this trajectory.
#' @return A list with the following three components:
#'   \item{probs}{a data frame with the fraction of random forest votes for each cell. Columns represent the target clusters. Column names are given by a concatenation of \code{t} and target cluster number.}
#'   \item{votes}{a data frame with the number of random forest votes for each cell. Columns represent the target clusters. Column names are given by a concatenation of \code{t} and target cluster number.}
#'   \item{tr}{list of vectors. Each component contains the IDs of all cells on the trajectory to a given target cluster. Component names are given by a concatenation of \code{t} and target cluster number.}
#'   \item{rfl}{list of randomForest objects for each iteration of the classification.}
#'   \item{trall}{vector of cell ids ordered by the random forest iteration in which they have been classified into one of the target clusters.}
#' @examples
#' x <- intestine$x
#' y <- intestine$y
#' tar <- c(6,9,13)
#' fb <- fateBias(x,y,tar,minnr=5,minnrh=20,adapt=TRUE,confidence=0.75,nbfactor=5)
#' head(fb$probs)
#' @importFrom stats binom.test cor median 
#' @importFrom utils head
#' @export
fateBias <- function(x,y,tar,z=NULL,minnr=NULL,minnrh=NULL,adapt=TRUE,confidence=0.75,nbfactor=5,use.dist=FALSE,seed=NULL,nbtree=NULL,...){
  if ( is.null(minnrh) ){
    minnrh <- max(round(median(aggregate(rep(1,sum(y %in% tar)),by=list(y[y %in% tar]),FUN=sum)[,2])/2,0),20)
  }
  if ( is.null(minnr) ){
    minnr  <- min(round(median(aggregate(rep(1,sum(y %in% tar)),by=list(y[y %in% tar]),FUN=sum)[,2])/2,0),20)
  }
  cat("minnr:",minnr,"\n")
  cat("minnrh:",minnrh,"\n")

  if (!is.null(seed) ) set.seed(seed)

  if (!is.null(z) ) z <- as.matrix(z)
  names(y) <- colnames(x)
  tr    <- list()
  trall <- c()
  for ( k in tar){
    tr[[paste("t",k,sep="")]] <- names(y)[y == k]
    trall <- append(trall, names(y)[y == k])
  }

  # create cell-to-cell distances if not supplied
  di <- if ( is.null(z) ) as.data.frame( 1 - cor(x) ) else as.data.frame(as.matrix(z))

  # order of columns/rows in z has to be the same as order of columns in x
  colnames(di) <- rownames(di) <- colnames(x)

  # assign whether distances or expression values should be used for random forest classification
  if ( use.dist ){
    d <- di
  }else{
    d <- x
  }

  # initialize matrix with probability of every cell to be assigned to each of the trajectories
  # cells of the target clusters are assigned to the target trajectory with probability 1
  probs <- as.data.frame(matrix(rep(1/length(tar),ncol(d)*length(tar)),ncol=length(tar)))
  rownames(probs) <- colnames(d)
  colnames(probs)    <- paste("t",tar,sep="")
  
  for ( k in tar){
    probs[colnames(d)[y == k],paste("t",k,sep="")] <- 1
    for ( j in tar){
      if ( j != k ) probs[colnames(d)[y == k],paste("t",j,sep="")] <- 0
    }
  }
  # initialize matrix of votes
  votes <- probs

  # iterate random forest classification until all cells are assigned
  i <- 0
  rfl <- list()
  weights <- rep(1,length(tar))
  names(weights) <- paste("t",tar,sep="")
  
  while ( length(trall) < ncol(d) ){
    i <- i + 1
    cat("test set size iteration",i,":",weights*minnr,"\n")
    # extract for each target trajectory the minnr cells from the pool of non-assigned cells closest to any of the cells on the trajectory
    u <- di[trall,colnames(d)[ ! colnames(d) %in% trall]]
    nl <- list()
    if ( sum(! colnames(d) %in% trall ) > 1 ){
      n <- c()
      for ( k in tar ){
        f <- trall %in% tr[[paste("t",k,sep="")]]
        v <- apply(u[f,],2,median)
        v <- v[order(v,decreasing=FALSE)]
        minnrL <- max(round(minnr * weights[paste("t",k,sep="")],0),1)
        nl[[k]] <- head(names(v),minnrL)
        n <- append(n,head(names(v),minnrL)) 
      }
      n <- unique(n)
    }else{
      n <- colnames(d)[ ! colnames(d) %in% trall]
      for ( k in tar ){
        nl[[k]] <- n
      }
    }

    # extract from each target trajectory the minnrh cells closest to any of the cells to be assigned in this iteration
    trt <- c()
    
    for ( k in tar){
      u <- di[nl[[k]],tr[[paste("t",k,sep="")]]]
      v <- apply(u,2,min)
      v <- v[order(v,decreasing=FALSE)]
      trt <- append(trt,head(names(v),minnrh))
    }
    trt <- unique(trt)  
    cat("randomforest iteration",i,"of", length(n),"cells\n")

    # assign feature matrix for random forest classification
    if ( use.dist ){
      # reference set
      xr <- di[trt,trt]
      # training set
      xt <- di[n,trt]
    }else{
      xr <- t(d[,trt])
      xt <- t(d[,n])
    }
    # assign partition for reference set
    pr <- apply(probs[trt,],1,function(x,y) y[which(x == max(x))][1],y=tar)
 
    if ( is.null(nbtree) ) nbtree = ncol(xr)*nbfactor
    rf <- randomForest(xr,as.factor(pr),xt,nbtree=nbtree,norm.votes=FALSE,importance=TRUE)
    rfl[[i]] <- rf
    # update probability matrix based on random forest votes for test set
    tv <- as.data.frame(rf$test$votes)
    names(tv) <- paste("t",colnames(tv),sep="")
    for ( np in colnames(probs) ){
      if ( ! np %in% colnames(tv) ){
        tv[,np] <- rep(0,nrow(tv))
      }
    }
    tvn <- tv
    tv  <- tv/apply(tv,1,sum)
    probs[n,] <- tv[,colnames(probs)]
    if ( i == 1 ) votes <- votes*max(apply(tvn,1,sum))
    votes[n,] <- tvn[,colnames(votes)]

    ## update weights
    if ( adapt ){
        for ( k in tar ){
            l <- paste("t",k,sep="")
            weights[l] <- sum( probs[n,l] > confidence ) + 1
            ##weights[l] <- sum( probs[n,l] > min(.9,2*apply(probs[n, colnames(probs) != l],1,max) ) ) + 1
            ##weights[l] <- max( sum( probs[nl[[k]],l] > 2*apply(probs[nl[[k]], colnames(probs) != l],1,max) ),1)/max(length(nl[[k]]),1)
        }
        weights <- weights/max(weights) 
    }
    # update trajectories with cells that exhibit significant bias
    for ( k in names(votes)){
      b <- bias(tvn)
      tr[[k]] <- append(tr[[k]],rownames(tvn)[b$bias[,k] > 1 & b$pv < .05])
    }
    trall <- append(trall, n)
  }
 
  return(list(probs=probs,votes=votes,tr=tr,rfl=rfl,trall=trall))
}

bias <- function(tvn){
  bias <- tvn/apply(tvn,1,function(x) x[order(x,decreasing=TRUE)][2] + 1e-3)
  pv <- apply(tvn,1,function(x){ h <- x[order(x,decreasing=TRUE)][1];  l <- x[order(x,decreasing=TRUE)][2]; binom.test( c(h, l), alternative="t" )$p.value})
  return(list(bias=bias,pv=pv))
}

#' @title Computation of dimensional reduction representations
#'
#' @description This function computes dimensional reduction representations to a specified number of dimensions using a number of different algorithms: t-SNE, cmd, lle, diffusion maps 
#' @param x expression data frame with genes as rows and cells as columns. Gene IDs should be given as row names and cell IDs should be given as column names. This can be a reduced expression table only including the features (genes) to be used in the analysis.
#' @param z Matrix containing cell-to-cell distances to be used in the fate bias computation. Default is \code{NULL}. In this case, a correlation-based distance is computed from \code{x} by \code{1 - cor(x)}
#' @param m a vector of dimensional reduction representations to be computed. By default, the following representations are computed: \code{lle} (locally-linear embedding), \code{cmd} (classical multidimensional scaling), \code{dm} (diffusion map), \code{tsne} (t-SNE map). The default value of m is \code{c("lle","cmd","dm","tsne")}. Any subset of these methods can be selected.
#' @param k vector of integers representing the dimensions for which the dimensional reduction representations will be computed. Default value is \code{c(2,3)}.
#' @param lle.n integer number for the number of neighbours used in the \code{lle} algorithm. Default value is 30.
#' @param dm.sigma parameter for the computation of the diffusion maps with the destiny package. See \url{https://bioconductor.org/packages/devel/bioc/html/destiny.html}. Default value is 1000.
#' @param dm.distance parameter for the computation of the diffusion maps with the destiny package. See \url{https://bioconductor.org/packages/devel/bioc/html/destiny.html}. Default value is \code{euclidean}.
#' @param tsne.perplexity positive number. Perplexity used in the t-SNE computation. Default value is 30.
#' @param seed integer seed for initialization. If equal to \code{NULL} then each run will yield slightly different results due to the randomness of the random forest algorithm. Default is \code{NULL}
#' @return A two-dimensional list with the dimensional reduction representation stored as data frames as components. Component names for the first dimension are given by one of the following algorithms:
#'   \item{lle}{locally linear embedding calculated by the lle function from the \pkg{lle} package.}
#'   \item{cmd}{classical multidimensional scaling computed by the \code{cmdscale} function of the \pkg{stats} package.}
#'   \item{dm}{diffusion map computed by the \code{DiffusionMap} function of the \pkg{destiny} package.}
#'   \item{tsne}{t-SNE map computed by the \code{Rtsne} function of the \pkg{Rtsne} package.}
#'
#' Component names of the second dimension are a concatenation of a capital D and an integer number of the dimension. There is one component for each dimension in \code{k}.
#' @examples
#'
#' x <- intestine$x
#' dr <- compdr(x,z=NULL,m="cmd",k=2,lle.n=30,dm.sigma=1000,dm.distance="euclidean",tsne.perplexity=30)
#' plot(dr[["cmd"]][["D2"]],pch=20,col="grey")
#'
#' @importFrom stats cor cmdscale as.dist
#' @importFrom lle lle
#' @importFrom Rtsne Rtsne
#' @export
compdr <- function(x,z=NULL,m=c("tsne","cmd","dm","lle"),k=c(2,3),lle.n=30,dm.sigma=1000,dm.distance="euclidean",tsne.perplexity=30,seed=12345){
  if (!is.null(seed) ) set.seed(seed)

  dr <- list()
  d  <- t(x)

  # create cell-to-cell distances if not supplied
  di <- if ( is.null(z) )  1 - cor(x) else as.matrix(z)

 
  # calculate diffusion map
  if ( "dm" %in% m ) x <- destiny::DiffusionMap(d,sigma = dm.sigma, distance = dm.distance)

  # compute dimensional reduction to the set of k dimensions with various methods
  for ( j in m ) dr[[j]] <- list()
  for ( j in k){
    jn <- paste("D",j,sep="")
    if ( "lle" %in% m ) dr[["lle"]][[jn]] <- as.data.frame(lle(d,m=j,k=lle.n)$Y)
    if ( "cmd" %in% m ) dr[["cmd"]][[jn]]  <- as.data.frame(cmdscale(di,k=j))
    if ( "dm" %in% m ) dr[["dm"]][[jn]]   <- as.data.frame(x@eigenvectors[,1:j])
    if ( "tsne" %in% m ) dr[["tsne"]][[jn]] <- as.data.frame(Rtsne(as.dist(di),dims=j,initial_config=cmdscale(di,k=j),perplexity=tsne.perplexity)$Y)
  }

   return(dr)
}

#' @title Computation of a principal curve for a given dimensional reduction representation
#'
#' @description This function computes a principal curve for a given dimensional reduction representation which is specified by component names of an object returned by \code{compdr} using the \pkg{princurve} package.
#' @param y clustering partition. A vector with an integer cluster number for each cell. The order of the cells has to be the same as for the columns of x.
#' @param fb fateBias object returned by the function \code{fateBias}.
#' @param dr list of dimensional reduction representations returned by the function \code{compdr}.
#' @param k integer number for the dimension to be used. This dimension has to be present in \code{dr}. Default value is 2.
#' @param m name of the dimensional reduction algorithms to be used for the principal curve computation. One of \code{lle}, \code{cmd}, \code{dm}, \code{tsne}. Default value is \code{cmd}. Has to be a component of \code{dr}, i.e. previously computed by \code{compdr}.
#' @param trthr real value representing the threshold of the fraction of random forest votes required for the inclusion of a given cell for the computation of the principal curve. If \code{NULL} then only cells with a significant bias >1 are included for each trajectory. The bias is computed as the ratio of the number of votes for a trajectory and the number of votes for the trajectory with the second largest number of votes. By this means only the trajectory with the largest number of votes will receive a bias >1. The siginifcance is computed based on counting statistics on the difference in the number of votes. A significant bias requires a p-value < 0.05. Default value is \code{NULL}.
#' @param start integer number representing a specified starting cluster number for all trajectories, i. e. a common progenitor cluster. The argument is optional. Default value is \code{NULL}.
#' @param  ... additional arguments to be passed to the low level function \code{principal_curve}.
#' @details The function computes a principal curve for each differentiation trajectory by considering only cells that are assigned to the trajectory with a significant fate bias >1 or at least \code{trthr} of the random forest votes, respectively.
#' @details For simulateneous computation and plotting of the principal curve, see function \code{plotFateMap}.
#' @return A list of the following two components:
#'   \item{pr}{A list of principal curve objects produced by the \code{principal_curve} function from the \pkg{princurve} package. Each component corresponds to one differentiation trajectory giving rise to one of the target clusters from the \code{fb} object.}
#'   \item{trc}{A list of ordered cell IDs for each trajectory in \code{pr}.}
#' @examples
#'
#' x <- intestine$x
#' y <- intestine$y
#' tar <- c(6,9,13)
#' fb <- fateBias(x,y,tar,z=NULL,minnr=5,minnrh=10,nbfactor=5,use.dist=FALSE,seed=NULL,nbtree=NULL)
#' dr <- compdr(x,z=NULL,m="cmd",k=2,lle.n=30,dm.sigma=1000,dm.distance="euclidean",tsne.perplexity=30)
#' pr <- prcurve(y,fb,dr,k=2,m="cmd",trthr=0.25,start=NULL)
#'
#' @importFrom stats median
#' @importFrom princurve principal_curve
#' @export
prcurve <- function(y,fb,dr,k=2,m="cmd",trthr=NULL,start=NULL,...){
  pr <- list()
  
  # calculate principal curve for each trajectory
  for ( j in names(fb$votes) ){
    # infer starting curve if starting clusters are given (start)
    if ( ! is.null(start) ){
      if ( !is.null(trthr) ){
        gs <- names(y) %in% rownames(fb$probs)[fb$probs[,j] > trthr]
      }else{
        b <- bias(fb$votes)
        gs <- names(y) %in% rownames(fb$votes)[b$bias[,j] > 1 & b$pv < .05]
      }
      gs <- gs & y %in% c(start,as.numeric(sub("t","",j)))
      st <- principal_curve(as.matrix(dr[[m]][[paste("D",k,sep="")]][gs,]),plot_iterations=FALSE,...)
    }
    # infer final curve
    if ( !is.null(trthr) ){
      g <- names(y) %in% rownames(fb$probs)[fb$probs[,j] > trthr]
    }else{
      b <- bias(fb$votes)
      g <- names(y) %in% rownames(fb$votes)[b$bias[,j] > 1 & b$pv < .05]
    }
    xa <- dr[[m]][[paste("D",k,sep="")]][g,]

    # increase weight of anchor clusters by a factor of ten
    xo <- dr[[m]][[paste("D",k,sep="")]][g,][y[g] == as.numeric(sub("t","",j)),]    
    xv <- list()
    for ( i in 1:ncol(xo) ){
      xv[[i]] <- rep(xo[,i],10)
    }
    xv <- as.matrix(as.data.frame(xv))
    xt <- t(cbind(t(xa),t(xv)))
    f <- rep(FALSE,nrow(xt))
    f[1:nrow(xa)] <- TRUE

    if ( is.null(start) ){
      pr[[j]] <- principal_curve(as.matrix(xt),plot_iterations=FALSE,...)
    }else{
      pr[[j]] <- principal_curve(as.matrix(xt),plot_iterations=FALSE,start=st$s[st$ord,],...)
    }
    pr[[j]]$s <- pr[[j]]$s[f,]
    pr[[j]]$ord <- pr[[j]]$ord[pr[[j]]$ord <= sum(f)]
    rownames(pr[[j]]$s) <- names(y)[g]
  }

  trc <- list()
  for ( j in names(fb$tr) ){
    n <- names(y)[names(y) %in% rownames(pr[[j]]$s)][pr[[j]]$ord]
    trc[[j]] <- n
    if (length(unique(y[n])) > 1){
      if ( median((1:length(n))[y[n] == sub("t","",j)]) < median((1:length(n))[y[n] != sub("t","",j)]) ) trc[[j]] <- rev(trc[[j]])
    }
  }
  return(list(pr=pr,trc=trc))
}

#' @title Inferenence of diffusion pseudotime (DPT) for all cells on a differentiation trajectory
#'
#' @description This function computes a pseudo-temporal ordering of cells reflecting the differentiation projects by using the \code{DPT} function from the \code{destiny} package.
#' @param x expression data frame with genes as rows and cells as columns. Gene IDs should be given as row names and cell IDs should be given as column names. This can be a reduced expression table only including the features (genes) to be used in the analysis.
#' @param y clustering partition. A vector with an integer cluster number for each cell. The order of the cells has to be the same as for the columns of x.
#' @param fb fateBias object returned by the function \code{fateBias}.
#' @param trthr real value representing the threshold of the fraction of random forest votes required for the inclusion of a given cell for the derivation of the diffusion pseudotime. If \code{NULL} then only cells with a signifcant fate bias >1 are included for each trajectory. Default value is \code{NULL}.
#' @param distance parameter for the computation of the underlying diffusion map computed by the function \code{DiffusionMap} from the \pkg{destiny} package. Default value is \code{euclidean}.
#' @param sigma parameter for the computation of the underlying diffusion map computed by the function \code{DiffusionMap} from the \pkg{destiny} package. Default value is 1000.
#' @param  ... additional arguments to be passed to the low level function \code{DiffusionMap}.
#' @details The function orders all cells assigned to a differentiation trajectory with a significant fate bias >1 or a probability greater \code{trthr} for a trajectory, respectively, by diffusion pseudotime.
#' @return trc A list of ordered cell IDs for each trajectory giving rise to one of the targer clusters in \code{fb}
#' @examples
#'
#' \donttest{
#' x <- intestine$x
#' y <- intestine$y
#' tar <- c(6,9,13)
#' fb <- fateBias(x,y,tar,z=NULL,minnr=5,minnrh=10,nbfactor=5,use.dist=FALSE,seed=NULL,nbtree=NULL)
#' dr <- compdr(x,z=NULL,m="cmd",k=2,lle.n=30,dm.sigma=1000,dm.distance="euclidean",tsne.perplexity=30)
#' trc <- dptTraj(x,y,fb,trthr=.25,distance="euclidean",sigma=1000)
#' }
#' @importFrom stats median
#' @export
dptTraj <- function(x,y,fb,trthr=NULL,distance="euclidean",sigma=1000,...){
  trc <- list()
  for ( j in names(fb$probs) ){
    if ( ! is.null(trthr) ){
      probs <- fb$probs
      n  <- rownames(probs)[probs[,j] > trthr]
    }else{
      votes <- fb$votes
      b <- bias(votes)
      n  <- rownames(votes)[b$bias[,j] > 1 & b$pv < .05]
    }
    dm <- destiny::DiffusionMap(as.matrix(t(x[,n])),distance=distance,sigma=sigma,...)
    root_idx <- destiny::random_root(dm)
    pt <- destiny::DPT(dm, root_idx)
    pto <- pt[root_idx, ]
    
    #b <- pt@branch[, 1]
    #tip_idx <- which(b==1 & !is.na(b) & pt@tips[, 1])
    #pto <- pt[tip_idx, ]
    
    n <- n[order(pto,decreasing=FALSE)]

    #ts <- Transitions(as.matrix(t(x[,n])),distance=distance,sigma=sigma,...)
    #pt <- dpt::dpt(ts, branching = FALSE)
    #n <- n[order(pt$DPT,decreasing=FALSE)]
  
    if ( median((1:length(n))[y[n] == sub("t","",j)]) < median((1:length(n))[y[n] != sub("t","",j)]) ) n <- rev(n)
    trc[[j]] <- n
  }
  return(trc)
}

plot2dmap <-  function(d,x,y,g=NULL,n=NULL,col=NULL,tp=1,logsc=FALSE){
 

  if ( !is.null(g) ){
    if ( is.null(n) ) n <- g[1]
    l <- apply(x[g,],2,sum)
    if (logsc) {
      f <- l == 0
      l <- log2(l + .01)
      l[f] <- NA
    }
    mi <- min(l,na.rm=TRUE)
    ma <- max(l,na.rm=TRUE)
    ColorRamp <- colorRampPalette(rev(brewer.pal(n = 7,name = "RdYlBu")))(100)
    ColorLevels <- seq(mi, ma, length=length(ColorRamp))
    v <- round((l - mi)/(ma - mi)*99 + 1,0)
    pardefault <- par()
    layout(matrix(data=c(1,3,2,4), nrow=2, ncol=2), widths=c(5,1,5,1), heights=c(5,1,1,1))
    par(mar = c(3,5,2.5,2))
    plot(d,xlab="",ylab="",main=n,pch=20,cex=0,col="grey",axes=FALSE)
    kk <- order(v,decreasing=FALSE)
    points(d[kk,1],d[kk,2],col=ColorRamp[v[kk]],pch=20,cex=1.5)
    par(mar = c(10,2.5,2.5,4))
    image(1, ColorLevels,
          matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
          col=ColorRamp,
          xlab="",ylab="",
          xaxt="n")
    layout(1)
    par(mar=pardefault$mar)
  }else{
    set.seed(111111)
    if ( is.null(col) ) col <- sample(rainbow(max(y)))
    plot(d,xlab="",ylab="",pch=20,cex=1.5,col=adjustcolor("lightgrey",tp),axes=FALSE)
    for ( i in 1:max(y) ){
      if ( sum(y == i) > 0 ) text(d[y == i,1],d[y == i,2],i,col=adjustcolor(col[i],tp),cex=.75,font=4)
    }
  }
}

plot3dmap <- function(d,x,y,g=NULL,n=NULL,col=NULL,tp=1,logsc=FALSE){
  plot3d(d[,1], d[,2], d[,3], xlab = "", ylab = "", zlab = "", alpha = tp, col = "grey", pch="16", type="p", size = 8, point_antialias = TRUE)
  if ( is.null(g) ){
    set.seed(111111)
    if ( is.null(col) ) col <- sample(rainbow(max(y)))
    for ( i in sort(unique(y)) ){ f <- y == i; text3d(d[f,1], d[f,2], d[f,3], rep(i,sum(f)), font=10, alpha=tp, size=9, depth_test = "always", color=col[i])}
  }
  if ( !is.null(g) ){
    if ( is.null(n) ) n <- g[1]
    l <- apply(x[g,],2,sum)
    if (logsc) {
      f <- l == 0
      l <- log2(l + .01)
      l[f] <- NA
    }
    mi <- min(l,na.rm=TRUE)
    ma <- max(l,na.rm=TRUE)
    ColorRamp <- colorRampPalette(rev(brewer.pal(n = 7,name = "RdYlBu")))(100)
    ColorLevels <- seq(mi, ma, length=length(ColorRamp))
    v <- round((l - mi)/(ma - mi)*99 + 1,0)
    kk <- order(v,decreasing=FALSE)
    points3d(d[kk,1],d[kk,2],d[kk,3],col=ColorRamp[v[kk]],pch="16",size=8)
    
#    apply(cbind(d[kk,],ColorRamp[v[kk]]),1,function(x) points3d(x[1],x[2],x[3],col=x[4],pch="16",size=8))
    image(1, ColorLevels,
          matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
          col=ColorRamp,
          xlab="",ylab="",
          xaxt="n")
  }
}

#' @title Plot dimensional reduction representation of the expression data
#'
#' @description This function plots a dimensional reduction representation using the output of the \code{compdr} function as input. It allows display of the input clusters as well as color coding of fate bias probabilities and gene expression.
#' @param y clustering partition. A vector with an integer cluster number for each cell. The order of the cells has to be the same as for the columns of x.
#' @param dr list of dimensional reduction representations returned by the function \code{compdr}.
#' @param x expression data frame with genes as rows and cells as columns. Gene IDs should be given as row names and cell IDs should be given as column names. This can be a reduced expression table only including the features (genes) to be used in the analysis. This input has to be provided if \code{g} (see below) is given and corresponds to a valid gene ID, i. e. one of the rownames of \code{x}. The default value is \code{NULL}. In this case, cluster identities are highlighted in the plot. 
#' @param g either the name of one of the trajectories from \code{fb} or a gene ID corresponding to one of the rownames of \code{x}. In the latter case, the input argument \code{x} needs to be provided. A vector of gene IDs can also be provided. In this case, the aggregated expression across all gene IDs is plotted. If \code{g} equals E, then the entropy of fate bias is displayed. The default value is \code{NULL}. In this case, cluster identities are highlighted in the plot.
#' @param n optional character string. This argument corresponds to a title for 2-dimensional plots. Default value is \code{NULL}. If not provided, and \code{g} is given, then \code{n} will equal \code{g} or g[1], respectively, if g is a vector of gene IDs.
#' @param prc logical. If \code{TRUE}, then a principal curve is computed and returned. Default is \code{FALSE}.
#' @param logsc logical. If \code{TRUE}, then gene expression of fate bias probabilities are plotted on a log2 scale. Default value is \code{FALSE}.
#' @param k integer number for the dimension to be used. This dimension has to be present in \code{dr}. Default value is 2.
#' @param m name of the dimensional reduction algorithms to be used for the principal curve computation. One of \code{lle}, \code{cmd}, \code{dm}, \code{tsne}. Default value is \code{cmd}. Has to be a component of \code{dr}, i.e. previously computed by \code{compdr}.
#' @param kr integer vector. If \code{k}>3 then \code{kr} indicates the dimensions to be plotted (either two or three of all possible dimensions). Default value is \code{NULL}. In this case, \code{kr} is given by \code{1:min(k,3)}.
#' @param col optional vector of valid color names for all clusters in \code{y} ordered by increasing cluster number. Default value is \code{NULL}.
#' @param fb fateBias object returned by the function \code{fateBias}. If \code{fb} is provided, then a principal curve is computed and shown in the plot. Default value is \code{NULL}. The curve is only displayed if \code{g} equal \code{NULL}.
#' @param trthr real value representing the threshold of the fraction of random forest votes required for the inclusion of a given cell for the computation of the principal curve. If \code{NULL} then only cells with a significant bias >1 are included for each trajectory. The bias is computed as the ratio of the number of votes for a trajectory and the number of votes for the trajectory with the second largest number of votes. By this means only the trajectory with the largest number of votes will receive a bias >1. The siginifcance is computed based on counting statistics on the difference in the number of votes. A significant bias requires a p-value < 0.05. Default value is \code{NULL}.
#' @param start integer number representing a specified starting cluster number for all trajectories, i. e. a common progenitor cluster. The argument is optional. Default value is \code{NULL}.
#' @param tp Transparency of points in the plot to allow better visibility of the principal curves. Default value is 1, i. e. non-transparent.
#' @param  ... additional arguments to be passed to the low level function \code{principal_curve}.
#' @return If \code{fb} is provided as input argument and \code{prc} equals \code{TRUE} then the output corresponds to the output of \code{prcurve}. Otherwise, only ouput is generated is \code{g} equals E. In this case a vector of fate bias entropies for all cells is given.
#' @examples
#'
#' x <- intestine$x
#' y <- intestine$y
#' # v contains all genes (no feature selection like in x)
#' v <- intestine$v
#' fcol <- intestine$fcol
#' tar <- c(6,9,13)
#' fb <- fateBias(x,y,tar,z=NULL,minnr=5,minnrh=10,nbfactor=5,use.dist=FALSE,seed=NULL,nbtree=NULL)
#' dr <- compdr(x,z=NULL,m="cmd",k=2,lle.n=30,dm.sigma=1000,dm.distance="euclidean",tsne.perplexity=30)
#'
#' # plot principal curves
#' pr <- plotFateMap(y,dr,k=2,prc=TRUE,m="cmd",col=fcol,fb=fb,trthr=0.25,start=NULL,tp=.5)
#'
#' # plot expression of gene Apoa1__chr9
#' plotFateMap(y,dr,x=v,g="Apoa1__chr9",prc=FALSE,k=2,m="cmd",col=intestine$fcol)
#'
#' @importFrom grDevices rainbow colorRampPalette adjustcolor
#' @importFrom graphics layout plot points text image abline axis box legend lines par
#' @importFrom RColorBrewer brewer.pal
#' @importFrom rgl plot3d points3d lines3d text3d
#' @export
plotFateMap <- function(y,dr,x=NULL,g=NULL,n=NULL,prc=FALSE,logsc=FALSE,k=2,m="cmd",kr=NULL,col=NULL,fb=NULL,trthr=NULL,start=NULL,tp=1,...){
  if ( is.null(kr) ) kr <- 1:min(k,3)
  d <- dr[[m]][[paste("D",k,sep="")]][,kr]
  if ( ! is.null(fb) & prc ){
    prcu <- prcurve(y=y,fb=fb,dr=dr,k=k,m=m,trthr=trthr,start=start,...)
    pr <- prcu$pr
  }
  if ( !is.null(g) & !is.null(fb) ){
    probs <- as.data.frame(t(fb$probs))
    if ( g %in% rownames(probs) ){
      x <- probs
    }
    if ( g[1] == "E" ){
      xe <- fb$probs + 1e-3
      x <- as.data.frame(t(-apply(xe*log(xe)/log(3),1,sum)))
      rownames(x) <- "E"
    }
  }
  if ( length(kr) == 2 ){
    plot2dmap(d,x=x,y=y,g=g,n=n,col=col,tp=tp,logsc=logsc)
    if ( is.null(g) & ! is.null(fb) & prc){
      for ( j in 1:length(names(pr)) ){
        lines(pr[[j]]$s[pr[[j]]$ord,kr[1]],pr[[j]]$s[pr[[j]]$ord,kr[2]],lwd=2,lty=j)
      }
      legend("topleft",names(pr),lty=1:length(names(pr)))
    }
  }
  if ( length(kr) == 3 ){
    plot3dmap(d,x=x,y=y,g=g,n=n,col=col,tp=tp,logsc=logsc)
    if ( is.null(g) & ! is.null(fb) & prc ){
      for ( j in 1:length(names(pr)) ){
        lines3d(pr[[j]]$s[pr[[j]]$ord,kr[1]],pr[[j]]$s[pr[[j]]$ord,kr[2]],pr[[j]]$s[pr[[j]]$ord,kr[3]],lwd=2,lty=j)
      }
    }
  }
  
  if ( ! is.null(fb) & prc){
    return(prcu)
  }
  if ( !is.null(g) ){
    if ( g[1] == "E" ){
      return(t(x)[,1])
    }
  }
}

#' @title Comparative plot of the expression levels of two genes
#'
#' @description This function produces a scatter plot of the expression levels of two genes. It allows plotting cells of selected clusters and permits highlighting of the fate bias.
#' @param x expression data frame with genes as rows and cells as columns. Gene IDs should be given as row names and cell IDs should be given as column names. This can be a reduced expression table only including the features (genes) to be used in the analysis.
#' @param y clustering partition. A vector with an integer cluster number for each cell. The order of the cells has to be the same as for the columns of x.
#' @param g1 gene id corresponding to a valid row names of x. Expression of gene \code{g1} versus gene \code{g2} will be plotted.
#' @param g2 gene id corresponding to a valid row names of x. Expression of gene \code{g1} versus gene \code{g2} will be plotted.
#' @param clusters vector of valid cluster ids. Expression is displayed for cells in any of the clusters contained in \code{clusters}. If the argument is not given, cells of all clusters are displayed. Default value is \code{NULL}.
#' @param fb fateBias object returned by the function \code{fateBias}. Default value is \code{NULL}. Only if both \code{tn} and \code{fb} are provided as input, the fate bias will be colour coded.
#' @param tn name of a target cluster, i. e. concatenation of a \code{t} and the number of a target cluster. Has to correspond to a column name of \code{fb$probs}. The default value is \code{NULL}. Only if both \code{tn} and \code{fb} are provided as input, the fate bias will be colour coded.
#' @param col optional vector of valid color names for all clusters in \code{y} ordered by increasing cluster number. Default value is \code{NULL}.
#' @param tp Transparency of points in the plot. Default value is 1, i. e. non-transparent.
#' @param plotnum logical value. If \code{TRUE}, then cluster numbers are displayed on top of the data points. Default value is \code{TRUE}.
#' @return None
#' @examples
#' x <- intestine$x
#' y <- intestine$y
#' v <- intestine$v
#' 
#' tar <- c(6,9,13)
#' fb <- fateBias(x,y,tar,z=NULL,minnr=5,minnrh=10,nbfactor=5,use.dist=FALSE,seed=NULL,nbtree=NULL)
#' gene2gene(v,y,"Muc2__chr7","Apoa1__chr9")
#' gene2gene(v,y,"Muc2__chr7","Apoa1__chr9",fb=fb,tn="t6",plotnum=FALSE)
#'
#' @importFrom grDevices rainbow colorRampPalette adjustcolor
#' @importFrom graphics layout plot points text image
#' @importFrom RColorBrewer brewer.pal
#' @export
gene2gene <- function(x,y,g1,g2,clusters=NULL,fb=NULL,tn=NULL,col=NULL,tp=1,plotnum=TRUE){
  set.seed(111111)
  if ( is.null(col) ) col <- sample(rainbow(max(y)))
  if ( is.null(clusters) ) clusters <- sort(unique(y))
  clust_n <- y[y %in% clusters]
  mnd <- cbind( data.frame(clust_n,color=col[clust_n]), t(x[c(g1,g2),y %in% clusters]) )
  if ( !is.null(fb) & !is.null(tn) ){
    l <- fb$probs[rownames(mnd),tn]
    mi <- min(l,na.rm=TRUE)
    ma <- max(l,na.rm=TRUE)
    ColorRamp <- colorRampPalette(rev(brewer.pal(n = 7,name = "RdYlBu")))(100)
    ColorLevels <- seq(mi, ma, length=length(ColorRamp))
    v <- round((l - mi)/(ma - mi)*99 + 1,0)
    pardefault <- par()
    layout(matrix(data=c(1,3,2,4), nrow=2, ncol=2), widths=c(5,1,5,1), heights=c(5,1,1,1))
    #par(mar = c(3,5,2.5,2))
    plot(mnd[,g1],mnd[,g2],col="grey",pch=20,cex=0, xlab=g1, ylab=g2, xlim=c(0,max(mnd[,g1])*1.1),ylim=c(0,max(mnd[,g2])*1.1),bty="n",main=paste("fate bias:",tn,sep=" "))
    kk <- order(v,decreasing=FALSE)
    points(mnd[kk,g1],mnd[kk,g2],col=adjustcolor(ColorRamp[v[kk]],tp),pch=20,cex=1.5)

    #for ( k in kk ){
    #  points(mnd[k,g1],mnd[k,g2],col=adjustcolor(ColorRamp[v[k]],tp),pch=20,cex=1.5)
    #}
    if (plotnum == TRUE){
      text(mnd[,g1],mnd[,g2],mnd$clust_n,col=adjustcolor("black",tp),cex=0.75,font=4)
    }
 
    #par(mar = c(3,2.5,2.5,2))
    par(mar = c(10,2.5,2.5,4))
    image(1, ColorLevels,
          matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
          col=ColorRamp,
          xlab="",ylab="",
          xaxt="n")
    layout(1)
    par(mar=pardefault$mar)
  }else{
    plot(mnd[,g1],mnd[,g2],col=adjustcolor("lightgrey",tp),pch=20,cex=1.5, xlab=g1, ylab=g2, xlim=c(0,max(mnd[,g1])*1.1),ylim=c(0,max(mnd[,g2])*1.1),bty="n")
    if (plotnum == TRUE){
      text(mnd[,g1],mnd[,g2],mnd$clust_n,col=adjustcolor(as.vector(mnd$color),tp),cex=.75,font=4)
    }
  }
}

#' @title Function for differential expression analysis
#'
#' @description This function performs differential expression analysis between two sets of single cell transcriptomes. The inference is based on a noise model or relies on the \code{DESeq2} approach.
#' @param x expression data frame with genes as rows and cells as columns. Gene IDs should be given as row names and cell IDs should be given as column names. This can be a reduced expression table only including the features (genes) to be used in the analysis. This input has to be provided if \code{g} (see below) is given and corresponds to a valid gene ID, i. e. one of the rownames of \code{x}. The default value is \code{NULL}. In this case, cluster identities are highlighted in the plot.
#' @param A vector of cell IDs corresponding column names of \code{x}. Differential expression in set \code{A} versus set \code{B} will be evaluated.
#' @param B vector of cell IDs corresponding column names of \code{x}. Differential expression in set \code{A} versus set \code{B} will be evaluated.
#' @param DESeq logical value. If \code{TRUE}, then \pkg{DESeq2} is used for the inference of differentially expressed genes. In this case, it is recommended to provide non-normalized input data \code{x}. Default value is \code{FALSE}
#' @param method either "per-condition" or "pooled". If DESeq is not used, this parameter determines, if the noise model is fitted for each set separately ("per-condition") or for the pooled set comprising all cells in \code{A} and \code{B}. Default value is "pooled".
#' @param norm logical value. If \code{TRUE} then the total transcript count in each cell is normalized to the minimum number of transcripts across all cells in set \code{A} and \code{B}. Default value is \code{FALSE}.
#' @param vfit function describing the background noise model. Inference of differentially expressed genes can be performed with a user-specified noise model describing the expression variance as a function of the mean expression. Default value is \code{NULL}.
#' @param locreg logical value. If \code{FALSE} then regression of a second order polynomial is perfomed to determine the relation of variance and mean. If \code{TRUE} a local regression is performed instead. Default value is \code{FALSE}.
#' @param  ... additional arguments to be passed to the low level function \code{DESeqDataSetFromMatrix}.
#' @return If \code{DESeq} equals \code{TRUE}, the function returns the output of \pkg{DESeq2}. In this case list of the following two components is returned:
#' \item{cds}{object returned by the \pkg{DESeq2} function \code{DESeqDataSetFromMatrix}.}
#' \item{res}{data frame containing the results of the \pkg{DESeq2} analysis.}
#'Otherwise, a list of three components is returned:
#' \item{vf1}{a data frame of three columns, indicating the mean \code{m}, the variance \code{v} and the fitted variance \code{vm} for set \code{A}.}
#' \item{vf2}{a data frame of three columns, indicating the mean \code{m}, the variance \code{v} and the fitted variance \code{vm} for set \code{B}.}
#' \item{res}{a data frame with the results of the differential gene expression analysis with the structure of the \code{DESeq} output, displaying mean expression of the two sets, fold change and log2 fold change between the two sets, the p-value for differential expression (\code{pval}) and the Benjamini-Hochberg corrected false discovery rate (\code{padj}).} 
#' @examples
#'
#' x <- intestine$x
#' y <- intestine$y
#' v <- intestine$v
#' 
#' tar <- c(6,9,13)
#' fb <- fateBias(x,y,tar,z=NULL,minnr=5,minnrh=10,nbfactor=5,use.dist=FALSE,seed=NULL,nbtree=NULL)
#' 
#' thr <- .3
#'
#' A <- rownames(fb$probs)[fb$probs[,"t6"]  > .3]
#' B <- rownames(fb$probs)[fb$probs[,"t13"] > .3]
#' de <- diffexpnb(v,A=A,B=B)
#'
#' @importFrom stats var approxfun fitted lm coef dnbinom p.adjust
#' @importFrom locfit locfit
#' @export
diffexpnb <- function(x,A,B,DESeq=FALSE,method="pooled",norm=FALSE,vfit=NULL,locreg=FALSE,...){
  if ( ! method %in% c("per-condition","pooled") ) stop("invalid method: choose pooled or per-condition")
  x <- x[,c(A,B)]
  if ( DESeq ){
    # run on sc@expdata
    des <- data.frame( row.names = colnames(x), condition = factor(c( rep(1,length(A)), rep(2,length(B)) )), libType = rep("single-end", dim(x)[2]))
    cds <- DESeq2::DESeqDataSetFromMatrix(countData=round(x,0),colData=des,design =~ condition,...) 
    cds <- DESeq2::DESeq(cds,fitType='local')
    res <- DESeq2::results(cds)
    list(des=des,cds=cds,res=res)
  }else{
    if (norm) x <- as.data.frame( t(t(x)/apply(x,2,sum))*min(apply(x,2,sum,na.rm=TRUE)) )
    fit <- list()
    m   <- list()
    v   <- list()
    for ( i in 1:2 ){
      g <- if ( i == 1 ) A else B
      m[[i]] <- if ( length(g) > 1 ) apply(x[,g],1,mean) else x[,g]
      v[[i]] <- if ( length(g) > 1 ) apply(x[,g],1,var)  else apply(x,1,var)
      if ( method == "pooled"){
        mg <- apply(x,1,mean)
        vg <- apply(x,1,var)
        vl <- log2(vg)
        ml <- log2(mg)
      }else{
        vl <- log2(v[[i]])
        ml <- log2(m[[i]])
      }

      if ( locreg ){
        f <- order(ml,decreasing=FALSE)
        u <- 2**ml[f]
        y <- 2**vl[f]
        lf <- locfit(y~lp(u,nn=.7),family="gamma",maxk=500)
        fit[[i]] <- approxfun(u, fitted(lf), method = "const")
      }else{
        if ( is.null(vfit) ){
          f <- ml > -Inf & vl > -Inf
          ml <- ml[f]
          vl <- vl[f]
          mm <- -8
          repeat{
            fit[[i]] <- lm(vl ~ ml + I(ml^2)) 
            if( coef(fit[[i]])[3] >= 0 | mm >= -1){
              break
            }
            mm <- mm + .5
            f <- ml > mm
            ml <- ml[f]
            vl <- vl[f]
          }
        }else{
          fit[[i]] <- vfit
        }
      }
    }

    if ( locreg ){
      vf  <- function(x,i) fit[[i]](x)
    }else{
      vf  <- function(x,i) 2**(coef(fit[[i]])[1] + log2(x)*coef(fit[[i]])[2] + coef(fit[[i]])[3] * log2(x)**2)
    }
    sf  <- function(x,i) x**2/(max(x + 1e-6,vf(x,i)) - x)

    #psp <- 1e-99
    #pv <- apply(data.frame(m[[1]],m[[2]]),1,function(x){ p12 <- (dnbinom(0:round(x[1]*length(A) + x[2]*length(B),0),mu=mean(x)*length(A),size=length(A)*sf(mean(x),1)) + psp)*(dnbinom(round(x[1]*length(A) + x[2]*length(B),0):0,mu=mean(x)*length(B),size=length(B)*sf(mean(x),2)) + psp); sum(p12[p12 <= p12[round(x[1]*length(A),0) + 1]])/sum(p12)} )
    pv <- apply(data.frame(m[[1]],m[[2]]),1,function(x){ p12 <- (dnbinom(0:round(x[1]*length(A) + x[2]*length(B),0),mu=mean(x)*length(A),size=length(A)*sf(mean(x),1)))*(dnbinom(round(x[1]*length(A) + x[2]*length(B),0):0,mu=mean(x)*length(B),size=length(B)*sf(mean(x),2))); if ( sum(p12) == 0 ) 0 else sum(p12[p12 <= p12[round(x[1]*length(A),0) + 1]])/(sum(p12))} )
    
    res <- data.frame(baseMean=(m[[1]] + m[[2]])/2,baseMeanA=m[[1]],baseMeanB=m[[2]],foldChange=m[[2]]/m[[1]],log2FoldChange=log2(m[[2]]/m[[1]]),pval=pv,padj=p.adjust(pv,method="BH"))
    vf1 <- data.frame(m=m[[1]],v=v[[1]],vm=vf(m[[1]],1))
    vf2 <- data.frame(m=m[[2]],v=v[[2]],vm=vf(m[[2]],2))
    rownames(res) <- rownames(x)
    rownames(vf1) <- rownames(x)
    rownames(vf2) <- rownames(x)
    list(vf1=data.frame(m=m[[1]],v=v[[1]],vm=vf(m[[1]],1)),vf2=data.frame(m=m[[2]],v=v[[2]],vm=vf(m[[2]],2)),res=res)
  }
}

#' @title Function for plotting differentially expressed genes
#'
#' @description This is a plotting function for visualizing the output of the \code{diffexpnb} function.
#' @param x output of the function \code{diffexpnb}.
#' @param pthr real number between 0 and 1. This number represents the p-value cutoff applied for displaying differentially expressed genes. Default value is 0.05. The parameter \code{padj} (see below) determines if this cutoff is applied to the uncorrected p-value or to the Benjamini-Hochberg corrected false discovery rate.
#' @param padj logical value. If \code{TRUE}, then genes with a Benjamini-Hochberg corrected false discovery rate lower than \code{pthr} are displayed. If \code{FALSE}, then genes with a p-value lower than \code{pthr} are displayed.
#' @param lthr real number between 0 and Inf. Differentially expressed genes are displayed only for log2 fold-changes greater than \code{lthr}. Default value is 0.
#' @param mthr real number between -Inf and Inf. Differentially expressed genes are displayed only for log2 mean expression greater than \code{mthr}. Default value is -Inf.
#' @param Aname name of expression set \code{A}, which was used as input to \code{diffexpnb}. If provided, this name is used in the axis labels. Default value is \code{NULL}.
#' @param Bname name of expression set \code{B}, which was used as input to \code{diffexpnb}. If provided, this name is used in the axis labels. Default value is \code{NULL}.
#' @param show_names logical value. If \code{TRUE} then gene names displayed for differentially expressed genes. Default value is \code{FALSE}.
#' @return None
#' @examples
#'
#' x <- intestine$x
#' y <- intestine$y
#' v <- intestine$v
#' 
#' tar <- c(6,9,13)
#' fb <- fateBias(x,y,tar,z=NULL,minnr=5,minnrh=10,nbfactor=5,use.dist=FALSE,seed=NULL,nbtree=NULL)
#' 
#' thr <- .3
#'
#' A <- rownames(fb$probs)[fb$probs[,"t6"]  > .3]
#' B <- rownames(fb$probs)[fb$probs[,"t13"] > .3]
#' de <- diffexpnb(v,A=A,B=B)
#' plotdiffgenesnb(de,pthr=.05)
#'
#' @importFrom grDevices rainbow colorRampPalette adjustcolor
#' @export
plotdiffgenesnb <- function(x,pthr=.05,padj=TRUE,lthr=0,mthr=-Inf,Aname=NULL,Bname=NULL,show_names=TRUE){
  y <- as.data.frame(x$res)
  if ( is.null(Aname) ) Aname <- "baseMeanA"
  if ( is.null(Bname) ) Bname <- "baseMeanB"

  plot(log2(y$baseMean),y$log2FoldChange,pch=20,xlab=paste("log2 ( ( #mRNA[",Aname,"] + #mRNA[",Bname,"] )/2 )",sep=""),ylab=paste("log2 #mRNA[",Bname,"] - log2 #mRNA[",Aname,"]",sep=""),col="grey")
  abline(0,0)
  if ( ! is.null(pthr) ){
    if ( padj ) f <- y$padj < pthr else f <- y$pval < pthr
    points(log2(y$baseMean)[f],y$log2FoldChange[f],col="red",pch=20)
  }
  if ( !is.null(lthr) ) f <- f & abs( y$log2FoldChange ) > lthr
  if ( !is.null(mthr) ) f <- f & log2(y$baseMean) > mthr
  if ( show_names ){
    if ( sum(f) > 0 ) text(log2(y$baseMean)[f],y$log2FoldChange[f],rownames(y)[f],cex=.5)
  }
}

#' @title Function filtering expression data
#'
#' @description This function discards lowly expressed genes from the expression data frame.
#' @param x expression data frame with genes as rows and cells as columns. Gene IDs should be given as row names and cell IDs should be given as column names. 
#' @param n ordered vector of cell IDs to be included. Cell IDs need to be column names of \code{x}. If not provided, then all cell IDs are included in arbitray order. Default value is \code{NULL}.
#' @param minexpr  positive real number. This is the minimum expression required for at least \code{minnumber} cells. All genes that do not fulfill this criterion are removed. The default value is 2.
#' @param minnumber positive integer number. This is the minimum number of cells in which a gene needs to be expressed at least at a level of \code{minexpr}. All genes that do not fulfill this criterion are removed. The default value is 1.
#' @return Reduced expression data frame with genes as rows and cells as columns in the same order as in \code{n}.
#' @examples
#' \donttest{
#' x <- intestine$x
#' y <- intestine$y
#' v <- intestine$v
#'
#' tar <- c(6,9,13)
#' fb <- fateBias(x,y,tar,z=NULL,minnr=5,minnrh=10,nbfactor=5,use.dist=FALSE,seed=NULL,nbtree=NULL)
#' trc <- dptTraj(x,y,fb,trthr=.25,distance="euclidean",sigma=1000)
#' n <- trc[["t6"]]
#' fs  <- filterset(v,n,minexpr=2,minnumber=1)
#' }
#' @export
filterset <- function(x,n=NULL,minexpr=2,minnumber=1){
  if ( is.null(n) ) n <- colnames(x)
  x[apply(x[,n] >= minexpr,1,sum) >= minnumber,n]
}

#' @title Topological ordering of pseudo-temporal expression profiles
#'
#' @description This function computes a topological ordering of pseudo-temporal expression profiles of all genes by using 1-dimensional self-organizing maps.
#' @param x expression data frame with genes as rows and cells as columns. Gene IDs should be given as row names and cell IDs should be given as column names. The pseudo-temporal expression profile of each gene is defined by the order of cell IDs, i. e. columns, in \code{x}.
#' @param nb positive integer number. Number of nodes of the self-organizing map. Default value is 1000.
#' @param alpha positive real number. Pseudo-temporal expression profiles are derived by a local regression of expression values across the ordered cells using the function \code{loess} from the package \pkg{stats}. This is the parameter, which controls the degree of smoothing. Larger values return smoother profiles. Default value is 0.5.
#' @return A list of the following three components:
#' \item{som}{a \code{som} object returned by the function \code{som} of the package \pkg{som}}
#' \item{x}{pseudo-temporal expression profiles, i. e. the input expression data frame \code{x} after smoothing by running mean or local regression, respectivey, and normalization. The sum of smoothened gene expression values across all cells is normalized to 1.}
#' \item{zs}{data frame of z-score transformed pseudo-temporal expression profiles.}
#' @examples
#'
#' \donttest{
#' x <- intestine$x
#' y <- intestine$y
#' v <- intestine$v
#' 
#' tar <- c(6,9,13)
#' fb <- fateBias(x,y,tar,z=NULL,minnr=5,minnrh=10,nbfactor=5,use.dist=FALSE,seed=NULL,nbtree=NULL)
#' dr <- compdr(x,z=NULL,m="cmd",k=2,lle.n=30,dm.sigma=1000,dm.distance="euclidean",tsne.perplexity=30)
#' pr <- prcurve(y,fb,dr,k=2,m="cmd",trthr=0.4,start=NULL)
#' n <- pr$trc[["t6"]]
#' fs  <- filterset(v,n,minexpr=2,minnumber=1)
#' s1d <- getsom(fs,nb=1000,alpha=.5)
#' }
#'
#' @importFrom stats var predict loess
#' @importFrom som som
#' @export
getsom <- function(x,nb=1000,alpha=.5){
    n <- colnames(x)
    x <- t(apply(x,1,function(x,alpha){ v <- 1:length(x); predict(loess( x ~ v, span=alpha ))},alpha=alpha))
    x <- t(apply(x,1,function(x){ x[x<0] <- .1; x }))
    x <- x/apply(x,1,sum)
    zs <- ( x - apply(x,1,mean) )/sqrt ( apply(x,1,var) )
    colnames(zs) <- colnames(x) <- n
    return( list(som=som(zs,1,nb),x=x,z=zs) )
}

#' @title Processing of self-organizing maps for pseudo-temporal expression profiles
#'
#' @description This function processes the self-organizing maps produced by the function \code{getsom}.
#' @param s1d output of function \code{getsom}
#' @param corthr correlation threshold, i. e. a real number between 0 and 1. The z-score of the average normalized pseudo-temporal expression profiles within each node of the self-organizing map is computed, and the correlation of these z-scores between neighbouring nodes is computed. If the correlation is greater than \code{corthr}, neighbouring nodes are merged. Default value is 0.85.
#' @param minsom positive integer number. Nodes of the self-organizing map with less than \code{minsom} transcripts are discarded. Default value is 3.
#' @return A list of the following seven components:
#' \item{k}{vector of Pearson's correlation coefficient between node \code{i} and node \code{i+1} of the populated nodes of the self-organizing map.}
#' \item{nodes}{vector with assignment of genes to nodes of the final self-organizing map (after merging). Components are node numbers and component names are gene IDs.}
#' \item{nodes.e}{data frame with average normalized pseudo-temporal expression profile for each node, ordered by node number.}
#' \item{nodes.z}{data frame with z-score transformed average normalized pseudo-temporal expression profile for each node, ordered by node number.}
#' \item{all.e}{data frame with normalized pseudo-temporal expression profile for all genes in the self-organizing map, ordered by node number.}
#' \item{all.z}{data frame with z-score transformed normalized pseudo-temporal expression profile for all genes in the self-organizing map, ordered by node number.}
#' \item{all.b}{data frame with binarized pseudo-temporal expression profile for all genes in the self-organizing map, ordered by node number. Expression is 1 in cells with z-score > 1 and -1 in cells with z-score < -1, and 0 otherwise.}
#' @examples
#'
#' \donttest{
#' x <- intestine$x
#' y <- intestine$y
#' v <- intestine$v
#' 
#' tar <- c(6,9,13)
#' fb <- fateBias(x,y,tar,z=NULL,minnr=5,minnrh=10,nbfactor=5,use.dist=FALSE,seed=NULL,nbtree=NULL)
#' dr <- compdr(x,z=NULL,m="cmd",k=2,lle.n=30,dm.sigma=1000,dm.distance="euclidean",tsne.perplexity=30)
#' pr <- prcurve(y,fb,dr,k=2,m="cmd",trthr=0.4,start=NULL)
#' n <- pr$trc[["t6"]]
#' fs  <- filterset(v,n,minexpr=2,minnumber=1)
#' s1d <- getsom(fs,nb=1000,alpha=.5)
#' ps <- procsom(s1d,corthr=.85,minsom=3)
#' }
#'
#' @importFrom stats var cor aggregate
#' @export
procsom <- function(s1d,corthr=.85,minsom=3){
  f  <- order(s1d$som$visual$y,decreasing=FALSE)
  x  <- s1d$x[f,]
  
  y <- t(apply(x,1,function(x){ z <- (x - mean(x))/sqrt(var(x)); z[z < 1 & z > -1] <- 0; z[z > 1] <- 1; z[z < -1] <- -1; z })) 

  xz  <- (x - apply(x,1,mean))/sqrt(apply(x,1,var))
   
  ma <- aggregate(s1d$x,by=list(x=s1d$som$visual$y),mean)
  ma <- ma[order(ma$x,decreasing=FALSE),]
  ma <- ma[-1]/apply(ma[,-1],1,sum)
  maz <- ( ma - apply(ma,1,mean))/sqrt(apply(ma,1,var))
    


  k <- c(1)
  for ( i in 2:nrow(maz) ) k <- append(k, cor(t(maz[(i-1):i,]),method="spearman")[1,2] )

  kp <- 1
  for ( i in 2:length(k) ){
    if ( k[i] < corthr ){
      kp[i] <- kp[i-1] + 1
    }else{
      kp[i] <- kp[i-1]
    }
  }
  u <- sort(unique(s1d$som$visual$y))
  v <- s1d$som$visual$y
  for ( i in unique(kp) ){
    v[ s1d$som$visual$y %in% u[kp == i] ] <- i
  }

  
   
  ma <- aggregate(s1d$x[f,],by=list(x=v[f]),mean)
  ma <- ma[order(ma$x,decreasing=FALSE),]
  ma <- ma[-1]/apply(ma[,-1],1,sum)
  maz <- ( ma - apply(ma,1,mean))/sqrt(apply(ma,1,var))
     
  q <- v[f]
  names(q) <- rownames(s1d$som$data)[f]
  g <- aggregate(rep(1,length(q)),by=list(x=q),sum)
  h <- q %in% g[g[,2]>minsom,1]
  ah <- sort(unique(q)) %in% g[g[,2]>minsom,1]
  q <- q[h]
  for  ( i in max(q):1){if (sum(q==i)==0) q[q>i] <- q[q>i] - 1 }

  
  ma  <- ma[ah,]
  maz <- maz[ah,]
  x   <- x[h,]
  xz  <- xz[h,]
  y   <- y[h,]

  return(list(k=k,nodes=q,nodes.e=ma,nodes.z=maz,all.e=x,all.z=xz,all.b=y))
}

#' @title Heatmap of expression profiles
#'
#' @description This function allows plotting of normalized or z-score transformed pseudo-temporal expression profiles and permits highlighting of partitioning along the x-axis and the y-axis
#' @param x data frame with input data to show. Columns will be displayed on the x-axis and rows on the y-axis in the order given in \code{x}. For example, columns can correspond to cells in pseudo-temporal order and rows contain gene expression, i. e. rows can represent pseudo-temporal gene expression profiles.
#' @param xpart optional vector with integer values indicating partitioning of the data points along the x-axis. For instance, \code{xpart} can be a cluster assignment of cell IDs. The order of the components has to be the same as for the columns in \code{x}. Default value is \code{NULL}.
#' @param xcol optional vector with valid color names. The number of components has to be equal to the number of different values on \code{xpart}. If provided, these colors are used to highlight partitioning along the x-axis based on \code{xpart}. Default value is \code{NULL}.
#' @param xlab logical value. If \code{TRUE} then the average position is indicated for each partition value along the x-axis. Default value is \code{TRUE}.
#' @param xgrid logical value. If \code{TRUE} then the partitioning along the x-axis is indicated by vertical lines representing the boundaries of all positions with a given value in \code{xpart}.
#' @param ypart optional vector with integer values indicating partitioning of the data points along the y-axis. For instance, \code{ypart} can be the assignment of gene IDs to nodes of a sel-organizing map. The order of the components has to be the same as for the rows in \code{x}. Default value is \code{NULL}.
#' @param ycol optional vector with valid color names. The number of components has to be equal to the number of different values on \code{ypart}. If provided, these colors are used to highlight partitioning along the y-axis based on \code{ypart}. Default value is \code{NULL}.
#' @param ylab logical value. If \code{TRUE} then the average position is indicated for each partition value along the y-axis. Default value is \code{TRUE}.
#' @param ygrid logical value. If \code{TRUE} then the partitioning along the y-axis is indicated by horizontal lines representing the boundaries of all positions with a given value in \code{ypart}. 
#' @return None
#' @examples
#'
#' \donttest{
#' x <- intestine$x
#' y <- intestine$y
#' v <- intestine$v
#' fcol <- intestine$col
#' 
#' tar <- c(6,9,13)
#' fb <- fateBias(x,y,tar,z=NULL,minnr=5,minnrh=10,nbfactor=5,use.dist=FALSE,seed=NULL,nbtree=NULL)
#' dr <- compdr(x,z=NULL,m="cmd",k=2,lle.n=30,dm.sigma=1000,dm.distance="euclidean",tsne.perplexity=30)
#' pr <- prcurve(y,fb,dr,k=2,m="cmd",trthr=0.4,start=NULL)
#' n <- pr$trc[["t6"]]
#' fs  <- filterset(v,n,minexpr=2,minnumber=1)
#' s1d <- getsom(fs,nb=1000,alpha=.5)
#' ps <- procsom(s1d,corthr=.85,minsom=3)
#' plotheatmap(ps$all.e,xpart=y[n],xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)
#' }
#'
#' @importFrom grDevices rainbow colorRampPalette adjustcolor
#' @importFrom graphics layout plot points text image abline axis box legend lines par rect
#' @importFrom RColorBrewer brewer.pal
#' @export
plotheatmap <- function(x,xpart=NULL,xcol=NULL,xlab=TRUE,xgrid=FALSE,ypart=NULL,ycol=NULL,ylab=TRUE,ygrid=FALSE){
  mi  <- min(x,na.rm=TRUE)
  ma  <- max(x,na.rm=TRUE)
  pardefault <- par()
  
  layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(5,1), heights=c(5,1))
  ColorRamp   <- rev(colorRampPalette(brewer.pal(n = 7,name = "RdYlBu"))(100))
  ColorLevels <- seq(mi, ma, length=length(ColorRamp))
  if ( mi == ma ){
    ColorLevels <- seq(0.99*mi, 1.01*ma, length=length(ColorRamp))
  }

  par(mar = c(3,5,2.5,2))
  image(t(as.matrix(x)),col=ColorRamp,axes=FALSE,ylim=c(-.02,1))
  box()
  set.seed(20)
  if (!is.null(xpart)) {
      tmp <- c()
      width <- ( 1/length(xpart) )/2
      k <- (0:(length(xpart) - 1)/(length(xpart) - 1))
      rect(k - width,rep(-.02,length(xpart)),k + width,rep(-.005,length(xpart)),col=xcol[xpart],border=NA)
      for ( u in unique(xpart) ){
          ol <- (0:(length(xpart) - 1)/(length(xpart) - 1))[xpart == u]
          tmp <- append(tmp, mean(ol))
          delta <- 0.5/(length(xpart) - 1)
          if (xgrid & max(ol) < 1) abline(v = max(ol) + delta, col = "grey", lty = 2)
      }
      if (xlab) axis(1, at = tmp, labels = unique(xpart))
  }
  set.seed(20)
  if ( !is.null(ypart) ){
    tmp <- c()
    for ( u in unique(ypart) ){
      ol <- (0:(length(ypart) - 1)/(length(ypart) - 1))[ypart == u]
      if ( !is.null(ycol) ) points(rep(0,length(ol)),ol,col=ycol[u + 1],pch=15,cex=.75)
      tmp <- append(tmp,mean(ol))
      delta <- .5/(length(ypart) - 1)
      if ( ygrid & max(ol) < 1) abline(a=max(ol) + delta,b=0,col="grey",lty=2)
    }
    if ( ylab ) axis(2,at=tmp,labels=unique(ypart))
  }
  
  par(mar = c(10,2,2.5,2))

  image(1, ColorLevels,
        matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
        col=ColorRamp,
        xlab="",ylab="",
        xaxt="n")
  layout(1)
  par(mar=pardefault$mar)
}

#' @title Plotting of pseudo-temporal expression profiles
#'
#' @description This function allows plotting pseudo-temporal expression profiles for single genes or groups of genes.
#' @param x expression data frame with genes as rows and cells as columns. Gene IDs should be given as row names and cell IDs should be given as column names.
#' @param y clustering partition. A vector with an integer cluster number for each cell. The order of the cells has to be the same as for the columns of \code{x}.
#' @param g a gene ID corresponding to one of the rownames of \code{x}. In the latter case, the input argument \code{x} needs to be provided. A vector of gene IDs can also be provided. In this case, the aggregated expression across all gene IDs is plotted.
#' @param n ordered vector of cell IDs to be included. Cell IDs need to be column names of \code{x}.
#' @param col optional vector of valid color names for all clusters in \code{y} ordered by increasing cluster number. Default value is \code{NULL}.
#' @param name optional character string. This argument corresponds to a title for the plot. Default value is \code{NULL}. If not provided, and \code{g} is given, then \code{name} will equal \code{g} or \code{g[1]}, respectively, if \code{g} is a vector of gene IDs.
#' @param cluster logical value. If \code{TRUE} then the partitioning along the x-axis is indicated be vertical lines representing the boundaries of all positions with a given value in \code{y}. The average position across all cells in a cluster will be indicated on the x-axis.
#' @param alpha positive real number. Pseudo-temporal expression profiles are derived by a local regression of expression values across the ordered cells using the function \code{loess} from the package \pkg{stats}. This is the parameter, which controls the degree of smoothing. Larger values return smoother profiles. Default value is 0.5.
#' @param types optional vector with IDs for different subsets of cells in \code{y}, e. g. different batches. All cells with the same ID will be displayed by the same symbol and color. Default value is \code{NULL}
#' @param cex size of data points. Default value is 3.
#' @param ylim vector of two numerical values: lower and upper limit of values shown on the y-axis. Default value is \code{NULL} and the whole range is shown.
#' @param map logical. If \code{TRUE} then data points are shown. Default value is \code{TRUE}. 
#' @param leg logical. If \code{TRUE} then axes and labels are shown. Default value is \code{TRUE}.
#' @return None
#' @examples
#'
#' \donttest{
#' x <- intestine$x
#' y <- intestine$y
#' v <- intestine$v
#' fcol <- intestine$col
#' tar <- c(6,9,13)
#' fb <- fateBias(x,y,tar,z=NULL,minnr=5,minnrh=10,nbfactor=5,use.dist=FALSE,seed=NULL,nbtree=NULL)
#' dr <- compdr(x,z=NULL,m="cmd",k=2,lle.n=30,dm.sigma=1000,dm.distance="euclidean",tsne.perplexity=30)
#' pr <- prcurve(y,fb,dr,k=2,m="cmd",trthr=0.4,start=NULL)
#' n <- pr$trc[["t6"]]
#' fs  <- filterset(v,n,minexpr=2,minnumber=1)
#' s1d <- getsom(fs,nb=1000,alpha=.5)
#' ps <- procsom(s1d,corthr=.85,minsom=3)
#' # plot average profile of all genes of node 1 in the self-organizing map
#' g <- names(ps$nodes)[ps$nodes == 1]
#' plotexpression(v,y,g,n,col=fcol,name="Node 1",cluster=FALSE,alpha=.5,types=NULL)
#' }
#'
#' @importFrom grDevices rainbow colorRampPalette adjustcolor
#' @importFrom graphics layout plot points text image abline axis box legend lines par
#' @importFrom stats loess predict
#' @export
plotexpression <- function (x, y, g, n, col = NULL, name = NULL, cluster = FALSE, alpha = 0.5, types = NULL, cex = 3, ylim = NULL, map = TRUE, leg = TRUE){
    
    cl <- unique(y[n])
    set.seed(111111)
    if (is.null(col)) 
        col <- sample(rainbow(max(y)))
    xlim <- c(1, length(n))
    if (!is.null(types)) 
        xlim[1] <- 1.25 * xlim[1]
    z <- if (length(g) == 1) 
        x[g, n]
    else t(apply(x[g, n], 2, sum))
    if (is.null(name)) 
        name <- g[1]
    if ( leg ){ ylab = "Expression" } else { ylab = NA }
    if ( leg ){ main = name } else { main = NA }
    if ( is.null(ylim) ){
        plot(c(1,length(n)), c(min(z),max(z)), cex = 0, axes = FALSE, xlab = "", 
             ylab = ylab, main = main, xlim = xlim)
    }else{
        plot(c(1,length(n)), c(min(z),max(z)), cex = 0, axes = FALSE, xlab = "", 
             ylab = ylab, main = main, xlim = xlim, ylim = ylim)
    }

    if ( map ){
        if (!is.null(types)) {
            coloc <- rainbow(length(unique(types)))
            syms <- c()
            for (i in 1:length(unique(types))) {
                f <- types == sort(unique(types))[i]
                syms <- append(syms, ((i - 1)%%25) + 1)
                points((1:length(n))[f], t(z)[f], col = coloc[i], 
                       pch = ((i - 1)%%25) + 1, cex = 1)
            }
        }
        else {
                                        #points((1:length(n)), t(z), col = "grey", pch = 20, cex = cex)
            for (i in 1:length(cl)) {
                f <- y[n] == cl[i]
                points((1:length(n))[f], t(z)[f], pch = 20, cex = cex, 
                       col = col[cl[i]])
            }
            
        }
    }
    if ( leg ){
        for (i in 1:length(cl)) {
            f <- y[n] == cl[i]
                                        #if (is.null(types)) {
                                        #    text((1:length(n))[f], t(z)[f], cl[i], font = 4, 
                                        #        col = col[cl[i]])
                                        #}
            zc <- if (i == 1) 
                      sum(f)
                  else append(zc, zc[i - 1] + sum(f))
            xc <- if (i == 1) 
                      sum(f)/2
                  else append(xc, zc[i - 1] + sum(f)/2)
            if (cluster) 
                abline(v = zc[i], col = "grey", lty = 2)
        }
        u <- 1:length(n)
        v <- as.vector(t(z))
        zc <- predict(loess(v ~ u, span = alpha))
        zc[zc < 0] <- 0.1
        lines(u, zc)
        if (!is.null(types)) 
            legend("topleft", legend = sort(unique(types)), col = coloc, 
                   pch = syms)
        axis(2)
       
        if (cluster) 
            axis(1, at = xc, labels = cl)
    }
    if ( !leg ) box(col="white") else box()
    
}

#' @title Extract genes with high importance values for random forest classification
#'
#' @description This function extracts all genes with an importance value for classifying cells into a given target cluster exceeding a given threshold for at least one of the random forest iterationns.
#' @param fb fateBias object returned by the function \code{fateBias}. If \code{fb} is provided, then a principal curve is computed and shown in the plot. Default value is \code{NULL}. The curve is only displayed if \code{g} equal \code{NULL}.
#' @param tn name of a target cluster, i. e. concatenation of a \code{t} and the number of a target cluster. Has to correspond to a column name of \code{fb$probs}. 
#' @param ithr positive real number. Threshold for the required importance measure (mean decrease in accuracy of classification upon removal, see \pkg{randomForest}) to include a gene into the output as important feature for classying cells in \code{tn}. Default value is 0.02.
#' @param zthr positive real number. Threshold for the required z-score of the importance measure (importance divided by the standard deviation of importance) to include a gene into the output as important feature for classying cells in \code{tn}. Default value is 2.
#' @return The function returns a list of two elements.
#' \item{d}{a data frame with mean importance values for all genes surviving the filtering by \code{ithr} and \code{zthr}. Columns correspond to random forest iterations, starting from the initial target cluster.}
#' \item{d}{a data frame with the standard deviation of importance values for all genes surviving the filtering by \code{ithr} and \code{zthr}. Columns correspond to random forest iterations, starting from the initial target cluster.}
#' The function produces a heatmap of \code{d} with hierarchical clustering of the rows using the function \code{pheatmap} from the \pkg{pheatmap} package.
#' @examples
#' x <- intestine$x
#' y <- intestine$y
#' tar <- c(6,9,13)
#' fb <- fateBias(x,y,tar,z=NULL,minnr=5,minnrh=10,nbfactor=5,use.dist=FALSE,seed=NULL,nbtree=NULL)
#' k <- impGenes(fb,"t6",ithr=.02,zthr=2)
#' @importFrom pheatmap pheatmap
#' @export
impGenes <- function(fb,tn,ithr=.02,zthr=2){
  n <- c()
  for ( i in 1:length(fb$rfl) ){
    k <- fb$rfl[[i]]$importance[,sub("t","",tn)]
    l <- fb$rfl[[i]]$importanceSD[,sub("t","",tn)]
    w <- k/(l + 1e-10)
    n <- unique(c(n,names(w)[w>zthr & k>ithr]))
  }
  
  for ( i in 1:length(fb$rfl) ){
    if ( i == 1 ) d  <- data.frame(fb$rfl[[i]]$importance[n,sub("t","",tn)])   else d[,i] <- fb$rfl[[i]]$importance[n,sub("t","",tn)]
    if ( i == 1 ) de <- data.frame(fb$rfl[[i]]$importanceSD[n,sub("t","",tn)]) else de[,i] <- fb$rfl[[i]]$importanceSD[n,sub("t","",tn)]
  }
  names(d) <- names(de) <- 1:ncol(d)
  pheatmap(d,cluster_cols=FALSE)
  return(list(d=d,de=de))
}
