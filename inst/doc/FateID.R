## ----echo=FALSE---------------------------------------------------------------
knitr::opts_chunk$set(fig.width=8, fig.height=8) 

## ----eval = FALSE-------------------------------------------------------------
#  install.packages("FateID")

## -----------------------------------------------------------------------------
library(FateID)

## -----------------------------------------------------------------------------
data(intestine)

## -----------------------------------------------------------------------------
x <- intestine$x
head(x[,1:5])

## -----------------------------------------------------------------------------
y <- intestine$y
head(y)

## -----------------------------------------------------------------------------
tar <- c(6,9,13)

## -----------------------------------------------------------------------------
FMarker <- list(c("Defa20__chr8","Defa24__chr8"), "Clca3__chr3", "Alpi__chr1")
xf <- getPart(x,FMarker,fthr=NULL,n=5)
head(xf$part)
head(xf$tar)
tar <- xf$tar
y <- xf$part

## -----------------------------------------------------------------------------
rc <- reclassify(x, y, tar, clthr=.75, nbfactor=5, use.dist=FALSE, seed=12345, nbtree=NULL, q=0.9)
y  <- rc$part

## -----------------------------------------------------------------------------
v <- intestine$v
rc <- reclassify(v, y, tar, clthr=.75, nbfactor=5, use.dist=FALSE, seed=12345, nbtree=NULL, q=0.9)
y  <- rc$part

## -----------------------------------------------------------------------------
x  <- rc$xf

## -----------------------------------------------------------------------------
x <- getFeat(v,y,tar,fpv=0.01)

## -----------------------------------------------------------------------------
tar <- c(6,9,13)
x <- intestine$x
y <- intestine$y
fb  <- fateBias(x, y, tar, z=NULL, minnr=5, minnrh=10, adapt=TRUE, confidence=0.75, nbfactor=5, use.dist=FALSE, seed=12345, nbtree=NULL)

## -----------------------------------------------------------------------------
dr  <- compdr(x, z=NULL, m=c("tsne","cmd","umap"), k=2, lle.n=30, tsne.perplexity=30, seed=12345)

## -----------------------------------------------------------------------------
plotFateMap(y,dr,k=2,m="umap")

## ----eval=FALSE---------------------------------------------------------------
#  plotFateMap(y,dr,k=3,m="umap")

## -----------------------------------------------------------------------------
plotFateMap(y,dr,k=2,m="umap",fb=fb,g="t6")

## -----------------------------------------------------------------------------
pr <- plotFateMap(y,dr,k=2,m="umap",trthr=.33,fb=fb,prc=TRUE)

## -----------------------------------------------------------------------------
v <- intestine$v
pr <-plotFateMap(y, dr, k=2, m="umap", g=c("Defa20__chr8", "Defa24__chr8"), n="Defa", x=v)

## -----------------------------------------------------------------------------
E <- plotFateMap(y,dr,k=2,m="umap",g="E",fb=fb)
head(E)

## -----------------------------------------------------------------------------
pr  <- prcurve(y,fb,dr,k=2,m="umap",trthr=0.33,start=3)

## -----------------------------------------------------------------------------
n <- pr$trc[["t6"]]

## -----------------------------------------------------------------------------
v <- intestine$v
fs  <- filterset(v,n=n,minexpr=2,minnumber=1)

## -----------------------------------------------------------------------------
s1d <- getsom(fs,nb=50,alpha=.5)

## -----------------------------------------------------------------------------
ps  <- procsom(s1d,corthr=.85,minsom=3)

## -----------------------------------------------------------------------------
fcol <- sample(rainbow(max(y)))

## -----------------------------------------------------------------------------
plotheatmap(ps$nodes.z, xpart=y[n], xcol=fcol, ypart=unique(ps$nodes), xgrid=FALSE, ygrid=TRUE, xlab=FALSE)

## -----------------------------------------------------------------------------
plotheatmap(ps$all.z, xpart=y[n], xcol=fcol, ypart=ps$nodes, xgrid=FALSE, ygrid=TRUE, xlab=FALSE)

## -----------------------------------------------------------------------------
plotheatmap(ps$all.e, xpart=y[n], xcol=fcol, ypart=ps$nodes, xgrid=FALSE, ygrid=TRUE, xlab=FALSE)

## -----------------------------------------------------------------------------
plotheatmap(ps$all.b, xpart=y[n], xcol=fcol, ypart=ps$nodes, xgrid=FALSE, ygrid=TRUE, xlab=FALSE)

## -----------------------------------------------------------------------------
g <- names(ps$nodes)[ps$nodes == 1]

## -----------------------------------------------------------------------------
plotexpression(fs, y, g, n, col=fcol, name="Node 1", cluster=FALSE, alpha=.5, types=NULL)

## -----------------------------------------------------------------------------
plotexpression(fs, y, "Clca4__chr3", n, col=fcol, cluster=FALSE, alpha=.5, types=NULL)

## -----------------------------------------------------------------------------
plotexpression(fs, y, g, n, col=fcol, name="Node 1", cluster=FALSE, types=sub("\\_\\d+","",n))

## -----------------------------------------------------------------------------
group <- head(g,6)
plotexpressionProfile(fs, y, group, n, name="Node 1", cluster=FALSE)

## -----------------------------------------------------------------------------
thr <- .5
a   <- "t13"
b   <- "t6"
cl  <- c(3,4,5)
A <- rownames(fb$probs)[fb$probs[,a] > thr]
A <- A[y[A] %in% cl]
B <- rownames(fb$probs)[fb$probs[,b] > thr]
B <- B[y[B] %in% cl]
de <- diffexpnb(v,A=A,B=B,DESeq=FALSE,norm=FALSE,vfit=NULL,locreg=FALSE)

## -----------------------------------------------------------------------------
plotdiffgenesnb(de,mthr=-4,lthr=0,Aname=a,Bname=b,padj=FALSE)

## -----------------------------------------------------------------------------
gene2gene(intestine$v,intestine$y,"Muc2__chr7","Apoa1__chr9")

## -----------------------------------------------------------------------------
gene2gene(intestine$v, intestine$y, "Muc2__chr7", "Apoa1__chr9", fb=fb, tn="t6", plotnum=FALSE)

## -----------------------------------------------------------------------------
k <- impGenes(fb,"t6")

