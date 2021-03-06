% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FateID_functions.R
\name{plotdiffgenesnb}
\alias{plotdiffgenesnb}
\title{Function for plotting differentially expressed genes}
\usage{
plotdiffgenesnb(
  x,
  pthr = 0.05,
  padj = TRUE,
  lthr = 0,
  mthr = -Inf,
  Aname = NULL,
  Bname = NULL,
  show_names = TRUE,
  ...
)
}
\arguments{
\item{x}{output of the function \code{diffexpnb}.}

\item{pthr}{real number between 0 and 1. This number represents the p-value cutoff applied for displaying differentially expressed genes. Default value is 0.05. The parameter \code{padj} (see below) determines if this cutoff is applied to the uncorrected p-value or to the Benjamini-Hochberg corrected false discovery rate.}

\item{padj}{logical value. If \code{TRUE}, then genes with a Benjamini-Hochberg corrected false discovery rate lower than \code{pthr} are displayed. If \code{FALSE}, then genes with a p-value lower than \code{pthr} are displayed.}

\item{lthr}{real number between 0 and Inf. Differentially expressed genes are displayed only for log2 fold-changes greater than \code{lthr}. Default value is 0.}

\item{mthr}{real number between -Inf and Inf. Differentially expressed genes are displayed only for log2 mean expression greater than \code{mthr}. Default value is -Inf.}

\item{Aname}{name of expression set \code{A}, which was used as input to \code{diffexpnb}. If provided, this name is used in the axis labels. Default value is \code{NULL}.}

\item{Bname}{name of expression set \code{B}, which was used as input to \code{diffexpnb}. If provided, this name is used in the axis labels. Default value is \code{NULL}.}

\item{show_names}{logical value. If \code{TRUE} then gene names displayed for differentially expressed genes. Default value is \code{FALSE}.}

\item{...}{Additional arguments for function \code{plot}.}
}
\value{
None
}
\description{
This is a plotting function for visualizing the output of the \code{diffexpnb} function as MA plot.
}
\examples{

x <- intestine$x
y <- intestine$y
v <- intestine$v

tar <- c(6,9,13)
fb <- fateBias(x,y,tar,z=NULL,minnr=5,minnrh=10,nbfactor=5,use.dist=FALSE,seed=NULL,nbtree=NULL)

thr <- .3

A <- rownames(fb$probs)[fb$probs[,"t6"]  > .3]
B <- rownames(fb$probs)[fb$probs[,"t13"] > .3]
de <- diffexpnb(v,A=A,B=B)
plotdiffgenesnb(de,pthr=.05)

}
