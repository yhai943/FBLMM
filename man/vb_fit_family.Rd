% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vb_fit_family.R
\name{vb_fit_family}
\alias{vb_fit_family}
\title{Fit Variational Bayes model for family data}
\usage{
vb_fit_family(
  y,
  genotype0,
  max_iter = 10000,
  weight.type = NULL,
  maf.filter.snp = 0.01,
  epsilon_conv = 1e-04,
  Bsigma2beta = 1,
  theta_beta = 0.1,
  theta_u = 0.1,
  verbose = TRUE,
  kernel = "Lin",
  kin.mat = kin.mat,
  dummy.fam = dummy.fam
)
}
\arguments{
\item{y}{a phenotype vector of length n}

\item{genotype0}{a list of genotype matrices of the target sequence}

\item{max_iter}{maximum number of iteration}

\item{weight.type}{type of weight function}

\item{maf.filter.snp}{a filtering threshold of minor allele frequency for the isolated predictors}

\item{epsilon_conv}{a numeric value giving the interval endpoint}

\item{Bsigma2beta}{a numeric value for sigma beta}

\item{theta_beta}{probability of causal variants}

\item{theta_u}{probability of causal regions}

\item{verbose}{informative iteration messages}

\item{kernel}{kernel type for covariance matrix}

\item{kin.mat}{kinship covariance matrix}

\item{dummy.fam}{family indicator matrix}
}
\description{
Take in genotype, phenotype and analyses the target sequence by using BLMM.
}
\details{
A hybrid model that includes a sparsity regression model and a LMM with multiple random effects.
The sparsity regression model part is designed to capture the strong predictive effects from isolated
predictors, whereas the LMM part is to capture the effects from a group of predictors located in
nearby genetic regions.
}
\examples{
datafile = system.file("extdata",package = "BLMM")
data <- read_file(gene_file = paste(datafile,"snp.bed",sep="/"),y_file = paste(datafile,"y.txt",sep="/"),kin.mat=T)
y = data$y
gene = data$genematrix
# Extract gene intervals from BED files
genotype0 <- list(as.matrix(gene[,1:200]), as.matrix(gene[,300:450]))
dummy.fam = data$dummy.fam
kin.mat = data$kin.mat
# choose model type: 1. "uw" for BLMM-UW; 2. "beta" for BLMM-BETA; 3. "wss" for BLMM-WSS;
fit <- vb_fit_family(y = y, genotype0 = genotype0,dummy.fam = dummy.fam,kin.mat=kin.mat)
}
