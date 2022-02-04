#' Preliminary processing of data
#'
#' @param gene_file File of bed to be imported
#' @param y phenotype file (need to include NA)
#' @param kin.mat whether to compute the kinship covariance matrix
#' 
#' @description
#' Normalize the data before BLMM analysis
#' @details
#' Write the genetic file as a matrix.
#' Generate a family indicator matrix. Compute the kinship covariance matrix
#' Phenotypic data needs to contain missing values
#' @examples
#' datafile = system.file("extdata",package = "BLMM")
#' data <- read_file(gene_file = paste(datafile,"snp.bed",sep="/"),y_file = paste(datafile,"y.txt",sep="/"),kin.mat=T)
#' y = data$y
#' gene = data$genematrix
#' # Extract gene intervals from BED files
#' genotype0 <- list(as.matrix(gene[,1:200]), as.matrix(gene[,300:450]))
#' dummy.fam = data$dummy.fam
#' kin.mat = data$kin.mat
#'
#' @importFrom BEDMatrix BEDMatrix
#' @export
read_file <- function(gene_file,y_file,kin.mat=F){
  #read file
  gene <- BEDMatrix(gene_file)
  #Y
  y <- read.table(y_file, header = F)[,3]
  #dummy.fam
  name <- read.table(y_file,header = F)[,1]
  sample.id.all <- as.character(name)
  n <- length(sample.id.all)
  famid.all <- sapply(1:n,function(i) strsplit(sample.id.all, split = "_")[[i]][1])
  n.fam <- length(unique(famid.all))
  famid.name <- unique(famid.all)
  dummy.fam <- matrix(0, nrow = n, ncol = n.fam)
  
  # give names
  rownames(dummy.fam) <- sample.id.all
  colnames(dummy.fam) <- famid.name
  
  # identify the family id (indicator)
  ind.fam <- match(famid.all,famid.name)
  for (i in 1:n) {
    dummy.fam[i, ind.fam[i]] <- 1
  }
 if (kin.mat==F){
  return(list(genematrix=gene,dummy.fam=dummy.fam,y=y))
  }else{
  M = (as.matrix(gene)-1)
  p_lower <- (apply(M,2,sum)+nrow(M))/(nrow(M)*2)
  p_upper <- 2*(p_lower-0.5)
  P <- matrix(p_upper,byrow=T,nrow=nrow(M),ncol=ncol(M))
  Z = M-P
  d <- 2*sum(p_lower*(1 - p_lower))
  G <- Z %*% t(Z) / d
  return(list(genematrix=gene,dummy.fam=dummy.fam,y=y,kin.mat=G))
  }
}