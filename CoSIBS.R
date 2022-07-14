###############################
##
## Project: CoSIBS 2022
##
## Purpose: Useful functions for spectral clustering
##
## Author: Charlie Carpenter
## Email: charles.carpenter@cuanschutz.edu
##
## Date Created: 2022-07-07
##
## ---------------------------
## Notes:
##   
##
## ---------------------------


# Kernels -----------------------------------------------------------------

Gaussian_kernel <- function(Z, rho){
  exp(-(1/rho)*as.matrix(dist(Z, method = "euclidean", upper = T)^2))
}

Zhang_kernel <- function(Z, p){
  d <- as.matrix( dist(Z) )
  
  si <- apply(d, 2, function(s) sort(s)[p+1] )
  db <- si %o% si
  dt <- -d^2/db
  
  exp(dt)
}

Spectrum_kernel <- function(Z, NN=3, NN2=7){
  tz <- t(Z)
  dtz <- as.data.frame(tz)
  Spectrum::CNN_kernel(mat=dtz, NN=NN, NN2=NN2)
}

knnGrf <- function(A, k, mutual=FALSE){
  knn <- function(aa, k){
    aa[order(aa)[(k+2):length(aa)]] <- 0
    aa[aa > 0] <- 1
    aa
  }
  
  akn <- apply(A,2,knn,k=k)
  
  if(mutual){
    for(i in 1:nrow(akn)){
      for(j in 1:ncol(akn)){
        akn[i,j] <- ifelse(akn[i,j]==akn[j,i], 
                           akn[i,j],
                           0)
      }
    }
  } else{
    for(i in 1:nrow(akn)){
      for(j in 1:ncol(akn)){
        akn[i,j] <- ifelse(akn[i,j]!=akn[j,i], 
                           max(akn[i,j],akn[j,i]),
                           akn[i,j])
      }
    }
  }
  akn
}

# Laplacian ---------------------------------------------------------------

## Laplacian 
Laplacian <- function(dat, rho=NULL, kernel='Gaussian',
                      lap.type = 'sym',
                      grf.type = 'full', k=5, p=5,
                      binary.grf = FALSE,
                      epsilon=NULL, mutual=FALSE){
                        
  if(is.null(rho)) rho <- median(dist(dat))
  if(kernel == 'Gaussian') A <- Gaussian_kernel(dat, rho = rho)
  
  if(kernel == 'Zhang') A <- Zhang_kernel(dat, p)
  if(kernel == 'Spectrum'){
    td <- t(dat)
    dtd <- as.data.frame(td)
    
    A <- Spectrum::CNN_kernel(dtd)
  } 
  
  if(kernel == 'Linear') A <- dat%*%t(dat)
  if(kernel == 'Cor') A <- cor(t(dat))
  
  diag(A) <- 0 ## Adjacency matrices need diag of 0
  
  if(grf.type == 'e-graph'){
    A[A < epsilon] <- 0
  }
  
  if(grf.type == 'knn'){
    A <- knnGrf(A, k, mutual=mutual)
  }
  
  if(binary.grf) A[A>0] <- 1
  
  deg <- rowSums(A)
  
  if(lap.type == 'sym'){
    ds <- ifelse(deg>0, 1/sqrt(deg), 0)
    # L = I - D^(-1/2) A D^(-1/2)
    L <- diag(nrow = nrow(A)) - diag(ds) %*% A %*% diag(ds)
  }
  
  ## 'Shifted' Laplacian
  if(lap.type == "shift"){
    ds <- ifelse(deg>0, 1/sqrt(deg), 0)
    # L = I + D^(-1/2) A D^(-1/2)
    L <- diag(nrow = nrow(A)) + diag(ds) %*% A %*% diag(ds)
  }
  
  ## Laplacian used by Ng (2002) 'On Spectral Clustering'
  if(lap.type == "Ng"){
    ds <- ifelse(deg>0, 1/sqrt(deg), 0)
    # L = D^(-1/2) A D^(-1/2)
    L <- diag(ds) %*% A %*% diag(ds)
  }
  
  ## Random walk Laplacian (best one)
  if(lap.type == 'rw'){
    ds <- ifelse(deg>0, 1/deg, 0)
    # L = I - D^-1 A
    L <- diag(nrow = nrow(A)) - diag(ds) %*% A
  } 
  
  L
}

# Flag Mean ---------------------------------------------------------------

## Check if matrix is orthogonal
orthoCheck <- function(x){
  xtx <- crossprod(x)
  I <- diag(nrow = nrow(xtx))
  
  # Same tolerance as 'isSymmetric'
  sum((I-xtx)^2) < 100*.Machine$double.eps
}
  
## Calculate flag mean subspace
## Marrinan, et al. "Finding the Subspace Mean or Median to Fit Your Need"
flagMean <- function(x){
  
  ## Check orthogonal
  # if(!all(sapply(x, orthoCheck))) stop("matrix not orthogonal")
  
  r <- max(sapply(x, ncol))
  
  ## Concatenating matrices
  X <- do.call(cbind, x)
  svd(X)
}



