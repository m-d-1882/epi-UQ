# © Copyright University of Exeter. This code is open source and released under the MIT Licence

#' Formatting data
#'
#' Formats data so that it is in the correct form for use in other functions, and calculates the (weighted) SVD basis of the ensemble
#'
#' @param data a matrix containing individual fields in the columns (i.e. the matrix has dimension lxn)
#' @param weightinv the inverse of lxl positive definite weight matrix W. If NULL, the identity matrix is used
#' @param RemoveMean if TRUE, centres the data prior to calculating the basis
#' @param StoreEigen if TRUE, stores Q, lambda from eigendecomposition of W (in order to make later calculations more efficient)
#'
#' @return \item{tBasis}{The (weighted) SVD basis of the centred ensemble if RemoveMean = TRUE, of the original data otherwise}
#' \item{CentredField}{The centred data if RemoveMean = TRUE, the original data otherwise.}
#' \item{EnsembleMean}{The mean across the columns of the data. A zero vector if RemoveMean = FALSE}
#' \item{}
#'
#' @export
MakeDataBasis <- function(data, weightinv = NULL, W = NULL, RemoveMean = TRUE, StoreEigen = TRUE){
  if (RemoveMean == TRUE){
    EnsembleMean <- apply(data, 1, mean)
    CentredField <- 0*data
    for (i in 1:dim(data)[2]){
      CentredField[,i] <- data[,i] - EnsembleMean
    }
  }
  else {
    EnsembleMean <- c(rep(0, dim(data)[1]))
    CentredField <- data
  }
  #if (is.null(weightinv)){
  #  weightinv <- diag(dim(data)[1])
  #}
  if (is.null(W)){
    tSVD <- wsvd(t(CentredField), weightinv = weightinv)
    tBasis <- tSVD$v
    if (StoreEigen == TRUE){
      Q <- tSVD$Q
      Lambda <- tSVD$Lambda
      return(list(tBasis = tBasis, CentredField = CentredField, EnsembleMean = EnsembleMean, Q = Q, Lambda = Lambda))
    }
    else {
      return(list(tBasis = tBasis, CentredField = CentredField, EnsembleMean = EnsembleMean))
    }
  }
  else if (!is.null(W) & is.null(weightinv)){
    eig <- eigen(W)
    Q <- eig$vectors
    Lambda <- 1 / eig$values
    Winv <- Q %*% diag(Lambda) %*% t(Q)
    attr(Winv, 'diagonal') <- FALSE
    attr(Winv, 'identity') <- FALSE
    tSVD <- wsvd(t(CentredField), weightinv = Winv, Q = Q, Lambda = Lambda)
    tBasis <- tSVD$v
    if (StoreEigen == TRUE){
      return(list(tBasis = tBasis, CentredField = CentredField, EnsembleMean = EnsembleMean, Q = Q, Lambda = Lambda, Winv = Winv))
    }
    else {
      return(list(tBasis = tBasis, CentredField = CentredField, EnsembleMean = EnsembleMean, Winv = Winv))
    }
  }
}


#' Weighted singular value decomposition
#'
#' Calculates the SVD basis across the output, given the inverse of W.
#'
#' @param data n x l matrix to calculate basis from (i.e. rows are output fields).
#' @param weightinv l x l inverse of W. If NULL, calculates standard SVD.
#' @param Q l x l matrix from eigen decomposition of W^{-1}, if already have this then speeds up calculation of basis
#' @param Lambda vector from eigen decomposition of W^{-1}, if already have this then speeds up calculation of basis
#'
#' @return The weighted SVD of the data.
#'
wsvd <- function(data, weightinv = NULL, Q = NULL, Lambda = NULL){
  if (is.null(weightinv)){
    svd_output <- svd(data)
  }
  else {
    stopifnot(dim(data)[2] == dim(weightinv)[1])
    if (is.null(Q) & attributes(weightinv)$diagonal == FALSE){
      eig <- eigen(weightinv)
      Q <- eig$vectors
      Lambda <- eig$values
      data_w <- data %*% Q %*% diag(sqrt(Lambda)) %*% t(Q)
      svd_output <- svd(data_w)
      svd_output$v <- t(t(svd_output$v) %*% Q %*% diag(1 / sqrt(Lambda)) %*% t(Q))
      svd_output$Q <- Q
      svd_output$Lambda <- Lambda
    }
    else if (is.null(Q) & attributes(weightinv)$diagonal == TRUE){
      diag_values <- diag(weightinv)
      data_w <- data %*% diag(sqrt(diag_values))
      svd_output <- svd(data_w)
      svd_output$v <- t(t(svd_output$v) %*% diag(1 / sqrt(diag_values)))
    }
    else if (!is.null(Q)){
      data_w <- data %*% Q %*% diag(sqrt(Lambda)) %*% t(Q)
      svd_output <- svd(data_w)
      svd_output$v <- t(t(svd_output$v) %*% Q %*% diag(1 / sqrt(Lambda)) %*% t(Q))
      svd_output$Q <- Q
      svd_output$Lambda <- Lambda
    }
  }
  return(svd_output)
}

#' Matrix inversion via cholesky decomposition
#'
#' Inverts matrix W, assigning attributes for whether W is diagonal, to speed up other calculations.
#'
#' @param W square positive definite variance matrix
#'
#' @return Inverse of W, with attributes 'identity' and 'diagonal', used by other functions in the package to make calculations more efficient.
#'
#' @examples Winv <- GetInverse(diag(100))
#' attributes(Winv) # diagonal = TRUE, identity = TRUE
#'
#' Winv2 <- GetInverse(runif(100,0.1,1)*diag(100))
#' attributes(Winv2) # diagonal = TRUE, identity = FALSE
#'
#' Winv3 <- GetInverse(seq(0.1,1,length=100) %*% t(seq(0.1,1,length=100)) + 0.1*diag(100))
#' attributes(Winv3) # diagonal = FALSE, identity = FALSE
#'
#' @export
GetInverse <- function(W){
  diagmat <- all(W[lower.tri(W)] == 0, W[upper.tri(W)] == 0)
  if (diagmat == TRUE){
    InvW <- diag(1 / diag(W))
  }
  else {
    Q <- chol(W)
    y <- backsolve(Q, diag(dim(W)[1]), transpose = TRUE)
    InvW <- crossprod(y, y)
  }
  attr(InvW, 'diagonal') <- diagmat
  if (all(diag(InvW) == 1) & diagmat == TRUE){
    attr(InvW, 'identity') <- TRUE
  }
  else {
    attr(InvW, 'identity') <- FALSE
  }
  return(InvW)
}





#' Field reconstructions from coefficients
#'
#' Given a vector of coefficients for a basis, calculates the field
#'
#' @param coeffs Coefficient vector
#' @param basis Basis matrix
#' @return Reconstructed field.
#'
#' @export
Recon <- function(coeffs, basis){
  if (is.null(dim(basis)[2])){
    q <- 1
  }
  else {
    q <- dim(basis)[2]
  }
  stopifnot(length(coeffs) == q)
  if (is.null(dim(basis)[2])){
    reconstruction <- basis*as.numeric(coeffs)
  }
  else {
    reconstruction <- basis%*%as.numeric(coeffs)
  }
  return(reconstruction)
}

#' Calculate a set of vectors from weights and a basis
#'
#' Give basis vectors from linear combinations
#'
#' @param weights A vector of weights
#' @param basis Original basis
#'
#' @return New basis vector(s)
#'
#' @export
ReconBasis <- function(weights, basis){
  n <- dim(basis)[2]
  q <- length(weights) / n
  if (q == 1){
    new.basis <- as.vector(tensor(basis, weights, 2, 1))
  }
  else {
    dim(weights) <- c(n, q)
    new.basis <- tensor(basis, weights, 2, 1)
  }
  return(new.basis)
}



#' Project and reconstruct a given field
#'
#' Gives the reconstruction of a field using a basis, by projecting and back-projecting on this basis.
#'
#' @param obs Vector over original field
#' @param basis Basis matrix
#'
#' @return Reconstruction of the original field.
#'
#' @examples
#'
#' @export
ReconObs <- function(obs, basis, ...){
  nb <- is.null(dim(basis))
  if(!nb)
    basis1 <- basis[,1]
  else
    basis1 <- basis
  obs <- c(obs)
  mask <- which(is.na(obs-basis1))
  if(length(mask)>0){
    recons <- rep(NA, length(obs))
    obs <- obs[-mask]
    if(nb)
      basis <- basis[-mask]
    else
      basis <- basis[-mask,]
    proj <- CalcScores(obs, basis, ...)
    recons.partial <- Recon(proj, basis)
    recons[-mask] <- recons.partial
  }
  else{
    proj <- CalcScores(obs, basis, ...)
    recons <- Recon(proj, basis)
  }
  return(recons)
}





#' Projection onto a basis
#'
#' Calculates the coefficients given by projecting data onto a basis
#'
#' @param data Data matrix to be projected, where each column is a representation on the original field
#' @param basis Basis matrix
#' @param weightinv If NULL, uses standard SVD projection. Otherwise, uses weighted projection.
#'
#' @return Matrix of basis coefficients
#'
#' @examples # First generate some data
#'
#' l <- 100 # dimension of output
#' n <- 10 # number of runs
#' DataBasis <- MakeDataBasis(data = matrix(runif(l*n), nrow=l, ncol=n), RemoveMean = TRUE) # data is 100x10
#'
#' # Project the (centred) ensemble onto the first 3 vectors of the SVD basis
#'
#' Coefficients <- CalcScores(data = DataBasis$CentredField, basis = DataBasis$tBasis[,1:3])
#'
#' # Instead of projecting using W = I, define a W with varying diagonal
#'
#' W <- runif(l, 1, 5) * diag(l) # 100x100 diagonal matrix
#' W_inv <- GetInverse(W) # inverse needed for projection
#' Coefficients_weighted <- CalcScores(data = DataBasis$CentredField, basis = DataBasis$tBasis[,1:3], weightinv = W_inv)
#'
#' @export
CalcScores <- function(data, basis, weightinv = NULL){
  d <- dim(data)[2]
  if (is.null(d)){
    d <- 1
  }
  p <- dim(basis)[2]
  l <- dim(basis)[1]
  if (is.null(p)){
    p <- 1
  }
  if (d == 1){
    data <- as.vector(data)
  }
  if (is.null(weightinv)){
    weightinv <- 0 # just need to set as something that isn't NULL so can give attribute
    attr(weightinv, 'diagonal') <- attr(weightinv, 'identity') <- TRUE
  }
  if (attributes(weightinv)$identity == TRUE){
    V <- t(basis) %*% basis
    Q <- chol(V)
    y <- backsolve(Q, diag(p), transpose = TRUE)
    x <- backsolve(Q, t(basis) %*% data, transpose = TRUE)
    scores <- crossprod(y, x)
  }
  else if (attributes(weightinv)$diagonal == TRUE) {
    V <- t(basis) %*% (diag(weightinv) * basis)
    Q <- chol(V)
    y <- backsolve(Q, diag(p), transpose = TRUE)
    tmp <- t(basis) %*% (diag(weightinv) * data)
    x <- backsolve(Q, tmp, transpose = TRUE)
    scores <- crossprod(y, x)
  }
  else {
    V <- t(basis) %*% weightinv %*% basis
    Q <- chol(V)
    y <- backsolve(Q, diag(p), transpose = TRUE)
    x <- backsolve(Q, t(basis) %*% weightinv %*% data, transpose = TRUE)
    scores <- crossprod(y, x)
  }
  return(t(scores))
}

#### Rename ####
#' Number of basis vectors required to explain proportion of data
#'
#' Finds the truncated basis that explains a set proportion of the variaiblity in the data.
#'
#' @param basis Basis matrix
#' @param data Data matrix
#' @param vtot The total proportion of variability in the data to be explained by the truncated basis
#' @param weightinv The inverse of W
#'
#' @return The number of basis vectors required to explain vtot of the data.
#'
#' @export
ExplainT <- function(DataBasis, vtot = 0.95, weightinv = NULL){
  v <- 0
  q <- 0
  while (v < vtot & q < dim(DataBasis$tBasis)[2]){
    v <- VarExplained(DataBasis$tBasis[,1:(q+1)], DataBasis$CentredField, weightinv)
    q <- q + 1
  }
  return(q)
}

#' Calculating the proportion of data explained by a basis
#'
#' Calculates the proportion of the data that is explained by projection onto a basis.
#'
#' @param basis The basis
#' @param data The data to be explained
#' @param weightinv Inverse of W (identity if NULL)
#' @param total_sum The total sum of squares of the data with respect to W
#' @param psi t(original_basis) %*% weightinv %*% original_basis, where the new basis is a linear combination of some original basis
#' @param basis_lincom Vector of linear combinations (if new basis is a linear combination of some original basis)
#'
#' @return The proportion of variability in the data that is explained by the basis
#'
#' @export
VarExplained <- function(basis, data, weightinv = NULL, total_sum = NULL, psi = NULL, basis_lincom = NULL){
  coeffs <- t(CalcScores(data, basis, weightinv))
  recon <- basis %*% coeffs
  if (is.null(weightinv)){
    explained <- crossprod(c(recon))/crossprod(c(data))
  }
  else {
    if (is.null(psi)){
      if (attributes(weightinv)$diagonal == TRUE){
        explained_num <- sum(t(recon)^2 %*% diag(weightinv))
      }
      else {
        explained_num <- sum(diag(t(recon) %*% weightinv %*% recon))
      }
    }
    else {
      stopifnot(!is.null(basis_lincom))
      explained_num <- t(coeffs) %*% t(basis_lincom) %*% psi %*%
        basis_lincom %*% coeffs
      explained_num <- sum(diag(explained_num))
    }
    #explained_num <- 0
    #for (i in 1:dim(data)[2]){
    #  explained_num <- explained_num + t(recon[,i]) %*% weightinv %*% recon[,i]
    #}
    #explained_den <- 0
    #for (i in 1:dim(data)[2]){
    #  explained_den <- explained_den + t(data[,i]) %*% weightinv %*% data[,i]
    #}
    if (is.null(total_sum)){
      if (attributes(weightinv)$diagonal == TRUE){
        explained_den <- sum(t(data)^2 %*% diag(weightinv))
      }
      else {
        explained_den <- sum(diag(t(data) %*% weightinv %*% data))
      }
    }
    else {
      explained_den <- total_sum
    }
    explained <- explained_num / explained_den
  }
  return(explained)
}


#' Reconstruction error
#'
#' Calculates the reconstruction error, R_W(basis, obs), of the observations given a basis and W.
#'
#' @param obs The observations
#' @param basis Basis to project and reconstruct the observations with
#' @param weightinv Inverse of weight matrix W. If NULL (default), calculates the mean squared error
#' @param scale If TRUE, scales by the dimension (so analogous to mean squared error)
#'
#' @return The reconstruction error
#'
#' @export
ReconError <- function(obs, basis, weightinv = NULL, scale = TRUE){
  if (is.null(weightinv)){
    weightinv <- 0
    attr(weightinv, 'diagonal') <- attr(weightinv, 'identity') <- TRUE
  }
  field <- ReconObs(obs, basis, weightinv)
  A <- c(obs) - field
  mask <- which(is.na(A))
  if(length(mask)>0){
    A <- A[-mask]
  }
  if (scale == TRUE){
    s <- length(c(obs))-length(mask)
  }
  else {
    s <- 1
  }
  if (attributes(weightinv)$diagonal == FALSE){
    if(length(mask)>0){
      warning("Implicit assumption that weight specified on the full field even though applying a mask to missing obs/ensemble grid boxes")
      weightinv <- weightinv[-mask,-mask]
    }
    wmse <- (t(A) %*% weightinv %*% A)/ s
  }
  else {
    if (attributes(weightinv)$identity == TRUE){
      wmse <- crossprod(A)/ s
    }
    else {
      wmse <- crossprod(A/(1/diag(weightinv)), A)/ s
    }
  }
  return(wmse)
}




#### More flexibility in specification, e.g. v, time allowed, make prior clearer, remove months etc. ####
#' Finding a calibration-optimal basis rotation
#'
#' Given a DataBasis object, observations, matrix W, vector v, applies the optimal rotation algorithm to find a basis more suitable for calibration.
#'
#' @param DataBasis An object containing lxn ensemble data and the lxn basis that will be rotated
#' @param obs A vector of length l with observations on the same scale as the ensemble
#' @param kmax Maximum number of iterations allowed (defaults to 5)
#' @param weightinv Inverse of positive definite weight matrix W. If set = NULL, uses the identity.
#' @param v Vector of minimum proportion of the ensemble data to be explained by the corresponding rotated basis vector
#' @param vtot Minimum proportion of ensemble variability to be explained by the truncated basis
#' @param MaxTime Maximum time (in seconds) to run the optimiser at each iteration.
#'
#' @return \item{tBasis}{Full rotated basis}
#' \item{CentredField}{The ensemble that was passed into the function}
#' \item{EnsembleMean}{The ensemble mean}
#' \item{scaling}{Initial scaling applied to the data}
#' \item{RW}{The reconstruction error after each iteration}
#' \item{VarExp}{The variance explained by each rotated basis vector}
#' \item{Opt}{Linear combination of the basis that gave rotated basis vectors} #### only true for first really ####
#'
#'@examples # First run an ensemble of idealised function fn
#'
#' n <- 60
#' sample <- as.data.frame(2*maximinLHS(n,6) - 1)
#' colnames(sample) <- c("x1","x2","x3","x4","x5","x6")
#' output <- array(c(rep(0,100*n)), dim=c(10,10,n))
#' for (i in 1:n){
#'   output[,,i] <- fn(as.numeric(sample[i,]))
#' }
#' dim(output) <- c(100, n)
#'
#' DataBasis <- MakeDataBasis(data = output, RemoveMean = TRUE)
#'
#' # Define the observations as a known value of x, plus some noise
#'
#' obs <- c(fn(c(0.7,0.01,0.01,0.25,0.8,-0.9)) + rnorm(100, mean = 0, sd = 0.1))
#' obs <- obs - DataBasis$EnsembleMean # centre observations so that comparable to the data
#'
#' # Look at the VarMSEplot for the SVD basis
#'
#' vSVD <- VarMSEplot(DataBasis = DataBasis, obs = obs)
#'
#' # Perform a rotation
#'
#' RotatedBasis <- RotateBasis(DataBasis = DataBasis, obs = obs, kmax = 3)
#'
#' # Editing the variance constraints so that the first vector explains at least 40% of the ensemble, at least 10% for later vectors
#'
#' RotatedBasis <- RotateBasis(DataBasis = DataBasis, obs = obs, kmax = 3, v = c(0.4,0.1,0.1))
#'
#' # So far assumed that W is the identity. Now add structure to W
#'
#'
#'
#' @export
RotateBasis <- function(DataBasis, obs, kmax = 5, weightinv = NULL, v = c(rep(0.1,5)), vtot = 0.95, prior = NULL,
                        StoppingRule = TRUE, MaxTime = 60,random.init = FALSE, ...){
  data <- DataBasis$CentredField
  basis <- matrix(DataBasis$tBasis,ncol = max(1,dim(DataBasis$tBasis)[2]))
  if (!dim(data)[1] == dim(basis)[1]){
    stop("Dimension of ensemble and basis (l) are not the same")
  }
  obs <- c(obs)
  if (!length(obs) == dim(basis)[1]){
    stop("Observations not the correct dimension (l)")
  }
  l <- dim(basis)[1]
  n <- dim(data)[2]
  minRw <- ReconError(obs, basis, weightinv)
  if (is.null(prior)){
    prior <- c(1:dim(basis)[2])
  }
  basis <- matrix(basis[,prior],ncol = length(prior))
  mse <- var <- numeric(kmax)
  x <- NULL
  new.basis <- NULL
  if (is.null(weightinv)){
    var_sum <- crossprod(c(data))
  }
  else {
    if (attributes(weightinv)$diagonal == TRUE){
      var_sum <- sum(t(data)^2 %*% diag(weightinv))
    }
    else {
      var_sum <- sum(diag(t(data) %*% weightinv %*% data))
    }
  }
  if (is.null(DataBasis$Q)){
    Q <- NULL
    Lambda <- NULL
  }
  else {
    Q <- DataBasis$Q
    Lambda <- DataBasis$Lambda
  }
  for (i in 1:kmax){
    p <- dim(basis)[2]
    if (is.null(weightinv)){
      psi <- t(basis) %*% diag(l) %*% basis
    }
    else {
      psi <- t(basis) %*% weightinv %*% basis
    }
    
    opt <- GenSA(c(1, rep(0, p-1)), WeightOptim,  lower = rep(-1, p*1),
                 upper = rep(1, p*1), basis = basis, obs = obs, data = data, weightinv = weightinv,
                 v = v[i], newvectors = new.basis, total_sum = var_sum, psi = psi, control = list(max.time = MaxTime), ...)
    best.patterns <- cbind(new.basis, ReconBasis(opt$par, basis))
    rank <- min(n, l)
    basis <- ResidBasis(best.patterns, data, weightinv, Q, Lambda)[,1:rank]
    x <- c(x, opt$par)
    q <- ExplainT(DataBasis, vtot, weightinv)
    mse[i] <- ReconError(obs, basis[,1:q], weightinv)
    var[i] <- VarExplained(basis[,i], data, weightinv)
    new.basis <- cbind(new.basis, basis[,i])
    basis <- basis[,-(1:i)]
    if (round(mse[i],4) == round(minRw,4)) break #### INSERT STOPPING RULE
  }
  new.basis <- cbind(new.basis, basis)[,1:rank]
  return(list(tBasis = new.basis, CentredField = DataBasis$CentredField,
              EnsembleMean = DataBasis$EnsembleMean, scaling = DataBasis$scaling,
              RW = mse, VarExp = var, Opt = x))
}


#' Find the residual basis
#'
#' Given basis vectors and data, calculate the residual basis to complete the basis for the data
#'
#' @param basisvectors Basis vector(s) to project the data onto
#' @param data Data matrix for to find a basis for
#'
#' @return The full basis for the data, i.e. the basis vectors passed into the function, with the residual basis vectors appended
#'
#' @export
ResidBasis <- function(basisvectors, data, weightinv = NULL, ...){
  if (is.null(weightinv)){
    basisvectors <- orthonormalization(basisvectors,basis=FALSE,norm=TRUE)
  }
  else {
    newvector <- basisvectors[,dim(basisvectors)[2]]
    basisvectors[,dim(basisvectors)[2]] <- newvector / as.numeric(sqrt(t(newvector)%*%weightinv %*% newvector))
  }
  l <- dim(data)[1]
  n <- dim(data)[2]
  recons <- matrix(numeric(l*n), nrow=l)
  for (i in 1:n){
    recons[,i] <- ReconObs(data[,i],basisvectors, weightinv)
  }
  resids <- data - recons
  if (is.null(weightinv)){
    svd.resid <- svd(t(resids))
  }
  else {
    svd.resid <- wsvd(t(resids), weightinv = weightinv, ...)
  }
  new.basis <- cbind(basisvectors, svd.resid$v)[,1:min(l,n)]
  return(new.basis)
}


#' Matrix projection
#'
#' Projects a variance matrix onto a given basis
#'
#' @param mat A square matrix to be projected onto the basis
#' @param basis The basis to project with
#' @param weightinv The inverse of positive definite matrix W. If NULL, uses the standard projection, otherwise projects in the norm given by W.
#'
#' @return The projection of the original matrix on the basis.
#'
#' @export
VarProj <- function(mat, basis, weightinv = NULL){
  if (is.null(weightinv)){
    proj <- t(basis) %*% mat %*% basis
  }
  else {
    V <- t(basis) %*% weightinv %*% basis
    Q <- chol(V)
    y <- backsolve(Q, diag(dim(basis)[2]), transpose = TRUE)
    x <- backsolve(Q, t(basis) %*% weightinv, transpose = TRUE)
    comp <- crossprod(y, x)
    proj <- comp %*% mat %*% t(comp)
  }
  return(proj)
}

#' Function used within optimiser for minimisation of the reconstruction error for rotated basis, subject to constraints
#'
#' Given a vector of weights, gives the new basis vector as a linear combination of the original basis, and calculates the reconstruction error, subject to a variability constraint
#'
#' @param x Vector giving the linear combination of the basis to use
#' @param basis The basis that is being rotated
#' @param obs Observation vector
#' @param data Ensemble data
#' @param weightinv Inverse of matrix W
#' @param v The proportion of variability to be explained by the basis vector
#' @param total_sum Common denominator used to calculate VarExplained
#' @param psi As new basis is linear combinarion of original, if pass psi = t(basis) %*% weightinv %*% basis adds efficiency
#' @param newvectors If the reconstruction error should account for any previous basis vectors
#'
#' @return The reconstruction error
#'
#' @export
WeightOptim <- function(x, basis, obs, data, weightinv, v = 0.1, total_sum = NULL, psi = NULL, newvectors = NULL){
  new.basis <- as.vector(tensor(basis, x, 2, 1))
  if (is.null(newvectors) == FALSE){
    new.basis <- cbind(newvectors, new.basis)
  }
  if (is.null(newvectors) == TRUE){
    v_new <- VarExplained(new.basis, data, weightinv, total_sum, psi, basis_lincom = x)
  }
  else {
    v_new <- VarExplained(new.basis[,dim(new.basis)[2]], data, weightinv, total_sum, psi, basis_lincom = x)
  }
  if (v_new < v){
    y <- 999999999
  }
  else {
    y <- ReconError(obs, new.basis, weightinv)
  }
  return(y)
}
