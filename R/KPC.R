#' \eqn{T_n} with geometric graphs
#'
#' Calculate \eqn{T_n} using directed Knn graph or minimum spanning tree (MST).
#'
#' \eqn{T_n} is an estimate of \eqn{E[E[k(Y_1,Y_1')|X]]}, with \eqn{Y_1}, \eqn{Y_1'} drawn iid from \eqn{Y|X}, given \eqn{X}.
#' For Knn graph, ties will be broken at random. Algorithm finding the MST is implemented the package \code{emstreeR}.
#'
#' @param X a matrix of predictors (n by dx)
#' @param Y a matrix of response (n by dy)
#' @param k a function \eqn{k(y, y')} of class \code{kernel}. It can be the kernel implemented in \code{kernlab} e.g. Gaussian kernel: \code{rbfdot(sigma = 1)}, linear kernel: \code{vanilladot()}. In practice, Gaussian kernel with empirical bandwidth \code{kernlab::rbfdot(1/(2*median(dist(Y))^2))} may be a good choice.
#' @param Knn the number of K-nearest neighbor to use; or "MST".
#' @export
#' @return Tn: a numeric number
TnKnn = function(Y,X,k,Knn=1) {
  if (!is.matrix(Y)) Y = as.matrix(Y)
  if (!is.matrix(X)) X = as.matrix(X)

  if (Knn == "MST") return(TnMST(Y,X,k))
  n = dim(Y)[1]

  # row i is the indices of KNN for Xi
  nn_index_X = get_neighbors(X,Knn)

  if (Knn == 1) {
    node_calculator = function(j) return(k(Y[j,],Y[nn_index_X[j,],]))
  }
  else {
    node_calculator = function(j) return(mean(kernelMatrix(k,Y[j,,drop=F],Y[nn_index_X[j,],,drop=F])))
  }

  return(mean(sapply(1:n, node_calculator)))
}

# Obtain KNN with ties broken at random
#
# @param X Matrix (n by dx)
# @param Knn number of nearest neighbors
# @return an n by Knn matrix showing the indices of KNN
get_neighbors = function(X,Knn) {
  if (!is.matrix(X)) X = as.matrix(X)
  # compute the nearest neighbor of X
  nn_X = RANN::nn2(X, query = X, k = Knn + 2)
  nn_index_X = nn_X$nn.idx[, 2:(Knn+1), drop=F]

  # find all data points that are not unique
  repeat_data = which(nn_X$nn.dists[, 2] == 0)
  if (length(repeat_data) > 0) {
    gp_size = id = NULL
    df_X = data.table::data.table(id = repeat_data, group = nn_X$nn.idx[repeat_data, 1])
    df_X[, gp_size := length(id), by = "group"]

    for (i in 1:length(repeat_data)) {
      if (df_X$gp_size[i] > Knn) {
        # The number of repeated data is more than Knn
        group_indices = df_X$id[df_X$group==df_X$group[i]]
        if (Knn == 1 & length(group_indices) == 2) {
          nn_index_X[df_X$id[i],] = setdiff(group_indices, df_X$id[i]) # avoid length 1 vector in sample function
        }
        else {
          nn_index_X[df_X$id[i],] = sample(setdiff(group_indices, df_X$id[i]),Knn)
        }
      }
      else {
        if (nn_X$nn.dists[df_X$id[i], Knn+1] < nn_X$nn.dists[df_X$id[i], Knn+2]) {
          # The number of repeated data is less than Knn
          # but there is no tie at the KNN
          nn_index_X[df_X$id[i],] = setdiff(nn_X$nn.idx[df_X$id[i], 1:(Knn+1)], df_X$id[i])
        }
        else {
          # The number of repeated data is less than Knn
          # There are ties at the kNN
          distances <- proxy::dist(matrix(X[df_X$id[i], ], ncol = ncol(X)), matrix(X[-df_X$id[i], ], ncol = ncol(X)))
          tie_dist <- sort(distances, partial = Knn)[Knn]
          id_small <- which(distances < tie_dist)
          id_small = id_small + (id_small >= df_X$id[i])
          nn_index_X[df_X$id[i],1:length(id_small)] = id_small
          id_equal = sample(which(distances == tie_dist),Knn-length(id_small))
          id_equal = id_equal + (id_equal >= df_X$id[i])
          nn_index_X[df_X$id[i],(1+length(id_small)):Knn] = id_equal
        }
      }
    }
  }
  ties = which(nn_X$nn.dists[, Knn+1] == nn_X$nn.dists[, Knn+2])
  ties = setdiff(ties, repeat_data)
  if (length(ties) > 0) {
    for (i in ties) {
      distances <- proxy::dist(matrix(X[i, ], ncol = ncol(X)), matrix(X[-i, ], ncol = ncol(X)))
      tie_dist <- sort(distances, partial = Knn)[Knn]
      id_small <- which(distances < tie_dist)
      if (length(id_small) > 0) {
        id_small = id_small + (id_small >= i)
        nn_index_X[i,1:length(id_small)] = id_small
      }
      id_equal = sample(which(distances == tie_dist),Knn-length(id_small))
      id_equal = id_equal + (id_equal >= i)
      nn_index_X[i,(1+length(id_small)):Knn] = id_equal
    }
  }
  return(nn_index_X)
}

# Calculate Tn using minimum spanning tree (MST).
TnMST = function(Y,X,k) {
  if (!is.matrix(Y)) Y = as.matrix(Y)
  n = dim(Y)[1]

  if (dim(X)[2] == 1) {
    Y = Y[order(X),,drop=F]

    node_calculator = function(j) {
      return(k(Y[j,],Y[j-1,]) + k(Y[j,],Y[j+1,]))
    }
    return((sum(sapply(2:(n-1), node_calculator))/2 + k(Y[1,],Y[2,]) + k(Y[n-1,],Y[n,]))/n)
  }

  if (!is.data.frame(X)) X = as.data.frame(X)

  out=emstreeR::ComputeMST(X,verbose = FALSE)
  tmp = matrix(0,n,2)
  # the first column is the degree of node i
  # the second column is the sum of k(xi,x_{N(i)})
  for (i in 1:(length(out$from)-1)) {
    tmp[out$from[i],1] = tmp[out$from[i],1] + 1
    tmp[out$to[i],1] = tmp[out$to[i],1] + 1
    tmp[out$from[i],2] = tmp[out$from[i],2] + k(Y[out$to[i],],Y[out$from[i],])
    tmp[out$to[i],2] = tmp[out$to[i],2] + k(Y[out$to[i],],Y[out$from[i],])
  }
  return(mean(tmp[,2]/tmp[,1]))
}





#' Kernel partial correlation with geometric graphs
#'
#' Calculate the kernel partial correlation (KPC) coefficient with directed Knn (nearest neighbor) graph or minimum spanning tree (MST).
#'
#' The kernel partial correlation squared (KPC) measures the conditional dependence
#' between \eqn{Y} and \eqn{Z} given \eqn{X}, based on an i.i.d. sample of \eqn{(Y, Z, X)}.
#' It converges to the population quantity which is between 0 and 1.
#' A small value indicates low conditional dependence between \eqn{Y} and \eqn{Z} given \eqn{X}, and
#' a large value indicates stronger conditional dependence.
#' If \code{X = NULL}, it returns the \code{Kmac(Y,Z,k,Knn)}, which measures the unconditional dependence between \eqn{Y} and \eqn{Z}.
#'
#' @param Y a matrix (n by dy)
#' @param X a matrix (n by dx) or \code{NULL} if \eqn{X} is empty
#' @param Z a matrix (n by dz)
#' @param k a function \eqn{k(y, y')} of class \code{kernel}. It can be the kernel implemented in \code{kernlab} e.g. Gaussian kernel: \code{rbfdot(sigma = 1)}, linear kernel: \code{vanilladot()}. In practice, Gaussian kernel with empirical bandwidth \code{kernlab::rbfdot(1/(2*median(dist(Y))^2))} may be a good choice.
#' @param Knn number of nearest neighbor to use; or "MST"
#' @param trans_inv TRUE or FALSE. Is \eqn{k(y, y)} free of \eqn{y}?
#'
#' @import data.table
#' @export
#'
#' @seealso \code{\link{Kmac}}, \code{\link{Klin}}
#'
#' @examples
#' library(kernlab)
#' n = 2000
#' x = rnorm(n)
#' z = rnorm(n)
#' y = x + z + rnorm(n,1,1)
#' KPCgraph(y,x,z,vanilladot(),Knn=1,trans_inv=FALSE)
#'
#' n = 1000
#' x = runif(n)
#' z = runif(n)
#' y = (x + z) %% 1
#' KPCgraph(y,x,z,rbfdot(5),Knn="MST",trans_inv=TRUE)
KPCgraph = function(Y,X,Z,k,Knn = 1,trans_inv=FALSE) {
  if (is.null(X)) return(Kmac(Y,Z,k,Knn))
  if (!is.matrix(Y)) Y = as.matrix(Y)
  if (!is.matrix(X)) X = as.matrix(X)
  if (!is.matrix(Z)) Z = as.matrix(Z)
  if ((nrow(Y) != nrow(X)) || (nrow(Y) != nrow(Z))) stop("Number of rows of the inputs should be equal.")
  if (Knn != "MST") {
    if ((floor(Knn) != Knn) || (Knn <= 0)) stop("Knn should be a positive integer or the string MST.")
    if (Knn + 2 > nrow(X)) stop("n should be greater than Knn + 1")
  }

  Tn_XZ = TnKnn(Y,cbind(X,Z),k,Knn)
  Tn_X = TnKnn(Y,X,k,Knn)

  if (trans_inv) {
    return((Tn_XZ - Tn_X)/(k(Y[1,],Y[1,])-Tn_X))
  }
  else {
    node_calculator = function(j) {
      return(k(Y[j,],Y[j,]))
    }

    return((Tn_XZ - Tn_X)/(mean(sapply(1:nrow(Y), node_calculator))-Tn_X))
  }
}

# Double-centering
#
# Double-centering a squared matrix
#
# Given a square matrix A, the function returns HAH, where H = diag(n) - 1/n is the centering matrix.
#
# @param M Matrix (n by n)
double_center = function(M){
  return(M - rowMeans(M) - matrix(colMeans(M),nrow(M),ncol(M),byrow = T) + mean(M))
}

#' Kernel partial correlation with RKHS method
#'
#' Compute estimate of Kernel partial correlation (KPC) coefficient using conditional mean embeddings in the reproducing kernel Hilbert spaces (RKHS).
#'
#' The kernel partial correlation (KPC) coefficient measures the conditional dependence
#' between \eqn{Y} and \eqn{Z} given \eqn{X}, based on an i.i.d. sample of \eqn{(Y, Z, X)}.
#' It converges to the population quantity which is between 0 and 1.
#' A small value indicates low conditional dependence between \eqn{Y} and \eqn{Z} given \eqn{X}, and
#' a large value indicates stronger conditional dependence.
#' If \code{X = NULL}, it measures the unconditional dependence between \eqn{Y} and \eqn{Z}.
#'
#' @param Y a matrix (n by dy)
#' @param X a matrix (n by dx) or \code{NULL} if \eqn{X} is empty
#' @param Z a matrix (n by dz)
#' @param ky a function \eqn{k(y, y')} of class \code{kernel}. It can be the kernel implemented in \code{kernlab} e.g. Gaussian kernel: \code{rbfdot(sigma = 1)}, linear kernel: \code{vanilladot()}. In practice, Gaussian kernel with empirical bandwidth \code{kernlab::rbfdot(1/(2*median(dist(Y))^2))} may be a good choice.
#' @param kx the kernel function for \eqn{X}
#' @param kxz the kernel function for \eqn{(X, Z)} or for \eqn{Z} if \eqn{X} is empty
#' @param eps a small positive regularization parameter for inverting the empirical cross-covariance operator
#' @param appro whether to use incomplete Cholesky decomposition for approximation
#' @param tol tolerance used for incomplete Cholesky decomposition (implemented by the function \code{inchol} in the package \code{kernlab})
#' @import kernlab
#' @export
#' @examples
#' n = 500
#' set.seed(1)
#' x = rnorm(n)
#' z = rnorm(n)
#' y = x + z + rnorm(n,1,1)
#' library(kernlab)
#' k = vanilladot()
#' KPCRKHS(y, x, z, k, k, k, 1e-5/n^(0.4), appro = FALSE)
#' # 0.4855134 (Population quantity = 0.5)
#' KPCRKHS(y, x, z, k, k, k, 1e-5/n^(0.4), appro = TRUE, tol = 1e-5)
#' # 0.4855134 (Population quantity = 0.5)
KPCRKHS = function(Y, X = NULL, Z, ky, kx, kxz, eps, appro = FALSE, tol = 1e-5) {
  if (!is.matrix(Y)) Y = as.matrix(Y)
  if (!is.null(X)) {
    if (!is.matrix(X)) X = as.matrix(X)
    if ((nrow(Y) != nrow(X))) stop("Number of rows of the inputs should be equal.")
  }
  if (!is.matrix(Z)) Z = as.matrix(Z)
  if ((nrow(Y) != nrow(Z))) stop("Number of rows of the inputs should be equal.")

  n = dim(Y)[1]
  if (!appro) {
    # exact computation
    tilde_Ky = double_center(kernlab::kernelMatrix(ky,Y))
    if (is.null(X)) {
      M = diag(n) - n*eps*solve(double_center(kernlab::kernelMatrix(kxz,Z))+n*eps*diag(n))
      numerator = sum(tilde_Ky * base::crossprod(M))
      denominator = sum(diag(tilde_Ky))
      return(numerator/denominator)
    }
    else {
      N = solve(double_center(kernlab::kernelMatrix(kx,X)) + n*eps*diag(n))
      denominator = sum(tilde_Ky * base::crossprod(N))
      numerator = sum(tilde_Ky * base::crossprod(solve(double_center(kernlab::kernelMatrix(kxz,cbind(X,Z))) + n*eps*diag(n)) - N))
      return(numerator/denominator)
    }
  }
  # Approximate computation with incomplete Cholesky decomposition
  if (is.null(X)) {
    Lz = inchol(Z, kxz, tol = tol)
    Lz = Lz - rep(colMeans(Lz), rep.int(n, ncol(Lz)))
    Ly = inchol(Y, ky, tol = tol)
    # a close examination of M shows we don't need to center Ly
    return(sum((t(Ly)%*%Lz%*%solve(dim(Y)[1]*eps*diag(dim(Lz)[2]) + t(Lz)%*%Lz)%*%t(Lz))^2)/sum(diag(double_center(kernlab::kernelMatrix(ky,Y)))))
  }
  L1 = inchol(X, kx, tol = tol)
  L2 = inchol(cbind(X,Z), kxz, tol = tol)
  L3 = inchol(Y, ky, tol = tol)
  L1 = L1 - rep(colMeans(L1), rep.int(n, ncol(L1)))
  L2 = L2 - rep(colMeans(L2), rep.int(n, ncol(L2)))
  L3 = L3 - rep(colMeans(L3), rep.int(n, ncol(L3)))
  N = diag(n) - L1%*%solve(n*eps*diag(dim(L1)[2]) + t(L1)%*%L1)%*%t(L1)
  denominator = sum((N%*%L3)^2)
  M = N - diag(n) + L2%*%solve(n*eps*diag(dim(L2)[2]) + t(L2)%*%L2)%*%t(L2)
  numerator = sum((M%*%L3)^2)
  return(numerator/denominator)
}

#' Kernel partial correlation with RKHS method using linear kernels
#'
#' Linear kernels are used for ky, kx, kxz in KPCRKHS, in which case the incomplete Cholesky decomposition is the data matrix itself and could speed up the computation to a great extent.
#'
#' Linear kernels are used for ky, kx, kxz in KPCRKHS.
#' The kernel partial correlation (KPC) coefficient measures the conditional dependence
#' between \eqn{Y} and \eqn{Z} given \eqn{X}, based on an i.i.d. sample of \eqn{(Y, Z, X)}.
#' It converges to the population quantity which is between 0 and 1.
#' A small value indicates low conditional dependence between \eqn{Y} and \eqn{Z} given \eqn{X}, and
#' a large value indicates stronger conditional dependence.
#' If \code{X = NULL}, it measures the unconditional dependence between \eqn{Y} and \eqn{Z}.
#'
#' @param Y a matrix (n by dy)
#' @param Z a matrix (n by dz)
#' @param X a matrix (n by dx) or \code{NULL} if \eqn{X} is empty
#' @param eps a small regularization parameter for inverting the empirical cross-covariance operator
#' @export
#' @examples
#' n = 2000
#' x = rnorm(n)
#' z = rnorm(n)
#' y = x + z + rnorm(n,1,1)
#' KPCRKHSlinear(y, x, z, 1e-5/n^(0.4))
KPCRKHSlinear = function(Y, X = NULL, Z, eps) {
  if (!is.matrix(Y)) Y = as.matrix(Y)
  if (!is.null(X)) {
    if (!is.matrix(X)) X = as.matrix(X)
    if ((nrow(Y) != nrow(X))) stop("Number of rows of the inputs should be equal.")
  }
  if (!is.matrix(Z)) Z = as.matrix(Z)
  if ((nrow(Y) != nrow(Z))) stop("Number of rows of the inputs should be equal.")

  n = dim(Y)[1]

  # Incomplete Cholesky decomposition is exact here
  if (is.null(X)) {
    Lz = Z
    Lz = Lz - rep(colMeans(Lz), rep.int(n, ncol(Lz)))
    Ly = Y
    tilde_Ly = Ly - rep(colMeans(Ly), rep.int(n, ncol(Ly)))
    norm_sq = function(x) return(sum(x^2))
    denominator = sum(apply(tilde_Ly, 1, norm_sq))
    return(sum((t(Ly)%*%Lz%*%solve(dim(Y)[1]*eps*diag(dim(Lz)[2]) + t(Lz)%*%Lz)%*%t(Lz))^2)/denominator)
  }
  L1 = X
  L2 = cbind(X,Z)
  L3 = Y
  L1 = L1 - rep(colMeans(L1), rep.int(n, ncol(L1)))
  L2 = L2 - rep(colMeans(L2), rep.int(n, ncol(L2)))
  L3 = L3 - rep(colMeans(L3), rep.int(n, ncol(L3)))
  N = diag(n) - L1%*%solve(n*eps*diag(dim(L1)[2]) + t(L1)%*%L1)%*%t(L1)
  denominator = sum((N%*%L3)^2)
  M = N - diag(n) + L2%*%solve(n*eps*diag(dim(L2)[2]) + t(L2)%*%L2)%*%t(L2)
  numerator = sum((M%*%L3)^2)
  return(numerator/denominator)
}



#' Kernel Feature Ordering by Conditional Independence
#'
#' Variable selection with KPC using directed Knn graph or minimum spanning tree (MST)
#'
#' A stepwise forward selection of variables using KPC. At each step the \eqn{Xj} maximizing \eqn{\hat{\rho^2}(Y,X_j | selected X_i)} is selected.
#' It is suggested to normalize the predictors before applying KFOCI.
#'
#' @param Y a matrix of responses (n by dy)
#' @param X a matrix of predictors (n by dx)
#' @param k a function \eqn{k(y, y')} of class \code{kernel}. It can be the kernel implemented in \code{kernlab} e.g. Gaussian kernel: \code{rbfdot(sigma = 1)}, linear kernel: \code{vanilladot()}. In practice, Gaussian kernel with empirical bandwidth \code{kernlab::rbfdot(1/(2*median(dist(Y))^2))} may be a good choice.
#' @param Knn the number of nearest neighbor; or "MST"
#' @param num_features the number of variables to be selected, cannot be larger than dx. The default value is NULL and in that
#'   case it will be set equal to dx. If stop == TRUE (see below), then num_features is the maximal number of variables to be selected.
#' @param stop If `stop == TRUE`, then the automatic stopping criterion (stops at the first instance of negative Tn, as mentioned in the paper) will be implemented and continued till `num_features` many variables are selected. If `stop == FALSE` then exactly `num_features` many variables are selected.
#' @param numCores number of cores that are going to be used for parallelizing the process.
#' @param verbose whether to print each selected variables during the forward stepwise algorithm
#' @export
#' @return a vector of the indices from 1,...,dx of the selected variables
#' @examples
#' n = 200
#' p = 10
#' X = matrix(rnorm(n * p), ncol = p)
#' Y = X[, 1] * X[, 2] + sin(X[, 1] * X[, 3])
#' KFOCI(Y, X, kernlab::rbfdot(1), Knn=1, numCores=1)
#'
#' ### install the package olsrr first
#' #surgical = olsrr::surgical
#' #for (i in 1:9) surgical[,i] = (surgical[,i] - mean(surgical[,i]))/sd(surgical[,i])
#' #ky = kernlab::rbfdot(1/(2*median(dist(surgical$y))^2))
#' #colnames(surgical)[KFOCI(surgical[,9],surgical[,1:8],ky,Knn=1)]
#' #### "enzyme_test" "pindex" "liver_test"  "alc_heavy"
#' \dontrun{
#' # This example may take several minutes on a personal computer.
#' n = 2000
#' p = 1000
#' set.seed(1)
#' X = matrix(rnorm(n * p), ncol = p)
#' Y = X[, 1] * X[, 2] + sin(X[, 1] * X[, 3])
#' KFOCI(Y, X, kernlab::rbfdot(1), Knn=1, numCores = 7, verbose=TRUE)
#' # 1 2 3
#' }
# code modified from Azadkia, M. and Chatterjee, S. (2019). A simple measure of conditional dependence.
KFOCI <- function(Y, X, k, Knn = 1, num_features = NULL, stop = TRUE, numCores = 1, verbose = FALSE){
  if (!is.matrix(X)) X = as.matrix(X)
  if (!is.matrix(Y)) Y = as.matrix(Y)
  if ((nrow(Y) != nrow(X))) stop("Number of rows of Y and X should be equal.")
  if (is.null(num_features)) num_features <- dim(X)[2]
  if (num_features > ncol(X)) stop("Number of features should not be larger than maximum number of original features.")
  if ((floor(num_features) != num_features) || (num_features <= 0)) stop("Number of features should be a positive integer.")
  if (Knn != "MST") {
    if ((floor(Knn) != Knn) || (Knn <= 0)) stop("Knn should be a positive integer or the string MST.")
    if (Knn + 2 > nrow(X)) stop("n should be greater than Knn + 1")
  }
  n = dim(Y)[1]
  p = ncol(X)
  Q = rep(0, num_features) # stores the values of Tn
  index_select = rep(0, num_features)
  # select the first variable
  estimateQFixedY <- function(id){
    return(TnKnn(Y, X[, id],k,Knn))
  }
  seq_Q = parallel::mclapply(seq(1, p), estimateQFixedY, mc.cores = numCores)
  seq_Q = unlist(seq_Q)


  Q[1] = max(seq_Q)
  if (Q[1] <= 0 & stop == TRUE) return(0)
  index_max = min(which(seq_Q == Q[1]))
  index_select[1] = index_max
  if (verbose) print(paste("Variable",index_max,"is selected"))
  count = 1

  # select rest of the variables
  while (count < num_features) {
    seq_Q = rep(0, p - count)
    # indices that have not been selected yet
    index_left = setdiff(seq(1, p), index_select[1:count])

    # find the next best feature
    estimateQFixedYandSubX <- function(id){
      return(TnKnn(Y, X[, c(index_select[1:count], id)],k,Knn))
    }

    if (length(index_left) == 1) {
      seq_Q = estimateQFixedYandSubX(index_left[1])
    } else {
      seq_Q = parallel::mclapply(index_left, estimateQFixedYandSubX, mc.cores = numCores)
      seq_Q = unlist(seq_Q)
    }
    Q[count + 1] = max(seq_Q)
    index_max = min(which(seq_Q == Q[count + 1]))
    if (Q[count + 1] <= Q[count] & stop == TRUE) break
    index_select[count + 1] = index_left[index_max]
    count = count + 1
    if (verbose) print(paste("Variable",index_select[count],"is selected"))
  }

  return(index_select[1:count])
}

#' Variable selection with RKHS estimator
#'
#' The algorithm performs a forward stepwise variable selection using RKHS estimators.
#'
#' A stepwise forward selection of variables using KPC. At each step the \eqn{Xj} maximizing \eqn{\tilde{\rho^2}(Y,X_j | selected X_i)} is selected.
#' It is suggested to normalize the features before applying the algorithm.
#'
#' @param Y a matrix of responses (n by dy)
#' @param X a matrix of predictors (n by dx)
#' @param ky a function \eqn{k(y, y')} of class \code{kernel}. It can be the kernel implemented in \code{kernlab} e.g. Gaussian kernel: \code{rbfdot(sigma = 1)}, linear kernel: \code{vanilladot()}. In practice, Gaussian kernel with empirical bandwidth \code{kernlab::rbfdot(1/(2*median(dist(Y))^2))} may be a good choice.
#' @param kS a function that takes X and a subset of indices S as inputs, and then outputs the kernel for X_S. The first argument of kS is X, and the second argument is a vector of positive integer.
#' @param num_features the number of variables to be selected, cannot be larger than dx.
#' @param eps a positive number; the regularization parameter for the RKHS estimator
#' @param appro whether to use incomplete Cholesky decomposition for approximation
#' @param tol tolerance used for incomplete Cholesky decomposition (\code{inchol} in package \code{kernlab})
#' @param numCores number of cores that are going to be used for parallelizing the process.
#' @param verbose whether to print each selected variables during the forward stepwise algorithm
#' @return a vector of the indices from \code{1,...,dx} of the selected variables
#' @export
#' @examples
#' n = 200
#' p = 10
#' X = matrix(rnorm(n * p), ncol = p)
#' Y = X[, 1] * X[, 2] + sin(X[, 1] * X[, 3])
#' library(kernlab)
#' kS = function(X,S) return(rbfdot(1/length(S)))
#' RKHS_select(Y, X, rbfdot(1), kS, 3, eps = 1e-3, appro = FALSE, numCores = 1)
#' kS = function(X,S) return(rbfdot(1/(2*stats::median(stats::dist(X[,S]))^2)))
#' RKHS_select(Y, X, rbfdot(1), kS, 3, eps = 1e-3, appro = FALSE, numCores = 1)
# code modified from Azadkia, M. and Chatterjee, S. (2019). A simple measure of conditional dependence.
RKHS_select <- function(Y, X, ky, kS, num_features, eps, appro = FALSE, tol = 1e-5, numCores = 1, verbose = FALSE){
  if (!is.matrix(X)) X = as.matrix(X)
  if (!is.matrix(Y)) Y = as.matrix(Y)
  if ((nrow(Y) != nrow(X))) stop("Number of rows of Y and X should be equal.")
  if (is.null(num_features)) num_features <- dim(X)[2]
  if (num_features > ncol(X)) stop("Number of features should not be larger than maximum number of original features.")
  if ((floor(num_features) != num_features) || (num_features <= 0)) stop("Number of features should be a positive integer.")

  n = dim(Y)[1]
  p = ncol(X)
  Q = rep(0, num_features) # stores the values of Tn
  index_select = rep(0, num_features)
  # select the first variable
  estimateQFixedY <- function(id){
    return(KPCRKHS_numerator(Y,NULL,X[,id],ky,NULL,kS(X,id),eps,appro,tol))
  }
  seq_Q = parallel::mclapply(seq(1, p), estimateQFixedY, mc.cores = numCores)
  seq_Q = unlist(seq_Q)


  Q[1] = max(seq_Q)
  index_max = min(which(seq_Q == Q[1]))
  index_select[1] = index_max
  if (verbose) print(paste("Variable",index_max,"is selected"))
  count = 1

  # select rest of the variables
  while (count < num_features) {
    seq_Q = rep(0, p - count)
    # indices that have not been selected yet
    index_left = setdiff(seq(1, p), index_select[1:count])

    # find the next best feature
    estimateQFixedYandSubX <- function(id){
      return(KPCRKHS_numerator(Y, X[,index_select[1:count]], X[, c(index_select[1:count], id)], ky, kS(X,index_select[1:count]), kS(X,c(index_select[1:count], id)), eps,appro,tol))
    }

    if (length(index_left) == 1) {
      seq_Q = estimateQFixedYandSubX(index_left[1])
    } else {
      seq_Q = parallel::mclapply(index_left, estimateQFixedYandSubX, mc.cores = numCores)
      seq_Q = unlist(seq_Q)
    }
    Q[count + 1] = max(seq_Q)
    index_max = min(which(seq_Q == Q[count + 1]))
    index_select[count + 1] = index_left[index_max]
    count = count + 1
    if (verbose) print(paste("Variable",index_select[count],"is selected"))
  }

  return(index_select[1:count])
}



# calculate the numerator of the RKHS estimator
# used for stepwise variable selection
KPCRKHS_numerator = function(Y, X = NULL, Z, ky, kx, kxz, eps, appro = FALSE, tol = 1e-5) {
  if (!is.matrix(Y)) Y = as.matrix(Y)
  if (!is.matrix(Z)) Z = as.matrix(Z)
  if (!is.null(X) & !is.matrix(X)) X = as.matrix(X)

  n = dim(Y)[1]
  if (!appro) {
    # exact computation
    tilde_Ky = double_center(kernlab::kernelMatrix(ky,Y))
    if (is.null(X)) {
      M = diag(n) - n*eps*solve(double_center(kernlab::kernelMatrix(kxz,Z))+n*eps*diag(n))
      numerator = sum(tilde_Ky * base::crossprod(M))
      return(numerator)
    }
    else {
      N = solve(double_center(kernlab::kernelMatrix(kx,X)) + n*eps*diag(n))
      numerator = sum(tilde_Ky * base::crossprod(solve(double_center(kernlab::kernelMatrix(kxz,cbind(X,Z))) + n*eps*diag(n)) - N))
      return(numerator)
    }
  }
  # Approximate computation with incomplete Cholesky decomposition
  if (is.null(X)) {
    Lz = inchol(Z, kxz, tol = tol)
    Lz = Lz - rep(colMeans(Lz), rep.int(n, ncol(Lz)))
    Ly = inchol(Y, ky, tol = tol)
    return(sum((t(Ly)%*%Lz%*%solve(dim(Y)[1]*eps*diag(dim(Lz)[2]) + t(Lz)%*%Lz)%*%t(Lz))^2))
  }
  L1 = inchol(X, kx, tol = tol)
  L2 = inchol(cbind(X,Z), kxz, tol = tol)
  L3 = inchol(Y, ky, tol = tol)
  L1 = L1 - rep(colMeans(L1), rep.int(n, ncol(L1)))
  L2 = L2 - rep(colMeans(L2), rep.int(n, ncol(L2)))
  L3 = L3 - rep(colMeans(L3), rep.int(n, ncol(L3)))
  M = - L1%*%solve(n*eps*diag(dim(L1)[2]) + t(L1)%*%L1)%*%t(L1) + L2%*%solve(n*eps*diag(dim(L2)[2]) + t(L2)%*%L2)%*%t(L2)
  numerator = sum((M%*%L3)^2)
  return(numerator)
}
# KPCRKHS_numerator = function(Y, X = NULL, Z, ky, eps, appro = FALSE, tol = 1e-5) {
#   if (!is.matrix(Y)) Y = as.matrix(Y)
#   if (!is.matrix(Z)) Z = as.matrix(Z)
#   if (!is.null(X) & !is.matrix(X)) X = as.matrix(X)
#
#   if (is.null(X)) {
#     kxz = kernlab::rbfdot(1/(2*stats::median(stats::dist(Z))^2))
#   }
#   else {
#     kx = kernlab::rbfdot(1/(2*stats::median(stats::dist(X))^2))
#     kxz = kernlab::rbfdot(1/(2*stats::median(stats::dist(cbind(X,Z)))^2))
#   }
#
#   n = dim(Y)[1]
#   if (!appro) {
#     # exact computation
#     tilde_Ky = double_center(kernlab::kernelMatrix(ky,Y))
#     if (is.null(X)) {
#       M = diag(n) - n*eps*solve(double_center(kernlab::kernelMatrix(kxz,Z))+n*eps*diag(n))
#       numerator = sum(tilde_Ky * base::crossprod(M))
#       return(numerator)
#     }
#     else {
#       N = solve(double_center(kernlab::kernelMatrix(kx,X)) + n*eps*diag(n))
#       numerator = sum(tilde_Ky * base::crossprod(solve(double_center(kernlab::kernelMatrix(kxz,cbind(X,Z))) + n*eps*diag(n)) - N))
#       return(numerator)
#     }
#   }
#   # Approximate computation with incomplete Cholesky decomposition
#   if (is.null(X)) {
#     Lz = inchol(Z, kxz, tol = tol)
#     Lz = Lz - rep(colMeans(Lz), rep.int(n, ncol(Lz)))
#     Ly = inchol(Y, ky, tol = tol)
#     return(sum((t(Ly)%*%Lz%*%solve(dim(Y)[1]*eps*diag(dim(Lz)[2]) + t(Lz)%*%Lz)%*%t(Lz))^2))
#   }
#   L1 = inchol(X, kx, tol = tol)
#   L2 = inchol(cbind(X,Z), kxz, tol = tol)
#   L3 = inchol(Y, ky, tol = tol)
#   L1 = L1 - rep(colMeans(L1), rep.int(n, ncol(L1)))
#   L2 = L2 - rep(colMeans(L2), rep.int(n, ncol(L2)))
#   L3 = L3 - rep(colMeans(L3), rep.int(n, ncol(L3)))
#   M = - L1%*%solve(n*eps*diag(dim(L1)[2]) + t(L1)%*%L1)%*%t(L1) + L2%*%solve(n*eps*diag(dim(L2)[2]) + t(L2)%*%L2)%*%t(L2)
#   numerator = sum((M%*%L3)^2)
#   return(numerator)
# }



#' \eqn{\hat{\eta}_n} (the unconditional version of graph-based KPC) with geometric graphs.
#'
#' Calculate \eqn{\hat{\eta}_n} (the unconditional version of graph-based KPC) using directed Knn graph or minimum spanning tree (MST).
#'
#' \eqn{\hat{\eta}_n} is an estimate of the population kernel measure of association, based on data \eqn{(X_1,Y_1),\ldots ,(X_n,Y_n)\sim \mu}.
#' For Knn graph, ties will be broken at random. MST is found using package \code{emstreeR}.
#' In particular,
#' \deqn{\hat{\eta}_n:=\frac{n^{-1}\sum_{i=1}^n d_i^{-1}\sum_{j:(i,j)\in\mathcal{E}(G_n)} K(Y_i,Y_j)-(n(n-1))^{-1}\sum_{i\neq j}K(Y_i,Y_j)}{n^{-1}\sum_{i=1}^n K(Y_i,Y_i)-(n(n-1))^{-1}\sum_{i\neq j}K(Y_i,Y_j)}},
#' where \eqn{G_n} denotes a MST or Knn graph on \eqn{X_1,\ldots , X_n}, \eqn{\mathcal{E}(G_n)} denotes the set of edges of \eqn{G_n} and
#' \eqn{(i,j)\in\mathcal{E}(G_n)} implies that there is an edge from \eqn{X_i} to \eqn{X_j} in \eqn{G_n}.
#'
#' @param Y a matrix of response (n by dy)
#' @param X a matrix of predictors (n by dx)
#' @param k a function \eqn{k(y, y')} of class \code{kernel}. It can be the kernel implemented in \code{kernlab} e.g. Gaussian kernel: \code{rbfdot(sigma = 1)}, linear kernel: \code{vanilladot()}. In practice, Gaussian kernel with empirical bandwidth \code{kernlab::rbfdot(1/(2*median(dist(Y))^2))} may be a good choice.
#' @param Knn the number of K-nearest neighbor to use; or "MST".
#' @return kmac: the empirical kernel measure of association
#' @export
#' @references Deb, N., P. Ghosal, and B. Sen (2020), “Measuring association on topological spaces using kernels and geometric graphs” <arXiv:2010.01768>.
#' @examples
#' library(kernlab)
#' Kmac(Y = rnorm(100), X = rnorm(100), k = rbfdot(1), Knn = 1)
Kmac = function(Y,X,k,Knn=1) {
  if (!is.matrix(Y)) Y = as.matrix(Y)
  if (!is.matrix(X)) X = as.matrix(X)
  if ((nrow(Y) != nrow(X))) stop("Number of rows of the inputs should be equal.")
  if (Knn != "MST") {
    if ((floor(Knn) != Knn) || (Knn <= 0)) stop("Knn should be a positive integer or the string MST.")
    if (Knn + 2 > nrow(X)) stop("n should be greater than Knn + 1")
  }
  kernelm = kernelMatrix(k,Y)
  dirsum=sum(diag(kernelm))
  crosssum=sum(kernelMatrix(k,Y))-dirsum
  n = dim(Y)[1]
  return((TnKnn(Y,X,k,Knn)-crosssum/(n*(n-1)))/(dirsum/n-crosssum/(n*(n-1))))
}

#' \eqn{\hat{\eta}_n^{\mbox{lin}}} with geometric graphs.
#'
#' Calculate \eqn{\hat{\eta}_n^{\mbox{lin}}} (the unconditional version of graph-based KPC) using directed Knn graph or minimum spanning tree (MST).
#' The computational complexity is O(nlog(n))
#'
#' \eqn{\hat{\eta}_n} is an estimate of the population kernel measure of association, based on data \eqn{(X_1,Y_1),\ldots ,(X_n,Y_n)\sim \mu}.
#' For Knn graph, \eqn{\hat{\eta}_n} can be computed in near linear time (in \eqn{n}).
#' In particular,
#' \deqn{\hat{\eta}_n^{\mbox{lin}}:=\frac{n^{-1}\sum_{i=1}^n d_i^{-1}\sum_{j:(i,j)\in\mathcal{E}(G_n)} K(Y_i,Y_j)-(n-1)^{-1}\sum_{i=1}^{n-1} K(Y_i,Y_{i+1})}{n^{-1}\sum_{i=1}^n K(Y_i,Y_i)-(n-1)^{-1}\sum_{i=1}^{n-1} K(Y_i,Y_{i+1})}},
#' where all symbols have their usual meanings as in the definition of \eqn{\hat{\eta}_n}.
#'
#' @param Y a matrix of response (n by dy)
#' @param X a matrix of predictors (n by dx)
#' @param k a function \eqn{k(y, y')} of class \code{kernel}. It can be the kernel implemented in \code{kernlab} e.g. \code{rbfdot(sigma = 1)}, \code{vanilladot()}
#' @param Knn the number of K-nearest neighbor to use; or "MST".
#' @return klin: an empirical kernel measure of association which can be computed in near linear time when Knn graphs are used.
#' @export
#' @references Deb, N., P. Ghosal, and B. Sen (2020), “Measuring association on topological spaces using kernels and geometric graphs” <arXiv:2010.01768>.
#' @examples
#' library(kernlab)
#' Klin(Y = rnorm(100), X = rnorm(100), k = rbfdot(1), Knn = 1)
Klin = function(Y,X,k,Knn=1) {
  if (!is.matrix(Y)) Y = as.matrix(Y)
  if (!is.matrix(X)) X = as.matrix(X)
  if ((nrow(Y) != nrow(X))) stop("Number of rows of the inputs should be equal.")
  if (Knn != "MST") {
    if ((floor(Knn) != Knn) || (Knn <= 0)) stop("Knn should be a positive integer or the string MST.")
    if (Knn + 2 > nrow(X)) stop("n should be greater than Knn + 1")
  }
  n = dim(Y)[1]
  kernelm = kernelMatrix(k,Y)

  node_calculator = function(j) return(k(Y[j,],Y[j,]))
  dirsum = sum(sapply(1:n, node_calculator))

  node_calculator = function(j) return(k(Y[j,],Y[j+1,]))
  crosssum = sum(sapply(1:(n-1), node_calculator))

  return((TnKnn(Y,X,k,Knn)-crosssum/(n-1))/(dirsum/n-crosssum/(n-1)))
}
