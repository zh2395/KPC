#' \eqn{T_n} with geometric graphs
#'
#' Calculate \eqn{T_n} using directed Knn graph or minimum spanning tree (MST).
#'
#' \eqn{T_n} is an estimate of \eqn{E[E[k(Y_1,Y_1')|X]]}, with \eqn{Y_1}, \eqn{Y_1'} drawn iid from \eqn{Y|X}, given \eqn{X}.
#' For Knn graph, ties will be broken at random. MST is found using package \code{emstreeR}.
#'
#' @param X a matrix of predictors (n by dx)
#' @param Y a matrix of response (n by dy)
#' @param k a function \eqn{k(y, y')} of class \code{kernel}. It can be the kernel implemented in \code{kernlab} e.g. \code{rbfdot(sigma = 1)}, \code{vanilladot()}
#' @param Knn the number of K-nearest neighbor to use; or "MST".
#' @return Tn: a numeric number
TnKnn = function(Y,X,k,Knn=1) {
  if (Knn == "MST") return(TnMST(Y,X,k))

  if(!is.matrix(Y)) Y = as.matrix(Y)

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
  # need Knn + 1 <= dx
  nn_X = RANN::nn2(X, query = X, k = Knn + 2)
  nn_index_X = nn_X$nn.idx[, 2:(Knn+1), drop=F]

  # find all data points that are not unique
  repeat_data = which(nn_X$nn.dists[, 2] == 0)
  if (length(repeat_data) > 0) {
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
          # but there is no tie at the kNN
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
  if(!is.data.frame(X)) X = as.data.frame(X)
  if(!is.matrix(Y)) Y = as.matrix(Y)

  n = dim(Y)[1]

  if (dim(X)[2] == 1) {
    Y = Y[order(X),,drop=F]

    node_calculator = function(j) {
      return(k(Y[j,],Y[j-1,]) + k(Y[j,],Y[j+1,]))
    }
    return((sum(sapply(2:(n-1), node_calculator))/2 + k(Y[1,],Y[2,]) + k(Y[n-1,],Y[n,]))/n)
  }

  out=emstreeR::ComputeMST(X,verbose = F)
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
#' Calculate the kernel partial correlation coefficient (KPC) with directed Knn graph or minimum spanning tree.
#'
#' The kernel partial correlation squared (KPC) measures the conditional dependency
#' between \eqn{Y} and \eqn{Z} given \eqn{X}, based on an i.i.d. sample of \eqn{(Y, Z, X)}.
#' It converges to the population quantity which is between 0 and 1.
#' A small value indicates low conditional dependency between \eqn{Y} and \eqn{Z} given \eqn{X}, and
#' a high value indicates stronger conditional dependence.
#'
#' @param Y a matrix (n by dy)
#' @param X a matrix (n by dx)
#' @param Z a matrix (n by dz)
#' @param k a function \eqn{k(y, y')} of class \code{kernel}. It can be the kernel implemented in \code{kernlab} e.g. \code{rbfdot(sigma = 1)}, \code{vanilladot()}
#' @param Knn number of nearest neighbor to use; or "MST"
#' @param trans_inv TRUE or FALSE. Is \eqn{k(y, y)} free of \eqn{y}?
#'
#' @examples
#' library(kernlab)
#' data("med")
#' KPCgraph(med$D,med$C,med$U,rbfdot(1/(2*median(dist(med$D))^2)),trans_inv=T)
#'
#'
#' n = 2000
#' x = rnorm(n)
#' z = rnorm(n)
#' y = x + z + rnorm(n,1,1)
#' KPCgraph(y,x,z,vanilladot(),Knn=1,trans_inv=F)
#'
#' n = 1000
#' set.seed(1)
#' x = runif(n)
#' z = runif(n)
#' y = (x + z) %% 1
#' KPCgraph(y,x,z,rbfdot(5),Knn=1,trans_inv=T)
KPCgraph = function(Y,X,Z,k,Knn = 1,trans_inv=FALSE) {
  if(!is.matrix(Y)) Y = as.matrix(Y)

  Tn_XZ = TnKnn(Y,cbind(X,Z),k,Knn)
  Tn_X = TnKnn(Y,X,k,Knn)

  if (trans_inv) {
    return((Tn_XZ - Tn_X)/(k(Y[1,],Y[1,])-Tn_X))
  }
  else {
    d = function(j) {
      return(k(Y[j,],Y[j,]))
    }

    return((Tn_XZ - Tn_X)/(mean(sapply(1:n, d))-Tn_X))
  }
}

# Double-centering
#
# Double-centering a squared matrix
#
# Given a square matrix A, the function return HAH, where H = diag(n) - 1/n is the centering matrix.
#
# @param M Matrix (n by n)
double_center = function(M){
  return(M - rowMeans(M) - matrix(colMeans(M),nrow(M),ncol(M),byrow = T) + mean(M[]))
}

#' Kernel partial correlation with CME method
#'
#' Calculate Kernel partial correlation coefficient (KPC) with cconditional mean embedding method.
#'
#' The kernel partial correlation squared (KPC) measures the conditional dependency
#' between \eqn{Y} and \eqn{Z} given \eqn{X}, based on an i.i.d. sample of \eqn{(Y, Z, X)}.
#' It converges to the population quantity which is between 0 and 1.
#' A small value indicates low conditional dependency between \eqn{Y} and \eqn{Z} given \eqn{X}, and
#' a large value indicates stronger conditional dependence.
#' If \code{X = NULL}, it measures the unconditional dependency between \eqn{Y} and \eqn{Z}.
#'
#' @param Y a matrix (n by dy)
#' @param X a matrix (n by dx) or \code{NULL} if \eqn{X} is empty
#' @param Z a matrix (n by dz)
#' @param ky a function \eqn{k(y, y')} of class \code{kernel}. It can be the kernel implemented in \code{kernlab} e.g. \code{rbfdot(sigma = 1)}, \code{vanilladot()}
#' @param kx the kernel function for \eqn{X}
#' @param kxz the kernel function for \eqn{(X, Z)} or for \eqn{Z} if \eqn{X} is empty
#' @param eps a small positive regularization parameter for inverting the empirical cross-covariance operator
#' @param appro whether to use incomplete Cholesky decomposition for approximation
#' @param tol tolerance used for incomplete Cholesky decomposition (\code{inchol} in package \code{kernlab})
#'
#' @examples
#' library(kernlab)
#' data("med")
#' KPCCME(med$D,med$C,med$U,rbfdot(1/(2*median(dist(med$D))^2)),rbfdot(1/(2*median(dist(med$C))^2)), rbfdot(1/(2*median(dist(cbind(med$C,med$U)))^2)), 1e-3,F) # 0.0003175344 D is independent of U given C
#' KPCCME(med$D,med$U,med$C,rbfdot(1/(2*median(dist(med$D))^2)),rbfdot(1/(2*median(dist(med$U))^2)), rbfdot(1/(2*median(dist(cbind(med$C,med$U)))^2)), 1e-3,F) # 0.6834605 D is associated with C controlling U
#' n = 2000
#' x = runif(n)
#' z = runif(n)
#' y = (x + z) %% 1
#' KPCCME(y, x, z, rbfdot(5), rbfdot(5), rbfdot(2), 1e-3/n^(0.49), appro = F)
#' KPCCME(y, x, z, rbfdot(5), rbfdot(5), rbfdot(2), 1e-4/n^(0.4), appro = F)
#' KPCCME(y, x, z, rbfdot(5), rbfdot(5), rbfdot(2), 1e-3/n^(0.49), appro = T, tol = 1e-5)
KPCCME = function(Y, X = NULL, Z, ky, kx, kxz, eps, appro = FALSE, tol = 1e-3) {
  if(!is.matrix(Y)) Y = as.matrix(Y)
  if(!is.matrix(Z)) Z = as.matrix(Z)
  if(!is.null(X) & !is.matrix(X)) X = as.matrix(X)

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

#' Kernel partial correlation with CME method using linear kernels
#'
#' Linear kernels are used for ky, kx, kxz in KPCCME.
#'
#' Linear kernels are used for ky, kx, kxz in KPCCME.
#' The kernel partial correlation squared (KPC) measures the conditional dependency
#' between \eqn{Y} and \eqn{Z} given \eqn{X}, based on an i.i.d. sample of \eqn{(Y, Z, X)}.
#' It converges to the population quantity which is between 0 and 1.
#' A small value indicates low conditional dependency between \eqn{Y} and \eqn{Z} given \eqn{X}, and
#' a large value indicates stronger conditional dependence.
#' If \code{X = NULL}, it measures the unconditional dependency between \eqn{Y} and \eqn{Z}.
#'
#' @param Y a matrix (n by dy)
#' @param Z a matrix (n by dz)
#' @param X a matrix (n by dx) or \code{NULL} if \eqn{X} is empty
#' @param eps a small regularization parameter for inverting the empirical cross-covariance operator
#' @examples
#' n = 2000
#' x = rnorm(n)
#' z = rnorm(n)
#' y = x + z + rnorm(n,1,1)
#' KPCCMElinear(y, x, z, 1e-5/n^(0.4))
KPCCMElinear = function(Y, X = NULL, Z, eps) {
  if(!is.matrix(Y)) Y = as.matrix(Y)
  if(!is.matrix(Z)) Z = as.matrix(Z)
  if(!is.null(X) & !is.matrix(X)) X = as.matrix(X)

  n = dim(Y)[1]

  # Incomplete Cholesky decomposition is exact here
  if (is.null(X)) {
    Lz = Z
    Lz = Lz - rep(colMeans(Lz), rep.int(n, ncol(Lz)))
    Ly = Y
    return(sum((t(Ly)%*%Lz%*%solve(dim(Y)[1]*eps*diag(dim(Lz)[2]) + t(Lz)%*%Lz)%*%t(Lz))^2)/sum(diag(double_center(kernlab::kernelMatrix(ky,Y)))))
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
#' Variable selection with KPC using directed Knn graph or minimum spanning tree
#'
#' A stepwise forward selection of variables using KPC. At each step the \eqn{Xj} maximizing \eqn{\hat{rho^2}(Y,Xj | selected Xi)} is selected.
#'
#' @param X a matrix of predictors (n by dx)
#' @param Y a matrix of responses (n by dy)
#' @param num_features the number of variables to be selected, cannot be larger than dx. The default value is NULL and in that
#'   case it will be set equal to dx. If stop == TRUE (see below), then num_features is the maximal number of variables to be selected.
#' @param stop Stops at the first instance of negative Tn, if TRUE.
#' @param numCores number of cores that are going to be used for parallelizing the process.
#' @param k a function \eqn{k(y, y')} of class \code{kernel} on the space of \eqn{Y}. It can be the kernel implemented in \code{kernlab} e.g. \code{rbfdot(sigma = 1)}, \code{vanilladot()}
#' @param Knn the number of nearest neighbor; or "MST"
#' @param verbose whether to print each selected variables during the forward stepwise algorithm
#' @return a vector of the indices from 1,...,dx of the selected variables
#' @examples
#' n = 200
#' p = 100
#' X = matrix(rnorm(n * p), ncol = p)
#' Y = X[, 1] * X[, 2] + sin(X[, 1] * X[, 3])
#' KFOCI(Y, X, kernlab::rbfdot(1), Knn=1, numCores = 1)
#' KFOCI(Y, X, kernlab::rbfdot(1), Knn=1, num_features = 2, numCores = 1)
#' colnames(X) = paste0(rep("X", p), seq(1, p))
#' FOCI::foci(Y, X, numCores = 1)$selectedVar$index
#'
#' surgical = olsrr::surgical
#' for (i in 1:9) surgical[,i] = (surgical[,i] - mean(surgical[,i]))/sd(surgical[,i])
#' colnames(surgical)[KFOCI(surgical[,9],surgical[,1:8],kernlab::rbfdot(1/(2*median(dist(surgical$y))^2)),Knn=1)]
#' #### "enzyme_test" "pindex" "liver_test"  "alc_heavy"
#' \dontrun{
#' # This example is computationally consuming. It may take 20 mins with 7 cores on a laptop.
#' n = 2000
#' p = 1000
#' set.seed(1)
#' X = matrix(rnorm(n * p), ncol = p)
#' Y = X[, 1] * X[, 2] + sin(X[, 1] * X[, 3]) + rnorm(n)
#' KFOCI(Y, X, kernlab::rbfdot(1), Knn=20, numCores = 7)
#' }
# code modified from Azadkia, M. and Chatterjee, S. (2019). A simple measure of conditional dependence.
KFOCI <- function(Y, X, k, Knn, num_features = NULL, stop = TRUE, numCores = 1, verbose = F){
  if(!is.matrix(X)) {
    X = as.matrix(X)
  }
  if(!is.matrix(Y)) {
    Y = as.matrix(Y)
  }

  if (is.null(num_features)) num_features = dim(X)[2]
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

#' Calculating the mean squared error
#'
#' Return the mean squared error between \eqn{k(Y_i,\cdot)} and the empirical CME \eqn{\hat{\mu}_{Y|Xi}}.
#'
#' The MSE provides a way to select the regularization parameter \code{eps} in the CME estimator.
#' \code{eps} should be chosen as small as possible while keeping the MSE reasonably small.
#'
#' @param X a matrix of predictors (n by dx)
#' @param Y a matrix of responses (n by dy)
#' @param ky a function \eqn{k(y, y')} of class \code{kernel} on the space of \eqn{Y}. It can be the kernel implemented in \code{kernlab} e.g. \code{rbfdot(sigma = 1)}, \code{vanilladot()}
#' @param kx the kernel function for \eqn{X}
#' @param eps the regularization parameter for the CME estimator
#' @examples
#' n = 200
#' p = 100
#' X = matrix(rnorm(n * p), ncol = p)
#' Y = X[, 1] * X[, 2] + sin(X[, 1] * X[, 3]) + rnorm(n)
#' for (eps in 10^(-(1:9))) print(MSE(X[,1],Y,kernlab::rbfdot(1),kernlab::rbfdot(1),eps))
#' for (eps in 10^(-(1:9))) print(MSE(X[,1:2],Y,kernlab::rbfdot(1/2),kernlab::rbfdot(1),eps))
#' for (eps in 10^(-(1:9))) print(MSE(X[,1:3],Y,kernlab::rbfdot(1/3),kernlab::rbfdot(1),eps))

MSE = function(X,Y,kx,ky,eps) {
  if(!is.matrix(X)) X = as.matrix(X)
  if(!is.matrix(Y)) Y = as.matrix(Y)
  n = dim(X)[1]
  A = solve(double_center(kernelMatrix(kx,X)) + n*eps*diag(n))
  return(sum(double_center(kernelMatrix(ky,Y)) * base::tcrossprod(A))*n*eps^2)
}


#' Variable selection with CEM estimator
#'
#' The algorithm performs a forward stepwise variable selection using CME estimators.
#'
#' A stepwise forward selection of variables using KPC. At each step the \eqn{Xj} maximizing \eqn{\tilde{rho^2}(Y,Xj | selected Xi)} is selected.
#'
#' @param X a matrix of predictors (n by dx)
#' @param Y a matrix of responses (n by dy)
#' @param num_features the number of variables to be selected, cannot be larger than dx.
#' @param numCores number of cores that are going to be used for parallelizing the process.
#' @param ky a function \eqn{k(y, y')} of class \code{kernel} on the space of \eqn{Y}. It can be the kernel implemented in \code{kernlab} e.g. \code{rbfdot(sigma = 1)}, \code{vanilladot()}
#' @param kx a list of length \code{num_features}; \code{kx[[k]]} is the kernel used for \code{(Xj1,...,Xjk)}, the first \code{k} selected variables.
#' @param eps a positive number; the regularization parameter for CME estimator
#' @param appro whether to use incomplete Cholesky decomposition for approximation
#' @param tol tolerance used for incomplete Cholesky decomposition (\code{inchol} in package \code{kernlab})
#' @param verbose whether to print each selected variables during the forward stepwise algorithm
#' @return a vector of the indices from \code{1,...,dx} of the selected variables
#' @examples
#' n = 200
#' p = 100
#' X = matrix(rnorm(n * p), ncol = p)
#' Y = X[, 1] * X[, 2] + sin(X[, 1] * X[, 3]) + rnorm(n)*0.5
#' kx = c(kernlab::rbfdot(1),kernlab::rbfdot(1/2),kernlab::rbfdot(1/3))
#' CME_select(Y, X, rbfdot(1), kx, 3, eps = 1e-3, appro = F, numCores = 1)
# code modified from Azadkia, M. and Chatterjee, S. (2019). A simple measure of conditional dependence.
CME_select <- function(Y, X, ky, kx, num_features, eps, appro = F, tol = 1e-3, numCores = 1, verbose = F){
  if(!is.matrix(X)) X = as.matrix(X)
  if(!is.matrix(Y)) Y = as.matrix(Y)

  n = dim(Y)[1]
  p = ncol(X)
  Q = rep(0, num_features) # stores the values of Tn
  index_select = rep(0, num_features)
  # select the first variable
  estimateQFixedY <- function(id){
    return(KPCCME_numerator(Y,NULL,X[,id],ky,NULL,kx[[1]],eps,appro,tol))
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
      return(KPCCME_numerator(Y, X[,index_select[1:count]], X[, c(index_select[1:count], id)], ky, kx[[count]], kx[[count+1]], eps,appro,tol))
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


# calculate the numerator of the CME estimator
# used for stepwise variable selection
KPCCME_numerator = function(Y, X = NULL, Z, ky, kx, kxz, eps, appro = FALSE, tol = 1e-3) {
  if(!is.matrix(Y)) Y = as.matrix(Y)
  if(!is.matrix(Z)) Z = as.matrix(Z)
  if(!is.null(X) & !is.matrix(X)) X = as.matrix(X)

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

