#' Tn with geometric graphs
#'
#' Calculate Tn using directed Knn graph or minimum spanning tree.
#'
#' Tn is an estimate of E[E[k(Y1,Y1')|X]], with Y1, Y1' drawn iid from Y|X, given X.
#' Kernels provided in the package "kernlab" are recommended, e.g. rbfdot(sigma), vanilladot().
#' For continuous data, "breaktie" can be set to FALSE for faster computation.
#'
#' @param X Matrix of predictors (n by dx)
#' @param Y Matrix of response (n by dy)
#' @param k a kernel function k(y,y')
#' @param Knn The number of K-nearest neighbor to use; or "MST".
#' @return Tn: a numeric number
TnKnn = function(Y,X,k,Knn=1) {
  if (Knn == "MST") {
    if(!is.data.frame(X)) {
      X = as.data.frame(X)
    }
    if(!is.matrix(Y)) {
      Y = as.matrix(Y)
    }
    n = dim(Y)[1]

    if (dim(X)[2] == 1) {
      Y = Y[order(X),]
      if(!is.matrix(Y)) {
        Y = as.matrix(Y)
      }
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

  id <- group <- NULL
  if(!is.matrix(X)) {
    X = as.matrix(X)
  }
  if(!is.matrix(Y)) {
    Y = as.matrix(Y)
  }

  n = dim(Y)[1]

  # compute the nearest neighbor of X
  # need to check Knn < p
  nn_X = RANN::nn2(X, query = X, k = Knn + 2)
  nn_index_X = nn_X$nn.idx[, 2:(Knn+1)]
  if(!is.matrix(nn_index_X)) {
    nn_index_X = as.matrix(nn_index_X)
  }

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
          nn_index_X[df_X$id[i],] = setdiff(group_indices, df_X$id[i]) # avoid length 1 in sample function
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

  # faster calculation for Gaussian kernel and linear kernel
  if (class(k) == "rbfkernel") {
    sig = kernlab::kpar(k)$sigma
    if (Knn == 1) {
      node_calculator = function(j) {
        return(exp(-sig * sum((Y[j,]-Y[nn_index_X[j,],])^2)))
      }
    }
    else {
      node_calculator = function(j) {
        return(sum(exp(-sig * (outer(sum(Y[j,]^2),base::rowSums(as.matrix(Y[nn_index_X[j,],]^2)),"+") - 2*base::tcrossprod(t(Y[j,]),as.matrix(Y[nn_index_X[j,],]))))))
      }
    }
  }
  else if (class(k) == "vanilladot") {
    if (Knn == 1){
      node_calculator = function(j) {
        return(sum(Y[j,]*Y[nn_index_X[j,],]))
      }
    }
    else {
      node_calculator = function(j) {
        return(sum(base::tcrossprod(t(Y[j,]),Y[nn_index_X[j,],])))
      }
    }
  }
  else {
    # calculation for general kernel
    node_calculator = function(j) {
      sum = 0
      for (i in 1:Knn) {
        sum = sum + k(Y[j,],Y[nn_index_X[j,i],])
      }
      return(sum)
    }
  }

  return(mean(sapply(1:n, node_calculator))/Knn)
}

#' KPC (with Knn graph)
#'
#' Calculate the kernel partial correlation squared (KPC) with directed Knn graph
#'
#' The kernel partial correlation squared (KPC) measures the conditional dependency
#' between Y and Z given X, based on an i.i.d. sample of (Y, Z, X).
#' The coefficient is asymptotically guaranteed to be between 0 and 1.
#'
#' @param Y Matrix (n by dy)
#' @param Z Matrix (n by dz)
#' @param X Matrix (n by dx)
#' @param k a kernel function k(y,y')
#' @param Knn number of nearest neighbor to use; or "MST"
#' @param trans_inv TRUE or FALSE. Is kernel k translation invariant?
#' @details The value returned by KPC can be positive or negative. Asymptotically, it is guaranteed
#'   to be between 0 and 1. A small value indicates low conditional dependence between Y and Z given X, and
#'   a high value indicates strong conditional dependence.
#'
#' @examples
#' data("medical")
#' KPCKnn(med$D,med$C,med$U,kernlab::rbfdot(1/(2*median(dist(med$D))^2)),trans_inv = T)
KPCKnn = function(Y,X,Z,k,Knn = 1,trans_inv=FALSE) {
  if(!is.matrix(Y)) {
    Y = as.matrix(Y)
  }
  Tn_XZ = TnKnn(Y,cbind(X,Z),k,Knn)
  Tn_X = TnKnn(Y,X,k,Knn)

  if (trans_inv) {
    return((Tn_XZ - Tn_X)/(k(Y[1,],Y[1,])-Tn_X))
  }
  else {
    d = function(j) {
      return(k(Y[j,],Y[j,]))
    }

    return(Tn_XZ - Tn_X)/(mean(sapply(1:n, d))-Tn_X)
  }
}

#' Double-centering
#'
#' Double-centering a squared matrix
#'
#' Given a square matrix A, the function return HAH, where H = diag(n) - 1/n is the centering matrix.
#'
#' @param M Matrix (n by n)
double_center = function(M){
  return(M - rowMeans(M) - matrix(colMeans(M),nrow(M),ncol(M),byrow = T) + mean(M[]))
}

#' KPCCME
#'
#' Calculate KPC rho^2(Y,Z|X) with CME method
#'
#' kx,ky,kxz should be kernels provided in the R package "kernlab", so that kernelMatrix can be applied
#'
#' @param Y Matrix (n by dy)
#' @param Z Matrix (n by dz)
#' @param X Matrix (n by dx) or NULL if X is empty
#' @param ky a kernel function k(y,y')
#' @param kx kernel function for X
#' @param kxz kernel function for (X,Z) or for Z if X is empty
#' @param eps A small regularization parameter for inverting the empirical cross-covariance operator
#' @param appro whether to use incomplete Cholesky decomposition for approximation
#' @param tol Tolerance used for incomplete Cholesky decomposition
#'
#' @examples
#' data("medical")
#' KPCCME(med$D,med$C,med$U,rbfdot(1/(2*median(dist(med$D))^2)),rbfdot(1/(2*median(dist(med$C))^2)), rbfdot(1/(2*median(dist(cbind(med$C,med$U)))^2)), 1e-3,F) # 0.0003175344 D is independent of U given C
#' KPCCME(med$D,med$U,med$C,rbfdot(1/(2*median(dist(med$D))^2)),rbfdot(1/(2*median(dist(med$U))^2)), rbfdot(1/(2*median(dist(cbind(med$C,med$U)))^2)), 1e-3,F) # 0.6834605 D is associated with C controlling U
KPCCME = function(Y, X = NULL, Z, ky, kx, kxz, eps, appro = TRUE, tol = 1e-3) {
  if(!is.matrix(Y)) {
    Y = as.matrix(Y)
  }
  if(!is.matrix(Z)) {
    Z = as.matrix(Z)
  }
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
  X = as.matrix(X)
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


#' Variable selection with KPC
#'
#' Variable selection with KPC using directed Knn graph
#'
#' A stepwise forward selection of variables using KPC. At each step the Xj maximizing rho^2(Y,Xj | selected Xi) is selected.
#'
#' @param X Matrix of predictors (n by dx)
#' @param Y Vector of responses (n by dy)
#' @param num_features Number of variables to be selected, cannot be larger than dx. The default value is NULL and in that
#'   case it will be set equal to dx. If stop == TRUE (see below), then num_features is irrelevant.
#' @param stop Stops at the first instance of negative Tn, if TRUE.
#' @param numCores number of cores that are going to be used for parallelizing the process.
#' @param k the kernel on Y: k(y,y')
#' @param Knn the number of nearest neighbor; or "MST"
#' @return A vector of the indices from 1,...,dx of the selected variables
#' @examples
#' n = 200
#' p = 100
#' set.seed(1)
#' X = matrix(rnorm(n * p), ncol = p)
#' Y = X[, 1] * X[, 2] + sin(X[, 1] * X[, 3]) + rnorm(n)*0.5
#' Varselect_Knn(Y, X, kernlab::rbfdot(1), Knn=10)
## It selects X1, X2, X3. FOCI selects X81
# colnames(X) = paste0(rep("X", p), seq(1, p))
# FOCI::foci(Y, X, numCores = 1)$selectedVar$index
#
# Example:
#surgical = olsrr::surgical
#for (i in 1:9) {
#  surgical[,i] = (surgical[,i] - mean(surgical[,i]))/sd(surgical[,i])
#}
#colnames(surgical)[Varselect_Knn(surgical[,9],surgical[,1:8],kernlab::rbfdot(1/(2*median(dist(surgical$y))^2)),1)]
#### "enzyme_test" "pindex" "liver_test"  "alc_heavy"
# code modified from Azadkia, M. and Chatterjee, S. (2019). A simple measure of conditional dependence.
Varselect_Knn <- function(Y, X, k, Knn, num_features = NULL, stop = TRUE, numCores = 1){
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
  }

  return(index_select[1:count])
}

