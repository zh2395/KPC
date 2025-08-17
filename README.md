# KPC
Kernel partial correlation (KPC) coefficient measures the strength of conditional association between Y and Z given X, with X, Y, Z being random variables taking values in general topological spaces. As the name suggests, KPC is defined in terms of kernels on reproducing kernel Hilbert spaces (RKHSs). The population KPC is a deterministic number between 0 and 1; it is 0 if and only if Y is conditionally independent of Z given X, and it is 1 if and only if Y is a measurable function of Z and X. This R package provides implementations of two empirical versions of KPC when X, Z are Euclidean, and Y, possibly taking values in a general topological space, can be stored as a vector.

One empirical KPC estimator is based on geometric graphs, such as K-nearest neighbor graphs (K-NN) and minimum spanning trees (MST), and is consistent under very weak conditions. The other empirical estimator, defined using conditional mean embeddings (CMEs) as used in the RKHS literature, is also consistent under suitable conditions. Users are allowed to define arbitrary kernels (see examples below).

Using KPC we also provide a stepwise forward variable selection algorithm KFOCI (using the graph based estimator of KPC), as well as a similar stepwise forward selection algorithm based on the RKHS based estimator. For more details on KPC, its empirical estimators and its application on variable selection, see https://www.jmlr.org/papers/v23/21-493.html.

When X is empty, KPC measures the unconditional dependence between Y and Z, which is also implemented in the package. The unconditional graph-based estimator has been described in *Deb, N., P. Ghosal, and B. Sen (2020). Measuring association on topological spaces using kernels and geometric graphs*. It is implemented in the functions `KMAc` and `Klin` in this package. The latter can be computed in near linear time. 

## Installation

This package depends on R (>= 4.0.0). You can install the package KPC by (the package `devtools` needs to be installed first):

``` r
devtools::install_github("zh2395/KPC")
```

You can uninstall the package by:
```r
remove.packages("KPC")
```

## Usage of the functions
Here we briefly introduce the functions in this package.
See the documentation (help page) of the R package for more details.

`KPCgraph` implements the KPC estimator based on geometric graphs.
The inputs are:
`Y`: a matrix of n rows;
`X`: a matrix of n rows, or NULL if X is empty, in which case it will return `KMAc(Y,Z,k,Knn)`, which measures the unconditional dependence between Y and Z.
`Z`: a matrix of n rows;
`k`: a function of class kernel. It can be the kernel implemented in `kernlab` e.g., Gaussian kernel `rbfdot(sigma = 1)`, linear kernel `vanilladot()`; In practice, Gaussian kernel with empirical bandwidth `kernlab::rbfdot(1/(2*stats::median(stats::dist(Y))^2))` may be a good choice.
`Knn`: a positive integer indicating the number of nearest neighbor to use; or "MST". A small Knn (e.g., Knn=1) is recommended for an accurate estimate of the population KPC.
`trans_inv`: whether k(y, y) is free of y.

``` r
library(kernlab)
library(KPC)
n = 1000
set.seed(1)
x = runif(n)
z = runif(n)
y = (x + z) %% 1
KPCgraph(Y = y, X = x, Z = z, k = rbfdot(5), Knn = 1, trans_inv = T)
# 0.9725613
# Theoretical KPC is 1 since y is a measurable function of x and z
```

`KPCRKHS` implements the KPC estimator based on RKHS method using CME formula.
The inputs are `Y`: a matrix of n rows; `X`: a matrix of n rows, or NULL if X is empty, in which case the coefficient
measures the unconditional association between Y and Z; `Z`: a matrix of n rows; `ky`, `kx`, `kxz`: the kernels used for the space of Y, X, (X,Z) respectively;
`eps`: a small positive regularization parameter for inverting the empirical cross-covariance operator;
`appro`: whether to use incomplete Cholesky decomposition for approximation;
`tol`: tolerance used for incomplete Cholesky decomposition (implemented by the function `inchol` in the package `kernlab`).

``` r
n = 1000
set.seed(1)
x = rnorm(n)
z = rnorm(n)
y = x + z + rnorm(n,1,1)
library(kernlab)
k = vanilladot() # linear kernel
KPCRKHS(y, x, z, k, k, k, 1e-3/n^(0.4), appro = F)
# 0.5158324 (Population quantity = 0.5)
KPCRKHS(y, x, z, k, k, k, 1e-3/n^(0.4), appro = T, tol = 1e-5)
# 0.5158324 (Population quantity = 0.5)
```

`KFOCI` implements variable selection with KPC using directed K-NN graph or MST.
The inputs are 
`Y` : a matrix of responses (n by dy);
`X`: a matrix of predictors (n by dx); 
`Z`: an integer vector of column indices in `X` to pre-condition on. These variables are always included in the conditioning set and are not re-selected. Formally, the goal is then to find $S \subset \lbrace 1, \dotsc, dx\rbrace\setminus Z$ such that $Y \perp X_{S^c}\mid (X_Z, X_S)$. The default `NULL` corresponds to no pre-conditioning;
`k`: the kernel function used for Y;
`Knn`: a positive integer indicating the number of nearest neighbor; or "MST". The suggested choice of Knn is 0.05n for samples up to a few hundred observations. For large n, the suggested Knn is sublinear in n. That is, it may grow slower than any linear function of n. The computing time is approximately linear in Knn. A smaller Knn takes less time.
`num_features`: the number of variables to be selected from the non-pre-conditioned set (which cannot be larger than $dx - |Z|$). The default value of `num_features` is `NULL` and in that
case it will be set equal to $dx - |Z|$. If `stop == TRUE` (see below), then `num_features` is the maximal number of variables to be selected (selection may stop earlier);
`stop`: If `stop == TRUE`, then the automatic stopping criterion (stops at the first instance of negative Tn, as mentioned in the paper) will be implemented and continued till `num_features` many variables are selected. If `stop == FALSE` then exactly `num_features` many variables are selected; 
`numCores`: number of cores that are going to be used for parallelizing the process;
`verbose`: whether to print each selected variables during the forward stepwise algorithm (default `FALSE`);

It is suggested to normalize the predictors before applying `KFOCI`. `KFOCI` returns a vector of the indices, from 1,...,dx, from the non-pre-conditioned set of the selected variables in the same order that they were selected. The variables at the front are expected to be more informative in predicting Y.

``` r
n = 200
p = 100
set.seed(1)
X = matrix(rnorm(n * p), ncol = p)
Y = X[, 1] * X[, 2] + sin(X[, 1] * X[, 3])
KFOCI(Y, X, k=kernlab::rbfdot(1), Knn=1, numCores = 1)
# 1 2 3

# Install package 'olsrr'
surgical = olsrr::surgical
for (i in 1:9) surgical[,i] = (surgical[,i] - mean(surgical[,i]))/sd(surgical[,i])
colnames(surgical)[KFOCI(surgical[,9],surgical[,1:8],k=kernlab::rbfdot(1/(2*median(dist(surgical$y))^2)),Knn=1)]
# "enzyme_test" "pindex" "liver_test"  "alc_heavy"
```


`KPCRKHS_VS` performs a forward stepwise variable selection using the RKHS based estimator of KPC.
One needs to pre-specify the number of variables to be selected.
The inputs are:
`Y`: a matrix of responses (n by dy);
`X`: a matrix of predictors (n by dx);
`num_features`: the number of variables to be selected, cannot be larger than dx;
`ky`: the kernel function for Y;
`kS`: a function that takes X and a subset of indices S as inputs, and then outputs the kernel for X_S. The first argument of kS is X, and the second argument is a vector of positive integer. If `kS == NULL`, Gaussian kernel with empitical bandwidth will be used, i.e., `kernlab::rbfdot(1/(2*stats::median(stats::dist(X[,S]))^2))`;
`eps`: a positive number, the regularization parameter for RKHS based KPC estimator;
`appro` whether to use incomplete Cholesky decomposition for approximation;
`tol`: tolerance used for incomplete Cholesky decomposition (implemented by `inchol` in package `kernlab`);
`numCores`: number of cores that are going to be used for parallelizing the process;
`verbose`: whether to print each selected variables during the forward stepwise algorithm

The algorithm returns a vector of the indices from 1,...,dx of the selected variables in the same order that they were selected. The variables at the front are expected to be more informative in predicting Y.
``` r
n = 200
p = 100
set.seed(1)
X = matrix(rnorm(n * p), ncol = p)
Y = X[, 1] * X[, 2] + sin(X[, 1] * X[, 3])
library(kernlab)
kS = function(X,S) return(rbfdot(1/length(S)))
KPCRKHS_VS(Y, X, num_features = 3, rbfdot(1), kS, eps = 1e-3, appro = FALSE, numCores = 1)
# 1 2 3
kS = function(X,S) return(rbfdot(1/(2*stats::median(stats::dist(X[,S]))^2)))
KPCRKHS_VS(Y, X, num_features = 3, rbfdot(1), kS, eps = 1e-3, appro = FALSE, numCores = 1)
# 1 2 3
```

## Real data example
The medical data included in the package are collected from *Edwards, D. (2012). Introduction to graphical modelling, Section 3.1.4, Springer Science & Business Media*.
``` r
# load medical data
data("med")
for (i in 1:3) med[,i] = (med[,i] - mean(med[,i]))/sd(med[,i]) # normalization
library(kernlab)
KPCgraph(med$D,med$C,med$U,rbfdot(1/(2*median(dist(med$D))^2)),Knn=1,trans_inv=T)
# 0.04069334
# Theoretical KPC is 0 since D is independent of U given C
set.seed(1) # There is randomness in breaking the ties
KPCgraph(med$D,med$U,med$C,rbfdot(1/(2*median(dist(med$D))^2)),Knn=1,trans_inv=T)
# 0.3255831 
# D is associated with C controlling U

# RKHS estimator
KPCRKHS(med$D,med$C,med$U,rbfdot(1/(2*median(dist(med$D))^2)),rbfdot(1/(2*median(dist(med$C))^2)), rbfdot(1/(2*median(dist(cbind(med$C,med$U)))^2)), 1e-2,F)
# 0.1502744
KPCRKHS(med$D,med$U,med$C,rbfdot(1/(2*median(dist(med$D))^2)),rbfdot(1/(2*median(dist(med$U))^2)), rbfdot(1/(2*median(dist(cbind(med$C,med$U)))^2)), 1e-2,F)
# 0.3852009
```


## Case study on general spaces
KPC can be defined for variables X, Y, Z taking values on general spaces, and if the variables can be stored as vectors, then the above functions can also be used.
Here we consider the special orthogonal group SO(3), which consists of 3 by 3 orthogonal matrices with determinant 1.
We store the matrices by concatenating the columns, and define a kernel on the vectorized SO(3) by
``` r
SO3ker = function(vecA,vecB){
  A = matrix(vecA,3,3)
  B = matrix(vecB,3,3)
  theta = base::acos(min(1,max((sum(diag(t(B)%*%A))-1)/2,0)))
  if (theta == 0 | theta == pi){
    return(pi)
  }
  return(theta*(pi-theta)/sin(theta))
}
class(SO3ker) <- "kernel"
```
Let `R1`, `R3` be the rotation matrices around x-axis and z-axis.
``` r
R1 = function(x) {
  cos_x = cos(x)
  sin_x = sin(x)
  return(matrix(c(1,0,0,0,cos_x,sin_x,0,-sin_x,cos_x),3,3))
}
R3 = function(z) {
  cos_z = cos(z)
  sin_z = sin(z)
  return(matrix(c(cos_z,sin_z,0,-sin_z,cos_z,0,0,0,1),3,3))
}

n = 1000
set.seed(1)
x = rnorm(n)
z = rnorm(n)
y1 = y2 = matrix(0,n,9)
for (i in 1:n) {
  y1[i,] = as.numeric(R1(x[i])%*%R3(z[i]))
  y2[i,] = as.numeric(R1(x[i])%*%R3(rnorm(1)))
}
KPCgraph(y1,x,z,SO3ker,Knn = 1,trans_inv=T)
# 0.8547098
# y1 is a function of x and z
KPCgraph(y2,x,z,SO3ker,Knn = 1,trans_inv=T)
# 0.00914022
# y2 is conditionally independent of z given x

KPCRKHS(y1, x, z, SO3ker, rbfdot(1), rbfdot(0.5), 1e-5, appro = F)
# 0.6198004
KPCRKHS(y2, x, z, SO3ker, rbfdot(1), rbfdot(0.5), 1e-5, appro = F)
# 0.05227157

# Variable selection
n = 100
p = 500
set.seed(1)
X = matrix(rnorm(n * p), ncol = p)
y = matrix(0,n,9)
for (i in 1:n) {
  y[i,] = as.numeric(R1(X[i,1])%*%R3(X[i,2]))
}
KFOCI(y, X, k=SO3ker, Knn=1, numCores = 1)
# 2 1
```

The 2017 Korea presidential election data, collected from https://github.com/OhmyNews/2017-Election,
consists of the voting results earned by the top five candidates from 250 electoral districts in Korea.
The top three candidates from three major parties representing progressivism, conservatism and centrism earned most of the votes,
so we will focus on the proportion of votes earned by each of these three candidates among them (Y), which can be viewed as a histogram-valued response.
The demographic information average age (X1), average years of education (X2), average housing price per square meter (X3) and average paid national health
insurance premium (X4) are available for each electoral district.
``` r
data(ElecData)
n = dim(ElecData)[1]/5
Y = matrix(0, n, 3)
# Each row of Y is the proportion of votes earned by each of the top three candidates among them
for (i in 1:n) {
  num_vote = ElecData$NumVote[(1+(i-1)*5):(3+(i-1)*5)]
  Y[i,] = num_vote/sum(num_vote)
}
X = ElecData[5*(1:n),4:7]
for (i in 1:4) X[,i] = (X[,i] - mean(X[,i]))/sd(X[,i]) # normalize the data
set.seed(1)
X[,5:8] = rnorm(n*4)

library(kernlab)
KFOCI(Y, X, k=rbfdot(1/(2*median(dist(Y))^2)), Knn=1, numCores = 1)
# 1 2 3 4

# define two kernels on histograms
k1 = function(a,b) return(1/prod(a+b+1))
k2 = function(a,b) return(exp(-sum(sqrt(a+b))))
class(k1) = class(k2) = "kernel"

KFOCI(Y, X, k=k1, Knn=4, numCores = 1)
# 1 2 4 3
KFOCI(Y, X, k=k2, Knn=4, numCores = 1)
# 1 2 4 3

KPCgraph(Y,X[,c(2,3,4)],X[,1],rbfdot(1/(2*median(dist(Y))^2)),Knn = 2,trans_inv=TRUE)
# 0.1543532
KPCRKHS(Y,X[,c(2,3,4)],X[,1],rbfdot(1/(2*median(dist(Y))^2)),rbfdot(1/(2*median(dist(X[,c(2,3,4)]))^2)), rbfdot(1/(2*median(dist(X[,1:4]))^2)), eps=1e-4, appro=F)
# 0.1473899

KPCgraph(Y,X[,c(1,2,4)],X[,3],rbfdot(1/(2*median(dist(Y))^2)),Knn = 2,trans_inv=TRUE)
# 0.05542749
KPCRKHS(Y,X[,c(1,2,4)],X[,3],rbfdot(1/(2*median(dist(Y))^2)),rbfdot(1/(2*median(dist(X[,c(1,2,4)]))^2)), rbfdot(1/(2*median(dist(X[,1:4]))^2)), eps=1e-4, appro=F)
# 0.06338199
# X3 and X4 are both measures of richness. The conditional association between X3 and Y given all other variables is weaker than that of X1 and Y.
```








