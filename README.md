# KPC
Kernel partial correlation coefficient (KPC) measures the strength of conditional association between Y and Z given X,
with X, Y, Z being random variables in topological spaces.
The population KPC is a deterministic number between 0 and 1.
It is 0 if and only if Y is conditionally independent of Z given X. It is 1 if and only if Y is a measurable function of Z and X.
This R package provides implementations of two empirical versions of KPC when X, Y, Z can be stored in vectors.
One is based on the geometric graph such as K-nearest neighbor graph (KNN) and minimum spanning tree (MST), and is consistent under very weak conditions.
The other is based on the conditional mean embedding (CME) formula in the kernel literature, which is also consistent under suitable conditions.
Users are free to define arbitrary kernels.
A stepwise forward variable selection algorithm KFOCI is given, as well as a similar stepwise forward selection algorithm based on CME.
For more details on KPC, its empirical estimators and its application on variable selection, see (link to the paper).

## Installation

You can install the package KPC by:

``` r
devtools::install_github("zh2395/KPC")
```

## Usage of the functions

`KPCgraph` implements the KPC estimator based on geometric graphs.
The inputs are `Y`, `X`, `Z`: matrices of n rows; `k`: a function of class kernel. It can be the kernel implemented in `kernlab` e.g. `rbfdot(sigma = 1)`, `vanilladot()`;
`Knn` the number of nearest neighbor to use, or "MST"; `trans_inv`: whether k(y, y) is free of y.

``` r
library(kernlab)
n = 1000
set.seed(1)
x = runif(n)
z = runif(n)
y = (x + z) %% 1
KPCgraph(Y = y, X = x, Z = z, k = rbfdot(5), Knn = 1, trans_inv = T)
# 0.9725613
# Theoretical KPC is 1 since y is a measurable function of x and z

# load medical data
data("med")
KPCgraph(med$D,med$C,med$U,rbfdot(1/(2*median(dist(med$D))^2)),trans_inv=T)
# 0.009254471
# Theoretical KPC is 0 since D is independent of U given C
```

`KPCCME` implements the KPC estimator based on CME formula.
The inputs are `Y`, `Z`: matrices of n rows; `X`: a matrix of n rows, or NULL if X is empty, in which case the coefficient
measures the unconditional association between Y and Z. `ky`, `kx`, `kxz`: the kernels used for the space of Y, X, (X,Z) respectively;
`eps`: a small positive regularization parameter for inverting the empirical cross-covariance operator;
`appro`: whether to use incomplete Cholesky decomposition for approximation;
`tol`: tolerance used for incomplete Cholesky decomposition (implemented in `inchol` in the package `kernlab`).

``` r
library(kernlab)
n = 1000
set.seed(1)
x = runif(n)
z = runif(n)
y = (x + z) %% 1
KPCCME(Y = y, X = x, Z = z, ky = rbfdot(5), kx = rbfdot(5), kxz = rbfdot(2), eps = 1e-3/n^(0.49), appro = F)
# 0.6854751
KPCCME(y, x, z, rbfdot(5), rbfdot(5), rbfdot(2), 1e-3/n^(0.49), appro = T, tol = 1e-5)
# 0.6854615

# load medical data
data("med")
KPCCME(med$D,med$C,med$U,rbfdot(1/(2*median(dist(med$D))^2)),rbfdot(1/(2*median(dist(med$C))^2)), rbfdot(1/(2*median(dist(cbind(med$C,med$U)))^2)), 1e-3,F)
# 0.0003175344
# D is independent of U given C
KPCCME(med$D,med$U,med$C,rbfdot(1/(2*median(dist(med$D))^2)),rbfdot(1/(2*median(dist(med$U))^2)), rbfdot(1/(2*median(dist(cbind(med$C,med$U)))^2)), 1e-3,F)
# 0.6834605
# D is associated with C controlling U
```

`KFOCI` implements variable selection with KPC using directed Knn graph or minimum spanning tree.
The inputs are `X` a matrix of predictors (n by dx); `Y` a matrix of responses (n by dy);
`num_features` the number of variables to be selected, cannot be larger than dx. The default value is `NULL` and in that
case it will be set equal to dx. If `stop == TRUE`, then `num_features` is then num_features is the maximal number of variables to be selected;
`stop` whether to stops at the first instance of negative Tn;
`numCores` number of cores that are going to be used for parallelizing the process;
`k` the kernel function used for Y;
`Knn` the number of nearest neighbor, or "MST".
`KFOCI` returns a vector of the indices from 1,...,dx of the selected variables.

``` r
n = 200
p = 100
set.seed(1)
X = matrix(rnorm(n * p), ncol = p)
Y = X[, 1] * X[, 2] + sin(X[, 1] * X[, 3])
KFOCI(Y, X, kernlab::rbfdot(1), Knn=1, numCores = 1)
# 1 2 3
KFOCI(Y, X, kernlab::rbfdot(1), Knn=1, num_features = 2, numCores = 1)
# 1 2

surgical = olsrr::surgical
for (i in 1:9) surgical[,i] = (surgical[,i] - mean(surgical[,i]))/sd(surgical[,i])
colnames(surgical)[KFOCI(surgical[,9],surgical[,1:8],kernlab::rbfdot(1/(2*median(dist(surgical$y))^2)),Knn=1)]
# "enzyme_test" "pindex" "liver_test"  "alc_heavy"
```


`MSE` returns the mean squared error between k(Y_i,) and the empirical CME \hat{\mu}_{Y|Xi}.
It provides a way to select the regularization parameter `eps` in the CME estimator.
`eps` should be chosen as small as possible while keeping the MSE reasonably small.
The inputs are `X` a matrix of predictors (n by dx);
`Y` a matrix of responses (n by dy);
`ky` the kernel function for Y;
`kx` the kernel function for X;
`eps` the regularization parameter for the CME estimator.
``` r
n = 200
p = 100
set.seed(1)
X = matrix(rnorm(n * p), ncol = p)
Y = X[, 1] * X[, 2] + sin(X[, 1] * X[, 3]) + rnorm(n)
for (eps in 10^(-(1:9))) print(MSE(X[,1],Y,kernlab::rbfdot(1),kernlab::rbfdot(1),eps))
# 0.313 0.285 0.290 0.340 0.637 4.258 135.981 1730.174 11931.17
for (eps in 10^(-(1:9))) print(MSE(X[,1:2],Y,kernlab::rbfdot(1/2),kernlab::rbfdot(1),eps))
# 0.307 0.245 0.216 0.348 1.388 7.630 57.920 472.576 2061.338
for (eps in 10^(-(1:9))) print(MSE(X[,1:3],Y,kernlab::rbfdot(1/3),kernlab::rbfdot(1),eps))
# 0.305 0.254 0.251 0.276 0.540 6.499 55.720 321.408 694.576
```

`CME_select` performs a forward stepwise variable selection using CME estimators.
One needs to pre-specify the number of variables to be selected.
The inputs are `X` a matrix of predictors (n by dx);
`Y` a matrix of responses (n by dy);
`num_features` the number of variables to be selected, cannot be larger than dx;
`numCores` number of cores that are going to be used for parallelizing the process;
`ky` the kernel function for Y.
`kx` a list of length `num_features`, where `kx[[k]]` is the kernel used for (Xj1,...,Xjk), the first k selected variables;
`eps` a positive number, the regularization parameter for CME estimator;
`appro` whether to use incomplete Cholesky decomposition for approximation;
`tol` tolerance used for incomplete Cholesky decomposition (implemented by `inchol` in package `kernlab`).
The algorithm returns a vector of the indices from \code{1,...,dx} of the selected variables
``` r
n = 200
p = 100
set.seed(2)
X = matrix(rnorm(n * p), ncol = p)
Y = X[, 1] * X[, 2] + sin(X[, 1] * X[, 3]) + rnorm(n)
kx = c(kernlab::rbfdot(1),kernlab::rbfdot(1/2),kernlab::rbfdot(1/3))
CME_select(Y, X, rbfdot(1), kx, 3, eps = 1e-3, appro = F, numCores = 1)
# 1 2 3
```

`KPCCMElinear` is the CME estimator when `ky`, `kx`, `kxz` are all linear kernels,
in which case the incomplete Cholesky decomposition is the data matrix itself and could speed up the computation to a great extent.
It is included here for reproducing the results in the paper.
``` r
n = 2000
set.seed(1)
x = rnorm(n)
z = rnorm(n)
y = x + z + rnorm(n,1,1)
KPCCMElinear(y, x, z, 1e-5/n^(0.4))
# 0.4859424
# Theoretical KPC is 0.5
# compared with classical partial correlation squared:
# ppcor::pcor.test(y, z, x)$estimate^2 
# 0.4859428
```


## Case study on general spaces
KPC can be defined on general spaces, and if the objects in the general spaces can be stored in vectors, then the above functions can also be used.
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
Let `R1`, `R3` be the rotation around x-axis and z-axis.
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
# y1 is conditionally independent of z given x

for (eps in 10^(-(1:9))) print(MSE(x,y1,kernlab::rbfdot(1),SO3ker,eps))
# 2.507 2.502 2.509 2.523 2.607 3.567 14.860 124.848 955.870
for (eps in 10^(-(1:9))) print(MSE(cbind(x,z),y2,kernlab::rbfdot(1/2),SO3ker,eps))
# 2.510 2.505 2.509 2.527 2.639 3.301 8.307 52.744 473.358
# Similarly we can also examine MSE((x,z),y1) and MSE(x,y2)
# Choose the smallest possible eps while controlling the MSE
# Let eps = 1e-5

KPCCME(y1, x, z, SO3ker, rbfdot(1), rbfdot(0.5), 1e-5, appro = F)
# 0.6198004
KPCCME(y2, x, z, SO3ker, rbfdot(1), rbfdot(0.5), 1e-5, appro = F)
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
KFOCI(y, X, SO3ker, Knn=1, numCores = 1)
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
X = ElecData[5*(1:n),9:12]
for (i in 1:4) X[,i] = (X[,i] - mean(X[,i]))/sd(X[,i]) # normalize the data

library(kernlab)
KFOCI(Y, X, rbfdot(1/(2*median(dist(Y)^2))), Knn=1, numCores = 1)
# 1 2 3 4

# define two kernels on histograms
k1 = function(a,b) return(1/prod(a+b+1))
k2 = function(a,b) return(exp(-sum(sqrt(a+b))))
class(k1) = class(k2) = "kernel"

KFOCI(Y, X, k1, Knn=3, numCores = 1)
# 2 1 4 3
KFOCI(Y, X, k2, Knn=5, numCores = 1)
# 1 2 4

KPCgraph(Y,X[,c(2,3,4)],X[,1],rbfdot(1/(2*median(dist(Y)^2))),Knn = 2,trans_inv=TRUE)
# 0.1543532
KPCCME(Y,X[,c(2,3,4)],X[,1],rbfdot(1/(2*median(dist(Y)^2))),rbfdot(1/(2*median(dist(X[,c(2,3,4)])^2))), rbfdot(1/(2*median(dist(X)^2))), eps=1e-4, appro=F)
# 0.1473899
KPCgraph(Y,X[,c(1,2,4)],X[,3],rbfdot(1/(2*median(dist(Y)^2))),Knn = 2,trans_inv=TRUE)
# 0.05542749
KPCCME(Y,X[,c(1,2,4)],X[,3],rbfdot(1/(2*median(dist(Y)^2))),rbfdot(1/(2*median(dist(X[,c(1,2,4)])^2))), rbfdot(1/(2*median(dist(X)^2))), eps=1e-4, appro=F)
# 0.06338199
# X3 and X4 are both measures of richness. The conditional association between X3 and Y given all other variables is weaker than that of X1 and Y.

# Check MSE
for (eps in 10^(-(1:9))) print(MSE(X,Y,rbfdot(1/(2*median(dist(X))^2)),rbfdot(1/(2*median(dist(Y))^2)),eps))
# 0.580 0.605 0.575 0.606 0.692 2.253 32.855 339.731 731.645
```








