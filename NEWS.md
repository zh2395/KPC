# KPC v0.1.2
In this version:
* The default choice of K for the K-nearest neighbor (K-NN) graph in the variable selection algorithm KFOCI() is changed to 0.05n. More details of our suggestions on the choice of K are added to the documentation.
* KFOCI() and KPCRKHS_VS() will run in parallel by default. 

# KPC v0.1.1
This is an update of the orignial version 0.1.0 of the package "KPC".
* The package "RcppMLPACK" used by KPC v0.1.0 previously to compute the Euclidean mimimum spanning tree (MST) is scheduled for archival on 2021-12-20, and we will use the package "mlpack" in v0.1.1 to compute the MST instead.
* The variable selection algorithm KPCRKHS_VS() uses the median of pairwise distances as the bandwidth of Gaussian kernel by default. It is likely that this median is 0. In this version we will raise a warning if the median happens to be 0 and use the mean instead. If the mean is also 0, then the algorithm raises an error, suggesting that there exists a feature of constants.