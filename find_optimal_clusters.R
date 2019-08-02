# Description:
#   A small wrapper function around the R package NbClust to find the optimal number of
#   clusters for a single cell dataset. Input should be a matrix where columns are
#   single cells and rows are genes.
#
# Arguments:
#   x = data matrix
#   min.nc = minimum number of clusters to test (must be > 1)
#   max.nc = maximum number of clusters to test
#   ncores = number of cores to use (one core per index)
#   method = the clustering method to evaluate, kmeans as default
#   runPCA = run principal component analysis on the data before the test
#
# Example:
#   ret <- find_opt_clus(x, min.nc=2, max.nc=10, ncores=30)
#
# Author:
#   Oscar Franzen <p.oscar.franzen@gmail.com>, July 2019
#
# References:
#   [1] Charrad et al (2014), Journal of Statistical Software
#       http://www.jstatsoft.org/v61/i06/
#   [2] https://cran.r-project.org/web/packages/NbClust

if (!require("NbClust")) install.packages("NbClust")
if (!require("parallel")) install.packages("parallel")

find_opt_clus <- function(x, min.nc=2, max.nc, ncores=2, method="kmeans",
                             runPCA=FALSE, comp=50) {
    if (missing(x)) stop("Input data is missing.")
    if (missing(max.nc)) stop("Maximum number of clusters (max.nc) must be specified.")
    
    indices = c("kl", "ch", "hartigan", "ccc", "scott", "marriot", "trcovw",
                "tracew", "friedman", "rubin", "cindex", "db", "silhouette",
                "duda", "pseudot2", "beale", "ratkowsky", "ball",
                "ptbiserial", "gap", "frey", "mcclain", "gamma", "gplus",
                "tau", "dunn", "hubert", "sdindex", "dindex", "sdbw")
    m <- t(x)
    if (runPCA) {
        pr <- prcomp(m)
        m <- pr$x[, 1:comp]
    }
    help <- function(ind) {
        print(paste0("Working on index ", ind))
        NbClust(m, min.nc=min.nc, max.nc=max.nc, method=method, index=ind)
    }
    ret <- mclapply(indices, FUN=help, mc.cores=ncores)
    sort(table(unlist(lapply(ret, function(x) unlist(x['Best.nc'])[1]))))
}
