# delete all existing objects in the global environment
rm(list=ls())

# Tsai (2010, Analysis of Financial Time Series)
# Chapter 9, Section 9.5, Statistical Factor Models
# Examples 9.2 and 9.4

#############################
# Example 9.2
# read.table("URL", header=T)
rtn <- read.table("https://faculty.chicagobooth.edu/-/media/faculty/ruey-s-tsay/teaching/fts3/m-5clog-9008.txt", header=T)

rtn_IBM <- rtn$IBM

# format: time series
rtn_IBM <- ts(rtn_IBM, start=c(1990,1,30), frequency = 12) 
# time plot
ts.plot(rtn_IBM)

# variance-covariance matrix
cov_rtn <- cov(rtn)
# correlation matrix
cor_rtn <- cor(rtn)

# The eigenvalues and eigenvectors are then computed from 
# the covariance matrix with the eigen() or svd() function.
cov_rtn_eigen <- eigen(cov_rtn)
cov_rtn_svd <- svd(cov_rtn)

#The eigenvalues and eigenvectors are then computed from 
# the correlation matrix with the eigen() or svd() function.
cor_rtn_eigen <- eigen(cor_rtn)
cor_rtn_svd <- svd(cor_rtn)

# factor loadings?

###############################################################
# computing the principal components via the function: prcomp()
# If center=T, scale=T, it's equivalent to find the
# eigenvectors of the sample correlation matrix.
pca_rtn <- prcomp(rtn, center = TRUE, scale. = TRUE)

names(pca_rtn)
# "sdev"   "rotation"   "center"   "scale"    "x"   

# PC1 (the egienvector corresponding to largest value of eigenvalues)
pca_rtn$rotation

# scree plot for eigenvalues
# type = c("barplot", "lines")
screeplot(pca_rtn, type="lines")

# R-square-type interpretation
summary(pca_rtn)


#############################
# Example 9.4
# read.table("URL", header=T)
rtn_10 <- read.table("https://faculty.chicagobooth.edu/-/media/faculty/ruey-s-tsay/teaching/fts3/m-barra-9003.txt", header=T)

# info about factanal()
# https://www.geo.fu-berlin.de/en/v/soga/Geodata-analysis/factor-analysis/A-simple-example-of-FA/index.html

stat_fac <- factanal(rtn_10, factors=2, method="mle")
stat_fac

stat_fac <- factanal(rtn_10, factors=3, method="mle")
stat_fac

loadings_factor_1 <- sort(stat_fac$loadings[,1], decreasing = TRUE)
barplot(loadings_factor_1, main="Factor 1,  loadings")

loadings_factor_2 <- sort(stat_fac$loadings[,2], decreasing = TRUE)
barplot(loadings_factor_2, main="Factor 2,  loadings")

loadings_factor_3 <- sort(stat_fac$loadings[,3], decreasing = TRUE)
barplot(loadings_factor_3, main="Factor 3,  loadings")


########################################################################
########################################################################
# Approximate/Large factor models
# determining the number of factors
# Bai and Ng (2002 Econometrica) information criteria 

# resolve the version conflicts
# install.packages("pspline")
# install.packages("https://cran.r-project.org/src/contrib/Archive/phtt/phtt_3.1.2.tar.gz", repos=NULL, type="source")

library(phtt)
# the package requries numeric maxtrix
# convert dataframe into numeric matrix
data_str_nume <- data.matrix(rtn_10, rownames.force = NA)

# number of factors: r
phtt::OptDim(data_str_nume, d.max=10, standardize=TRUE, level=0.05,
             criteria = c("IC1", "IC2", "IC2", "PC1", "PC2", "PC3"))



