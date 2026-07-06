########################################################################
# Factor analysis, principal components
# Jau-er Chen 
########################################################################

rm(list=ls(all=TRUE))
# set working directory
#setwd("/Users/jauer/Dropbox/workshop2019_factor/prg/str_pca")
setwd('/Users/chenshengyuan/Desktop/NTU_ECON/Course/Advanced econometric topic/factor_models/hands-on_R_factors/str_pca')

########################################################################
# data  sciencer package
# install.packages("tidyverse")
library(tidyverse)
# pipe %>%; %<>%
# assign a data.frame to data argument %$%
library(magrittr)

########################################################################

# import csv dataset 
ca_school <- read_csv(file = "ca_school.csv")

# remove test score related variables: testscore, elarts_score, math_score
# remove string variables: countyname, districtname, schoolname
# remove 0-1 binary variables: charter_s, unified_d, la_unified_d, sd_unified_d, yrcal_s

# We will do pca on the resulting dataset (data_str)
data_str <- dplyr::select(ca_school, -c(testscore, elarts_score, math_score, countyname, districtname, schoolname, charter_s, unified_d, la_unified_d, sd_unified_d, yrcal_s))

########################################################################
# computing the principal components via the function: prcomp()
pca_data_str <- prcomp(data_str, center = TRUE, scale. = TRUE)
summary(pca_data_str)

names(pca_data_str)
# "sdev"   "rotation"   "center"   "scale"    "x"   

# The "center" and "scale" components correspond to the means
# and standard deviations of the variables that were used for
# scaling prior to implementing PCA

# The "rotation" matrix provides the principal component loadings;
# each column of pca_data_str$rotation contains the corresponding 
# principal component loading vector, that is, 
# when we matrix-multiply the X matrix by pca_data_str$rotation,
# it gives us the coordinates of the data in the rotated coordinate system.
# These coordinates are the principal component scores.

# the pca object contains the relationship between the initial variables and 
# the principal components ($rotation): i.e. variable loadings: w_{j}
# PC_j = w_{j1} x_1 + ... + w_{jp} x_p

top2.pca.eigenvector <- pca_data_str$rotation[, 1:2] # the first two columns of the rotation matrix
top2.pca.eigenvector

first.pca <- top2.pca.eigenvector[, 1]   # PC_1
second.pca <- top2.pca.eigenvector[, 2]  # PC_2

# sorting values and plot dotchart
# PC_1
first.pca[order(first.pca, decreasing=FALSE)] 

dotchart(first.pca[order(first.pca, decreasing=FALSE)][79:99]  ,   
         main="Loading Plot for PC1",                      
         xlab="Variable Loadings",                        
         col="red")    

# ploting PC_1 and PC_2. 
# (name PC_1, for example: income, family background: white, housing owner)

library(ggfortify) # use autoplot

autoplot(pca_data_str, data = data_str, colour = 'med_income_z')
autoplot(pca_data_str, data = data_str, colour = 're_wht_frac_s')


########################################################################
# PC_2
second.pca[order(second.pca, decreasing=FALSE)]  

dotchart(second.pca[order(second.pca, decreasing=FALSE)][79:99]  ,  
         main="Loading Plot for PC2",                       
         xlab="Variable Loadings",                          
         col="blue")   

# ploting PC_1 and PC_2. 
# (name PC_2, for example: str, teacher's salary, school quality)

autoplot(pca_data_str, data = data_str, colour = 'str_s')
autoplot(pca_data_str, data = data_str, colour = 'te_salary_avg_d')

########################################################################
# $x has its columns the principal component score vectors; that is,
# the k-th column is the k-th principal component score vector.
dim(pca_data_str$x)

# biplot: plot the first two principal components
# we have many variables here; the biplot is messy though.
biplot(pca_data_str, scale=0)
# the scale=0 argument to biplot() ensures that
# the arrows are scaled to represent the loadings.


########################################################################
# rule of thumb: number of PCs

plot(pca_data_str,         
     type="line", 
     ylim=c(0, 20),
     main="Scree Plot for PCs") 
abline(h=1, col="red") # Kaiser eigenvalue-greater-than-one rule

########################################################################
########################################################################
# Approximate/Large factor models
# determining the number of factors
# Bai and Ng (2002 Econometrica) information criteria 

library(phtt)
# the package requries numeric maxtrix
# convert dataframe into numeric matrix
data_str_nume <- data.matrix(data_str, rownames.force = NA)

# number of factors: r
phtt::OptDim(data_str_nume, d.max=10, criteria = c("IC1", "IC2", "IC2",
                                                   "PC1", "PC2", "PC3"))

########################################################################
########################################################################
# Prediction
# Determining the number of factors, r, via cross validation

# Principal Component Regression ( or Factor Augmented Regression)
# involes constructing the first r PCs, and then using these components
# as predictors in a linear regression model that is fit using OLS.

# we need the following library and function pls::pcr()
# we now apply this to the data_str, in order to predict "testscore".

# construct testscore, and add this "outcome varialbe" to the dataset
testscore <- (ca_school$elarts_score + ca_school$math_score)/2
data_str$testscore <- testscore

# we then construct the training set and the test set:
x <- model.matrix(testscore~., data_str)[, -dim(data_str)[2]] 
y <- data_str$testscore
set.seed(1)
train <- sample(1:nrow(x), nrow(x)/2)
test <- (-train)
y.test <- y[test]

# training the model
library(pls)
f.augmented.reg <- pls::pcr(testscore~., data=data_str, subset=train, scale=TRUE,
                            validation="CV")
# "subset": an optional vector specifying a subset of observations 
# to be used in the fitting process.

# setting "scale=TRUE" has the effect of standardizing each predictor,
# prior to generating the principal components
# setting "CV" computes the 10-fold cross-validation error 
# for each possible "r", # of PCs used.

summary(f.augmented.reg) # CV score for each possible r (r=0 onwards).
# much easier to pin down "r" through the polot of cross-validation scores
# using the validationplot()
pls::validationplot(f.augmented.reg, val.type="MSEP")
# "MSEP" will cause the cross-validation MSE to be plotted.
# we find that the lowest cross-validation error occurs when r = 79
# training model done !

#######################################
# now evaluate its test set performance 

#f.augmented.reg.pred <- predict(f.augmented.reg, x[test, ], ncomp=79)
# compute the test MSE as follows
#mean((f.augmented.reg.pred - y.test)^2)

#######################################

# Finally, we fit the factor augmented reg on the full data set,
# using r = 79, the number of components identified by cross-validation







