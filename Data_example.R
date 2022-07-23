

list.of.packages<-c("dplyr", "mvtnorm","MASS","Matrix","matrixcalc", "numDeriv" )

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

require(mvtnorm)
require(MASS)
require(Matrix)
require(matrixcalc)
require(numDeriv)

#setwd("") #put the appropriate directory where Data_generation.R, CMPLE.R,and Data_example.R has been saved  

p=2
q=4

source("Data_generation.R") # This simulates a data set with p=2 and q=4 for a 
# sample size of n=500. 
source("CMPLE.R")

results<-CMPLE(mydata,p=2,q=4,epsilon = 0.001)
# The results is a list of objects, named out1, out2, out3. results$out3 is data frame
# containing the estimates, standard errors, p-value and 95% confidence intervals. 
