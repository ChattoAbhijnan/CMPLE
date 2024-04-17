list.of.packages<-c("dplyr","mvtnorm","MASS","Matrix","matrixcalc","numDeriv" )

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

require(mvtnorm)
require(MASS)
require(Matrix)
require(matrixcalc)
require(numDeriv)

#setwd("") #put the appropriate directory where Data_generation1.R,
#Data_generation2.R , CMPLE.R,and Data_example.R has been saved  


#Example 1 which simulates a data set with p=2 and q=4 for a 
# sample size of n=600. 

p=2
q=4

source("Data_generation1.R") 
source("CMPLE.R")

results_1<-CMPLE(mydata,p=2,q=4,epsilon = 0.001)
# The results is a list of objects, named out1, out2, out3. Notably, out3 is data frame
# containing the estimates, standard errors, p-value and 95% confidence intervals. 
# Here, the alpha parameters denote regression parameters concerning standard deviation
# and delta parameters denotes regression parameters concerning pairwise correlations
# A demo output can be found at https://github.com/ChattoAbhijnan/CMPLE/blob/main/Demo_output.txt

#############################################################################################

#Example 2 which simulates a data set with p=10 and q=4 for a 
# sample size of n=1000. 

p=10
q=4

source("Data_generation2.R") 
source("CMPLE.R")

results_2<-CMPLE(mydata,p=10,q=4,epsilon = 0.001)
# The results is a list of objects, named out1, out2, out3. Notably, out3 is data frame
# containing the estimates, standard errors, p-value and 95% confidence intervals. 
# Here, the alpha parameters denote regression parameters concerning standard deviation
# and delta parameters denotes regression parameters concerning pairwise correlations
# A demo output can be found at https://github.com/ChattoAbhijnan/CMPLE/blob/main/Demo_output.txt
