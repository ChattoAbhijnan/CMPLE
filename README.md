# CMPLE
Correlation Modeling under Pairwise Likelihood Estimation using the MM algorithm.

The code is used to calculate estimates and the standard errors of the correlation
and standard deviation parameters using the pairwise composite likelihood method.  The method is abbreviated as CMPLE for Correlation Modeling
under Pairwise Likelihood Estimation (CMPLE).


This folder contains the following files besides this readme.txt. 

Data_generation1.R: R file (code) for simulating data on four responses (phenotypes) and two covariates (SNPs).

Data_generation2.R: R file (code) for simulating data on four responses (phenotypes) and ten covariates (SNPs).

CMPLE.R: R file (code) to perform estimation, standard error calculation, and obtain p-values. 

Demo_examples.R: R file (code) shows how to apply the CMPLE.R function to analyze two simulated datasets. Practitioners should save 
all four files in the same folder. Then execute each line of code of the file Demo_example.R. 

The following example illustrates how the code can be used.

Input arguments for the function CMPLE.R:

q: Number of responses.

p: Number of covariates.

mydata: Data frame containing the covariates followed by the responses in this particular order.

epsilon: Threshold for the MM optimization scheme (sum of the absolute relative difference of the parameter estimates 
between subsequent iterations). We suggest using 0.001.


Output of function CMPLE.R:


out1= Parameter estimates (first q*(p+1) entries are \alpha parameters, followed by the \delta parameters). The first (p+1) parameters of \alpha represent the 
regression parameters of the standard deviation for the first response, and so on. \alpha can be identified from the labeling of the 
regression parameters of the standard deviation of the response variables The first (p+1) parameters of \delta represent the 
regression parameters of the pairwise correlation of the response variables.  \delta can be identified from the labeling of the 
regression parameters of the pairwise correlation of the response variables

out2= Standard error of parameter estimates (In the same order as out1).

out3= Data frame containing parameter names, estimates, standard error, p-value, lower and upper limits of the 95\% confidence interval. 

One demo output can be found at the file: new_readme.txt
