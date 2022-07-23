The code is for calculating estimates and the standard errors of the parameters of the correlation 
and standard deviation using the pairwise composite likelihood method.  The method is abbreviated as CMPLE for Correlation Modeling
under Pairwise Likelihood Estimation (CMPLE).

This folder contains the following files besides this readme.txt. 

Data_generation.R: R file (code) for simulating data on four responses (phenotypes) and two covariates (SNPs).

CMPLE.R: R file (code) to perform estimation, standard error calculation, and obtain p-values. 

Demo_example.R: R file (code) shows how to apply the CMPLE.R function to analyze simulated data. Practitioners should save 
all three files in the same folder. Then execute each line of code of the file Demo_example.R. 

The following example illustrates how the code can be used.

Input arguments for the function CMPLE.R:

q: Number of responses.

p: Number of covariates.

mydata: Data frame containing the covariates followed by the responses, in this particular order.

epsilon: Threshold for the MM optimization scheme (sum of the absolute relative difference of the parameter estimates 
between subsequent iterations). We suggest using 0.001.


Output of function CMPLE.R:


out1= Parameter estimates (first q*(p+1) entries are \alpha parameters, followed by the \delta parameters). The first (p+1) parameters of \alpha represent the 
regression parameters of the standard deviation for the first response, and so on. \alpha can be identified from the labeling of the 
regression parameters of the standard deviation of the response variables The first (p+1) parameters of \delta represents the 
regression parameters of the pairwise correlation of the response variables.  \delta can be identified from the labeling of the 
regression parameters of the pairwise correlation of the response variables

out2= Standard error of parameter estimates (In the same order as out1).

out3= Data frame containing parameter names, estimates, standard error, p-value, lower and upper limits of the 95\% confidence interval. 

Example:
>source("Data_generation.R")
 # This code generates data for two predictors (p=2) and four responses (q=4). 
># However, our CMPLE function works for general p and q. 
>source("CMPLE.R")
>p=2
>q=4
>Covariates<- as.matrix(mydata[,1:p])
>Responses<- mydata[,(p+1):(p+q)]
>X<- model.matrix(~ Covariates)
>Y<- as.matrix(Responses)  
>result<- CMPLE(mydata,p=2,q=4,epsilon = 0.001)


>result$out3
     Parameter Estimates Standard_error P.value Lower_conf Upper_conf
1    alpha_1,0     -1.03           0.05    0.00      -1.13      -0.93
2    alpha_1,1      1.01           0.06    0.00       0.89       1.13
3    alpha_1,2      0.27           0.06    0.00       0.15       0.39
4    alpha_2,0      2.06           0.05    0.00       1.96       2.16
5    alpha_2,1      0.08           0.06    0.15      -0.04       0.20
6    alpha_2,2     -0.51           0.06    0.00      -0.63      -0.39
7    alpha_3,0      0.95           0.06    0.00       0.83       1.07
8    alpha_3,1      0.43           0.07    0.00       0.29       0.57
9    alpha_3,2      0.12           0.07    0.07      -0.02       0.26
10   alpha_4,0     -1.01           0.05    0.00      -1.11      -0.91
11   alpha_4,1     -0.44           0.06    0.00      -0.56      -0.32
12   alpha_4,2      1.02           0.06    0.00       0.90       1.14
13 delta_1,2,0      0.36           0.13    0.01       0.11       0.61
14 delta_1,2,1      0.49           0.17    0.00       0.16       0.82
15 delta_1,2,2      0.82           0.17    0.00       0.49       1.15
16 delta_1,3,0      0.33           0.17    0.05       0.00       0.66
17 delta_1,3,1      0.36           0.18    0.05       0.01       0.71
18 delta_1,3,2      0.50           0.18    0.01       0.15       0.85
19 delta_1,4,0      0.20           0.16    0.22      -0.11       0.51
20 delta_1,4,1      0.58           0.17    0.00       0.25       0.91
21 delta_1,4,2     -0.52           0.17    0.00      -0.85      -0.19
22 delta_2,3,0     -0.05           0.16    0.74      -0.36       0.26
23 delta_2,3,1      0.92           0.19    0.00       0.55       1.29
24 delta_2,3,2     -0.20           0.18    0.26      -0.55       0.15
25 delta_2,4,0      0.41           0.13    0.00       0.16       0.66
26 delta_2,4,1      0.57           0.16    0.00       0.26       0.88
27 delta_2,4,2     -0.99           0.16    0.00      -1.30      -0.68
28 delta_3,4,0     -0.28           0.14    0.05      -0.55      -0.01
29 delta_3,4,1      0.67           0.17    0.00       0.34       1.00
30 delta_3,4,2     -0.40           0.17    0.02      -0.73      -0.07