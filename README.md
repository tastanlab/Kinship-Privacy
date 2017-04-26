## Finding families in the Genomic Database
hierclustering_opensnp.r applies hierchical clustering on opensnp data (VCF format). (See [SNPRelate](http://corearray.sourceforge.net/tutorials/SNPRelate/) )
## A Utility Maximizing and Privacy Preserving Approach for Protecting Kinhsip in Genomic Databases
hide_snps code generates kinship and outlier constraints. It can solve the optimization model if outlier constraints are relaxed.  If kinship constraints are relaxed, those constraints become non-linear. Unfortunately, Rcplex cannot solve non-linear constrainted integer programming problems, so this R code only generates the constraints and writes to txt files and then ga_lp2.m code is used to find minimum phi value that satisfies given non-linear constraints. We replace the optimal phi value into non_linear kinship constraints and make them linear. When the problem becomes linear with found phi value, we can solve original problem in Rcplex.
These cases are denoted as choice variable and it is the first input paremeter to run the code:
#### If choice = 0:
Solve the problem where outlier constraints are relaxed and kinship constraints are satisfied. (This is linear integer programming problem). Also minimize slack variables (quantity of being outliers.) 
When choice = 0, run the code with following input parameters:  
Rscript hide_snps.r 0 27300 27454 15019  
The last three parameters denote outlier values o10,o11,o12 found in the population.  
#### If choice = 1:
Solution is based on relaxing kinship constraints and satisfying outlier constraints. It generates constraints to find optimal phi value. (This is a non-linear constrained integer problem). These constraints are used in Matlab (See [Solution for choice = 1](https://github.com/tastanlab/Kinship-Privacy/blob/master/README.md#solution-for-choice--1) )
No extra parameters are needed to run this case:
Rscript hide_snps.r 1
#### If choice = 2:
Solve the original problem by replacing optimal phi value found in the previous step. (This is linear integer programming problem). 
This case requires only one external parameter; Phi. i.e. 0.07:
Rscript hide_snps.r 0 2 0.07
### Solution for choice = 1
ga_lp2.m solves the optimization model with non-linear constraints. The objective function is phi (the minimum kinship value that protects privacy of a family by satisfying outlier constraints.) 
