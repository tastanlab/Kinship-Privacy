
# A Utility Maximizing and Privacy Preserving Approach for Protecting Kinhsip in Genomic Databases
Assuming that family members arrive to the database in a sequential time order, a family f wants to hide their familial relationships while sharing their genomic data at the same time. We achieve this through minimizing privacy risks we have identified. There are two types of privacy risks that might reveal kinship information:
1. Being outliers in the population  
2. Having a kinship coefficient > 0 
  
_hide_snps.r_ code transform the problem into an optimization model by generating kinship and outlier constraints. However, in general there is not a solution that satisfies two types privacy risks. Therefore, we have developed a solution while satisfying one type of constraints and relaxing the other one. To achieve sharing maximum amount of data, the objective function is the minimum number of snsp to be hidden.  
  
hide_snps.r code solves the optimization model if outlier constraints are relaxed. Unfortunately, if kinship constraints are relaxed, those constraints become non-linear and Rcplex cannot solve non-linear constrainted integer programming problems, when the solution is based on relaxing kinship constraints, hide_snps.r only generates the constraints and writes to txt files. Then, _ga_lp2.m_ Matlab code is used to find minimum phi value that satisfies given non-linear constraints. After replacing the optimal phi value into non_linear kinship constraints, the model becomes non-linear and it is available to be solved by Rcplex.   
These cases are denoted as choice variable and it is the first input paremeter to run the code:
### If choice = 0:
Solve the problem where outlier constraints are relaxed and kinship constraints are satisfied. (This is linear integer programming problem). Also it minimizes slack variables (quantity of being outliers.) 
When choice = 0, run the code with following input parameters:  
```shell
Rscript hide_snps.r 0 o10 o11 o12
```
Note that the initial outlier values found in the population are: o10 = 27300, o11 = 27454, o12 = 15019   
### If choice = 1:
Solution is based on relaxing kinship constraints and satisfying outlier constraints. It generates constraints to find optimal phi value. (This is a non-linear constrained integer problem). These constraints are used in Matlab (See [Solution for choice = 1](https://github.com/tastanlab/Kinship-Privacy/blob/master/README.md#solution-for-choice--1) )
No extra parameters are needed to run this case:  
Rscript hide_snps.r 1  
### If choice = 2:
Solve the original problem by replacing optimal phi value found in the previous step. (This is linear integer programming problem). 
This case requires only one external parameter; Phi. i.e. 0.07:  
```shell
Rscript hide_snps.r 0 2 0.07  
```
### Solution for choice = 1
_ga_lp2.m_ solves the optimization model with non-linear constraints. The objective function is phi (the minimum kinship value that protects privacy of a family by satisfying outlier constraints.) 
## Finding families in the Genomic Database
_hierclustering_opensnp.r_ applies hierchical clustering on opensnp data (VCF format). (See [SNPRelate](http://corearray.sourceforge.net/tutorials/SNPRelate/) )
