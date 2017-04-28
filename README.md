
# A Utility Maximizing and Privacy Preserving Approach for Protecting Kinhsip in Genomic Databases
Assuming that family members arrive to the database in a sequential time order, a family f wants to hide their familial relationships while sharing their genomic data at the same time. We achieve this through minimizing privacy risks we have identified. There are two types of privacy risks that might reveal kinship information:
1. Being outliers in the population  
2. Having a kinship coefficient > 0 
  
_hide_snps.r_ code transforms the problem into an optimization model by generating kinship and outlier constraints. However, in general there is not a solution that satisfies two types privacy risks (constraints). Therefore, we have developed a solution while satisfying one type of constraints and relaxing the other one. To achieve maximum amount of data to be shared in public databases, we defined the objective function as the minimum number of SNPs to be hidden.
  
hide_snps.r code solves the optimization model if outlier constraints are relaxed. Unfortunately, if kinship constraints are relaxed, those constraints become non-linear and Rcplex cannot solve non-linear constrainted integer programming problems, when the solution is based on relaxing kinship constraints, hide_snps.r only generates the constraints and writes to txt files. Then, _optimal_phi.m_ Matlab code is used to find minimum phi value that satisfies given non-linear constraints. After replacing the optimal phi value into non_linear kinship constraints, the model becomes non-linear and it is available to be solved by Rcplex.   
These cases are denoted as _choice_ variable and it is the first input paremeter to run the code:
### If choice = 0:
Solve the problem where outlier constraints are relaxed and kinship constraints are satisfied. (This is linear integer programming problem). Also it minimizes slack variables (quantity of being outliers.) 
When choice = 0, there are 3 more parameters needed:
```shell
Rscript hide_snps.r 0 o10 o11 o12
```
Note that the initial outlier values found in the population are: o10 = 27300, o11 = 27454, o12 = 15019   
### If choice = 1:
Solution is based on relaxing kinship constraints and satisfying outlier constraints. It generates constraints to find optimal phi value. (This is a non-linear constrained integer problem). These constraints are used in Matlab (See [Solution for choice = 1](https://github.com/tastanlab/Kinship-Privacy/blob/master/README.md#solution-for-choice--1) )
No extra parameters are needed to run this case:  
```shell
Rscript hide_snps.r 1  
```
which outputs _phi.txt_ , _size_constraints.txt_ ,  _kin_constraints.txt_. These txt files will be read as inputs in Matlab. Be sure that all files are in the same folder.    
Now run _optimal_phi.m_ Matlab code to solve the optimization model with non-linear constraints. Here, the objective function is phi and the aim is to find the minimum kinship value that protects privacy of a family by satisfying outlier constraints. After finding optimal phi, run _hide_snps.r_ as _choice_ = 2.

### If choice = 2:
Solve the original problem by replacing optimal phi value found in the previous step. (This is linear integer programming problem). 
This case requires only one extra parameter; Phi. i.e. 0.07:  
```shell
Rscript hide_snps.r 2 0.07  
```
### You should have:
[Rcplex](https://cran.r-project.org/web/packages/Rcplex/index.html), [Bioconductor](https://www.bioconductor.org) packages gdsfmt and [SNPRelate](http://corearray.sourceforge.net/tutorials/SNPRelate/) under R and
[Global optimization toolbox](https://www.mathworks.com/products/global-optimization.html) under Matlab


## Finding families in the Genomic Database
We found the families in openSNP by applying hierarchical clustering as clusters. Code is available in _hierclustering_opensnp.r_ and we used ([SNPRelate](http://corearray.sourceforge.net/tutorials/SNPRelate/)).
