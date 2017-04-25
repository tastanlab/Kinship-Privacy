library(methods)
source('functions.r')
source("relationship_matrix_generator.r")
library(gdsfmt)
library(SNPRelate)
library(Rcplex)


#Input parameters. The first parameter is choice:
### choice = 0 : Solution for the model where outlier constraints are relaxed and kinship constraints are satisfied.
### choice = 1 : Creates outlier and kinship constraints to solve model where kinship constraints are relaxed. (This constraint matrix is used by Matlab to find minimum Phi value that satisfies created constraints. 
#Matlab is used because Rcplex do not solve nonlinear constrainted integer programming problems.)
### choice = 2 : Solution for the model where kinship constraints are relaxed. It replaces phi value obtained from the original problem. (Phi is solved by Matlab using parameters that is obtained from choice = 1).

args = commandArgs(trailingOnly=TRUE)
choice <- as.numeric(args[1])
if(choice==0) 
{
  #reduced outlier values in the population
  o10 <- as.numeric(args[2])
  o11 <- as.numeric(args[3])
  o12 <- as.numeric(args[4])
  #which family is evaluated:
  family <- as.numeric(args[5])
  
}else if(choice == 1)
{
  phi <- 0 
  #initial outlier values found in the population before hiding any SNPs.
  o10<-27300
  o11<-27454
  o12 <- 15019
  family <- as.numeric(args[1])
}else
{
  phi<-read.table("phi.txt")
  phi<-phi$V1
  o10<-27300
  o11<-27454
  o12 <- 15019
  family <- as.numeric(args[1])
}

flag_10 <- FALSE
flag_11 <- FALSE
flag_12 <- FALSE


## Define relationship matrices in the family by reading VCF files (of existing and incoming members in the family). ###
genotype <- relationship_matrix()



