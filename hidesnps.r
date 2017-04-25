library(methods)
source('functions.r')
args = commandArgs(trailingOnly=TRUE)
choice <- as.numeric(args[1])
if(choice==0) #take outlier values from user
{
  o10 <- as.numeric(args[2])
  o11 <- as.numeric(args[3])
  o12 <- as.numeric(args[4])
  
}else if(choice == 1)
{
  phi <- 0
  o10<-27300
  o11<-27454
  o12 <- 15019
}else
{
  phi<-read.table("phi.txt")
  phi<-phi$V1
  print(phi)
  o10<-27300
  o11<-27454
  o12 <- 15019
}
