library(methods)
source("functions.r")
source("relationship_matrix_generator.r")
library(gdsfmt)
library(SNPRelate)
library(Rcplex)

#Input parameters. The first parameter is choice:
### choice = 0 : Solution for the model where outlier constraints are relaxed and kinship constraints are satisfied.
### choice = 1 : Creates outlier and kinship constraints to solve model where kinship constraints are relaxed. (This constraint matrix is used by Matlab to find minimum Phi value that satisfies created constraints. 
#Matlab is used because Rcplex do not solve nonlinear constrainted integer programming problems.)
### choice = 2 : Solution for the model where kinship constraints are relaxed. It replaces phi value obtained to the original problem. (Phi is solved by Matlab using parameters that is obtained from choice = 1).

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
#read the genotype information of family members
files<-list.files(path=getwd(),pattern=".vcf")
genof <- snpgdsOpen("c.gds")
g <- read.gdsn(index.gdsn(genof, "genotype"))
snp.id <- read.gdsn(index.gdsn(genof, "snp.id"))
snp.position <- read.gdsn(index.gdsn(genof, "snp.position"))
snp.chromosome <- read.gdsn(index.gdsn(genof, "snp.chromosome"))

#Calculate Kinship Coefficients
ibd.robust <- snpgdsIBDKING(genof, sample.id=1:length(files), family.id=NULL)
snp.id<-ibd.robust$snp.id
genotype<-g[snp.id,]
snp.pos <- snp.position[snp.id]
print(ibd.robust$kinship)

#create genotype matrix of family members where reverse arrival order to the database.
liste <- relationship_matrix()
genotype <- liste$genotype
relations <- liste$relations

flag_10 <- FALSE
flag_11 <- FALSE
flag_12 <- FALSE

# If there are more than 2 relationships:
if(dim(relations)[1]>=2)
{ 
  removal_index <- 1 # snps will be hidden from always first person (latest added family member).
  arr <- which(relations == removal_index, arr.ind = TRUE)
  
  #Relatives of latest added member
  removal_rel <- relations[arr[,1],]
  
  # Create X matrix
  number_of_people<-length(files) # number of people in family
  indices<-1:number_of_people
  X_general <- c()
  
  if (is.null(nrow(removal_rel)))
  {
    removal_rel <- t(removal_rel)
  }
  
  #Fill X matrix for every pair of relationship i.e. for a three membered family create all combinations where only pairwise relationships are considered:
  #i.e. assume that A,B,C are the three membered family members. A-B and B-C are related. A is the latest arrived member. Only latest arrived members snps will be hidden.
  #For A-B relationship only x11* is considered, so x1*1=[x101, x111, x121]
  #similarly for B-C relationship:
  #x*11 is considered which is [x011,x111,x211]. However since snps are hidden only from the latest arrived member, only x1..'s are considered.
  #As a result X_general matrix consists of [x101,x111,x112]
  #If A and C would be related as well, then X_general will consist [x101,x111,x112,x101,x121]. 
  for (z in 1:nrow(removal_rel))
  {
    index1<-removal_rel[z,1] 
    index2<-removal_rel[z,2]
    new_indices<-indices[-removal_rel[z,]] # other people's indices
    
    X <- matrix( nrow = (3^(number_of_people - 2)), ncol = number_of_people)
    col1<-rep(c(0,1,2), times=3^(number_of_people - 2)/3)
    X[,new_indices[1]]<-col1
    if(number_of_people ==4) 
    {
      col2<-rep(c(0,1,2), each=3^(number_of_people - 2)/3)
      X[,new_indices[2]]<-col2
    }
    
    X[,index1] <- rep(1,((3^(number_of_people - 2)))) # inceledigin ikiliye 1 ata
    X[,index2] <- rep(1,((3^(number_of_people - 2)))) # "
    X_general<- rbind(X_general,X)
  }
  
  X_general <- unique(X_general)
  
  #### Generate kinship constraints for optimization model  ####
  Amat<-c() #kinship constraint coefficients
  bvec<-c()#RHS of kinship constraints
  phi_rhs <-c()
  phi_x_array<- c()
  #calculate kinship constraints for every pair of relationships
  for (r in 1:nrow(relations))
  {
    #r<-3
    rel <- relations[r,]
    i<-rel[1]
    j<-rel[2]
    g1<-length(which(genotype[,i]==1,arr.ind=TRUE))
    g2<-length(which(genotype[,j]==1, arr.ind=TRUE))
    d <- length(which(genotype[,i]==1 & genotype[,j]==1, arr.ind=TRUE))
    c <- length(which(genotype[,i]==2 & genotype[,j]==0, arr.ind=TRUE)) + length(which(genotype[,i]==0 & genotype[,j]==2, arr.ind=TRUE))
    
    if(g1 > g2)
    {
      b<-g2
      a<-g1
    }else
    {
      b<-g1
      a<-g2
    } # if choice = 0 or 1 Phi=0.
    
    if(choice <= 1)
    {
      phi<-0
    }
    
    RHS <- 2*d - 4*c -a + b - 4*b*phi;
    if(choice == 1)
    {
      phi_rhs <- c(phi_rhs,4*b)
    }
    
    bvec<- c(bvec,RHS)
    related_X_indices <- get_relatedX(X_general, rel)
    related_X_indices <- sort(related_X_indices)
    related_x <- X_general[related_X_indices,]
    obj_arr<-matrix(0,nrow = nrow(X_general), ncol = 1, byrow = TRUE)
    phi_arr <- matrix(0,nrow = nrow(X_general), ncol = 1, byrow = TRUE)  
    #these are calculated based on KING coefficient:
    x_11 <- which(related_x[,i] == 1 & related_x[,j]==1 , arr.ind = TRUE)
    x_20 <- c( which(related_x[,i] == 2 & related_x[,j]==0 , arr.ind = TRUE), which(related_x[,i] == 0 & related_x[,j]==2, arr.ind = TRUE))
    p <- which(compMat2(related_x[x_11,],X_general))
    q <- which(compMat2(related_x[x_20,],X_general))
    obj_arr[p] <- obj_arr[p] + 2
    obj_arr[q] <- obj_arr[q] - 4
    
    #kinship formula is not symmetric:
    #i > j
    if(g1 > g2)
    {
      x_1 <- which(related_x[,i] == 1 , arr.ind = TRUE)  # a
      x__1 <- which(related_x[,j] == 1  , arr.ind = TRUE)  #b
      s <- which(compMat2(related_x[x__1,],X_general))
      t <- which(compMat2(related_x[x_1,],X_general))
      
    }else
    {
      x_1 <-which(related_x[,j] == 1, arr.ind = TRUE) #a
      x__1 <-which(related_x[,i] == 1, arr.ind = TRUE) #b
      s <- which(compMat2(related_x[x__1,],X_general))
      t <- which(compMat2(related_x[x_1,],X_general))
    }
    
    obj_arr[t] <- obj_arr[t] - 1;
    obj_arr[s] <- obj_arr[s] + (1 - 4*phi); 
    phi_arr[s]<- -4
    phi_x_array<-rbind(phi_x_array, t(phi_arr))
    Amat <- rbind(Amat,t(obj_arr)) #constraint matrix
  }
  
  if(choice == 1)
  {
    df <-data.frame( ncol=2*nrow(X_general) +2 )
    df <- cbind(Amat,phi_x_array, t(t(bvec)), t(t(phi_rhs)))
    colnames(df) <- c(rep("x",nrow(X_general)), rep("phi_x",nrow(X_general)), "RHS", "phi_RHS")
  }
  
  #Generate outlier constraints
  RHS2<-c()
  other_matrix <- c()
  for(a in 1:nrow(X_general))
  {
    x <- X_general[a,]
    str <- "";
    for(num in 1: number_of_people)
    { 
      if(num!=number_of_people)
      {
        str<-paste0(str,"genotype[,",num, "]==x[", num,"] & ") 
      }else
      {
        str<-paste0(str,"genotype[,",num, "]==x[", num,"]")  
      }  
    }
    str<-paste0("length(which(",str,",arr.ind=TRUE))" )
    n_x<-eval(parse(text=str))
    if(choice==0)
    {
      o_row<-rep(0,nrow(X_general) + 3)
    }else
    {
      o_row<-rep(0,nrow(X_general))
    }
    o_row[a] <- 1
    other_matrix <- rbind(other_matrix, o_row)
    RHS2 <- rbind(RHS2,n_x)
    
    if(choice==0)
    {
      for(index1 in 1:(number_of_people-1))
      {
        for(index2 in (index1+1):number_of_people)
        {
          if(!is.na(rowmatch(relations,matrix(, data = c(index1,index2), ncol = 2))))
          {
            if(x[index1]==1 & x[index2]==1)
            {
              n11 <- length(which(genotype[,index1]==1 & genotype[,index2]==1, arr.ind = TRUE))
              o_row<-rep(0,nrow(X_general) + 3)
              o_row[a] <- 1
              o_row[nrow(X_general) + 2] <- (-1)
              other_matrix <- rbind(other_matrix, o_row)
              
              new_rhs<-n11-o11
              if(new_rhs < 0)
              {
                new_rhs <- 0
              }
              RHS2 <-rbind(RHS2,new_rhs)
            }
            if((x[index1]==1 & x[index2]==0) || (x[index1]==0 & x[index2]==1))
            {
              n10 <- length(which(genotype[,index1]==x[index1] & genotype[,index2]==x[index2], arr.ind = TRUE))
              o_row<-rep(0,nrow(X_general) + 3)
              o_row[a] <- 1
              o_row[nrow(X_general) + 1] <- (-1)
              other_matrix <- rbind(other_matrix, o_row)
              new_rhs<-n10-o10
              if(new_rhs < 0)
              {
                new_rhs <- 0
              }
              RHS2 <-rbind(RHS2,new_rhs)
            }
            if((x[index1]==1 & x[index2]==2) || (x[index1]==2 & x[index2]==1))
            {
              n12 <- length(which(genotype[,index1]==x[index1] & genotype[,index2]==x[index2], arr.ind = TRUE))
              o_row<-rep(0,nrow(X_general) + 3)
              o_row[a] <- 1
              o_row[nrow(X_general) + 3] <- (-1)
              other_matrix <- rbind(other_matrix, o_row)
              new_rhs<-n12-o12
              if(new_rhs < 0)
              {
                new_rhs <- 0
              }
              
              RHS2 <-rbind(RHS2,new_rhs)
            }
          }
        }
      }
    }else
    {
      for(index1 in 1:(number_of_people-1))
      {
        for(index2 in (index1+1):number_of_people)
        {
          if(!is.na(rowmatch2(relations,matrix(, data = c(index1,index2), ncol = 2))))
          {
            if(x[index1]==1 & x[index2]==1)
            {
              n11 <- length(which(genotype[,index1]==1 & genotype[,index2]==1, arr.ind = TRUE))
              o_row<-rep(0,nrow(X_general))
              o_row[a] <- 1
              other_matrix <- rbind(other_matrix, o_row)
              new_rhs<-n11-o11
              new_rhs <- (n11-o11)
              if(new_rhs < 0)
              {
                print(c(index1,index2))
                print("n11 < o11")
                print(n11 < o11)
                flag_11 <- TRUE
                new_rhs <- 0
              }
              RHS2 <- rbind(RHS2,new_rhs) 
            }
            #constraints where index1 and index2 user have 1 or 0 snp values
            if((x[index1]==1 & x[index2]==0) || (x[index1]==0 & x[index2]==1))
            {
              n10 <- length(which(genotype[,index1]==x[index1] & genotype[,index2]==x[index2], arr.ind = TRUE))
              o_row<-rep(0,nrow(X_general))
              o_row[a] <- 1
              other_matrix <- rbind(other_matrix, o_row)
              new_rhs<-n10-o10
              if(new_rhs < 0)
              {
                flag_10 <- TRUE
                print("n10 < o10")
                new_rhs <- 0
              }
              RHS2 <- rbind(RHS2,(new_rhs))
            }
            #constraints for where index1 and index2 user have 1 or 2 snp values
            if((x[index1]==1 & x[index2]==2) || (x[index1]==2 & x[index2]==1))
            {
              n12 <- length(which(genotype[,index1]==x[index1] & genotype[,index2]==x[index2], arr.ind = TRUE))
              o_row<-rep(0,nrow(X_general))
              o_row[a] <- 1
              other_matrix <- rbind(other_matrix, o_row)
              new_rhs <- n12-o12
              if(new_rhs < 0)
              {
                flag12 <- TRUE
                print("n12 < o12")
                new_rhs <- 0
              }
              RHS2 <- rbind(RHS2,new_rhs) 
            }
          }
        }
      }
    }
  }
}else # Only one kinship relationship exist no need to calculate X_general. Only x11 is considered.
{
  i<-relations[1]
  j<-relations[2]
  g1<-length(which(genotype[,i]==1,arr.ind=TRUE))
  g2<-length(which(genotype[,j]==1, arr.ind=TRUE))
  d <- length(which(genotype[,i]==1 & genotype[,j]==1, arr.ind=TRUE))
  c <- length(which(genotype[,i]==2 & genotype[,j]==0, arr.ind=TRUE)) + length(which(genotype[,i]==0 & genotype[,j]==2, arr.ind=TRUE))
  
  if(g1 > g2)
  {
    b<-g2
    a<-g1
  }  else
  {
    b<-g1
    a<-g2
  }
  #solution for choice 0
  if(choice == 0)
  {
    phi<-0
    RHS <- (2*d - 4*c -a + b - 4*b*phi)/2;
    gen_11 <- which(genotype[,1]==1 & genotype[,2]==1)
    gen_10 <- which(genotype[,1]==1 & genotype[,2]==0)
    gen_12 <- which(genotype[,1]==1 & genotype[,2]==2)
    #find minimum number of snps to be hidden 
    Amat <- rbind(c(1,0),c(1,-1)) #constraints
    bvec <- c(RHS, length(gen_11)-o11) #Right handside
    cvec <- c(1,0) #objective function
    res<-Rcplex(cvec, Amat, bvec, Qmat = NULL, control=list(solnpoolintensity=4),
                objsense = c("min"), sense = c("G", "L"), vtype = rep("I",length(cvec)), n=1000)
    #find minimum number of slack variables
    Amat2 <- rbind(c(1,0),c(1,-1),c(1,0)) #additional constraint is taken from the previous model's solution: x11 = res[[1]]$obj
    bvec2 <- c(RHS, length(gen_11)-o11,res[[1]]$obj) 
    cvec2 <- c(0,1)
    res2<-Rcplex(cvec2, Amat2, bvec2, Qmat = NULL, control=list(solnpoolintensity=4),
                 objsense = c("min"), sense = c("G", "L","E"), vtype = rep("I",length(cvec)), n=1000)
    
    e10 <- 0
    e12 <- 0
    if(length(gen_10)-o10 < 0)
    {
      e10 <- o10 - length(gen_10)
    }
    if(length(gen_12)-o12 < 0)
    {
      e12 <- o12 - length(gen_12)
    }
    
    #slack variable for o11
    e11 <- res2[[1]]$obj
    x11<- res[[1]]$obj
    output <- data.frame(t(c(e10, e11, e12, x11)))
    colnames(output) <- c("e10", "e11", "e12", "x11")
    write.table(output,"out.txt", row.names = FALSE)
    
    #select randomly x11 SNPs
    snp_pos11 <- snp.pos[ which(genotype[,1]==1 & genotype[,2]==1) ]
    snp_pos11 <- sample(snp_pos11)
    snp_pos11 <- snp_pos11[1:x11]
    write.table(snp_pos11, "hide_snps.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, append = FALSE)
  }else  ## kinship is relaxed
  {
    n10 <- length(which(genotype[,1]==1 & genotype[,2]==0, arr.ind = TRUE))
    n11 <- length(which(genotype[,1]==1 & genotype[,2]==1, arr.ind = TRUE))
    n12 <- length(which(genotype[,1]==1 & genotype[,2]==2, arr.ind = TRUE))
    if( n10 < o10)
    {
      flag_10 <- n10 - o10
      print("n10 < o10")
    }
 
    if(n11 < o11)
    {
      flag_11 <-n11 - o11
      print("n11 < o11")
    }
    
    if(n12 < o12)
    {
      flag_12 <- n12 - o12
      print("n12 < o12")
    }
    
    if(choice == 2  & flag_10 == FALSE & flag_11 == FALSE & flag_12 == FALSE)
    {     
      RHS_x11 <- (2*d-4*c-a+b-4*b*phi)/(2*(1-2*phi))
      Amat<-matrix(c(1,1), nrow = 2)
      bvec <- c(RHS_x11, n11-o11)
      res<-Rcplex(c(1), Amat, bvec, Qmat = NULL, control=list(solnpoolintensity=4),
                  objsense = c("min"), sense = c("G","L"), vtype = c("I"), n=500) 
   
      #take the positions of snps to be hidden
      pos_11 <- snp.pos[sample(which(genotype[,1]==1 & genotype[,2] ==1 ))]
      snp.pos_11 <-pos_11[1:res[[1]]$obj[1]]  
      write.table(snp.pos_11,"hide_snps.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, append = FALSE) 
    }else if(choice == 1 ) #only generate constraints for matlab to solve nonlinear model.
    {
      if ( flag_10 < 0 | flag_11 < 0 | flag_12 < 0 )
      {
        write.table(matrix(c("n10 < o10", flag_10,"n11 < o11", flag_11, "n12 < o12", flag_12), ncol = 2, byrow = TRUE),"flag_constraints.txt", sep = "\t", col.names = FALSE, row.names = FALSE) 
      }
      
      ## kinship_constraints coefficients
      RHS <-  2*d - 4*c -a + b;
      phi_rhs <- -4*b
      x_phi <- -4
      x_11 <- 2
      df_2 <- data.frame(t(c(x_11, x_phi, RHS, phi_rhs)))
      write.table(df_2, "kin_constraints.txt", sep = " ",row.names = FALSE, col.names = FALSE)
      ## size_constraints 
      write.table(c(n11-o11), "size_constraints.txt", sep = " ",row.names = FALSE, col.names = FALSE)
      
    }else
    {
      print("no solution exist! n1. < o1.")
    }
  } 
}


#write the constraints to a txt file so matlab code can use it.
if(choice==1 & length(files) > 2)
{
  new_size_matrix <- c()
  for(index in 1:dim(X_general)[1])
  {
    new_size_matrix<- c(new_size_matrix,min(RHS2[which(other_matrix[,index]==1)]))
  }
  write.table(new_size_matrix,"size_constraints.txt",col.names = FALSE, row.names = FALSE)  
  write.table(df,"kin_constraints.txt",col.names = FALSE, row.names = FALSE)
  if(flag_10 < 0 | flag_11 < 0 | flag_12 < 0)
  {
    write.table(matrix(c("n10 < o10", flag_10,"n11 < o11", flag_11, "n12 < o12", flag_12), ncol = 2, byrow = TRUE),"flag_constraints.txt", sep = "\t", col.names = FALSE, row.names = FALSE) 
  }
}


# Solve the model which minimizes the number of hidden snps where optimal phi value is already found by matlab and calculated constraints previously. 
if(choice==2 & length(files) > 2)
{
  
  new_size_matrix <- c()
  for(index in 1:dim(X_general)[1])
  {
    new_size_matrix<- c(new_size_matrix,min(RHS2[which(other_matrix[,index]==1)]))
  }
  
  Amat <- rbind(Amat, diag(1,dim(X_general)[1],dim(X_general)[1]))
  bvec<-c(bvec,new_size_matrix)
  
  sense_vec<- c(rep("G", dim(relations)[1]), rep("L",length(new_size_matrix)))
  flag_minus <- FALSE
  minus_ind <- c()
  for( abc in 1:length(bvec))
  {
    if(bvec[abc] < 0)
    {
      flag_minus <- TRUE
      minus_ind <- c(abc, minus_ind)
    }
  }
  
  if(flag_minus == TRUE)
  {
    Amat[minus_ind,] <- -(Amat[minus_ind,])
    bvec[minus_ind] <- -bvec[minus_ind]
    sense_vec[minus_ind] <- "L"
  }
  
  cvec <- c(rep(1,nrow(X_general)))
  res<-Rcplex(cvec, Amat, bvec, Qmat = NULL, control=list(solnpoolintensity=4),
              objsense = c("min"), sense = sense_vec, vtype = rep("I",length(cvec)), n=500)
  print(res[[1]])
  
  
  remove_snps<-c()
  if(length(res)>0)
  {
    for(row_ind in 1:dim(X_general)[1])
    {
      snp_position <- snp.pos[which(apply(genotype,1, function(x) all.equal(x,X_general[row_ind,])) == "TRUE")]
      snp_position <- sample(snp_position)
      snp_position <- snp_position[1:(res[[1]]$xopt[row_ind])]
      remove_snps <- c(remove_snps, snp_position)
    }
    write.table(remove_snps, "hide_snps.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, append = FALSE)
  }else
  {
    print("there is no solution!")  
  }
  
}
                                          
if (choice ==0 & length(files) > 2) {
  vec<-duplicated(other_matrix)
  duplic_rows <- other_matrix[vec,]
  
  if(!(is.null(duplic_rows)))
  {
    row_ind <- rowmatch2(duplic_rows, other_matrix )
    unik <- na.omit(unique(row_ind))
    
    for (k in 1:length(unik))
    {
      row_ind <- rowmatch2(duplic_rows, other_matrix )
      unique_indices <- na.omit(unique(row_ind))
      ind_array <- which(row_ind==unique_indices[k])
      new_rhs_value <- min(RHS2[(ind_array)])
      other_matrix <- other_matrix[-ind_array[2:length(ind_array)],]
      RHS2 <- RHS2[-ind_array[2:length(ind_array)]]
      RHS2[ind_array[1]] <- new_rhs_value
    }
  } 
  Amat<-cbind(Amat,matrix(0,nrow = nrow(Amat), ncol=3))
  Amat <-rbind(Amat,other_matrix)
  
  bvec<-c(bvec,t(RHS2))
  
  #minimize the number of snps to be hidden 
  cvec <- c(rep(1,nrow(X_general)),0,0,0) #objective
  res<-Rcplex(cvec, Amat, bvec, Qmat = NULL, control=list(solnpoolintensity = 4),
              objsense = c("min"), sense = c(rep("G", nrow(relations) ), rep("L",dim(Amat)[1]-nrow(relations))), vtype = rep("I",length(cvec)), n =100)

  #take the first solution and write it as another constraint i.e. x101+x111+x121 = 4000. Then solve the problem to find minimum number of e1+e2+e3. (Find minimum outlier leakege)
  
  #minimize e1 + e2 + e3 slack variables
  Amat_yeni <- rbind(Amat, c(rep(1,nrow(X_general)),0,0,0))
  bvec_yeni <- c(bvec, res[[1]]$obj)
  cvec_yeni <- c(rep(0,nrow(X_general)),1,1,1)
  res2 <- Rcplex(cvec_yeni,Amat_yeni, bvec_yeni, Qmat = NULL, control=list(solnpoolintensity=4),
                 objsense = c("min"), sense = c(rep("G", nrow(relations) ), rep("L",length(RHS2)), "E"), vtype = rep("I",length(cvec)), n=100)
  
  if(length(res)==0 | length(res2)==0)
  {
    print("no solut??on ex??st")
  }else
  {
    write.table(res2[[1]]$xopt, "out.txt")
    #removal of snps
    remove_snps <-c()
    for(row_ind in 1:dim(X_general)[1])
    {
      snp_position <- snp.pos[which(apply(genotype,1, function(x) all.equal(x,X_general[row_ind,])) == "TRUE")]
      snp_position <- sample(snp_position)
      snp_position <- snp_position[1:(res2[[1]]$xopt[row_ind])]
      remove_snps <- c(remove_snps, snp_position)
    }    
    write.table(remove_snps,"hide_snps.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, append = FALSE) 
  }
}
