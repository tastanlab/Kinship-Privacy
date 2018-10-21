library(methods)
#setwd("/Users/apple/Dropbox/kinship_opt2/phi2/teyze-cocuk/vcf/")

#o10<-27300
#o11<-13037
#o12<-15019

#o10<-27300
#o11<-27454
#o12 <- 15019
#choice <-0
get_relatedX <- function(X_general , rel) 
{
  pairs<-matrix(c(1,2,1,0,0,1,2,1,0,2,2,0,1,1), ncol = 2, byrow = TRUE) # kinship calculation is affected by these pairs
  related_X_indices<-c()
  for( p in 1:nrow(pairs))
  {
    p1<-pairs[p,1]
    p2<-pairs[p,2]
    ind <- which(X_general[,rel[1]]==p1 & X_general[,rel[2]]==p2, arr.ind=TRUE);
    related_X_indices<- c(ind,related_X_indices);
    
  }
  
  return(related_X_indices)
}

compMat2 <- function(A, B) {  # rows of B present in A 
  B0 <- B[!duplicated(B), ] 
  na <- nrow(A); 
  nb <- nrow(B0) ;
  if( is.null(na))
    na <- 1
  if(is.null(B))
    nb <- 1;
  
  AB <- rbind(A, B0) 
  ab <- duplicated(AB)[(na+1):(na+nb)] 
  
  return(ab) 
} 


rowmatch <- function(x, want) {
  isTRUE(all.equal(x, want))
}



rowmatch2 <- function(A,B) { 
  # Rows in A that match the rows in B
  f <- function(...) paste(..., sep=":")
  if(!is.matrix(B)) B <- matrix(B, 1, length(B))
  a <- do.call("f", as.data.frame(A))
  b <- do.call("f", as.data.frame(B))
  match(b, a)
}



args = commandArgs(trailingOnly=TRUE)
choice <- as.numeric(args[1])
if(choice==0)
{
  o10 <- as.numeric(args[2])
  o11 <- as.numeric(args[3])
  o12 <- as.numeric(args[4])
  
}else
{
  phi<- as.numeric(args[2])
  o10<-27300
  o11<-27454
  o12 <- 15019
}


library(gdsfmt)
library(SNPRelate)
library(Rcplex)


files<-list.files(path=getwd(),pattern=".vcf")


genof <- snpgdsOpen("c.gds")
g <- read.gdsn(index.gdsn(genof, "genotype"))
snp.id <- read.gdsn(index.gdsn(genof, "snp.id"))
snp.position <- read.gdsn(index.gdsn(genof, "snp.position"))
snp.chromosome <- read.gdsn(index.gdsn(genof, "snp.chromosome"))


ibd.robust <- snpgdsIBDKING(genof, sample.id=1:length(files), family.id=NULL)
snp.id<-ibd.robust$snp.id
genotype<-g[snp.id,]
snp.pos <- snp.position[snp.id]
print(ibd.robust$kinship)





## Define relationship matirx ###
relations<- matrix(c(3,5, 0,4,5,0, 3,4,0, 1,2,0 ,1,3,0,2,3,0), ncol = 3 ,byrow = TRUE)
genotype_index<-matrix(,length(files),2)

for(i in 1:length(files))
{
  if(isTRUE(grep("1860",files[i]) == 1)) {
    genotype_index[i,1]<-1
    genotype_index[i,2]<-"baba"
    if(any(relations[,1]==1)){
      rows = relations[,1]==1
      relations[rows,3] = relations[rows,3] + 1
    }
    if(any(relations[,2]==1)){
      rows = relations[,2]==1
      relations[rows,3] = relations[rows,3] + 1
    }
  } else if(isTRUE(grep("1861",files[i]) == 1)){
    genotype_index[i,1] <- 2
    genotype_index[i,2] <-"amca"
    if(any(relations[,2]==2))
    {
      rows = relations[,2]==2
      relations[rows,3] = relations[rows,3] + 1
    }
    if(any(relations[,1]==2))
    {
      rows = relations[,1]==2
      relations[rows,3] = relations[rows,3] + 1
    } 
    
  } else if(isTRUE(grep("1862",files[i]) == 1)){
    genotype_index[i,1]<-3
    genotype_index[i,2]<-"cocuk"
    if(any(relations[,2]==3))
    {
      rows = relations[,2]==3
      relations[rows,3] = relations[rows,3] + 1
    }
    if(any(relations[,1]==3))
    {
      rows = relations[,1]==3
      relations[rows,3] = relations[rows,3] + 1
    } 
    
  } else if(isTRUE(grep("1863",files[i]) == 1)){
    genotype_index[i,1] <- 4
    genotype_index[i,2] <- "anne"
    
    if(any(relations[,2]==4))
    {
      rows = relations[,2]==4
      relations[rows,3] = relations[rows,3] + 1
    }
    if(any(relations[,1]==4))
    {
      rows = relations[,1]==4
      relations[rows,3] = relations[rows,3] + 1
    } 
  }else  {
    genotype_index[i,1]<-5
    genotype_index[i,2]<-"teyze"
    if(any(relations[,2]==5))
    {
      rows = relations[,2]==5
      relations[rows,3] = relations[rows,3] + 1
    }
    if(any(relations[,1]==5))
    {
      rows = relations[,1]==5
      relations[rows,3] = relations[rows,3] + 1
    }
  }
}

relations = relations[relations[,3]==2,1:2]
if (is.null(dim(relations)))
{
  relations <- t(relations)
}  

pos = gregexpr('/', getwd()) 
foldername = substr(getwd(),pos[[1]][length(pos[[1]])-1]+1,pos[[1]][length(pos[[1]])]-1)
family_order <- sapply(strsplit(foldername, "-"), "[")

new_relations<-matrix(,nrow = dim(relations)[1], ncol = dim(relations)[2])
for (i in 1:length(family_order))
{
  if(family_order[i]=="teyze")  
  {
    indices <- which(relations == 5)
  }else if(family_order[i]=="cocuk")  
  {
    indices <- which(relations == 3)
  }else if(family_order[i]=="baba")  
  {
    indices <- which(relations == 1)
  }else if(family_order[i]=="amca")  
  {
    indices <- which(relations == 2)
  } else
  {
    indices <- which(relations == 4)
  }
  print("i:")
  print(i)
  
  new_relations[indices] = length(family_order)-i+1
  print(length(family_order)-i+1)
}
relations <- t(apply(new_relations,1,sort))


new_genotype<-matrix(,nrow=nrow(genotype),ncol=ncol(genotype))
rev_fam_order <-rev(family_order)
for(i in 1:length(rev_fam_order))
{
  new_ind <- which(genotype_index[,2]==rev_fam_order[i])
  new_genotype[,i] <- genotype[,new_ind]
  print(new_ind)
}

genotype<-new_genotype

### end relationship matrix

#choice <- 0 # : outlierlar degisken
#choice <- 2 # : outlierlar sabit, phi = 2. problemin coz??m??nden bulunuyor. i.e. phi=0.76
#choice <- 1 # : phi kac olursa cozum olur icin constraint matrixi yaratma k??sm??
flag_10 <- FALSE
flag_11 <- FALSE
flag_12 <- FALSE


if(dim(relations)[1]>=2)
{ 
  removal_index <- 1 # kimden silineceginin indexi her zaman 1
  arr <- which(relations == removal_index, arr.ind = TRUE)
  
  #constraints  
  removal_rel <- relations[arr[,1],]
  
  # Create X matrix
  
  number_of_people<-length(files) # number of people in family
  indices<-1:number_of_people
  X_general <- c()
  
  
  if (is.null(nrow(removal_rel)))
  {
    removal_rel <- t(removal_rel)
  }
  
  #X matrixini doldurmak:
  for (z in 1:nrow(removal_rel))
  {
    index1<-removal_rel[z,1] 
    index2<-removal_rel[z,2]
    new_indices<-indices[-removal_rel[z,]] # ??b??r ki??ilerin indexi
    
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
  
  ## create objective matrix for optimization 
  Amat<-c()
  ## create RHS matrix for optimization
  bvec<-c()
  phi_rhs <-c()
  phi_x_array<- c()
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
    }  else
    {
      b<-g1
      a<-g2
    }
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
    
    #objective variables:
    related_X_indices <- get_relatedX(X_general, rel)
    related_X_indices <- sort(related_X_indices)
    related_x <- X_general[related_X_indices,]
    
    obj_arr<-matrix(0,nrow = nrow(X_general), ncol = 1, byrow = TRUE)
    phi_arr <- matrix(0,nrow = nrow(X_general), ncol = 1, byrow = TRUE)
    
    x_11 <- which(related_x[,i] == 1 & related_x[,j]==1 , arr.ind = TRUE)
    x_20 <- c( which(related_x[,i] == 2 & related_x[,j]==0 , arr.ind = TRUE), which(related_x[,i] == 0 & related_x[,j]==2, arr.ind = TRUE))
    
    p <- which(compMat2(related_x[x_11,],X_general))
    q <- which(compMat2(related_x[x_20,],X_general))
    
    obj_arr[p] <- obj_arr[p] + 2
    obj_arr[q] <- obj_arr[q] - 4
    
    
    #i > j
    if(g1 > g2)
    {
      x_1 <- which(related_x[,i] == 1 , arr.ind = TRUE)  # a
      x__1 <- which(related_x[,j] == 1  , arr.ind = TRUE)  #b
      s <- which(compMat2(related_x[x__1,],X_general))
      t <- which(compMat2(related_x[x_1,],X_general))
      
    } else
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
    #print(obj_arr)
    Amat <- rbind(Amat,t(obj_arr))
    
  }
  
  if(choice ==1)
  {
    df <-data.frame( ncol=2*nrow(X_general) +2 )
    df <- cbind(Amat,phi_x_array, t(t(bvec)), t(t(phi_rhs)))
    colnames(df) <- c(rep("x",nrow(X_general)), rep("phi_x",nrow(X_general)), "RHS", "phi_RHS")
  }
  
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
    
    # slack variablelar varm???? gibi ????z
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
              #o_row[17]<- (-1)
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
                print("index1 and index2:")
                print(c(index1,index2))
                print("n11 < o11")
                print(n11 < o11)
                flag_11 <- TRUE
                new_rhs <- 0
              }
              RHS2 <- rbind(RHS2,new_rhs) 
              
            }
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
            if((x[index1]==1 & x[index2]==2) || (x[index1]==2 & x[index2]==1))
            {
              n12 <- length(which(genotype[,index1]==x[index1] & genotype[,index2]==x[index2], arr.ind = TRUE))
              #o_row<-rep(0,18)
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
}else # regular removal, only one relationship exist
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
  
  if(choice == 0)
  {
    phi<-0
    RHS <- (2*d - 4*c -a + b - 4*b*phi)/2;
    
    gen_11 <- which(genotype[,1]==1 & genotype[,2]==1)
    gen_10 <- which(genotype[,1]==1 & genotype[,2]==0)
    gen_12 <- which(genotype[,1]==1 & genotype[,2]==2)
    Amat <- rbind(c(1,0),c(1,-1))
    bvec <- c(RHS, length(gen_11)-o11)
    cvec <- c(1,0)
    res<-Rcplex(cvec, Amat, bvec, Qmat = NULL, control=list(solnpoolintensity=4),
                objsense = c("min"), sense = c("G", "L"), vtype = rep("I",length(cvec)), n=1000)
    
    Amat2 <- rbind(c(1,0),c(1,-1),c(1,0))
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
    
    
    e11 <- res2[[1]]$obj
    x11<- res[[1]]$obj
    output <- data.frame(t(c(e10, e11, e12, x11)))
    colnames(output) <- c("e10", "e11", "e12", "x11")
    write.table(output,"out.txt", row.names = FALSE)
    
    #select randomly x11 SNPs
    
    snp_pos11 <- snp.pos[ which(genotype[,1]==1 & genotype[,2]==1) ]
    snp_pos11 <- sample(snp_pos11)
    snp_pos11 <- snp_pos11[1:x11]
    
    write.table(snp_pos11, "remove_snps.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, append = FALSE)
    
    
  }else  ## kinship degisken
  {
    n10 <- length(which(genotype[,1]==1 & genotype[,2]==0, arr.ind = TRUE))
    n11 <- length(which(genotype[,1]==1 & genotype[,2]==1, arr.ind = TRUE))
    n12 <- length(which(genotype[,1]==1 & genotype[,2]==2, arr.ind = TRUE))
    if( n10 < o10)
    {
      flag_10 <- TRUE
      print("n10 < o10")
    }
    
    if(n11 < o11)
    {
      flag_11 <-TRUE
      print("n11 < o11")
    }
    
    if(n12 < o12)
    {
      flag_12 <- TRUE
      print("n12 < o12")
    }
    
    if(choice == 2 & flag_10 == FALSE & flag_11 == FALSE & flag_12 == FALSE)
    {
        RHS_x11 <- (2*d-4*c-a+b-4*b*phi)/(2*(1-2*phi))
        Amat<-matrix(c(1,1), nrow = 2)
        bvec <- c(RHS_x11, n11-o11)
        res<-Rcplex(c(1), Amat, bvec, Qmat = NULL, control=list(solnpoolintensity=4),
                    objsense = c("min"), sense = c("G","L"), vtype = c("I"), n=500)
    
        

        pos_11 <- snp.pos[sample(which(genotype[,1]==1 & genotype[,2] ==1 ))]

        snp.pos_11 <-pos_11[1:res[[1]]$obj[1]]

        
        write.table(snp.pos_11,"remove_snps.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, append = FALSE) 
        
    }else if(choice == 1 & flag_10 == FALSE & flag_11 == FALSE & flag_12 == FALSE)
    {
      
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


#print the coefficients run the matlab code
if(choice==1 & length(files) > 2)
{
  new_size_matrix <- c()
  for(index in 1:dim(X_general)[1])
  {
    new_size_matrix<- c(new_size_matrix,min(RHS2[which(other_matrix[,index]==1)]))
  }
  write.table(new_size_matrix,"size_constraints.txt",col.names = FALSE, row.names = FALSE)  
  write.table(df,"kin_constraints.txt",col.names = FALSE, row.names = FALSE)
  if(flag_10 == TRUE | flag_11 == TRUE | flag_12 == TRUE)
  {
    write.table(matrix(c("n10 < o10", flag_10,"n11 < o11", flag_11, "n12 < o12", flag_12), ncol = 3, byrow = TRUE),"flag_constraints.txt", sep = "\t", col.names = FALSE, row.names = FALSE) 
  }
}


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
    write.table(remove_snps, "remove_snps.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, append = FALSE)
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
  
  #minimize x1 + x2 + x3 slack variables
  cvec <- c(rep(1,nrow(X_general)),0,0,0)
  res<-Rcplex(cvec, Amat, bvec, Qmat = NULL, control=list(),
              objsense = c("min"), sense = c(rep("G", nrow(relations) ), rep("L",length(RHS2))), vtype = rep("I",length(cvec)))
  
  print(res)
  #minimize e1 + e2 + e3 slack variables
  Amat_yeni <- rbind(Amat, c(rep(1,nrow(X_general)),0,0,0))
  bvec_yeni <- c(bvec, res$obj)
  cvec_yeni <- c(rep(0,nrow(X_general)),1,1,1)
  res2 <- Rcplex(cvec_yeni,Amat_yeni, bvec_yeni, Qmat = NULL, control=list(solnpoolintensity=4),
                 objsense = c("min"), sense = c(rep("G", nrow(relations) ), rep("L",length(RHS2)), "E"), vtype = rep("I",length(cvec)), n=5000)
  
  print(res2)
  if(length(res)==0 | length(res2)==0)
  {
    print("no solut??on ex??st")
  }else
  {
    write.table(res2[[1]]$xopt, "out.txt")
    #removal of snps
    if(dim(X_general)[1] == 3)
    {
      pos_110 <- snp.pos[sample(which(genotype[,1]==1 & genotype[,2] ==1 & genotype[,3]==0))]
      pos_111 <- snp.pos[sample(which(genotype[,1]==1 & genotype[,2] ==1 & genotype[,3]==1))]
      pos_112 <- snp.pos[sample(which(genotype[,1]==1 & genotype[,2] ==1 & genotype[,3]==2))]
      snp.pos_110 <-pos_110[1:res2[[2]]$xopt[1]]
      snp.pos_111 <-pos_111[1:res2[[2]]$xopt[2]]
      snp.pos_112 <-pos_112[1:res2[[2]]$xopt[3]]
      remove_snps <- c(snp.pos_110,snp.pos_111, snp.pos_112)
      write.table(remove_snps,"remove_snps.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, append = FALSE) 
    }
    else{
      print("4. kisi ekleneillir, kodu modifiye et.")
    }
  }
  
}


### END ###

