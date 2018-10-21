library(gdsfmt)
library(SNPRelate)


files<-list.files(path=getwd(),pattern=".vcf")

allGenos <- c()
count = 0
for (f in files)
{
  vcf<-read.table(f)
  scores<-vcf$V10
  y <- vector(mode="numeric", length=length(scores))
  for (i in 1:length(scores))
  {
    str<-scores[i]
    str1<-substring(str, 1, 1)
    str2<-substring(str,3,3)
    if(str1==0 && str2==0)
      score=0
    if(str1==0 && str2==1)
      score=1
    if(str1==1 && str2==0)
      score=1
    if(str1==1 && str2==1)
      score=2
    if(str1==".")
      score=3
    
    y[i]=score
    
    
  }
  count = count + 1
  genotype=data.matrix(y)
  allGenos <- cbind(allGenos, genotype)
}


snp.chromosome = vcf$V1
snp.position=vcf$V2
sample.id = 1:dim(allGenos)[2]
snp.id=1:length(snp.chromosome)


snpgdsCreateGeno("c.gds", genmat = allGenos, sample.id=sample.id, snp.id = snp.id , snp.chromosome =snp.chromosome, snp.position = snp.position, snp.allele =NULL, snpfirstdim=TRUE)


genof <- snpgdsOpen("c.gds")
#genof <- snpgdsOpen("10mb.gds")
ibs<-snpgdsIBS(genof)
ibs$ibs

genotype <- read.gdsn(index.gdsn(genof, "genotype"))
ibd.robust <- snpgdsIBDKING(genof, sample.id=1:length(files), family.id=NULL)
kinship=ibd.robust$kinship
kinship
gen11<-which(genotype[,1]==1 & genotype[,2]==1,arr.ind = TRUE)
print("[1 1]")


snp.id <- read.gdsn(index.gdsn(genof, "snp.id"))
ibd.robust <- snpgdsIBDKING(genof, sample.id=1:length(files), family.id=NULL)
snp.id<-ibd.robust$snp.id

genotype<-genotype[snp.id,]
gen11<-which(genotype[,1]==1 & genotype[,2]==1,arr.ind = TRUE)
print("[1 1]")
print(length(gen11))