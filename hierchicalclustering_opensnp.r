library(gdsfmt)
library(SNPRelate)


files<-list.files(path=getwd(),pattern=".vcf")

#create users x snps x genotype tensor
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

#create gds file
snpgdsCreateGeno("c.gds", genmat = allGenos, sample.id=sample.id, snp.id = snp.id , snp.chromosome =snp.chromosome, snp.position = snp.position, snp.allele =NULL, snpfirstdim=TRUE)
genof <- snpgdsOpen("c.gds")

#hierchical clustering on OpenSNP data
genotype <- read.gdsn(index.gdsn(genof, "genotype"))
diss <- snpgdsDiss(genof)
hc <- snpgdsHCluster(diss)
rv <- snpgdsCutTree(hc,z.threshold=15,outlier.n=10)
rv$clust.count

plot(rv$dendrogram, leaflab="none", main="Kinship")
