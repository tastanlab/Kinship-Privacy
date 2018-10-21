
files<-list.files(path=".",pattern=".23andme.txt")

file1=files[1]
f1=read.table(file1,fill=TRUE)
kesisim=as.character(f1$V1)
count=0
kesisim2=kesisim
for (f in files)
{
    count=count+1
    data1<-read.table(f,fill=TRUE)
	  rsid=as.character(data1$V1)
	  kesisim2=kesisim
		kesisim=intersect(kesisim,rsid)
    f2=f
}


write.table(kesisim, "kesisim.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, append=FALSE)
