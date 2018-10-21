 #!/bin/bash   


for f in ` ls kesisim/*`;
do
  filename=$(basename "$f")
  filename="${filename%.*}"
  perl 23andme2vcfyeni.pl $f "vcf/$filename.vcf"
done


