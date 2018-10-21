#!/bin/bash

mv kesisim.r ilkhalleri
mv sed.sh ilkhalleri
(cd ilkhalleri; sh sed.sh
Rscript kesisim.r
mv kesisim.txt ..) 

mkdir kesisim
./grep_script.sh
mv *.txt kesisim
mv kesisim/kesisim.txt ..
mkdir vcf
chmod +x generatevcf.sh
./generatevcf.sh
mv clustering.r vcf
(cd vcf; Rscript clustering.r)
