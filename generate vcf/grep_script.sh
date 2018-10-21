#!/bin/bash

for f in ` ls ilkhalleri/*`;
do
  filename=$(basename "$f")
  #filename="${filename%.*}"  
  grep -wFf kesisim.txt $f > $filename
  #grep "MT" $f > $filename
done
