#!/bin/bash

for f in *;
do
  sed --in-place '/D/d' $f
  sed --in-place '/--/d' $f
  sed --in-place '/I/d' $f
   
done
