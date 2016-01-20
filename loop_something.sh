#!/bin/bash

i="20"

while [ $i -lt 900 ]
do
   #wget http://www.compsoc.com/~gmurphy/chombo_plt_cnt_01$i.hdf5.gz
	echo ./print.sh $i 
	./print.sh $i 3

i=$[$i+20]
done
