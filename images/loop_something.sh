#!/bin/bash

i="10"

while [ $i -lt 4400 ]
do
   #wget http://www.compsoc.com/~gmurphy/chombo_plt_cnt_01$i.hdf5.gz
#	echo ./print.sh $i 14
#	./print.sh $i 14
	export arg1=`printf "%04d" $i`
	#convert  file$arg1.png file$arg1.ppm
	convert  file$arg1.png file$arg1.jpeg

i=$[$i+10]
done
