#!/bin/bash



#!/bin/bash

i="20"

while [ $i -lt 1920 ]
do
   #wget http://www.compsoc.com/~gmurphy/chombo_plt_cnt_01$i.hdf5.gz
	echo ./print.sh $i 14
	./print.sh $i 14

	gs -dSAFER -dBATCH -dNOPAUSE -sDEVICE=png16m - -r300 -sOutputFile=pm3d.15.ghostscript.enormous.png pm3d.15.ps
i=$[$i+20]
done
