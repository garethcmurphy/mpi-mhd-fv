#!/bin/bash
if [ "$1" = "" ] ||  [ "$2" = "" ]
then
echo Usage: $0 filename element-name
file=out_0001
element=1
	
else
file=$1
export arg1=`printf "%04d" $1`
file=out_$arg1.dat
element=$2
fi

cat << EOF  > gnu.plt
#set terminal postscript color enhanced solid lw 2 "Times-Roman" 10
set terminal jpeg 
set output 'file$arg1.jpeg' 
#set output 'file$arg1.ps' 
#set size ratio 2
set size square
set pm3d map
set palette defined (-8 "blue", -6 "white", 1 "red")
#set contour
#set grid back
#set cntrparam levels 100
#set samples 100
#set contour base
splot '$file' using $element
#pause -1
EOF
gnuplot gnu.plt
