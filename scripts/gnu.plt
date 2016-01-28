#set terminal postscript enhanced color 
#set output 'file0001.ps' 
#set size ratio  2.0
set size  square
set pm3d map
#set cntrparam levels 35
#set contour base
#show contour
#set grid front
#set cntrparam levels 100
#set samples 100
#set contour base
splot 'out_0001.dat' using 15
pause -1
