#set terminal png transparent nocrop enhanced font arial 8 size 420,320 
#set output 'pm3d.15.png'
set terminal postscript color enhanced
set output 'pm3d.15.ps'
set border 4095 lt -1 lw 1.000
set view 130, 10, 1, 1
set samples 50, 50
set isosamples 50, 50
unset surface
set title "\"set pm3d scansbackward\" makes this as viewed from above" 0.000000,0.000000  font ""
set pm3d at s
set pm3d scansbackward flush begin noftriangles nohidden3d transparent implicit corners2color mean
splot sin(sqrt(x**2+y**2))/sqrt(x**2+y**2)

