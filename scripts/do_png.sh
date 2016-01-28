#!/bin/sh
gnuplot pm3d.15.gnu

convert pm3d.15.ps pm3d.15.imagemagick.png 

# Ghostscript is long-winded but high-quality. Note use of same postscript source each time.

gs -dSAFER -dBATCH -dNOPAUSE -sDEVICE=png16m -dTextAlphaBits=4 -dGraphicsAlphaBits=4 -r36 -sOutputFile=pm3d.15.ghostscript.small.png pm3d.15.ps
gs -dSAFER -dBATCH -dNOPAUSE -sDEVICE=png16m -dTextAlphaBits=4 -dGraphicsAlphaBits=4 -r72 -sOutputFile=pm3d.15.ghostscript.med.png pm3d.15.ps
gs -dSAFER -dBATCH -dNOPAUSE -sDEVICE=png16m -dTextAlphaBits=4 -dGraphicsAlphaBits=4 -r144 -sOutputFile=pm3d.15.ghostscript.large.png pm3d.15.ps
gs -dSAFER -dBATCH -dNOPAUSE -sDEVICE=png16m -dTextAlphaBits=4 -dGraphicsAlphaBits=4 -r300 -sOutputFile=pm3d.15.ghostscript.enormous.png pm3d.15.ps
