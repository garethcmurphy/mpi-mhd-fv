
#mencoder "mf://*.jpeg" -mf fps=10 -o test.avi -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=800 
 #mencoder mf://*.jpeg -mf w=800:h=600:fps=25:type=jpg -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o output.avi
 mencoder mf://*.jpeg -mf w=800:h=600:fps=25:type=jpg -ovc lavc -lavcopts vcodec=msmpeg4v2:mbd=2:trell -oac copy -o output.avi

