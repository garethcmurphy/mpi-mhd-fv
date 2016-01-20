
  pro test

  array = intarr(6)
  infile = 'test.txt'
  openr, inlun, infile, /get_lun
  readf, inlun, array
  close, inlun
  free_lun, inlun
  print, array



Device, Decomposed=0, Retain=2
levels = 64

  infile = 'density.dat'
  openr, inlun, infile, /get_lun

!x.style = 1
!y.style = 1

  readf, inlun, nx
  readf, inlun, ny
!x.range = [0.,nx-1]
!y.range = [0.,ny-1]

Window, 0, Title='Density',xsize=700,ysize=800


LoadCT, 33
white = GetColor('White', 1)
black = GetColor('Black', 2)
LoadCT, 33, NColors=64, Bottom=3


  array2 = fltarr(nx,ny)
  readf, inlun, array2
  close, inlun
  free_lun, inlun
  shade_Surf, array2 ,$
  AX = 090, AZ = 000, $
        Position=[ 0.1,0.1, 0.98, 0.8], $
		          font=6, Background=1, Color=black, shades = bytscl(array2)

		  
  ColorBar, NColors=64, Bottom=3, Divisions=6, $
     Range=[Min(array2), Max(array2)], Format='(G10.3)', $
	        Position=[ 0.1,0.9, 0.98, 0.98], $
			                  Charsize=2.0,$
									                Color=black
														 


  return
  end
