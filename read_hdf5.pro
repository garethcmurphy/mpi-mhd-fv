file = DIALOG_PICKFILE(/READ, FILTER = "*.h5")
;file = 'hdf5_out_2d_00050.h5'
file_id = H5F_OPEN(file)
print, file

 dataset_id1 = H5D_OPEN(file_id,'Density')
  space_id = H5D_GET_SPACE(dataset_id1)
   dims = H5S_GET_SIMPLE_EXTENT_DIMS(space_id)
	 density = H5D_READ(dataset_id1)
	density=transpose(density)

 dataset_id1 = H5D_OPEN(file_id,'Velx')
  space_id = H5D_GET_SPACE(dataset_id1)
   dims = H5S_GET_SIMPLE_EXTENT_DIMS(space_id)
	 velx = H5D_READ(dataset_id1)
	velx=transpose(velx)

 dataset_id1 = H5D_OPEN(file_id,'Vely')
  space_id = H5D_GET_SPACE(dataset_id1)
   dims = H5S_GET_SIMPLE_EXTENT_DIMS(space_id)
	 vely = H5D_READ(dataset_id1)
	vely=transpose(vely)


 dataset_id1 = H5D_OPEN(file_id,'Velz')
  space_id = H5D_GET_SPACE(dataset_id1)
   dims = H5S_GET_SIMPLE_EXTENT_DIMS(space_id)
	 velz = H5D_READ(dataset_id1)
	velz=transpose(velz)

 dataset_id1 = H5D_OPEN(file_id,'Energy')
  space_id = H5D_GET_SPACE(dataset_id1)
   dims = H5S_GET_SIMPLE_EXTENT_DIMS(space_id)
	 energy = H5D_READ(dataset_id1)
	energy=transpose(energy)


 dataset_id1 = H5D_OPEN(file_id,'Bx')
  space_id = H5D_GET_SPACE(dataset_id1)
   dims = H5S_GET_SIMPLE_EXTENT_DIMS(space_id)
	 b1 = H5D_READ(dataset_id1)
	b1=transpose(b1)

 dataset_id1 = H5D_OPEN(file_id,'By')
  space_id = H5D_GET_SPACE(dataset_id1)
   dims = H5S_GET_SIMPLE_EXTENT_DIMS(space_id)
	 b2 = H5D_READ(dataset_id1)
	b2=transpose(b2)

 dataset_id1 = H5D_OPEN(file_id,'Bz')
  space_id = H5D_GET_SPACE(dataset_id1)
   dims = H5S_GET_SIMPLE_EXTENT_DIMS(space_id)
	 b3 = H5D_READ(dataset_id1)
	b3=transpose(b3)

device, decomposed=0, retain=2
Window, 0, Title='Density Plot',xsize=800,ysize=400


rescol
LoadCT, 33

shade_surf, density, $
            shades=bytscl(density), $
				Background=1, $
				position=[0.1,0.1,0.8,0.9],$
				Color=black, $
            AX = 090, AZ = 000


;white = GetColor('White', 1)
;black = GetColor('Black', 2)
LoadCT, 33, NColors=64, Bottom=3
;LoadCT, 33, NColors=64

ColorBar, NColors=64, Bottom=3, Divisions=7, $
   Range=[Min(density), Max(density)], Format='(G10.2)', $
	position=[0.8,0.1,0.9,0.9],$
   Charsize=3., $
   Vertical=1, $
   Color=black


;.r resize_draw
;draw_resize, density, 'density'
;draw_resize, pressure, 'Pressure'
;draw_resize, velx, 'velx'
;draw_resize, vely, 'vely'
;draw_resize, energy, 'Energy'

END
