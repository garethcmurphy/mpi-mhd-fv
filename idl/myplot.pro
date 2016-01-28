file = DIALOG_PICKFILE(/READ, FILTER = "*.h5")
;file = 'hdf5_out_2d_00026.h5'
file_id = H5F_OPEN(file)
print, file
 dataset_id1 = H5D_OPEN(file_id,'Density')
  space_id = H5D_GET_SPACE(dataset_id1)
   dims = H5S_GET_SIMPLE_EXTENT_DIMS(space_id)
	 density = H5D_READ(dataset_id1)
;	 density = alog10(density)


 dataset_id1 = H5D_OPEN(file_id,'Pressure')
  space_id = H5D_GET_SPACE(dataset_id1)
   dims = H5S_GET_SIMPLE_EXTENT_DIMS(space_id)
	 temperature = H5D_READ(dataset_id1)

 dataset_id1 = H5D_OPEN(file_id,'Velx')
  space_id = H5D_GET_SPACE(dataset_id1)
   dims = H5S_GET_SIMPLE_EXTENT_DIMS(space_id)
	 velx = H5D_READ(dataset_id1)
	 velx=velx/density;

 dataset_id1 = H5D_OPEN(file_id,'Vely')
  space_id = H5D_GET_SPACE(dataset_id1)
   dims = H5S_GET_SIMPLE_EXTENT_DIMS(space_id)
	 vely = H5D_READ(dataset_id1)
	 vely=vely/density;

 dataset_id1 = H5D_OPEN(file_id,'Velz')
  space_id = H5D_GET_SPACE(dataset_id1)
   dims = H5S_GET_SIMPLE_EXTENT_DIMS(space_id)
	 velz = H5D_READ(dataset_id1)
	 velz=velz/density;

 dataset_id1 = H5D_OPEN(file_id,'Energy')
  space_id = H5D_GET_SPACE(dataset_id1)
   dims = H5S_GET_SIMPLE_EXTENT_DIMS(space_id)
	 energy = H5D_READ(dataset_id1)

 dataset_id1 = H5D_OPEN(file_id,'Bx')
  space_id = H5D_GET_SPACE(dataset_id1)
   dims = H5S_GET_SIMPLE_EXTENT_DIMS(space_id)
	 Bx = H5D_READ(dataset_id1)

 dataset_id1 = H5D_OPEN(file_id,'By')
  space_id = H5D_GET_SPACE(dataset_id1)
   dims = H5S_GET_SIMPLE_EXTENT_DIMS(space_id)
	 By = H5D_READ(dataset_id1)


 dataset_id1 = H5D_OPEN(file_id,'Bz')
  space_id = H5D_GET_SPACE(dataset_id1)
   dims = H5S_GET_SIMPLE_EXTENT_DIMS(space_id)
	 Bz = H5D_READ(dataset_id1)

H5D_OPEN(file_id)

device, decomposed=0, retain=2
Window, 0, Title='Density Plot',xsize=800,ysize=800


rescol



;white = GetColor('White', 1)
;black = GetColor('Black', 2)
;LoadCT, 33, NColors=64

!p.multi=[0,3,3]
plot, density(5,*), ytitle='density' ,psym=4
plot, temperature(5,*), ytitle='Pressure' ,psym=4
plot, velx(5,*), ytitle='Velx' ,psym=4
plot, vely(5,*), ytitle='Vely' ,psym=4
plot, velz(5,*), ytitle='Velz' ,psym=4
plot, energy(5,*), ytitle='Energy' ,psym=4
plot, Bx(5,*), ytitle='Bx' ,psym=4
plot, By(5,*), ytitle='By' ,psym=4
plot, Bz(5,*), ytitle='Bz' ,psym=4


END
