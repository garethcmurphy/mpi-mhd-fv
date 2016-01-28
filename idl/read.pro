file = DIALOG_PICKFILE(/READ, FILTER = "*.h5")
;file = 'hdf5_0000.h5'


;n=i
;num=strtrim(n,4)
;num=STRING(i, FORMAT='(I4.4)')
;file = 'hdf5_'+num+'.h5'

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

	v2= velx^2 +  vely^2 + velz^2
	bsquare= b1^2 + b2^2 + b3^2
	;
	pressure=energy - 0.5*density*(v2) - 0.5*(bsquare)
	pressure=pressure*3.0/5.0
	end
