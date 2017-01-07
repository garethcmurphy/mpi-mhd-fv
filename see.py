import h5py
f = h5py.File('hdf5_0010.h5','r')
import matplotlib 
%matplotlib inline
print f['Bx']
import matplotlib.pyplot as plt

bx=f['Bx']
subplot(1,2,2)
plt.subimage(f['Bx'])
plt.subimage(f['By'])
