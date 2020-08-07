from frac_cloud import frac_physics,readSHDOMRad
import matplotlib.pyplot as plt
import numpy as np

fc = frac_physics('1dfract',dpath=None)
#fc.read_frach5('profile_00001.hdf5')
fc.read_frach5('spot3D_fractal.hdf5')

Rad,VZA,VAA,SZA,SAA = readSHDOMRad("spot3D0p865CER12r.out")

ref = np.pi*Rad[0,:]/np.deg2rad(SZA)
plt.figure(1)
plt.plot(fc.x,ref)

plt.figure(2)
plt.plot(fc.x,fc.tau)

plt.figure(3)
plt.plot(fc.x,ref/fc.tau)
plt.show()

