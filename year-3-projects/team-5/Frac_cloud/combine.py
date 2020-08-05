from frac_cloud import frac_physics,readSHDOMRad
import matplotlib.pyplot as plt
import numpy as np
import h5py

for pid in range(4000,4001):
    fc = frac_physics('1dfract',dpath=None)
    fc.read_frach5('/home/kirana/cybertraining2020_team5/research/RT/profile/profile_%05d.hdf5' %pid)

    cer = fc.re[0];
    sfname1 = "profile_%05d_WVL_%01.2f_CER_%02.2f_SZA_%02.1f" %(pid,0.645,cer,60.)
    sfname2 = "profile_%05d_WVL_%01.2f_CER_%02.2f_SZA_%02.1f" %(pid,0.555,cer,60.)
    sfname3 = "profile_%05d_WVL_%01.2f_CER_%02.2f_SZA_%02.1f" %(pid,0.469,cer,60.)

    Rad1,VZA,VAA,SZA,SAA = readSHDOMRad("./TestOut/" + sfname1.replace(".","p") + "r.out")
    Rad2,VZA,VAA,SZA,SAA = readSHDOMRad("./TestOut/" + sfname2.replace(".","p") + "r.out")
    Rad3,VZA,VAA,SZA,SAA = readSHDOMRad("./TestOut/" + sfname3.replace(".","p") + "r.out")
    ref1 = np.pi*Rad1[0,:]/abs(np.cos(np.deg2rad(SZA)))
    ref2 = np.pi*Rad2[0,:]/abs(np.cos(np.deg2rad(SZA)))
    ref3 = np.pi*Rad3[0,:]/abs(np.cos(np.deg2rad(SZA)))

    f = h5py.File('/home/kirana/cybertraining2020_team5/research/DATA/Out1/profile_%05d.hdf5' %pid,'w')
    PCentry=f.create_dataset('tau',data=fc.tau)
    PCentry.dims[0].label='x'
    PCentry.attrs["long_name"]='Cloud_optical_thickness'
    
    PCentry=f.create_dataset('lwp',data=fc.lwp)
    PCentry.dims[0].label='x'
    PCentry.attrs["long_name"]='liquid_water_path'
    PCentry.attrs["unit"]='g/m^2'
    
    PCentry=f.create_dataset('re',data=fc.re)
    PCentry.dims[0].label='x'
    PCentry.attrs["long_name"]='effective_radius'
    PCentry.attrs["unit"]='microns'
    
    PCentry=f.create_dataset('red',data=ref1)
    PCentry.dims[0].label='x'
    PCentry.attrs["long_name"]='reflectance_red'
    
    PCentry=f.create_dataset('grn',data=ref2)
    PCentry.dims[0].label='x'
    PCentry.attrs["long_name"]='reflectance_green'
    
    PCentry=f.create_dataset('blu',data=ref3)
    PCentry.dims[0].label='x'
    PCentry.attrs["long_name"]='reflectance_blue'
    
    PCentry=f.create_dataset('x',data=fc.x)
    PCentry.attrs["long_name"]='x_dimension'
    PCentry.attrs["unit"]='km'
    
    f.create_dataset('y',data=fc.y)
    
    f.close()
