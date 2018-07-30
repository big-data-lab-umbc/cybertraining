#Import general functions
import h5py, os
import numpy as np
import sys

# Import class from fractal cloud generator
from fractal_cloud_generation import frac_physics

print('Number of arguments:', len(sys.argv), 'arguments.')
print('Argument List:', str(sys.argv))

re = float(sys.argv[1])
lwp = float(sys.argv[2])
xdist = float(sys.argv[3])
filenm = sys.argv[4]

if len(sys.argv)>5:
	rseed = int(sys.argv[5])
	print(rseed)


def fractal_cloud(re,lwp,xdist,filenm,rseed):
        
        print('Number of arguments:', len(sys.argv), 'arguments.')
        print('Argument List:', str(sys.argv))

        '''
        Generate 1D Fractal cloud - Cahalan et. al. (1994)
        re: re in mu
        lwp: LWP in g/m^2 , 90 for stratocumulus cloud
        f: Fractal parameter fn=fc^n
        xdist: horizontal distance of the cloud in km
        xorder: Order of the fractal cloud
        rseed: Random number generator seed. 
        filenm: File name to save as .nc
        '''
        # Define constants
        f = 0.5
        xorder = 12
        
	fcld = frac_physics(filenm)
	fcld.f = 0.37 #what is this?
       
	 # Change seed if you want to reproduce
        np.random.seed(rseed)
        
        # Define fractal cloud
        fcld.generate_fractal(re,lwp,xorder,xdist)
        
        #Save
        fcld.save_frach5(filenm)
        
        # return
        return 1
        
fractal_cloud(re,lwp,xdist,filenm,rseed)
