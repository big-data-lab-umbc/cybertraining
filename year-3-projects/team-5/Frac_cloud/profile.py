import numpy as np
import os
import matplotlib.pyplot as plt
from frac_cloud import frac_physics,readSHDOMRad,setProp


def wrap(re,lwp,fn,out_dir):
    fc = frac_physics(fn,dpath=None)
    fc.generate_fractal(re,lwp,xorder=12)
    fc.save_frach5("%s/%s" %(out_dir,fn)) 
    
if __name__ == "__main__":
  out_dir = "./profile" 
  
  sub_size = 10
  nre, nlwp = 20, 20
  re_0  = np.linspace(5,20,nre)
  lwp_0 = np.logspace(1,2,nlwp)
  
  #d_re  = re_0[1]/re_0[0]
  d_re = 1.06
  d_lwp = lwp_0[1]/lwp_0[0]
  
  
  for i in range(sub_size):
    for j in range(nre):
      for k in range(nlwp):  
        fname = "profile_%05d" %(i*nre*nlwp+ j*nlwp + k + 1)
        rn1 = np.random.uniform(1/d_re,d_re,1)
        rn2 = np.random.uniform(1/d_lwp,d_lwp,1)

        #print(re_0[j]*rn1,lwp_0[k]*rn2)
        wrap(re_0[j]*rn1,lwp_0[k]*rn2,fname,out_dir)    
    
    
    
