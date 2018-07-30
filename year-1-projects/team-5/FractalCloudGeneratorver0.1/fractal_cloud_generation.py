#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
********************************************
Created on Sun Mar 18 22:47:42 2018
by
Chamara Rajapakshe
(cpn.here@umbc.edu)
********************************************
Generating fractal clouds. 
"""
import h5py, os
import numpy as np


class frac_physics(object):
    def __init__(self,fname,dpath=None):
        if dpath==None:
            self.fname=fname
        else:
            self.fname=dpath+fname
    def generate_fractal(self,re,lwp,f=0.5,xdist=None,xorder=12):
        '''
        Generate 1D Fractal cloud - Cahalan et. al. (1994)
        re: re in mu
        lwp: LWP in g/m^2 , 90 for stratocumulus cloud
        f: Fractal parameter fn=fc^n
        xdist: horizontal distance of the cloud in km
        xorder: Order of the fractal cloud
        '''
        self.f=f
        self.orderX=xorder
        bc=lwp
        varReduction=0.8#c=2^{-1/3}=0.8
        var=self.f#fractal parameter fn=fc^n;
        def next_step_1d(bc,var):
            l=np.size(bc)
            bc2=np.zeros(l*2)
            sign=np.random.uniform(-1,1,size=(l))
            sign=sign/abs(sign)
            bc2[0:2*l:2]=bc+sign*var*bc
            bc2[1:2*l:2]=bc-sign*var*bc
            return bc2
            
        for i in np.arange(0,self.orderX):
            bc=next_step_1d(bc,var)
            var=var*varReduction
        self.lwp=bc
        self.tau=bc*3.0/2.0/re
        self.y=np.array([1.0])
        self.re=np.repeat(re,len(bc))
        if xdist==None:
            self.x=np.linspace(0,len(bc)/100.0,len(bc)+1)#let's say this is 0.01km resolution
        else:
            self.x=np.linspace(0,xdist,len(bc)+1)
        self.x=(self.x[0:-1]+self.x[1:])/2#assigning the middle
    def save_frach5(self,fname=None):
        '''
        Save fractal cloud properties as an hdf5 file.
        '''
        if not(fname==None):
            self.fname=fname
        
        if os.path.isfile(self.fname):
            rd=raw_input('File already exist. Replace?: ')
        else:
            rd='y'
        if rd=='y':
            f=h5py.File(self.fname,'w')
            PCentry=f.create_dataset('tau',data=self.tau)
            PCentry.dims[0].label='x'
            PCentry.attrs["long_name"]='Cloud_optical_thickness'
            
            PCentry=f.create_dataset('orderX',data=self.orderX)
            PCentry.attrs["long_name"]='x_dimension_order'
            
            PCentry=f.create_dataset('lwp',data=self.lwp)
            PCentry.dims[0].label='x'
            PCentry.attrs["long_name"]='liquid_water_path'
            PCentry.attrs["unit"]='g/m^2'
            
            PCentry=f.create_dataset('re',data=self.re)
            PCentry.dims[0].label='x'
            PCentry.attrs["long_name"]='effective_radius'
            PCentry.attrs["unit"]='microns'
            
            PCentry=f.create_dataset('x',data=self.x)
            PCentry.attrs["long_name"]='x_dimension'
            PCentry.attrs["unit"]='km'
            
            f.create_dataset('y',data=self.y)
            f.create_dataset('f',data=self.f)
            
            f.close()
            print(self.fname+' SAVED!')
    def read_frach5(self,fname=None):
        if not(fname==None):
            self.fname=fname
        f=h5py.File(self.fname,'r')
        self.re=f['re'][:]
        self.tau=f['tau'][:]
        self.f=f['f']
        self.orderX=f['orderX']
        self.lwp=f['lwp'][:]
        self.x=f['x'][:]
        self.y=f['y'][:]
        f.close()
        
    def plot_lwp_pdf(self,):
        lwp,lwpax=plt.subplots()
        weights = np.ones_like(self.lwp)/len(self.lwp)
        lwpax.hist(self.lwp,bins = 10 ** np.linspace(1, 3, 50),weights=weights,histtype=u'step')
        #val,bins=np.histogram(self.lwp,bins = 10 ** np.linspace(1, 3, 50),weights=weights)
        #width = np.diff(bins)
        #center = (bins[:-1] + bins[1:]) / 2
        #lwpax.bar(center,val,align='center',width=width)
        lwpax.set_xscale('log')
        #x=np.linspace(ss.lognorm.ppf(0.01,0.954),ss.lognorm.ppf(0.99,0.954),100)
        #lwpax.plot(x,ss.lognorm.pdf(x,0.954))
        lwp.show()
        return lwp,lwpax

    
    def plot_cld(self,fig=None,ax=None):
        
        if ax==None:
            fig,ax=plt.subplots()
        ax.plot(self.x,self.re,'r',label='re')
        ax.set_ylabel(r'Effective radius ($\mu$)')
        ax2=ax.twinx()
        ax2.set_ylabel(r'Optical depth')
        ax2.plot(self.x,self.tau,'g',label=r'$\tau$')
        ax.legend(loc='best',prop={'size':10})
        ax.set_xlabel('km')
        fig.show()
        
        return fig
 

#Uncomment below to run example 
#lwp=90#gm^{-2}
#re=12.0#microns
#fcld=frac_physics('fractal_cloud_example.nc')
#fcld.f=0.37#fractal parameter fn=fc^n
#fcld.generate_fractal(re,lwp,xorder=12,xdist=40.0)
#fcld.plot_lwp_pdf()
#fcld.plot_cld()
