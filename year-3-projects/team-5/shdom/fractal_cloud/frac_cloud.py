#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
********************************************
Created on Tue Mar 24 07:07:16 2020
by
Chamara Rajapakshe
(cpn.here@umbc.edu)
********************************************
- frac_physics(): Fractal cloud generation and etc.
- readSHDOMRad() reads SHDOM radiance output for fractal cloud 
- setProp() writes *.prp files for SHDOM

"""
import numpy as np
import matplotlib.pyplot as plt
import h5py,os
class frac_physics(object):
    '''
    To generate 1D and 2D fracta clouds based on Cahalan et. al. (1994)
    fname: filename
    dpath: path
    dims: 1 or 2 (1D or 2D cloud)
    '''
    def __init__(self,fname,dpath=None,dims=1):
        self.dims=dims
        if dpath==None:
            self.fname=fname
        else:
            self.fname=dpath+fname
    def read(self,):
        # Me Dange fractal cloud eka kiyawanta
        f=h5py.File(self.fname,'r')
        self.re=f['re'][0]
        self.tau=f['tau'][0]
        self.ve=f['ve'][0]
        self.x=f['x'][0]
        self.y=f['y'][0]
        f.close()
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
    def generate_fractal2D(self,re,lwp,f=0.5,xydist=None,xyorder=12):
        '''
        **** NOT TESTED !!! DON'T USE !! ****
        Generate 2D Fractal cloud. Extending same theory in generate_fractal()
        re: re in mu
        lwp: LWP in g/m^2 , 90 for stratocumulus cloud
        f: Fractal parameter fn=fc^n
        xydist: horizontal distance of the cloud in km
        xyorder: Order of the fractal cloud
        '''
        self.f=f
        self.orderXY=xyorder
        bc=lwp
        varReduction=0.8#c=2^{-1/3}=0.8
        var=self.f#fractal parameter fn=fc^n;
        def next_step_2D(bc,var):
            l=np.size(bc)
            xy=int(np.sqrt(l))
            bc2=np.zeros((xy*2,xy*2),dtype=float)
            bc3=np.zeros((xy*2,xy*2),dtype=float)
            #4 random signs
            sign1=np.random.uniform(-1,1,size=(xy,xy))
            sign1=sign1/abs(sign1)
            sign2=np.random.uniform(-1,1,size=(xy,xy))
            sign2=sign2/abs(sign2)
            sign3=np.random.uniform(-1,1,size=(xy,xy))
            sign3=sign3/abs(sign3)
            sign4=np.random.uniform(-1,1,size=(xy,xy))
            sign4=sign4/abs(sign4)
            
            sign4expand=np.zeros((xy,xy),dtype=int)
            for i in np.arange(0,xy):
                for j in np.arange(0,xy):
                    sign4expand[i,j]=sign4[i/2,j/2]

            bc2[0:2*l:2,0:2*l:2]=bc+sign1*var*bc+sign2*var*var*bc
            bc2[0:2*l:2,1:2*l:2]=bc+sign1*var*bc-sign2*var*var*bc
            bc2[1:2*l:2,0:2*l:2]=bc-sign1*var*bc+sign3*var*var*bc
            bc2[1:2*l:2,1:2*l:2]=bc-sign1*var*bc-sign3*var*var*bc

            bc3[0:2*l:2,0:2*l:2]=bc+sign1*var*bc+sign2*var*var*bc
            bc3[1:2*l:2,0:2*l:2]=bc+sign1*var*bc-sign2*var*var*bc
            bc3[0:2*l:2,1:2*l:2]=bc-sign1*var*bc+sign3*var*var*bc
            bc3[1:2*l:2,1:2*l:2]=bc-sign1*var*bc-sign3*var*var*bc

            bc2[np.where(sign4expand==-1)]=bc3[np.where(sign4expand==-1)]

            return bc2        
        for i in np.arange(0,self.orderXY):
            bc=next_step_2D(bc,var)
            var=var*varReduction
        self.lwp=bc
        self.tau=bc*3.0/2.0/re
        self.re=np.repeat(re,len(bc))
        if xydist==None:
            self.x=np.linspace(0,len(bc)/100.0,len(bc)+1)#let's say this is 0.01km resolution
        else:
            self.x=np.linspace(0,xydist,len(bc)+1)
        self.x=(self.x[0:-1]+self.x[1:])/2#assigning the middle
        self.y=self.x
    def save_frach5(self,fname=None):
        '''
        Save fractal cloud properties as an hdf5 file.
        '''
        if not(fname==None):
            self.fname=fname
        
        if os.path.isfile(self.fname+'.hdf5'):
            rd=input('File already exist. Replace?: ')
        else:
            rd='y'
        if rd=='y':
            f=h5py.File(self.fname+'.hdf5','w')
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
            print(self.fname+'.hdf5 SAVED!')
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
        '''
        Plot lwp PDF for both 1D and 2D case.
        '''
        lwp,lwpax=plt.subplots()
        if np.size(self.lwp.shape)>1:
            lwpval=self.lwp.reshape(self.lwp.shape[0]**2)
        else:
            lwpval=self.lwp
        weights = np.ones_like(lwpval)/len(lwpval)
        lwpax.hist(lwpval,bins = 10 ** np.linspace(1, 3, 50),weights=weights,histtype=u'step')
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
def readSHDOMRad(filename,hl=19,nx=4096):
    '''
    filename='cloudExr.out'
    hl = 19 #Header lines
    nx = 4096 #number of columns
    '''
    fl=open(filename,'r')
    sr = (fl.readlines()[11:12][0]).split('=',3)[2:]
    fl.close()
    musza = float(sr[0].split('   ',1)[0])
    SZA = np.rad2deg(np.arccos(musza))
    SAA = float(sr[1])
    muphi1 = np.genfromtxt(filename, skip_header=hl,max_rows=1)[1:3]
    rad1   = np.genfromtxt(filename, skip_header=hl+1,max_rows=nx)
    muphi2 = np.genfromtxt(filename, skip_header=hl+1+nx,max_rows=1)[1:3]
    rad2   = np.genfromtxt(filename, skip_header=hl+1+nx+1,max_rows=nx)
    VAA = np.array([muphi1[1],muphi2[1]])
    VZA = np.rad2deg(np.arccos([muphi1[0],muphi2[0]]))
    Rad = np.array([rad1[:,2],rad2[:,2]])
    return Rad,VZA,VAA,SZA,SAA
def setProp(scat_file, wavelength,fc_obj,nx=4096):
    '''
    scat_file = 'wcld_b0.86_re12.scat'
    fc_obj = frac_cloud.frac_physics
    wavelength = <in um>
    '''
    fc = fc_obj
    #Setting-up SHDOM prop file
    cot1 = fc.tau
    cer1 = fc.re
    cth1 = np.zeros_like(fc.re,dtype=float)+0.5  #convert into km
    
    #create lookup table
    ny=1
    print(nx, ny)
    #nz=20 #1.0/20=50meters
    nz=20 #1.0/20=50meters
    dx=1.0 
    dy=1.0 
    z0=0.0 
    z1=1.0
    temp=0.0
    #extinct1=10
    #albedo=1.0
    iphase=1
    
    wv = wavelength
    print("wv:", wv)
    extv=cot1
    cerv=cer1
    cthv=cth1 #in km
    #        cthv=cthv/cthv.max()*1.0 #scale to 1.0
    print(extv.max(),cerv.max(), extv.shape, cerv.shape)
    prp_name = scat_file.split('.',1)[0]
    f=open(prp_name.replace('.','p')+'.prp','w') 

    nphase=1 #number of phase function
    f.write("Tabulated phase function property file\n")
    f.write("%d %d %d\n" % (nx, ny, nz))
    f.write("%5.3f %5.3f " % (dx, dy))
    for i in range(nz):
        hi=z0+(z1-z0)/(nz-1)*i
        f.write("%5.3f " % hi)
    f.write("\n")
    print("nphase:", nphase)
    f.write("%d\n"% nphase)
    albedov=np.zeros(nphase)
    for n in range(nphase):
#                scat_file=scatpath+'scat_'+str(wv)+'_'+str(reff)
        print(scat_file)
        dat1=np.genfromtxt(scat_file, skip_header=6,max_rows=1)
        nlg=dat1[3]
        f.write("%d"% nlg)
        albedov[n]=dat1[2]
        dat1=np.genfromtxt(scat_file, skip_header=7,max_rows=1)
        dat2=np.genfromtxt(scat_file, skip_header=8,max_rows=1)
        dat3=np.genfromtxt(scat_file, skip_header=9,max_rows=1)
        dat0=np.concatenate([dat1,dat2,dat3])
        if(len(dat0)-1<nlg):
            dat4=np.genfromtxt(scat_file, skip_header=10,max_rows=1)
            #        dat0=np.concatenate([dat1,dat2,dat3,dat4])
#f.write(dat0.shape)
            for datmp in dat1[1:],dat2, dat3, dat4:
                #    datmp=dat1[1:]
                len1=len(datmp)
                for i in range(len1):
                    if(i<len1-1):
                        f.write("  %3.2f"% datmp[i])
                    else:
                        f.write("  %3.2f\n"% datmp[i])
        else:
            for datmp in dat1[1:],dat2, dat3:
                #    datmp=dat1[1:]
                len1=len(datmp)
                for i in range(len1):
                    if(i<len1-1):
                        f.write("  %3.2f"% datmp[i])
                    else:
                        f.write("  %3.2f\n"% datmp[i])

    print(extv.max(),cerv.max(), extv.shape, cerv.shape)                
    for i in range(nx):
        for k in range(nz):
            nk=int(cthv[i]/(z1-z0)*nz) #cloud grid vertical #1km->20, 0.5km->10
            if(k==nk):
                extinct=extv[i]/0.1 #/nk, no need vertical is 1km
                albd=albedov[iphase-1]
            else:
                extinct=0
                albd=1.0
            
            f.write("%d   1  %d  %5.0f  %5.3f %5.3f %d\n" % (i+1, k+1, temp, extinct, albd, iphase))
    f.close()    
    print(prp_name+'.prp SAVED !!')
if __name__=='__main__':
    #Generating Fractal Cloud
    re = 12 # um
    lwp = 90 # g/m^2
    fc = frac_physics('fractal_phys',dpath=None)
#    fc.generate_fractal(re,lwp)
    fc.read_frach5('fractal_phys.hdf5')
#    fc.plot_cld()
#    fc.plot_lwp_pdf()

