# Author: Meng Gao
# DATE: 2018
import numpy as np
import h5py

filename='./MOD06_L2.A2000306.2325._cot_0900_0250COTCHT.nc'
#read in modis data cloud field
f = h5py.File(filename, 'r')
cot1 = f['tau'][()].T
cer1 = f['CER'][()].T
cth1 = f['CHT'][()].T/1000 #convert into km
print(cot1.shape, cer1.shape, cth1.shape)

#
#create lookup table
nx=50
ny=50
print(nx, ny)
#nz=20 #1.0/20=50meters
nz=50 #1.0/20=50meters
dx=1.0 
dy=1.0 
z0=0.0 
z1=1.0
temp=0.0
#extinct1=10
#albedo=1.0
iphase=1

#wv=0.86 #2.1
#wv=2.1
#reff=10
#scatpath='/home/mgao/team5_mgao/3DRT/test/shdom/cloud_lookup/scat_table/'
scatpath='/home/mgao/team5_mgao/3DRT/test/shdom/cloud_2d_modis3/scat_table/'

for wv in [0.86,2.1]:
    print("wv:", wv)
    for index in [22]:
        extv=cot1
        cerv=cer1
        cthv=cth1 #in km
        #        cthv=cthv/cthv.max()*1.0 #scale to 1.0
        print(extv.max(),cerv.max(), extv.shape, cerv.shape)
        f=open('atmos'+str(wv)+'_index'+str(index)+'.prp','w') 

        nphase=1 #number of phase function
        f.write("Tabulated phase function property file\n")
        f.write("%d %d %d\n" % (nx, ny, nz))
        f.write("%5.3f %5.3f " % (dx, dy))
        for i in range(nz):
            hi=z0+(z1-z0)/(nz-1)*i
            f.write("%5.3f " % hi)
        f.write("\n")
        reffv=np.array(range(8,18))
        nphase=len(reffv) #number of phase function
        print("nphase:", nphase)
        f.write("%d\n"% nphase)
        albedov=np.zeros(nphase)
        for n in range(nphase):
            reff=reffv[n]
            scat_file=scatpath+'scat_'+str(wv)+'_'+str(reff)
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
            for j in range(ny):
                for k in range(nz):
                    nk=int(cthv[i,j]/(z1-z0)*nz) #cloud grid vertical #1km->20, 0.5km->10
                    if(k<=nk):
                        extinct=nz*extv[i, j]/nk/(z1-z0) #/nk, no need vertical is 1km
                        cert=int(cerv[i,j])
                    else:
                        extinct=0
                        cert=0
                    
                    if(cert<=8):
                        iphase=1
                    elif(cert>=17):
                        iphase=10
                    else:
                        iphase=cert-8+1
                    f.write("%d   %d  %d  %5.0f  %5.3f %5.3f %d\n" % (i+1, j+1, k+1, temp, extinct, albedov[iphase-1], iphase))
        f.close()
