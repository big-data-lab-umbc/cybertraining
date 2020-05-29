#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
********************************************
Created on Thu Mar  9 20:19:08 2017
by
Chamara Rajapakshe
(cpn.here@umbc.edu)
******************************************** 
MODIS library developed during Prof. Vanderlie's PHYS 722 class (an intermediate version).

SDSinfo: Generate *.info text file which contains the meta data of the SDS
drawmaps: Draws mapprojections for "cyl" and "ortho" projections.
rgbproj: Create projected RGB matrix
JulianDate_to_MMDDYYY: Convert to Julian date
rmvxtrms: remove extrieme values from a data set.

04/04/2017: Changed interpolation from "cubic" to "linear" to handle nan
    z_01 = griddata(xy1_igrid, z_igrid_01, (xi, yi), method='linear')
    z_02 = griddata(xy1_igrid, z_igrid_02, (xi, yi), method='linear')
    z_03 = griddata(xy1_igrid, z_igrid_03, (xi, yi), method='linear')
    
04/06/2017:
    drawrgb: To plot a rgb image.
    drawContourf: To draw a contourf figer of a reflectance/ radiance array
04/16/2017:
    Going to include level 2 products
    SDSinfo adapted to handdle different MODIS products
    Masking part corrected in dsetMOD021K
"""
from pyhdf.SD import SD, SDC
import numpy as np
import mpl_toolkits.basemap as bm
import matplotlib.pyplot as plt
import scipy.io as sio
from scipy.interpolate import griddata
import calendar
from matplotlib.patches import Polygon
from rebin import rebin
import numpy.ma as ma

class dsetMOD021K(object):
    def __init__(self,res):
        self.res=res# resolution 250/400/1KM
        self.bandnames=[]
        self.RefSB=[]
        self.RefSB_units=[]
        self.RefSB_rad=[]
        self.RefSB_rad_units=[]
        self.radscale=[]
        self.refscale=[]
        self.cscale=[]
        self.coff=[]
        self.radoff=[]
        self.fillval=[]
        self.swath_lon=[]
        self.swath_lat=[]
    def readres(self, hdf):
        #res=EV_250_Aggr1km
        RefSB=hdf.select(self.res)[:][:]
        self.radscale=np.array(hdf.select(self.res).attributes().get('radiance_scales'))
        self.refscale=np.array(hdf.select(self.res).attributes().get('reflectance_scales'))
        self.cscale=np.array(hdf.select(self.res).attributes().get('corrected_counts_scales'))
        self.coff=np.array(hdf.select(self.res).attributes().get('corrected_counts_offsets'))
        self.radoff=np.array(hdf.select(self.res).attributes().get('radiance_offsets'))
        self.fillval=hdf.select(self.res).attributes().get('_FillValue')
        mask=RefSB==self.fillval

        self.RefSB=np.empty_like(RefSB,dtype=float)
        self.RefSB_rad=np.empty_like(RefSB,dtype=float)
        for i in range(0,len(self.refscale)):
            self.RefSB[i,:,:]=(RefSB[i,:,:]-self.coff[i])*self.refscale[i]
            self.RefSB_rad[i,:,:]=(RefSB[i,:,:]-self.radoff[i])*self.radscale[i]
        bandnames=hdf.select(self.res).attributes().get('band_names')
        self.bandnames=bandnames.split(',')
        self.RefSB_units=hdf.select(self.res).attributes().get('corrected_counts_units')
        self.RefSB_rad_units=hdf.select(self.res).attributes().get('radiance_units')   
        
        lats=hdf.select('Latitude')[:][:]
        fval=hdf.select('Latitude').attributes().get('_FillValue')
        lats[lats==fval]=np.nan
        lons=hdf.select('Longitude')[:][:]
        fval=hdf.select('Longitude').attributes().get('_FillValue')
        lons[lons==fval]=np.nan
        cross_dim=np.shape(self.RefSB)[1]
        along_dim=np.shape(self.RefSB)[2]
        # lat lon interpolation
        self.swath_lon=np.zeros([cross_dim,along_dim])
        self.swath_lat=np.zeros([cross_dim,along_dim])
        xp=np.linspace(0,lons.shape[0]-1,lons.shape[0])
        x=np.linspace(0,lons.shape[0],cross_dim)
        
        for j in np.arange(0,along_dim):
            self.swath_lon[:,j]=np.interp(x,xp,lons[:,j/5])
            self.swath_lat[:,j]=np.interp(x,xp,lats[:,j/5])
    
        self.RefSB=ma.masked_array(self.RefSB,mask=mask)
        
    def setCldMask(self,ref470std,threshold=None):
        #spatial variability cloud mask
        #ref470 std of the reflectance in 470 band (676 by 451)
        shp=self.RefSB.shape
        RefSB=np.empty((shp[0],shp[1]/3,shp[2]/3),dtype=float)
        RefSB_rad=np.empty_like(RefSB)
        if threshold==None:
            threshold=0.01
        self.cldMask=ref470std>threshold
        for i in np.arange(0,shp[0]):
            RefSB[i,:,:]=rebin(self.RefSB[i,:,:],factor=(3,3))
            RefSB[i,self.cldMask]=np.nan
            RefSB_rad[i,:,:]=rebin(self.RefSB_rad[i,:,:],factor=(3,3))
            RefSB_rad[i,self.cldMask]=np.nan
        swath_lon=rebin(self.swath_lon,factor=(3,3))
        swath_lat=rebin(self.swath_lat,factor=(3,3))
        self.RefSB=RefSB
        self.RefSB_rad=RefSB_rad
        self.swath_lon=swath_lon
        self.swath_lat=swath_lat
        
            
def SDSinfo(data_path,fname,product=None):
    '''
    #product: MOD021KM,MYD04_L2,.. MODIS product type (just the output file name)
    #for *.hdf files
    '''
    file=data_path+fname
    hdf=SD(file)
    dsets=hdf.datasets()
    keys=dsets.keys()
    if product==None or product=='MOD021K':
        product=fname[0:7]
        f1=open(product+".info","w")
    else:
        f1=open(product+".info","w")
    f1.write(fname+'\n')
    f1.write('Number of variables: '+str(np.size(keys))+'\n')
    f1.write('Key list'+'\n')
    for i in np.arange(np.size(keys)):
        f1.write('['+str(i)+']'+keys[i]+'\n')
        
    for i in np.arange(np.size(keys)):
        d=hdf.select(keys[i])
        f1.write('==========================================================='+'\n')
        f1.write('['+str(i)+'] '+'Variable name:'+keys[i]+'\n')
        f1.write(str(d.dimensions())+'\n')
        f1.write('==========================================================='+'\n')
    
        atrb=d.attributes()
        for j in np.arange(np.size(atrb.keys())):
            f1.write('['+str(i)+'.'+str(j)+']'+(atrb.keys())[j]+'\n')
            f1.write('\t'+str(atrb.get((atrb.keys())[j]))+'\n')

        
        f1.write('\n\n\n')
    f1.close()
    print(product+".info SAVED!!")

def drawmaps(ax,lat_mn,lon_mn,lat_mx,lon_mx,swath_lon,swath_lat,proj=None,gvmap=False):
    #lat_mn, lat_mx : integers
    #proj : "cyl" and "ortho" 
    if proj==None:
        proj='ortho'
    if proj=='cyl':
        mapproj = bm.Basemap(ax=ax,projection='cyl',llcrnrlat= lat_mn, \
            llcrnrlon= lon_mn,urcrnrlat= lat_mx, urcrnrlon= lon_mx)
    elif proj=='ortho':
        lon_0=lon_mn+(lon_mx-lon_mn)/2
        lat_0=lat_mn+(lat_mx-lat_mn)/2
        #lower right and upper right
        m1 = bm.Basemap(projection='ortho',lon_0=lon_0,lat_0=lat_0,resolution=None)
#
        xp0, yp0 = m1(lon_0,lat_0) 
#    
        xp1, yp1 = m1(swath_lon[0,0],swath_lat[0,0]) 
        xp2, yp2 = m1(swath_lon[0,-1],swath_lat[0,-1]) 
        xp3, yp3 = m1(swath_lon[-1,-1],swath_lat[-1,-1])
        xp4, yp4 = m1(swath_lon[-1,0],swath_lat[-1,0])
            
        llx = min(xp1,xp2,xp3,xp4) - xp0  # lower left
        lly = min(yp1,yp2,yp3,yp4) - yp0

        urx = max(xp1,xp2,xp3,xp4) - xp0  # upper right
        ury = max(yp1,yp2,yp3,yp4) - yp0
        mapproj = bm.Basemap(ax=ax,projection='ortho',lon_0=lon_0,\
            lat_0=lat_0,resolution='l',llcrnrx=llx,llcrnry=lly,\
            urcrnrx=urx,urcrnry=ury)
    else:
        print('Only ortho or cyl projections!!')
    lonproj, latproj = mapproj(swath_lon, swath_lat)
    latlines = np.arange(lat_mn,lat_mx,5.0)
    lonlines = np.arange(lon_mn,lon_mx,10.0)
        
    mapproj.drawcoastlines()
    mapproj.drawparallels(latlines, labels=[1,0,0,0])
    mapproj.drawmeridians(lonlines, labels=[0,0,0,1])
    ax.patch.set_color('gray')
    
    if gvmap:
        return lonproj,latproj,mapproj
    else:
        return lonproj,latproj
    
def rgbproj(lon_mn,lon_mx,lat_mn,lat_mx,swath_lat,swath_lon,red,green,blue,ax):
    
    lon_0=lon_mn+(lon_mx-lon_mn)/2
    lat_0=lat_mn+(lat_mx-lat_mn)/2
    #lower right and upper right
    m1 = bm.Basemap(projection='ortho',lon_0=lon_0,lat_0=lat_0,resolution=None)
    xp0, yp0 = m1(lon_0,lat_0)  
    xp1, yp1 = m1(swath_lon[0,0],swath_lat[0,0]) 
    xp2, yp2 = m1(swath_lon[0,-1],swath_lat[0,-1]) 
    xp3, yp3 = m1(swath_lon[-1,-1],swath_lat[-1,-1])
    xp4, yp4 = m1(swath_lon[-1,0],swath_lat[-1,0])
    
    llx = min(xp1,xp2,xp3,xp4) - xp0  # lower left
    lly = min(yp1,yp2,yp3,yp4) - yp0

    urx = max(xp1,xp2,xp3,xp4) - xp0  # upper right
    ury = max(yp1,yp2,yp3,yp4) - yp0

#    m = bm.Basemap(projection='cyl',lon_0=lon_0,lat_0=lat_0,resolution='l',\
#        llcrnrlon=lon_mn,llcrnrlat=lat_mn,urcrnrlon=lon_mx,urcrnrlat=lat_mx,\
#        llcrnrx=llx,llcrnry=lly,urcrnrx=urx,urcrnry=ury)
    m = bm.Basemap(ax=ax,projection='ortho',lon_0=lon_0,lat_0=lat_0,resolution='l',\
        llcrnrx=llx,llcrnry=lly,urcrnrx=urx,urcrnry=ury)

    x_igrid, y_igrid = m(swath_lon,swath_lat)    
    x_igrid = x_igrid - xp0
    y_igrid = y_igrid - yp0

    x1_igrid = x_igrid.ravel()
    y1_igrid = y_igrid.ravel()
    z_igrid_01 = red.ravel()
    z_igrid_02 = green.ravel()
    z_igrid_03 = blue.ravel()
    
    xy1_igrid = np.vstack((x1_igrid, y1_igrid)).T
    xi, yi = np.mgrid[llx:urx:1000j, lly:ury:1000j]

    z_01 = griddata(xy1_igrid, z_igrid_01, (xi, yi), method='linear')
    z_02 = griddata(xy1_igrid, z_igrid_02, (xi, yi), method='linear')
    z_03 = griddata(xy1_igrid, z_igrid_03, (xi, yi), method='linear')
    
    rgb_projected = np.zeros((1000, 1000,3))
    rgb_projected[:,:,0] = z_01[:,:]
    rgb_projected[:,:,1] = z_02[:,:]
    rgb_projected[:,:,2] = z_03[:,:]
    rgb_projected=np.rot90(np.flipud(rgb_projected),k=3)
    
    return m,rgb_projected
    

def JulianDate_to_MMDDYYY(y,jd):
    month = 1
    day = 0
    while jd - calendar.monthrange(y,month)[1] > 0 and month <= 12:
        jd = jd - calendar.monthrange(y,month)[1]
        month = month + 1
    print(month,jd,y)
    
def rmvxtrms(dM):
    if dM.size!=0:
        q3=np.percentile(dM,75);q1=np.percentile(dM,25)
        dMmin=q1-1.5*(q3-q1);dMmax=q3+1.5*(q3-q1)
        dM=dM[dM>dMmin];dM=dM[dM<dMmax]
#    else:
#        print('Zero size dM!!!')
    return dM
    
def draw_screen_poly( lats, lons, m, ax):
    x, y = m( lons, lats )
    xy = zip(x,y)
    poly = Polygon( xy, edgecolor='red', facecolor='#999999',alpha=0.0,linewidth=1.0)
    ax.add_patch(poly)

def draw_box(rgn,m,clr=None):
#    rgn=np.array([-62,-57,-45,-40])#lon_min,lon_max,lat_min,lat_max
#    m mapproj
    if clr==None:
        clr='r'
    bx_lon=[rgn[0],rgn[1],rgn[1],rgn[0]]
    bx_lat=[rgn[2],rgn[2],rgn[3],rgn[3]]
    m.drawgreatcircle(bx_lon[0], bx_lat[0],bx_lon[1], bx_lat[1] ,linewidth=2,color=clr)
    m.drawgreatcircle(bx_lon[1], bx_lat[1],bx_lon[2], bx_lat[2] ,linewidth=2,color=clr)
    m.drawgreatcircle(bx_lon[2], bx_lat[2],bx_lon[3], bx_lat[3] ,linewidth=2,color=clr)
    m.drawgreatcircle(bx_lon[3], bx_lat[3],bx_lon[0], bx_lat[0] ,linewidth=2,color=clr)
    
def drawrgb(res250,res500,ax,bfac=None,rtn=False):
    if bfac==None:
        bfac=1.0
    swath_lat=res250.swath_lat
    swath_lon=res250.swath_lon
    lon_mn,lon_mx=int(res250.swath_lon.min()),int(res250.swath_lon.max())
    lat_mn,lat_mx=int(res250.swath_lat.min()),int(res250.swath_lat.max())
    red = res250.RefSB[0,:,:]/bfac
    green = res500.RefSB[1,:,:]/bfac
    blue = res500.RefSB[0,:,:]/bfac   
    m,rgb_projected=rgbproj(lon_mn,lon_mx,lat_mn,lat_mx,swath_lat,\
                                    swath_lon,red,green,blue,ax)
    
    rgb_projected[rgb_projected>1.0]=1.0
    rgb_projected[rgb_projected<0.0]=0.0
    cmap=plt.cm.jet
    cmap.set_bad('gray',alpha=0)
    m.imshow(rgb_projected, interpolation='nearest', origin='lower',ax=ax,cmap=cmap)

    print(str(np.nanmin(rgb_projected))+','+str(np.nanmax(rgb_projected)))
    
    m.drawcoastlines()
    latlines = np.arange(lat_mn,lat_mx,10.0)
    lonlines = np.arange(lon_mn,lon_mx,10.0)
    m.drawparallels(latlines, labels=[1,0,0,0])
    m.drawmeridians(lonlines, labels=[0,0,0,1])
    if rtn:
        return m

def drawContourf(resAr,resAr_obj,ax,cmap=None,v=None,gvmap=False):
    #resAr 2d array to plot, resAr_obj corresponding dsetMOD021K object
    if cmap==None:
        cmap=plt.cm.jet
    if v==None:
        v=np.linspace(resAr.min(),resAr.max(),100)
    swath_lat=resAr_obj.swath_lat
    swath_lon=resAr_obj.swath_lon
    lon_mn,lon_mx=int(resAr_obj.swath_lon.min()),int(resAr_obj.swath_lon.max())
    lat_mn,lat_mx=int(resAr_obj.swath_lat.min()),int(resAr_obj.swath_lat.max())
    [lonproj,latproj,mapproj]=drawmaps(ax,lat_mn,lon_mn,lat_mx,lon_mx,swath_lon,swath_lat,\
                            gvmap=True)
    ctf = ax.contourf(lonproj,latproj,resAr,v,cmap=cmap,extend='both')
    if gvmap:
        return ctf, mapproj
    else:
        return ctf

