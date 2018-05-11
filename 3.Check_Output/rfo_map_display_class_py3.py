import os.path
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as cls
import mpl_toolkits.basemap as bmp

class rfo_map_display(object):
    """
    Display map of Relative-Frequency-of-Occurrence (RFO)

    parameters
    ----------
    ncl: # of clusters
    nrow: # of row for figure
    ncol: # of column for figure
    lons2d, lats2d: meshgrid result of input data lons and lats
    lon_limit, lat_limit: lon and lat boundary list to show on the map

    Attributes
    ----------


    """

    def __init__(self,ncl=1):
        self.ncl=ncl

    def read_bin_data(self,fname,dtp=np.float32):
        """
        Open a binary file, and read data

        fname : file name
        dtp   : data type; np.float32 or np.float64, etc.
        """
        if not os.path.isfile(fname):
            print("File does not exist:"+fname)
            sys.exit()

        with open(fname,'rb') as fd:
            bin_mat = np.fromfile(file=fd,dtype=dtp)

        return bin_mat


    def calc_rfomap(self,crnums):
        """
        From the CR_num data, calculate mean rate of occurreance of each cluster
        """
        nday,nlat,nlon=crnums.shape
        rmap=np.empty([self.ncl,nlat,nlon])

        for k in range(1,self.ncl+1,1):
            idx= crnums==k
            rmap[k-1,:,:]=idx.sum(axis=0)/float(nday)

        return rmap

    def set_parameters(self,nrow,ncol,lons2d,lats2d,lon_limit=[0,360],lat_limit=[-90,90]):
        self.nrow=nrow
        self.ncol=ncol
        self.lons2d=lons2d
        self.lats2d=lats2d
        self.lon_limit=lon_limit
        self.lat_limit=lat_limit

    def map_display_main(self,rfomap,cf,page_size_inches=[6,4],suptit=""):

        self.page_size_inches=page_size_inches

        cm = plt.cm.get_cmap('CMRmap_r',100) #'CMRmap_r' 'YlOrBr' 'Accent' 'afmhot_r'
        cmnew = cm(np.arange(100))
        cmnew = cmnew[1:96,:]
        newcm = cls.LinearSegmentedColormap.from_list("newCMR",cmnew)
        newcm.set_under("white")
        self.cm=newcm

        fig = plt.figure()
        self.fig=fig
        if suptit!="":
            fig.suptitle(suptit,fontsize=18,y=1.)
        axes=self._setup_page()

        for i,ax in enumerate(axes):
            if i<self.ncl:
                tmpmap=rfomap[i,:,:]*100.; grfo=tmpmap.mean()
                pic1=self._map_show(ax,i,tmpmap)
                self._map_show_common(ax,i,cf[i],grfo)

#                if i==int((self.ncl-1)/self.ncol)*self.ncol:
#                    ax.set_xlabel('Optical Thickness',fontsize=13)
#                    ax.set_ylabel('Pressure (hPa)',fontsize=13,labelpad=0)
                if i==self.ncl-1:
                    ### Add Color Bar
                    tt=np.arange(0,96,10)
                    tt2=[str(x)+'%' for x in tt]
                    cb=self._add_colorbar_horizontal(ax,pic1,tt,tt2,ext='max')
            else:
                ax.set_visible(False)



    def _setup_page(self):
        """
        Axes Setup

        1. Set page size
        2. Return axes
        """

        fig=self.fig
        fig.set_size_inches(*self.page_size_inches)

        lf=0.06;rf=0.94
        bf=0.06;tf=0.92
        #fig.subplots_adjust(hspace=0.,wspace=0.,left=lf,right=rf,top=tf,bottom=bf)

        wdt=12;hgt=3.6;gap=1
        x_sum=(wdt*self.ncol+gap*(self.ncol-1))/(rf-lf)
        xw=wdt/float(x_sum); xgap=gap/float(x_sum)
        y_sum=(hgt*self.nrow+gap*(self.nrow-1))/(tf-bf)
        yh=hgt/float(y_sum); ygap=gap/float(y_sum)

        axes=[]
        iy=tf-yh-ygap
        for j in range(self.nrow):
            ix=lf
            for i in range(self.ncol):
                axes.append(fig.add_axes([ix,iy,xw,yh]))  # [Left,Bottom,Width,Height]
                ix=ix+xw+xgap
            iy=iy-yh-ygap

        self.fig=fig
        return axes

    def _add_colorbar_horizontal(self,ax1,pic1,tt,tt2=None,ext='both'):
        ### Add colorbar
        if tt2==None:
            tt2=tt
        pos1=ax1.get_position().bounds  ##<= (left,bottom,width,height)
        cb_ax = self.fig.add_axes([0.1,pos1[1]-0.1,0.8,0.026])
        #cb_ax = fig.add_axes([1.0,0.2,0.03,0.6])  ##<= (left,bottom,width,height)
        cb = self.fig.colorbar(pic1,cax=cb_ax,orientation='horizontal',ticks=tt,extend=ext)
        cb.ax.set_xticklabels(tt2,size=12,stretch='condensed')
        return cb

    def _map_show_common(self,ax1,i,cf,rfo):

        ### add a title.
        subtit='CR{} [CF={:.1f}%]: RFO={:.1f}%]'.format(i+1,cf,rfo)
        print(subtit)
        ax1.set_title(subtit,x=0.0,ha='left',fontsize=12,stretch='semi-condensed')

        ### Ticks
        ax1.tick_params(axis='both',which='major',labelsize=11)

        return

    def _map_show(self,ax1,i,rfomap):

        m = bmp.Basemap(llcrnrlon=self.lon_limit[0],llcrnrlat=self.lat_limit[0],urcrnrlon=self.lon_limit[1],urcrnrlat=self.lat_limit[1],projection='cyl', fix_aspect=False, ax=ax1)

        m.drawcoastlines(color='dimgray',linewidth=0.9)

        props = dict(norm=cls.Normalize(vmin=0,vmax=90),cmap=self.cm,alpha=0.9)
        pic1 = m.pcolormesh(self.lons2d,self.lats2d,rfomap,shading='gouraud',latlon=True,**props)

        # draw parallels and meridians, but don't bother labelling them.
        if i<self.ncl-self.ncol:
            parl_labm=[False,False,False,False]
        else:
            parl_labm=[False,False,False,True]

        if i%self.ncol==0:
            parl_labp=[True,False,False,False]
        elif (i+1)%self.ncol==0 or  i==self.ncl-1:
            parl_labp=[False,True,False,False]
        else:
            parl_labp=[False,False,False,False]

        m.drawparallels(np.arange(-15.,15.1,15.),color='silver',linewidth=0.6,dashes=[3,3],labels=parl_labp)
        m.drawmeridians(np.arange(0.,361.,60.),color='silver',linewidth=0.6,dashes=[3,3],labels=parl_labm)

        return pic1

    def savefig(self,fname):
        ### Show or Save
        #plt.show()
        self.fig.savefig(fname,bbox_inches='tight',dpi=175)
#        self.fig.savefig(fname,dpi=175)

        #if os.path.isfile(outdir+fnout) and not os.path.isfile('./Pics/'+fnout):
        #    call(["ln","-s",outdir+fnout,"./Pics/"])
