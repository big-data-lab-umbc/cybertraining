#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
********************************************
Created on Mon Apr  2 10:14:27 2018
by
Chamara Rajapakshe
(cpn.here@umbc.edu)
********************************************
Frequently used common python library
"""
import numpy as np
import matplotlib.pyplot as plt
import os,string
import scipy.signal as signal
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import pickle
import mpl_toolkits.basemap as bm

class pkl_classes(object):
    def __init__(self,):
        self.class_names=[]

def movingaverage (values, window):
    '''
    Return a new array with the moving averaged values of the given window.
    values: The array of the values
    window: moving average window
    '''
    weights = np.repeat(1.0, window)/window
    sma = np.convolve(values, weights, 'valid')
    return sma
def movingaverage2D(values,window):
    '''
    Return a new array with the moving averaged values of the given window.
    values: The 2D array of the values
    window: (intiger) Dimension of a side of the moving average window. (ex. 3 for 3by3 moving average)
    '''    
    weights = np.ones((window,window),dtype=float)/window/window
    sma = signal.convolve2d(values, weights, 'same')
    return sma
    
def find_CDF(data,bins=None):
    '''
    To find cumulative distribution function (CDF)
    data: 1D data array
    bins: bins for the histogram
    '''
    weights=np.ones_like(data)/len(data)
    if bins==None:
        val, base = np.histogram(data, weights=weights)
    else:
        val, base = np.histogram(data, bins=bins,weights=weights)
    return base[1:], np.cumsum(val)

def rmvxtrms(dM):
    '''
    remove extreim value from a data set
    '''
    if dM.size!=0:
        q3=np.percentile(dM,75);q1=np.percentile(dM,25)
        dMmin=q1-1.5*(q3-q1);dMmax=q3+1.5*(q3-q1)
        dM=dM[dM>dMmin];dM=dM[dM<dMmax]
#    else:
#        print('Zero size dM!!!')
    return dM
    
def save_obj(obj, name ):
    '''
    To temporally save object/dictonary
    File names will be OVERWRITTEN!!
    '''
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)
    print(name + '.pkl SAVED!')

def load_obj(name ):
    '''
    To load temporally saved object/dictionary
    '''
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)
'''
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FIGURES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'''
def setup_figures(plt):
    plt.rc('font', family='serif',weight='bold')
#    plt.rc('xtick', labelsize='x-small')
#    plt.rc('ytick', labelsize='x-small')
    plt.rc('xtick', labelsize=12)
    plt.rc('ytick', labelsize=12)
    plt.rc('lines',linewidth=2.0)
    
def sub_labels(ax,clr='k'):
    '''
    Add subfigure labels
    '''
    axs=ax.flat
    for n, ax in enumerate(axs):
        ax.text(0.01, 0.9, '('+string.ascii_lowercase[n]+')', transform=ax.transAxes,
                size=15,color=clr)
        
def add_common_cb(fig,ctf,ts=None,label=None):
    '''
    To add a common color bar for a subplot
    fig: Figure object matplotlib.figure.Figure
    ctf:contourf output instance
    '''
    fig.subplots_adjust(right=0.81)
    cbar_ax = fig.add_axes([0.84, 0.15, 0.01, 0.7])
    if not(ts==None) and not(label==None):
        fig.colorbar(ctf, cax=cbar_ax, ticks=ts,label=label)
    else:
        fig.colorbar(ctf,cax=cbar_ax)
def add_cb(fig,ctf,ax,ticks=None,orientation='horizontal',label='label',pad=0.2):
    divider = make_axes_locatable(ax)
    if orientation=='horizontal':
        cax = divider.append_axes("bottom", size="5%", pad=pad)
    elif orientation=='vertical':
        cax = divider.append_axes("right", size="5%", pad=pad)
    if ticks==None:    
        fig.colorbar(ctf, cax=cax,orientation=orientation,label=label)
    else:
        fig.colorbar(ctf, cax=cax,ticks=ticks,orientation=orientation,label=label)

def savefig(fig,fig_ttl,path=None):
    '''
    fig: Figure object matplotlib.figure.Figure
    fig_ttl: figure title string (some specila characterls will be removed from the file name)
    '''
    rp=False
    for ch in [' ','[',']']:
        if ch in fig_ttl:
            fig_ttl=fig_ttl.replace(ch,'_')
    fig_ttl=fig_ttl.replace('.','p')
    if path==None:
        filename=fig_ttl
    else:
        filename=path+fig_ttl
    if os.path.isfile(filename+'.png'):
        usr=raw_input('Replace existing file?: ')
        if usr=='y':
            rp=True
    else:
        rp=True
    if rp:
        fig.savefig(filename+'.png', format='png', dpi=200)
        print(filename+'.png SAVED.')

def add_subplot_axes(ax,rect,axisbg='w'):
    '''
    example1():
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111)
    rect = [0.2,0.2,0.7,0.7]
    ax1 = add_subplot_axes(ax,rect)
    ax2 = add_subplot_axes(ax1,rect)
    ax3 = add_subplot_axes(ax2,rect)
    plt.show()
    '''
    fig = plt.gcf()
    box = ax.get_position()
    width = box.width
    height = box.height
    inax_position  = ax.transAxes.transform(rect[0:2])
    transFigure = fig.transFigure.inverted()
    infig_position = transFigure.transform(inax_position)    
    x = infig_position[0]
    y = infig_position[1]
    width *= rect[2]
    height *= rect[3]  # <= Typo was here
    subax = fig.add_axes([x,y,width,height],axisbg=axisbg)
    x_labelsize = subax.get_xticklabels()[0].get_size()
    y_labelsize = subax.get_yticklabels()[0].get_size()
    x_labelsize *= rect[2]**0.5
    y_labelsize *= rect[3]**0.5
    subax.xaxis.set_tick_params(labelsize=x_labelsize)
    subax.yaxis.set_tick_params(labelsize=y_labelsize)
    return subax

def Corr_plot_frames(ax,lineS='b-',linewidth=1.0):
    '''
    Draws frames for a correlation plot. 
    '''
    xb=ax.get_xbound()
    yb=ax.get_ybound()
    xmin,xmax=float(xb[0]),float(xb[1])
    ymin,ymax=float(yb[0]),float(yb[1])
    bmx=np.max([xmax,ymax])
    bmn=np.min([xmin,ymin])    
    ax.plot([bmn,bmx],[bmn,bmx],lineS,linewidth=linewidth)
    ax.plot([0.0,0.0],[bmn,bmx],lineS,linewidth=linewidth)
    ax.plot([bmn,bmx],[0.0,0.0],lineS,linewidth=linewidth)
    ax.set_xlim(bmn,bmx)
    ax.set_ylim(bmn,bmx)

def mapProject(ax,lat={'mn':-90,'mx':90},lon={'mn':-180,'mx':180},line_step={'lon':45.0,'lat':30.0}):
    lat_int=line_step['lat']
    lon_int=line_step['lon']

    mapproj = bm.Basemap(ax=ax,projection='cyl',llcrnrlat= lat['mn'], llcrnrlon= lon['mn'],urcrnrlat= lat['mx'], urcrnrlon= lon['mx'])
    
    latlines = np.arange(lat['mn'],lat['mx'],lat_int)
    lonlines = np.arange(lon['mn'],lon['mx'],lon_int)
    
    mapproj.drawcoastlines()
    mapproj.drawparallels(latlines, labels=[1,0,0,0])
    mapproj.drawmeridians(lonlines, labels=[0,0,0,1])
    
    return mapproj
