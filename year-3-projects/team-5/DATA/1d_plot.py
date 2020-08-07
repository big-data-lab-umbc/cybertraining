#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 22:30:53 2020

@author: Jackyue
"""
import numpy as np
import h5py
import matplotlib.pyplot as plt

def load_var(fname):
  f=h5py.File(fname,'r')
  
  tau = f['tau'][:]  
  x   = f['x'][:]
  re  = f['red'][:]
  lwp = f['grn'][:]
  lwp = f['blu'][:]
  f.close()
  
  return x, tau, re, lwp


if __name__ == "__main__":
  n_profile = 4000
  out_dir = "./../Out1"

  fname = "./../Out1/profile_%05d.hdf5" %(1)
  hf = h5py.File(fname, 'r')
  x = np.array(hf.get('x'))
  x_size = np.size(x)
  hf.close()

  red = np.empty((n_profile, x_size), dtype=float)
  green = np.empty((n_profile, x_size), dtype=float)
  blue = np.empty((n_profile, x_size), dtype=float)
  tau = np.empty((n_profile, x_size), dtype=float)
  for i in range(0, 4000):
      fname = "./../Out1/profile_%05d.hdf5" %(i+1)
      hf = h5py.File(fname, 'r')
      red[i, :] = np.array(hf.get("red"))
      green[i, :] = np.array(hf.get("grn"))
      blue[i, :] = np.array(hf.get("blu"))
      tau[i, :] = np.array(hf.get("tau"))
      hf.close()

  print(red.shape)
  print(green.shape)
  print(blue.shape)
  print(tau.shape)
  print(x.shape)
  # plot
  fig, axs = plt.subplots(nrows=3,ncols=1,sharex=True)

  ax = axs[0]
  ax.plot(x,red[:100].T,'k',alpha=0.2,color='r')
  ax.set_title(r"Reflectance $(Y)$, sza $= 60\degree$")
  ax.set_ylabel(r"$0.645\quad[\mu m]$")
  
  ax = axs[1]
  ax.plot(x,green[:100].T,'k',alpha=0.2,color='g')
  ax.set_ylabel(r"$0.555\quad[\mu m]$")

  ax = axs[2]
  ax.plot(x,blue[:100].T,'k',alpha=0.2,color='b')
  ax.set_ylabel(r"$0.469\quad[g/m^2]$")
  ax.set_xlabel('X [km]')

  plt.savefig("1dplot.png",dpi=200,bbox_inches='tight')

  
  
