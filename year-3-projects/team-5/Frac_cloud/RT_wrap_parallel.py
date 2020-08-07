from mpi4py import MPI
from frac_cloud import setProp, readSHDOMRad, frac_physics
import numpy as np
import os

def RT_wrap(rank,wvl,cer,CEV,SZA,caseid,pro_dir,out_dir):
  # Inputs
  alpha  = "%0.1f"%(1.0/CEV-3)
  sfname = "profile_%05d_WVL_%01.2f_CER_%02.2f_SZA_%02.1f" %(caseid,wvl,cer,SZA)
  #sfname = "profile_%05d" %caseid+"_"+wvl+'CER'+cer+'SZA'+str(SZA)
  wvl = str(wvl) # Wavelength in micron
  cer = str(cer)  # Cloud Effective Radius in micron
  rank = str(rank)

  SAA = 0 # degree
  VZA = 0 
  
  #Creating Mie Scattering Table
  cr = [] #checking system commands
  sfname = sfname.replace('.','p')
  cm1 = "sed \'s/<WLN>/"+wvl+"/g\' run_mie > Tmp/"+rank+"tmp1"
  cm2 = "sed \'s/<CER>/"+cer+"/g\' Tmp/"+rank+"tmp1 > Tmp/"+rank+"tmp2"
  cm3 = "sed \'s/<NAME>/"+sfname+".scat/g\' Tmp/"+rank+"tmp2 > Tmp/"+rank+"tmp3"
  cm4 = "sed \'s/<ALPHA>/"+alpha+"/g\' Tmp/"+rank+"tmp3 > Tmp/"+rank+"tmp4"
  cr =[os.system(st) for st in [cm1,cm2,cm3,cm4]]
  cr= np.append(cr,os.system('chmod +x Tmp/'+rank+'tmp4'))
  cr= np.append(cr,os.system('./Tmp/'+rank+'tmp4'))
  cr= np.append(cr,os.system('mv Tmp/'+rank+'tmp4 run_mie_%s' %sfname))
  cr= np.append(cr,os.system('rm Tmp/'+rank+'tmp*'))
  scat_file = sfname.replace('.','p')+'.scat'
  print('Mie Scattering table ...')
  print(cr)     
  
  fc = frac_physics('1dfract',dpath=None)
  fc.read_frach5('%s/profile_%05d.hdf5' %(pro_dir,caseid))

  if sum(cr)==0:
      #Generating *.prp file for SHDOM
      setProp(scat_file, wvl,fc,dx=0.2)
      #setProp(scat_file, wvl,fc,dx=0.02)
      #Running SHDOM
      mu_SZA = np.cos(np.deg2rad(SZA))
      mu_VZA = np.cos(np.deg2rad(VZA))
      cm1 = "sed \'s/<PRPN>/"+sfname+".prp/g\' run_shdom > Tmp/"+rank+"tmp1"
      cm2 = "sed \'s/<BASN>/"+sfname+"/g\' Tmp/"+rank+"tmp1 > Tmp/"+rank+"tmp2"
      cm3 = "sed \'s/<WVL>/"+wvl+"/g\' Tmp/"+rank+"tmp2 > Tmp/"+rank+"tmp3"
      cm4 = "sed \'s/<MSZA>/%0.2f"%mu_SZA+"/g\' Tmp/"+rank+"tmp3 > Tmp/"+rank+"tmp4"
      cm5 = "sed \'s/<DSAA>/%0.1f"%(SAA)+"/g\' Tmp/"+rank+"tmp4 > Tmp/"+rank+"tmp5"
      cm6 = "sed \'s/<MVZA>/%0.1f"%(mu_VZA)+"/g\' Tmp/"+rank+"tmp5 > Tmp/"+rank+"tmp6"
      cr =[os.system(st) for st in [cm1,cm2,cm3,cm4,cm5,cm6]]
      cr =np.append(cr,os.system('chmod +x Tmp/'+rank+'tmp6'))
      cr =np.append(cr,os.system('./Tmp/'+rank+'tmp6'))
      cr =np.append(cr,os.system('mv Tmp/'+rank+'tmp6 run_shdom_%s' %sfname))
      cr =np.append(cr,os.system('rm Tmp/'+rank+'tmp*'))
      cr =np.append(cr,os.system('mv %sr.out %s/%sr.out' %(sfname,out_dir,sfname)))
      cr =np.append(cr,os.system('rm run_shdom_%s' %sfname ))
      cr =np.append(cr,os.system('rm run_mie_%s' %sfname ))
      cr =np.append(cr,os.system('rm %sm.out' %sfname ))
      cr =np.append(cr,os.system('rm %s.prp' %sfname ))
      cr =np.append(cr,os.system('rm %s.scat' %sfname ))
      print('SHDOM run ...')
      print(cr)

  else:
      print('Something went wrong in Mie computation')

  if sum(cr)>0:
      print('Error occured!!')


if __name__ == "__main__":
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    in_dir = "/home/kirana/cybertraining2020_team5/research/RT/profile"

    nprof = 5 # number of profiles to run
    l_nprof = int(nprof / size)
    remainder = nprof % size
    if (rank < remainder):
        start = (l_nprof + 1) * rank + 1
        end = l_nprof +1 + start
    else:
        start = (l_nprof + 1) * remainder + l_nprof * (rank - remainder) + 1
        end = l_nprof + start

    out_dir = "./LargeOut"
    wvl = [0.645, 0.555, 0.469]
    cev = 0.05
    sza = 60.

    for i in range(1, nprof + 1):
      fn = in_dir + "/profile_%05d.hdf5" %(i)
      fc = frac_physics('1dfract',dpath=None)
      fc.read_frach5(fn)
      re = fc.re[0]
      RT_wrap(rank,wvl[0],re,cev,sza,i,in_dir,out_dir)
      RT_wrap(rank,wvl[1],re,cev,sza,i,in_dir,out_dir)
      RT_wrap(rank,wvl[2],re,cev,sza,i,in_dir,out_dir)
