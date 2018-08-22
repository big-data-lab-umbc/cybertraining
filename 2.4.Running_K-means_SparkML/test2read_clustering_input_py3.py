import numpy as np
import sys
import os.path
import csv

def bin_file_read2mtx(fname,dtp=np.float32):
    """ Open a binary file, and read data
        fname : file name
        dtp   : data type; np.float32 or np.float64, etc. """

    if not os.path.isfile(fname):
        print("File does not exist:"+fname)
        sys.exit()

    fd = open(fname,'rb')
    bin_mat = np.fromfile(file=fd,dtype=dtp)
    fd.close()
    return bin_mat


indir = '/umbc/xfs1/cybertrn/cybertraining2018/team2/research/kmeans/'
infile = indir+'aquad3c6tvppcl.noMissing.20050101200512313445612x42.float32.dat'

nelem=42
chist=bin_file_read2mtx(infile)
### In case the size of input data is too large
#chist = np.memmap(infile,dtype=np.float32,mode='r')

n=chist.shape[0]
m=n/nelem
chist=chist.reshape([m,nelem])

print(chist.shape)
print(chist[000,:])

with open("kMeansData1",'w') as f:
	writer = csv.writer(f,lineterminator='\n')
	writer.writerows(chist)
