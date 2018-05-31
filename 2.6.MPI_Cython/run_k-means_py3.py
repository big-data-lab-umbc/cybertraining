from k_means_class_py3 import K_means

# indir = '/data1/djin1/Scratch/'
indir = '../'
infile = indir+'aqua_d3_c6_tvp_pcl.noMissing.20050101-20051231_3445612x42.float32.dat'

domain_size = [30,360]
nelem=42
nthreads=4

outdir = './'
outfilehead = outdir+'MODIS_Aqua_b42_TR'

### Define Object
km=K_means(domain_size=domain_size,nelem=nelem,epsilon=0.00001)

### Read input data for clustering
indata=km.read_bin_data(infile)
###** If file size is too large
###** import numpy as np; indata=np.memmap(infile,dtype=np.float32,mode='r')

### Initializing
indata=km.initialize(indata,num_threads=nthreads)
# print(indata.shape,indata.dtype)

### Loop for several K numbers and different initial coditions
# num_try=40
num_try=1
# for kk in range(4,6,1):
for kk in range(10,11,1):
    # for iid in range(1,num_try+1,1):
    for iid in range(1,2):
        km.set_knum_id(knum=kk,id_=iid)
        ictd=km.get_initial_ctd(indata)
        ctd=km.K_means_main(indata,ictd) #[nelem,knum]
        km.print('CF: ',ctd.sum(axis=1))

        km.write_centroid(outfilehead,ctd,ftype='b')
