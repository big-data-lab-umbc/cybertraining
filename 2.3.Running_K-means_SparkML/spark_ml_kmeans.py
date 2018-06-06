import os.path
import numpy as np
import pandas as pd
import sys
from datetime import datetime
from numpy import array
from pyspark.ml.clustering import KMeans
from pyspark.ml.evaluation import ClusteringEvaluator
from pyspark.sql import SparkSession
from pyspark.ml.linalg import Vectors
from pyspark.ml.feature import VectorAssembler

if __name__ == "__main__":
    startTime = datetime.now()
    print("starting time: ", startTime)

    k = int(sys.argv[1])
    sid = int(sys.argv[2])  # tials

    spark = SparkSession \
        .builder \
        .appName("SparkMLKMeans") \
        .config("spark.sql.execution.arrow.enabled", "true") \
        .getOrCreate()

def bin_file_read2mtx(fname, dtp=np.float32):
    """ Open a binary file, and read data
        fname : file name
        dtp   : data type; np.float32 or np.float64, etc. """

    if not os.path.isfile(fname):
        print("File does not exist:" + fname)
        sys.exit()

    fd = open(fname, 'rb')
    bin_mat = np.fromfile(file=fd, dtype=dtp)
    fd.close()
    return bin_mat

indir = '/umbc/xfs1/cybertrn/cybertraining2018/team2/research/kmeans/'
infile = indir + 'aquad3c6tvppcl.noMissing.20050101200512313445612x42.float32.dat'

nelem = 42
chist = bin_file_read2mtx(infile)
print(chist.shape)

n = chist.shape[0]
m = n / nelem
chist = chist.reshape([m, nelem])

print(chist.shape)
print(chist[000, :])

df = pd.DataFrame(data=chist)

sparkdf = spark.createDataFrame(df)

# dataFrame = spark.read.csv("/umbc/xfs1/cybertrn/cybertraining2018/team2/research/kmeans/kMeansData2008.csv",
#                            header=False, inferSchema=True)
# sparkdf.printSchema()
# dataFrame.head()

assembler = VectorAssembler(
    inputCols=["0", "1", "2", "3", "4", "5", "6", "7",
               "8", "9", "10", "11", "12", "13", "14",
               "15", "16", "17", "18", "19", "20", "21",
               "22", "23", "24", "25", "26", "27", "28",
               "29", "30", "31", "32", "33", "34", "35",
               "36", "37", "38", "39", "40", "41"],
    outputCol="features")

output = assembler.transform(sparkdf)
# print("Assembled columns to vector column 'features'")
output.select("features").show()

# kmeans = KMeans().setK(10).setSeed(sid)
kmeans = KMeans(k=k, maxIter=sid)
model = kmeans.fit(output)
# Make predictions
predictions = model.transform(output)
# Evaluate clustering by computing Silhouette score
evaluator = ClusteringEvaluator()
silhouette = evaluator.evaluate(predictions)
# print("Silhouette with squared euclidean distance = " + str(silhouette))

# Shows the result.
# k = []
centers = model.clusterCenters()
print("Silhouette with squared euclidean distance = " + str(silhouette))
# for center in centers:
#     # print(center)
#     # print(np.shape(center))
#     # print(center[0])
#     # print(center[41])
#     # for index in range(42):
#     #     print(center[index])
#     matrix = []
#     chunks = 6
#     for loopi in range(len(center)/chunks):
#         matrix.append(center[loopi*chunks:(loopi+1)*chunks])
#     print matrix
#     plt.matshow(matrix)
#     plt.colorbar()
#     plt.show()

ctd = centers
print("Cluster Centers: ")

ctd = array(centers)
# ctd = ctd.reshape([10,42])
print(ctd)

# def write_centroid(fname,ctd,ftype):
#     """
#     Sorting the centroid and then write to a file
#
#     ftype='b': binary
#     ftype='t': text
#
#     """
#     # ctd=ctd.T  #[knum,nelem]
#     # ctd=_sort_centroid(ctd)
#     print('Sorted_CF: ',ctd.sum(axis=1))
#
#     # fname = "SparkOutput"
#
#     if ftype=='b':
#         with open(fname+str(id)+'.float64_dat','wb') as fd:
#             ctd.tofile(fd)
#     elif ftype=='t':
#         np.savetxt(fname+'.txt',ctd,fmt='%.8f',delimiter=' ')
#
#     return

# def _sort_centroid(ctd):
#     """
#     Sort the centroid
#
#     Thick and high first, thin high second, and thin low last.
#     The lowest CF one (less than 50%) is at the end.
#     Input: centriod, dimension=[knum,nelem]
#     Output: sorted centroid
#     """
#
#     cf=ctd.sum(axis=1)
#     idx= cf<0.5
#     ctd2=ctd[~idx,:].reshape([-1,7,3,2]).sum(axis=3)
#     ctd2[:,0,:]=ctd2[:,0:3,:].sum(axis=1)
#     ctd2[:,1,:]=ctd2[:,3:5,:].sum(axis=1)
#     ctd2[:,2,:]=ctd2[:,5:7,:].sum(axis=1)
#     ctd2=ctd2[:,0:3,:].reshape([-1,9])
#
#     wt=np.arange(1,10,1).reshape([3,3])[::-1,:].reshape(-1)
#     wcf=np.average(ctd2,weights=wt,axis=1)
#     ctd0=ctd[~idx,:][np.argsort(wcf)[::-1],:]
#
#     if idx.sum()>0:
#         xx=np.argsort(cf[idx])[::-1]
#         ctd2=ctd[idx,:].reshape([-1,self.nelem])[xx,:]
#         ctd0=np.concatenate((ctd0,ctd2))
#     return ctd0

# write_centroid("sparkML",ctd,'b')

print(datetime.now() - startTime)
