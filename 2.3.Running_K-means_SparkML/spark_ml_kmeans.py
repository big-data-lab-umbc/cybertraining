from pyspark.ml.clustering import KMeans
from pyspark.ml.evaluation import ClusteringEvaluator
from pyspark.ml.linalg import Vectors
from pyspark.ml.feature import VectorAssembler
from pyspark.sql import SparkSession
from datetime import datetime
import numpy as np
from numpy  import array
import sys
import os.path
# import csv
from subprocess import call
# import matplotlib.colors as cls
# import matplotlib.pyplot as plt

if __name__ == "__main__":

    startTime = datetime.now()
    print("starting time: ", startTime)

    id = int(sys.argv[1])

    spark = SparkSession \
        .builder \
        .appName("SparkMLKMeans") \
        .getOrCreate()

    # def bin_file_read2mtx(fname, dtp=np.float32):
    #     """ Open a binary file, and read data
    #         fname : file name
    #         dtp   : data type; np.float32 or np.float64, etc. """
    #
    #     if not os.path.isfile(fname):
    #         print("File does not exist:" + fname)
    #         sys.exit()
    #
    #     fd = open(fname, 'rb')
    #     bin_mat = np.fromfile(file=fd, dtype=dtp)
    #     fd.close()
    #     return bin_mat
    #
    # indir = '/umbc/xfs1/cybertrn/cybertraining2018/team2/research/kmeans/'
    # infile = indir + 'aquad3c6tvppcl.noMissing.20050101200512313445612x42.float32.dat'
    #
    # nelem = 42
    # chist = bin_file_read2mtx(infile)
    #
    # n = chist.shape[0]
    # m = n / nelem
    # chist = chist.reshape([m, nelem])
    #
    # print(chist.shape)
    # print(chist[000, :])

    dataFrame = spark.read.csv("/umbc/xfs1/cybertrn/cybertraining2018/team2/research/kmeans/kMeansData1.csv",
                               header=False, inferSchema=True)
    dataFrame.printSchema()
    dataFrame.head()

    assembler = VectorAssembler(
        inputCols=["_c0", "_c1", "_c2", "_c3", "_c4", "_c5", "_c6", "_c7",
                   "_c8", "_c9", "_c10", "_c11", "_c12", "_c13", "_c14",
                   "_c15", "_c16", "_c17", "_c18", "_c19", "_c20", "_c21",
                   "_c22", "_c23", "_c24", "_c25", "_c26", "_c27", "_c28",
                   "_c29", "_c30", "_c31", "_c32", "_c33", "_c34", "_c35",
                   "_c36", "_c37", "_c38", "_c39", "_c40", "_c41"],
        outputCol="features")

    output = assembler.transform(dataFrame)
    print("Assembled columns to vector column 'features'")
    output.select("features").show()

    kmeans = KMeans().setK(10).setSeed(id)
    model = kmeans.fit(output)
    # Make predictions
    predictions = model.transform(output)
    # Evaluate clustering by computing Silhouette score
    evaluator = ClusteringEvaluator()
    silhouette = evaluator.evaluate(predictions)
    print("Silhouette with squared euclidean distance = " + str(silhouette))

    # Shows the result.
    # k = []
    centers = model.clusterCenters()
    print("Cluster Centers: ")
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

    # ctd = centers
    ctd = array(centers)
    # ctd = ctd.reshape([10,42])
    print(ctd)

    def write_centroid(fname,ctd,ftype):
        """
        Sorting the centroid and then write to a file

        ftype='b': binary
        ftype='t': text

        """
        # ctd=ctd.T  #[knum,nelem]
        # ctd=_sort_centroid(ctd)
        print('Sorted_CF: ',ctd.sum(axis=1))

        # fname = "SparkOutput"

        if ftype=='b':
            with open(fname+str(id)+'.float64_dat','wb') as fd:
                ctd.tofile(fd)
        elif ftype=='t':
            np.savetxt(fname+'.txt',ctd,fmt='%.8f',delimiter=' ')

        return

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

    write_centroid("sparkML",ctd,'b')

    print(datetime.now() - startTime)
