from pyspark.ml.linalg import Vectors
from pyspark.ml.stat import Correlation
# from pyspark import SparkContext, SparkConf
from pyspark.sql import SparkSession
import numpy as np
import sys
import os.path
# import csv

if __name__ == "__main__":

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
                               header=False, nferSchema=True)

    dataFrame.printSchema()




