from pyspark.ml.linalg import Vectors
from pyspark.ml.stat import Correlation
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
    # ### In case the size of input data is too large
    # # chist = np.memmap(infile,dtype=np.float32,mode='r')
    #
    # n = chist.shape[0]
    # m = n / nelem
    # chist = chist.reshape([m, nelem])
    #
    # print(chist.shape)
    # print(chist[000, :])

    sparkDf = spark.read.option("header","false").csv("/umbc/xfs1/cybertrn/cybertraining2018/team2/research/kmeans/kMeansData1.csv")

    sparkDf.first()

    data = [(Vectors.sparse(4, [(0, 1.0), (3, -2.0)]),),
            (Vectors.dense([4.0, 5.0, 0.0, 3.0]),),
            (Vectors.dense([6.0, 7.0, 0.0, 8.0]),),
            (Vectors.sparse(4, [(0, 9.0), (3, 1.0)]),)]

    df = spark.createDataFrame(data, ["features"])

    r1 = Correlation.corr(df, "features").head()
    print("Pearson correlation matrix:\n" + str(r1[0]))

    r2 = Correlation.corr(df, "features", "spearman").head()
    print("Spearman correlation matrix:\n" + str(r2[0]))

    # rdd1 = spark.parallelize(chist)
    # rdd2 = rdd1.map(lambda x: [int(i) for i in x])
    #
    # df2 = rdd2.toDF(["X1", "X2", "X3", "X4", "X5", "X6", "X7",
    #                  "X8", "X9", "X10", "X11", "X12", "X13", "X14",
    #                  "X15", "X16", "X17", "X18", "X19", "X20", "X21",
    #                  "X22", "X23", "X24", "X25", "X26", "X27", "X28",
    #                  "X29", "X30", "X31", "X32", "X33", "X34", "X35",
    #                  "X36", "X37", "X38", "X39", "X40", "X41", "X42"])
    # df2.first()
    # df.printSchema()
