from helpers import parApply, fullParApply, addLabels, distCalc
from multiprocessing_on_dill import cpu_count
from multiprocessing_on_dill import Pool as pool
from pandas import concat
from numpy import array_split
import numpy as np
def addMaggiLabels(detM, lowEThresh, upperEThresh):
    # lowEThresh = 0.05
    # upperEThresh = 2.7
    otherLabels = ['TD','TT','DtoT1','DtoT2','DtoT3','CO','WOD','WOT132','WOT213','WOT231','WOT312','WOT321']
    # Setup the singles
    sings = detM[detM['e2'].isnull()]
    sings = sings[(sings['e1']>=lowEThresh) & (sings['e1']<=upperEThresh)]
    singsData = np.zeros((sings.shape[0],12))
    s = np.zeros(12)
    s[5]  = 1
    s[5:] = np.nan
    print("Adding labels to singles...")
    f = lambda row: addLabels(row, labels=otherLabels, values=s)
    sings = fullParApply([sings], f)[0]
    # Setup the doubles
    doubs = detM[detM['e3'].isnull() & detM['e2'].notnull()]
    doubs = doubs[(doubs['e1']>=lowEThresh) & (doubs['e1']<=upperEThresh) & (doubs['e2']>=lowEThresh) & (doubs['e2']<=upperEThresh)]
    d = np.zeros(12)
    d[0], d[5] = (1, 1)
    d[7:] = np.nan
    print("Adding labels to doubles...")
    f = lambda row: addLabels(row, labels=otherLabels, values=d)
    doubs = fullParApply([doubs], f)[0]
    doubs = doubs.sample(frac=1).reset_index(drop=True)
    # Setup the triples
    trips = detM[detM['e3'].notnull() & detM['e2'].notnull()]
    trips = trips[(trips['e1']>=lowEThresh) & (trips['e1']<=upperEThresh) & (trips['e2']>=lowEThresh) & (trips['e2']<=upperEThresh) & (trips['e3']>=lowEThresh) & (trips['e3']<=upperEThresh)]
    t = np.zeros(12)
    t[1], t[5] = (1, 1)
    print("Adding labels to triples...")
    f = lambda row: addLabels(row, labels=otherLabels, values=t)
    trips = fullParApply([trips], f)[0]
    trips = trips.sample(frac=1).reset_index(drop=True)
    return sings, doubs, trips



def parApply(df, sz=None, axis=1):
    # def tempApply(row):
        # balancedData(row,size=tsize)
    return df.apply(lambda row: balanceData(row, size=sz), axis=axis)

def fullParApply(frames, func):
    # Balance out the shuffled data in parallel because it takes a long time
    # This technically works for an arbitrary function and frame
    split_frames = list(map(lambda frame: array_split(frame, cpu_count()), frames))
    # sTrips = array_split( trips, 36)
    full_frames = []
    with pool(cpu_count()) as p:
        for frame in split_frames:
            full_frames.append(concat(p.map(func, frame)))
            # tTriples = concat(p.map(parApply, sTrips))
    return full_frames

def addLabels(row, labels=None, values=None):
    for i in range(len(labels)):
        row[labels[i]] = values[i]
    return row
