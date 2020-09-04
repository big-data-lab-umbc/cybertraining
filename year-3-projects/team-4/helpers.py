from multiprocessing_on_dill import cpu_count
from multiprocessing_on_dill import Pool as pool
import numpy as np
def addLabels(row, labels=None, values=None):
    for i in range(len(labels)):
        row[labels[i]] = values[i]
    return row

def parApply(df, sz=None, axis=1):
    # def tempApply(row):
        # balancedData(row,size=tsize)
    return df.apply(lambda row: balanceData(row, size=sz), axis=axis)

def fullParApply(frames, func):
    # Balance out the shuffled data in parallel because it takes a long time
    # This technically works for an arbitrary function and frame
    split_frames = list(map(lambda frame: array_split(frame, 36), frames))
    # sTrips = array_split( trips, 36)
    full_frames = []
    with pool(cpu_count()) as p:
        for frame in split_frames:
            full_frames.append(concat(p.map(func, frame)))
            # tTriples = concat(p.map(parApply, sTrips))
    return full_frames

def distCalc(x1,y1,z1,x2,y2,z2):
    '''
    distVal = distCalc(x1,y1,z1,x2,y2,z2)
    
    returns the euclidian distance between two points
    '''
    distVal = np.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
    return distVal

def distCalcVec(x, y):
    return distCalc(x[0], x[1], x[2], y[0], y[1], y[2])
