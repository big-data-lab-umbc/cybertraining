# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 11:07:45 2019

@author: paul.maggi
"""
# For parallel apply
from multiprocessing_on_dill import Pool as pool
from numpy import array_split
from pandas import concat


import numpy as np
# import keras
# from keras.models import Sequential
# from keras.layers import Dense, Activation
# from keras import regularizers
from matplotlib.pyplot import *
import pandas as pd
import numpy.random as random
import scipy.sparse as sprs


#%%
def distCalc(x1,y1,z1,x2,y2,z2):
    '''
    distVal = distCalc(x1,y1,z1,x2,y2,z2)
    
    returns the euclidian distance between two points
    '''
    distVal = np.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
    return distVal

def distCalcVec(x, y):
    return distCalc(x[0], x[1], x[2], y[0], y[1], y[2])

def detPix(csvData,resXZ=0.86,resY=0.7):
    '''
    pixelized = detPix(csvData,resXZ,resY)
    
    a crude pixelation just based on data flooring. Accepts singles, doubles or triples in normal 4, 8, or 12 column format
    
    input:
        csvData : interaction data in normal (edep, x, y, z) format
        resXZ : pixelation resolution in X and Z; assumd symetric
        resY : y data 'resolution'. In reality this is incorrect; use this method to just approximate uncertainy
    output:
        pixelized : pixelated data.
    '''
    pixelized = csvData.copy()
    temp = np.floor(csvData[:,1]/resXZ)
    pixelized[:,1] = (temp*resXZ)
    temp = np.floor((csvData[:,2]/resY))
    pixelized[:,2] = (temp*resY)
    temp = np.floor((csvData[:,3]/resXZ))
    pixelized[:,3] = (temp*resXZ)
    if csvData.shape[1]>4:
        temp = np.floor((csvData[:,5]/resXZ))
        pixelized[:,5] = (temp*resXZ)
        temp = np.floor((csvData[:,6]/resY))
        pixelized[:,6] = (temp*resY)
        temp = np.floor((csvData[:,7]/resXZ))
        pixelized[:,7] = (temp*resXZ)
    if csvData.shape[1]>8:
        temp = np.floor((csvData[:,9]/resXZ))
        pixelized[:,9] = (temp*resXZ)
        temp = np.floor((csvData[:,10]/resY))
        pixelized[:,10] = (temp*resY)
        temp = np.floor((csvData[:,11]/resXZ))
        pixelized[:,11] = (temp*resXZ)

    return pixelized


def modListGen(loc):
    '''
    modList = modListGen(loc):
    
    looks at the 3 valued array from the input and translates it to a module number

    input:
        loc : a three valued array specificying x, y and z locations of a module
    outout:
        modList : module number from 0 to 15
    '''
    modList = np.zeros((loc.shape[0],))
    modList[np.sum((loc[:] == ((1,-1,2))),axis=1)==3] = 0
    modList[np.sum((loc[:] == ((1,-1,1))),axis=1)==3] = 1
    modList[np.sum((loc[:] == ((1,-1,-1))),axis=1)==3] = 2
    modList[np.sum((loc[:] == ((1,-1,-2))),axis=1)==3] = 3
    modList[np.sum((loc[:] == ((-1,-1,-2))),axis=1)==3] = 4
    modList[np.sum((loc[:] == ((-1,-1,-1))),axis=1)==3] = 5
    modList[np.sum((loc[:] == ((-1,-1,1))),axis=1)==3] = 6
    modList[np.sum((loc[:] == ((-1,-1,2))),axis=1)==3] = 7
    modList[np.sum((loc[:] == ((1,1,-2))),axis=1)==3] = 8
    modList[np.sum((loc[:] == ((1,1,-1))),axis=1)==3] = 9
    modList[np.sum((loc[:] == ((1,1,1))),axis=1)==3] = 10
    modList[np.sum((loc[:] == ((1,1,2))),axis=1)==3] = 11
    modList[np.sum((loc[:] == ((-1,1,2))),axis=1)==3] = 12
    modList[np.sum((loc[:] == ((-1,1,1))),axis=1)==3] = 13
    modList[np.sum((loc[:] == ((-1,1,-1))),axis=1)==3] = 14
    modList[np.sum((loc[:] == ((-1,1,-2))),axis=1)==3] = 15

    return modList

def getModSings(singlesRaw):
    '''
    mdls = getModSings(singlesRaw)
    
    this function takes in the data from a single set of interactions (edep,x,y,z) and returns the module number of each interaction.
    Assumes 16 modules. It can be passed any interaction quartet, e.g. (edep2,x2,y2,z2) from a triple.

    input: 
        singlesRaw : list of (edep,x,y,z) values
    output: 
        mdls : module numbers [0,15]
    '''
    singles = detPix(singlesRaw,0.1,0.1)
    yMean = np.unique(singles[:,2]).mean() #two options
    xMean = np.unique(singles[:,1]).mean() #two options over 4 crystal locs
    zMean = np.unique(singles[:,3]).mean() #four options over 8 crystal locs 
    zLeft = (singles[singles[:,3]<zMean,3].min()+singles[singles[:,3]<zMean,3].max())/2
    zRight = (singles[singles[:,3]>zMean,3].max()+singles[singles[:,3]>zMean,3].min())/2
    locs = np.zeros((singles.shape[0],3))
    locs[singles[:,2]<yMean,1] = -1
    locs[singles[:,2]>yMean,1] = 1
    locs[singles[:,1]<xMean,0] = -1
    locs[singles[:,1]>xMean,0] = 1
    locs[singles[:,3]<zLeft,2] = -2
    locs[np.logical_and(singles[:,3]>zLeft,singles[:,3]<zMean),2] = -1
    locs[np.logical_and(singles[:,3]>zMean,singles[:,3]<zRight),2] = 1
    locs[singles[:,3]>zRight,2] = 2
    
    mdls = modListGen(locs)
    
    return np.array((mdls)).astype(int)


def getModl2(doublesRaw):
    '''
    mdls = getModl2(doublesRaw)
    
    this function takes in the data from a single set of interactions (edep,x,y,z) and returns the module number of each interaction.
    Assumes 16 modules. It can be passed any interaction quartet, e.g. (edep2,x2,y2,z2) from a triple.

    input: 
        singlesRaw : list of (edep,x,y,z) values
    output: 
        mdls : module numbers [0,15]
    '''
    doubles = detPix(doublesRaw,0.1,0.1)
    yMean = np.unique(doubles[:,(2,6)]).mean() #two options
    xMean = np.unique(doubles[:,(1,5)]).mean() #two options over 4 crystal locs
    zMean = np.unique(doubles[:,(3,7)]).mean() #four options over 8 crystal locs
    zLeft = (doubles[doubles[:,3]<zMean,3].min()+doubles[doubles[:,3]<zMean,3].max())/2
    zRight = (doubles[doubles[:,3]>zMean,3].max()+doubles[doubles[:,3]>zMean,3].min())/2
    locs1 = np.zeros((doubles.shape[0],3))
    locs1[doubles[:,2]<yMean,1] = -1
    locs1[doubles[:,2]>yMean,1] = 1
    locs1[doubles[:,1]<xMean,0] = -1
    locs1[doubles[:,1]>xMean,0] = 1
    locs1[doubles[:,3]<zLeft,2] = -2
    locs1[np.logical_and(doubles[:,3]>zLeft,doubles[:,3]<zMean),2] = -1
    locs1[np.logical_and(doubles[:,3]>zMean,doubles[:,3]<zRight),2] = 1
    locs1[doubles[:,3]>zRight,2] = 2

    locs2 = np.zeros((doubles.shape[0],3))
    locs2[doubles[:,6]<yMean,1] = -1
    locs2[doubles[:,6]>yMean,1] = 1
    locs2[doubles[:,5]<xMean,0] = -1
    locs2[doubles[:,5]>xMean,0] = 1
    locs2[doubles[:,7]<zLeft,2] = -2
    locs2[np.logical_and(doubles[:,7]>zLeft,doubles[:,7]<zMean),2] = -1
    locs2[np.logical_and(doubles[:,7]>zMean,doubles[:,7]<zRight),2] = 1
    locs2[doubles[:,7]>zRight,2] = 2
    
    mdls1 = modListGen(locs1)
    mdls2 = modListGen(locs2)
    
    return mdls1.astype(int),mdls2.astype(int)



def genFalseDoub(singsDF,numD=100):
    '''
    dfOut = genFalseDoub(singsDF,numD)
    
    randomly combines 2 singles from the same module into false doubles. Ensures that they aren't too close (within 1/2 pixel) from each other
    
    inputs:
            singsDF : pandas dataframe of the singles, with the appropriate labels.
            numD : number of false doubles to generate
    outputs:
        dfOut : pandas dataframe of the false doubles, with appropriate labels
    '''
    singles = np.array(singsDF.loc[:,'e1':'z1'])
    sMax = singles.shape[0]
    sMods = getModSings(singles)
    modTots = np.zeros((16,))
    for iii in range(16):
        modTots[iii] = (sMods==iii).sum()

    numSingPerMod = np.ceil(numD * modTots/(modTots.sum()))
    doubOut = np.zeros((numSingPerMod.sum().astype(int),8))
    doubCount = 0
    for jjj in range(16):
        sList = singles[sMods==jjj,:]
        mT = modTots[jjj]
        for iii in range(numSingPerMod[jjj].astype(int)):
            ind1 = np.random.randint(mT)
            ind2 = np.random.randint(mT)
            while ind1 == ind2:
                ind2 = np.random.randint(mT)
            d1 = distCalc(sList[ind1,1],sList[ind1,2],sList[ind1,3],sList[ind2,1],sList[ind2,2],sList[ind2,3])
            while d1 <= 1.72:
                ind1 = np.random.randint(mT)
                ind2 = np.random.randint(mT)
                while ind1 == ind2:
                    ind2 = np.random.randint(mT)
                d1 = distCalc(sList[ind1,1],sList[ind1,2],sList[ind1,3],sList[ind2,1],sList[ind2,2],sList[ind2,3])            
            doubOut[doubCount,:4] = sList[ind1,:].copy()
            doubOut[doubCount,4:] = sList[ind2,:].copy()
            doubCount += 1
            
    random.shuffle(doubOut)
    random.shuffle(doubOut)
    doubOut = doubOut[:numD,:].copy()
    dfOut = pd.DataFrame(np.zeros((doubOut.shape[0],len(singsDF.columns))),columns=singsDF.columns)
    dfOut.loc[:,'e1':'z2']=doubOut.copy()
    dfOut['PS']=np.nan
    dfOut['FP']=np.nan
    dfOut.loc[:,'CO':'e0']=np.nan
       
    return dfOut


def genFalseTrip(singsDF,numD=100):
    '''
    dfOut = genFalseTrip(singsDF,numD)
    
    randomly combines 3 singles from the same module into false doubles. Ensures that they aren't too close (within 1/2 pixel) from each other
    
    inputs:
            singsDF : pandas dataframe of the singles, with the appropriate labels.
            numD : number of false doubles to generate
    outputs:
        dfOut : pandas dataframe of the false triples (didn't change variable names), with appropriate labels
    '''
    singles = np.array(singsDF.loc[:,'e1':'z1'])
    sMax = singles.shape[0]
    sMods = getModSings(singles)
    modTots = np.zeros((16,))
    for iii in range(16):
        modTots[iii] = (sMods==iii).sum()

    numSingPerMod = np.ceil(numD * modTots/(modTots.sum()))
    doubOut = np.zeros((numSingPerMod.sum().astype(int),12))
    doubCount = 0
    for jjj in range(16):
        sList = singles[sMods==jjj,:]
        mT = modTots[jjj]
        for iii in range(numSingPerMod[jjj].astype(int)):
            ind1 = np.random.randint(mT)
            ind2 = np.random.randint(mT)
            ind3 = np.random.randint(mT)
            while (ind1 == ind2) | (ind1 == ind3) | (ind2 == ind3):
                ind1 = np.random.randint(mT)
                ind2 = np.random.randint(mT)
                ind3 = np.random.randint(mT)
            d1 = distCalc(sList[ind1,1],sList[ind1,2],sList[ind1,3],sList[ind2,1],sList[ind2,2],sList[ind2,3])
            d2 = distCalc(sList[ind2,1],sList[ind2,2],sList[ind2,3],sList[ind3,1],sList[ind3,2],sList[ind3,3])
            d3 = distCalc(sList[ind1,1],sList[ind1,2],sList[ind1,3],sList[ind3,1],sList[ind3,2],sList[ind3,3])
            while (d1 <= 1.72) | (d2 <= 1.72) | (d3 <= 1.72):
                ind1 = np.random.randint(mT)
                ind2 = np.random.randint(mT)
                ind3 = np.random.randint(mT)
                while (ind1 == ind2) | (ind1 == ind3) | (ind2 == ind3):
                    ind1 = np.random.randint(mT)
                    ind2 = np.random.randint(mT)
                    ind3 = np.random.randint(mT)
                d1 = distCalc(sList[ind1,1],sList[ind1,2],sList[ind1,3],sList[ind2,1],sList[ind2,2],sList[ind2,3])
                d2 = distCalc(sList[ind2,1],sList[ind2,2],sList[ind2,3],sList[ind3,1],sList[ind3,2],sList[ind3,3])
                d3 = distCalc(sList[ind1,1],sList[ind1,2],sList[ind1,3],sList[ind3,1],sList[ind3,2],sList[ind3,3])         
            doubOut[doubCount,:4] = sList[ind1,:].copy()
            doubOut[doubCount,4:8] = sList[ind2,:].copy()
            doubOut[doubCount,8:] = sList[ind3,:].copy()
            doubCount += 1
            
    random.shuffle(doubOut)
    random.shuffle(doubOut)
    doubOut = doubOut[:numD,:].copy()
    dfOut = pd.DataFrame(np.zeros((doubOut.shape[0],len(singsDF.columns))),columns=singsDF.columns)
    dfOut.loc[:,'e1':'z3']=doubOut.copy()
    dfOut['PS']=np.nan
    dfOut['FP']=np.nan
    dfOut.loc[:,'CO':'e0']=np.nan
       
    return dfOut
# Currently in use inside genDtoT
def tooClose(record1, record2):
    # Computes the distance between the single and double 
    v1 = [record1[1], record1[2], record1[3]]
    v2 = [record2['x1'], record2['y1'], record2['z1']]
    v3 = [record2['x2'], record2['y2'], record2['z2']]
    d1 = distCalcVec(v1, v2)
    d2 = distCalcVec(v2, v3)
    d3 = distCalcVec(v1, v3)
    # Return True if too close together in any dimension
    return (d1 <= 1.72) or (d2 <= 1.72) or (d3 <= 1.72)

def findSingleForDouble_par(row, sList=None, dList=None, mTS=None, mTD=None):
    # Being used in parallel. Create a new numpy seed on each process...
    np.random.seed()
    return findSingleForDouble(row, sList=sList, dList=dList, mTS=mTS, mTD=mTD)

def findSingleForDouble(row, sList=None, dList=None, mTS=None, mTD=None):
    ind1, ind2 = 0, 0
    while ind1 == ind2 or tooClose(sList[ind1], dList.iloc[ind2]):
        # Picks two random indices
        ind1 = np.random.randint(mTS)
        ind2 = np.random.randint(mTD)
    row.loc['PS':] = dList.loc[ind2,'PS':]
    row['TD']  = 0
    row['WOT132':'WOT321'] = 0
    row['e1':'z1'] = dList.loc[ind2,'e1':'z1']
    row['e2':'z2'] = dList.loc[ind2,'e2':'z2']
    row['e3':'z3'] = sList[ind1,:]
    row['DtoT3']   = 1
    return row

def genDtoT(singDF,doubleDF,numD=100):
    '''
    doubOut = genDtoT(singDF,doubleDF,numD)
    
    randomly combines 1 single with 1 double from the same module into D-toT events. Ensures that they aren't too close (within 1/2 pixel) from each other
    
    inputs:
            singDF : pandas dataframe of the singles, with the appropriate labels.
            doubleDF : pandas dataframe of the doubles, with the appropriate labels.
            numD : number of false triple to generate
    outputs:
        doubOut : pandas dataframe of the D-to-T (didn't change variable names), with appropriate labels.
                the event ordering asssumings the double order doesn't change, for example
                a D-to-T-2 is where the single is in place 2, and interactions from the double are 1 and 3
    '''
    singles = np.array(singDF.loc[:,'e1':'z1'])
    doubles = doubleDF.copy()#np.array(doubleDF.loc[:,'e1':'z2'])
   
    sMax = np.min((singles.shape[0],doubles.shape[0]))
    sMods = getModSings(singles)
    dMods,_ = getModl2(np.array(doubles.loc[:,'e1':'z2']))
    
    modTotS = np.zeros((16,))
    modTotD = np.zeros((16,))
    
    # Totals the single and doubles per modulus
    for iii in range(16):
        modTotS[iii] = (sMods==iii).sum()
        modTotD[iii] = (dMods==iii).sum()
        
    modTots = np.max((modTotS,modTotD,modTotD),axis=0)
    # The number of singles requested is divided evenly among the modulus
    numSingPerMod = np.ceil(numD * modTots/(modTots.sum()))
    # Preallocate the suitable dataframe 
    # doubOut = pd.DataFrame(np.zeros((numSingPerMod.sum().astype(int),len(singDF.columns))),columns=singDF.columns)
    doubOut = []

    doubCount = 0
    
    # Loops over each modulus individually
    for jjj in range(16):
        # Add a new frame to the pile
        doubOut.append(
            pd.DataFrame(np.zeros((numSingPerMod[jjj].astype(int),len(singDF.columns))),columns=singDF.columns)
        )
        # Gathers all singles and doubles from that modulus
        sList = singles[sMods==jjj,:]
        # dList = np.array(doubles.loc[:,'e1':'z2'])[dMods==jjj,:]
        # dLabeled = doubles.loc[dMods==jjj,:].reset_index(drop=True)
        dList = doubles.loc[dMods==jjj,:].reset_index(drop=True)
        # Stores the number of singles and doubles for this modulus...
        mTS = modTotS[jjj]
        mTD = modTotD[jjj]
        # Iterates over the NUMBER of singles per modulus
        f = lambda row: findSingleForDouble_par(row=row, sList=sList, dList=dList, 
                mTS=mTS, mTD=mTD)
        g = lambda df: df.apply(f, axis=1)
        # doubOut[-1] = doubOut[-1].apply(f, axis=1)
        doubOut[jjj] = fullParApply([doubOut[jjj]], g)[0]
    doubOut  = concat(doubOut)
    doubtOut = doubOut.sample(frac=1).reset_index(drop=True)
    # print(doubOut.shape)
    doubOut = doubOut.iloc[:numD].copy()
    
#    random.shuffle(doubOut)
#    random.shuffle(doubOut)
#    doubOut = doubOut[:numD,:].copy()
#    dfOut = pd.DataFrame(np.zeros((doubOut.shape[0],len(singsDF.columns))),columns=singsDF.columns)
#    dfOut.loc[:,'e1':'z3']=doubOut.copy()
#    dfOut['PS']=np.nan
#    dfOut['FP']=np.nan
#    dfOut.loc[:,'CO':'e0']=np.nan
    
    return doubOut
    
    
    


def swapOrderDoubles(doubDF1,percentSwap=0.5):
    '''
    doubDF = swapOrderDoubles(doubDF1,percentSwap)
    
    this function swaps a percentage of the input doubles, and returns the complete, modified stack
    
    inputs:
        doubDF1 : doubles dataframe, with labels
        percentSwap : percentage of doubles to swap
    outputs:
        doubDF : modifed data-stack, with swaped and unswaped data, and modified labels
    
    '''
    doubDF = doubDF1.copy()
    randList = np.random.rand(doubDF.shape[0])
    swapInd = randList<=percentSwap
    temp = np.array(doubDF.loc[swapInd,'e2':'z2']).copy()
    doubDF.loc[swapInd,'e2':'z2'] = np.array(doubDF.loc[swapInd,'e1':'z1']).copy()
    doubDF.loc[swapInd,'e1':'z1'] = temp.copy()
    doubDF.loc[swapInd,'CO'] = 0
    doubDF.loc[swapInd,'WOD'] = 1
            
    return doubDF          


def swapOrderTriples(tripDF1,percentSwap=0.8333):
    '''
    tripDF = swapOrderTriples(tripDF1,percentSwap)
    
    this function swaps a percentage of the input triples, and returns the complete, modified stack
    
    inputs:
        tripDF1 : triples dataframe, with labels
        percentSwap : percentage of doubles to swap
    outputs:
        tripDF : modifed data-stack, with swaped and unswaped data, and modified labels
            note: for the swapped triples, it randomly selects between the 5 other orders,
            that is 1-3-2, 2-1-3, 2-3-1, 3-1-2, 3-2-1.
    
    '''
    tripDF = tripDF1.copy()
    randList = np.random.rand(tripDF.shape[0])

    swapInd132 = randList<=percentSwap/5
    swapInd213 = (randList>percentSwap/5) & (randList<=2*percentSwap/5)
    swapInd231 = (randList>2*percentSwap/5) & (randList<=3*percentSwap/5)
    swapInd312 = (randList>3*percentSwap/5) & (randList<=4*percentSwap/5)
    swapInd321 = (randList>4*percentSwap/5) & (randList<=percentSwap)
    
    temp1 = np.array(tripDF.loc[:,'e1':'z1']).copy()
    temp2 = np.array(tripDF.loc[:,'e2':'z2']).copy()
    temp3 = np.array(tripDF.loc[:,'e3':'z3']).copy()
    
    #132
    tripDF.loc[swapInd132,'e1':'z1'] = temp1[swapInd132,:].copy()
    tripDF.loc[swapInd132,'e2':'z2'] = temp3[swapInd132,:].copy()
    tripDF.loc[swapInd132,'e3':'z3'] = temp2[swapInd132,:].copy()
    tripDF.loc[swapInd132,'WOT132'] = 1
    
    #213
    tripDF.loc[swapInd213,'e1':'z1'] = temp2[swapInd213,:].copy()
    tripDF.loc[swapInd213,'e2':'z2'] = temp1[swapInd213,:].copy()
    tripDF.loc[swapInd213,'e3':'z3'] = temp3[swapInd213,:].copy()
    tripDF.loc[swapInd213,'WOT213'] = 1
    
    #231
    tripDF.loc[swapInd231,'e1':'z1'] = temp2[swapInd231,:].copy()
    tripDF.loc[swapInd231,'e2':'z2'] = temp3[swapInd231,:].copy()
    tripDF.loc[swapInd231,'e3':'z3'] = temp1[swapInd231,:].copy()
    tripDF.loc[swapInd231,'WOT231'] = 1
    
    #312
    tripDF.loc[swapInd312,'e1':'z1'] = temp3[swapInd312,:].copy()
    tripDF.loc[swapInd312,'e2':'z2'] = temp1[swapInd312,:].copy()
    tripDF.loc[swapInd312,'e3':'z3'] = temp2[swapInd312,:].copy()
    tripDF.loc[swapInd312,'WOT312'] = 1
    
    #321
    tripDF.loc[swapInd321,'e1':'z1'] = temp3[swapInd321,:].copy()
    tripDF.loc[swapInd321,'e2':'z2'] = temp2[swapInd321,:].copy()
    tripDF.loc[swapInd321,'e3':'z3'] = temp1[swapInd321,:].copy()
    tripDF.loc[swapInd321,'WOT321'] = 1
    
    tripDF.loc[randList<=percentSwap,'CO']=0
            
    return tripDF        

def balanceData(row, size=None):
    """
    
    Takes in a pandas row and changes it based on index in the dataframe.

    Does a lot of nice things. Auto balances Triples, Doubles, and DtoT for
    machine learning algorithm.
    """
    from makeNewLabels import rowToRow
    from numpy import isnan
    # Is this row a double or triple?
    if row["TD"] != 1 and not isnan(row['FP']): # DtoT, tTriple:
        # TT and DtoT we swap them in 1/6 portions for all possible cases
        toSwap = int(size*1/6)
        index = row.name
        order = int(index / toSwap)
        # Easier to have this then 5 if statements
        aOrders = {
                1: [("1", "3", "2"),"WOT132", "DtoT2"],
                2: [("2", "1", "3"),"WOT213", "DtoT3"],
                3: [("2", "3", "1"),"WOT231", "DtoT2"],
                4: [("3", "2", "1"),"WOT321", "DtoT1"],
                5: [("3", "1", "2"),"WOT312", "DtoT1"]}
        if order == 0:
            pass
        elif order == 6:
            print("ERROR order == 6?")
            print(index, toSwap, order, size)
        else: # Double to Triple
            # Easily swap the row values using python syntax magic
            for t in ["e", "x", "y", "z"]:
                # Swap row values
                row[t + '1'],                  row[t + '2'],                  row[t + '3'] =( 
                row[t + aOrders[order][0][0]], row[t + aOrders[order][0][1]], row[t + aOrders[order][0][2]] )
            if row["TT"] == 1:
                # Set WOTXXX value
                row[aOrders[order][1]] = 1
            else:
                # Set DtoTX value
                row["DtoT3"] = 0
                row[aOrders[order][2]] = 1
            row["CO"] = 0
        # Extra Distance the third camera event
        row['euc2'] = distCalc(row['x2'],row['y2'],row['z2'],row['x3'],row['y3'],row['z3'])
        row['euc3'] = distCalc(row['x3'],row['y3'],row['z3'],row['x1'],row['y1'],row['z1'])
        row['de2'] = row['e2'] - row['e3']
        # row['de3'] = 0
    elif row['TD'] == 1: # True double
        # TD we swap half of them
        toSwap = int(size*1/2)
        index = row.name
        if row.name < toSwap:
            # Do nothing, here for completeness
            pass
        else:
            row['e1'], row['e2'] = (row['e2'], row['e1'])
            row['x1'], row['x2'] = (row['x2'], row['x1'])
            row['y1'], row['y2'] = (row['y2'], row['y1'])
            row['z1'], row['z2'] = (row['z2'], row['z1'])
            row['CO'] = 0
            row['WOD'] = 1
        row['de2']  = 0
        row['euc2'] = 0
        row['euc3'] = 0
    else:
        # False Triples are the only special case
        if row['e3'] != 0:
            row['de2'] = row['e2'] - row['e3']
            row['euc2'] = distCalc(row['x2'],row['y2'],row['z2'],row['x3'],row['y3'],row['z3'])
            row['euc3'] = distCalc(row['x3'],row['y3'],row['z3'],row['x1'],row['y1'],row['z1'])
        else: # For false doubles we need to see all the same 0s as regular doubles
            row['de2']  = 0
            row['euc2'] = 0
            row['euc3'] = 0
    # Present for doubles, triples, false doubles, false triples
    row['de1']  = row['e1'] - row['e2']
    row['euc1'] = distCalc(row['x1'],row['y1'],row['z1'],row['x2'],row['y2'],row['z2'])
    # Also make new labels while you're here
    return rowToRow(row)

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
    with pool(36) as p:
        for frame in split_frames:
            full_frames.append(concat(p.map(func, frame)))
            # tTriples = concat(p.map(parApply, sTrips))
    return full_frames

def addLabels(row, labels=None, values=None):
    for i in range(len(labels)):
        row[labels[i]] = values[i]
    return row

def addMaggiLabels(detM):
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

if __name__ == "__main__":      
    from sys import argv
    outputDir = argv[argv.index("-o") + 1]
    inputFile = argv[argv.index("-i") + 1]
    from os import mkdir
    try:
        mkdir(outputDir)
    except Exception as e:
        e = str(e)
        if "exists" not in e:
            print(e)
            exit()
    lowEThresh = 0.05
    upperEThresh = 2.7
    #data labels from file: e1,x1,y1,z1,e2,x2,y2,z2,e3,x3,y3,z3,PS,FP,e0
    #PS = primary (1) or scatter (0)
    #FP = full (1) or partial (0)
    #e0 = initial energy
    otherLabels = ['TD','TT','DtoT1','DtoT2','DtoT3','CO','WOD','WOT132','WOT213','WOT231','WOT312','WOT321']
    #TD (data prep label 0) = true double (bool)
    #TT (data prep label 1) = true triple (bool)
    #DtoT1 (data prep label 2) = double-to-triple, 1 is wrong (bool)
    #DtoT2 (data prep label 3) = double-to-triple, 2 is wrong (bool)
    #DtoT3 (data prep label 4) = double-to-triple, 3 is wrong (bool)
    #CO (data prep label 5) = correct order? (bool)
    #WOD (data prep label 6) = wrong order double (bool)
    #WOT132 (data prep label 7) = wrong order triple, 1,3,2 (bool) 
    #WOT213 (data prep label 8) = wrong order triple, 2,1,3 (bool) 
    #WOT231 (data prep label 9) = wrong order triple, 2,3,1 (bool) 
    #WOT312 (data prep label 10) = wrong order triple, 3,1,2 (bool) 
    #WOT321 (data prep label 11) = wrong order triple, 3,2,1 (bool) 
    # detM = pd.read_csv(r'/umbc/xfs1/gobbert/common/research/oncology/maggi/data/1GEvents_150MeV_sp_E0_nulls.csv')
    # detM = pd.read_csv(r'/umbc/xfs1/gobbert/common/research/oncology/maggi/data/5GEvents_150MeV_sp_E0_nulls.csv')
    print("Loading initial data file...")
    # detM = pd.read_csv(r'/umbc/xfs1/gobbert/common/research/oncology/maggi/data/40GEvents_150MeV_sp_E0_nulls.csv')
    detM = pd.read_csv(inputFile)
    
    # Setup the singles
    sings = detM[detM['e2'].isnull()]
    sings = sings[(sings['e1']>=lowEThresh) & (sings['e1']<=upperEThresh)]
    singsData = np.zeros((sings.shape[0],12))
    # singsData[:,5] = 1
    # singsData[:,5:] = np.nan
    s = np.zeros(12)
    s[5]  = 1
    s[5:] = np.nan
    # for iii in range(len(otherLabels)):
        # sings.insert(sings.shape[1]-1,otherLabels[iii],singsData[:,iii])
    print("Adding labels to singles...")
    f = lambda row: addLabels(row, labels=otherLabels, values=s)
    sings = fullParApply([sings], f)[0]
    sings = sings.sample(frac=1).reset_index(drop=True)


    # Setup the doubles
    doubs = detM[detM['e3'].isnull() & detM['e2'].notnull()]
    doubs = doubs[(doubs['e1']>=lowEThresh) & (doubs['e1']<=upperEThresh) & (doubs['e2']>=lowEThresh) & (doubs['e2']<=upperEThresh)]
    # doubsData = np.zeros((doubs.shape[0],12))
    # doubsData[:,0] = 1
    # doubsData[:,5] = 1
    # doubsData[:,7:] = np.nan
    d = np.zeros(12)
    d[0], d[5] = (1, 1)
    d[7:] = np.nan
    # for iii in range(len(otherLabels)):
        # doubs.insert(doubs.shape[1]-1,otherLabels[iii],doubsData[:,iii])
    print("Adding labels to doubles...")
    f = lambda row: addLabels(row, labels=otherLabels, values=d)
    doubs = fullParApply([doubs], f)[0]
    doubs = doubs.sample(frac=1).reset_index(drop=True)

    # Setup the triples
    trips = detM[detM['e3'].notnull() & detM['e2'].notnull()]
    trips = trips[(trips['e1']>=lowEThresh) & (trips['e1']<=upperEThresh) & (trips['e2']>=lowEThresh) & (trips['e2']<=upperEThresh) & (trips['e3']>=lowEThresh) & (trips['e3']<=upperEThresh)]
    t = np.zeros(12)
    t[1], t[5] = (1, 1)
    # for iii in range(len(otherLabels)):
        # trips.insert(trips.shape[1]-1,otherLabels[iii],tripsData[:,iii])
    print("Adding labels to triples...")
    f = lambda row: addLabels(row, labels=otherLabels, values=t)
    trips = fullParApply([trips], f)[0]
    trips = trips.sample(frac=1).reset_index(drop=True)
        

        
    ### The full suite of preparation before one hot
    # First we need to split the triples in 1/6 portions
    # Get all data
    # aTriples = swapOrderTriples( trips, percentSwap=.50)
    print("Computing sundry values for downsizing data...")
    tsize = trips.shape[0]
    oneSixth = tsize // 6
    # We have a ton of doubles so just split them in half
    doubs1 = doubs.loc[doubs.index  < int(doubs.shape[0]/2)]
    doubs2 = doubs.loc[doubs.index >= int(doubs.shape[0]/2)]
    # Compute DtoT
    dsize = doubs1.shape[0]
    tsize = oneSixth*6
    print("Generating DtoT samples...")
    # Just trash the extra data
    aDtoT = genDtoT(sings, doubs2, numD=tsize)
    print("Removing extra samples...")
    trips  =  trips.loc[trips.index  < tsize]
    doubs1 = doubs1.loc[doubs1.index < 2*oneSixth]
    aDtoT  =   aDtoT.loc[aDtoT.index < tsize]


    # tTriples = concat(p.map(parApply, sTrips))
    # Add all of the labels
    parF = lambda df: parApply(df, sz=tsize, axis=1)
    print("Balancing triples, doubles to triples...")
    tTriples, aDtoT = fullParApply(
        [trips, aDtoT], parF)

    print("Balancing doubles...")
    parF = lambda df: parApply(df, sz=2*oneSixth, axis=1)
    tDoubles = fullParApply([doubs1], parF)[0]
        # [trips, doubs1, aDtoT], parApply)
    print("Generating false doubles and triples...")
    # Generate the correct number of false cases
    fDoubles = genFalseDoub(sings, oneSixth//2)
    fTriples = genFalseTrip(sings, oneSixth)
    print("Adding labels and balancing false doubles and false triples...")
    # Standard parApply does not work the way we need
    parF = lambda df: parApply(df, sz=oneSixth*2, axis=1)
    fTriples = fullParApply(
            [fTriples], parF)[0]
    parF = lambda df: parApply(df, sz=oneSixth//2, axis=1)
    fDoubles = fullParApply(
            [fDoubles], parF)[0]
    # Convert empty things into NaN
    tTriples = tTriples.replace(r'^\s+$', np.nan, regex=True)
    fTriples = fTriples.replace(r'^\s+$', np.nan, regex=True)
    tDoubles = tDoubles.replace(r'^\s+$', np.nan, regex=True)
    fDoubles = fDoubles.replace(r'^\s+$', np.nan, regex=True)
    aDtoT    = aDtoT.replace(r'^\s+$', np.nan, regex=True)
    # Add my new labels and write all of them outs
    # aDtoT    =    aDtoT.apply(rowToRow, axis=1)
    # Write out three separate data pieces to be processed later
    print("Saving all data")
    tTriples.to_csv("{}/tTriples_all.csv".format(outputDir), na_rep='NaN', index=False)
    tDoubles.to_csv("{}/tDoubles_all.csv".format(outputDir), na_rep='NaN', index=False)
    fDoubles.to_csv("{}/fDoubles_all.csv".format(outputDir), na_rep='NaN', index=False)
    fTriples.to_csv("{}/fTriples_all.csv".format(outputDir), na_rep='NaN', index=False)
    aDtoT.to_csv(      "{}/aDtoT_all.csv".format(outputDir), na_rep='NaN', index=False)

    # Write it out with only my intended labels

    allData = concat([tTriples, fTriples, aDtoT, tDoubles, fDoubles])
    # allData.to_csv("./balanced5/all_data_5G.csv", na_rep="NaN", index=False)
    # For easier one hot with doubles
    allData.to_csv("{}/all_data_40G.csv".format(outputDir), na_rep="0.0", index=False)
    allData.to_csv("{}/all_data_40G_new_labels.csv".format(outputDir), columns=[
         'de1','euc1','e1', 'x1','y1','z1',
         'de2','euc2','e2','x2','y2', 'z2', 
         'euc3','e3','x3','y3','z3',
         'e0', 'perm_index'],
            na_rep="NaN", index=False)
