# All labels will be converted into this format
# 1->1,1->2,1->3,2->1,2->2,2->3,3->1,3->2,3->1, PS, FP
# Take in a data frame and append these new columns to it
from itertools import combinations as combos
from multiprocessing_on_dill import Pool as pool
from numpy import array_split
from pandas import DataFrame
from numpy import isnan
from helpers import distCalc
import numpy as np
vals = []
# Adding the extra 1/0 allows for softmax classification of doesn't belong or
# "invalid entry"
Z = [0,0,0,1]
Z = "4"
O = [1,0,0,0]
O = "1"
T = [0,1,0,0]
T = "2"
R = [0,0,1,0]
R = "3"
# 123 124 132 134 142 143 213 214 231 234 241 243 312 314 321 324 341 342 412 413 421 423 431 432
#   0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22  23
index = {'123': '0', '124': '1', '132': '2', '134': '3', '142': '4', '143': '5', 
'213': '6', '214': '7', '231': '8', '234': '9', '241': '10', '243': '11', 
'312':'12', '314': '13', '321': '14', '324': '15', '341': '16', '342': '17', 
'412': '18', '413': '19', '421': '20', '423': '21', '431': '22', '432': '23', 
"444":"24"}

def rowToRow(row):
    """
    Give a data frame and an i,j label we interpret whether or not i->j is 0 or
    1 based on already existing columns of a single row
    """
    if row['TD'] == 1:
        if row['CO'] == 1:
            oneHot = "{}{}{}".format(O, T, Z)
        else:
            oneHot = "{}{}{}".format(T, O, Z)
    elif not isnan(row['FP']): # DtoT, tTriple
        # It is a true triple
        if row['TT'] == 1:
            if row['CO'] == 1:
                oneHot = "{}{}{}".format(O, T, R)
                # oneHot = 0
            elif row['WOT132'] == 1:
                oneHot = "{}{}{}".format(O, R, T)
                # oneHot = 2
            elif row['WOT213'] == 1:
                oneHot = "{}{}{}".format(T, O, R)
                # oneHot = 6
            elif row['WOT231'] == 1:
                oneHot = "{}{}{}".format(T, R, O)
            elif row['WOT312'] == 1:
                oneHot = "{}{}{}".format(R, O, T)
            elif row['WOT321'] == 1:
                oneHot = "{}{}{}".format(R, T, O)
        else: # It is a DtoT
            if row['DtoT1'] == 1:
                if row['CO'] == 1:
                    oneHot = "{}{}{}".format(Z, O, T)
                else:
                    oneHot = "{}{}{}".format(Z, T, O)
            elif row['DtoT2'] == 1:
                if row['CO'] == 1:
                    oneHot = "{}{}{}".format(O, Z, T)
                else:
                    oneHot = "{}{}{}".format(T, Z, O)
            else: #row['DtoT3'] == 0:
                if row['CO'] == 1:
                    oneHot = "{}{}{}".format(O, T, Z)
                else:
                    oneHot = "{}{}{}".format(T, O, Z)
    else:
        # False Event
        oneHot = "{}{}{}".format(Z, Z, Z)
    oneHot = index[oneHot]
    row['perm_index'] = oneHot
    return row

def computePermuteLabels(row, outputs=True):
    """
    
    Takes in a pandas row and changes it based on index in the dataframe.

    Does a lot of nice things. Auto balances Triples, Doubles, and DtoT for
    machine learning algorithm.
    """
    from numpy import isnan
    # Is this row a double or triple?
    if row["TD"] != 1 and not isnan(row['FP']): # DtoT, tTriple:
        # Extra Distance the third camera event
        row['euc2'] = distCalc(row['x2'],row['y2'],row['z2'],row['x3'],row['y3'],row['z3'])
        row['euc3'] = distCalc(row['x3'],row['y3'],row['z3'],row['x1'],row['y1'],row['z1'])
        row['de2'] = row['e2'] - row['e3']
        # row['de3'] = 0
    elif row['TD'] == 1: # True double
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
    if outputs:
        return rowToRow(row)
    else:
        return row
