from pandas import DataFrame as df
import pandas as pd
import numpy
from numpy import array
from numpy import zeros
from sys import argv
def loadFiles(files=None):
    files = list(map(lambda x: pd.read_csv(x), files))
        # ["tTriples_all.csv", "fTriples_all.csv", "aDtoT_all.csv"]))
    return files
def merge(files=None):
    data = pd.concat(files)
    return data

def toOneHotArray(data, outputs=True):
    x = []
    y = []
    for row in data.itertuples():
       # print(row)
        # exit()
        # Pandas(Index=12539, e1=0.14648599999999998, x1=-36.7706, y1=295.541,
                # z1=-17.8135, e2=0.0770949, x2=-32.7857, y2=295.127, z2=-21.3662,
                # e3=0.07502, x3=-36.7693, y3=295.462, z3=-18.0182, PS=1.0,
                # FP=0.0, TD=0.0, TT=1.0, DtoT1=0.0, DtoT2=0.0, DtoT3=0.0,
                # CO=0.0,^C WOD=0.0, WOT132=1.0, WOT213=0.0, WOT231=0.0,
                # WOT312=0.0, WOT321=0.0, e0=0.510999, euc2=5.214441289342512,
                # euc3=0.2194191878573996, euc1=5.354670979621437, _31=1.0,
                # _32=0.0, _33=0.0, _34=0.0, _35=0.0, _36=0.0, _37=1.0, _38=0.0,
                # _39=0.0, _40=1.0, _41=0.0, _42=0.0)
        x.append([row.de1, row.euc1, row.e1, row.x1, row.y1, row.z1,
                       row.de2, row.euc2, row.e2, row.x2, row.y2, row.z2,
                                row.euc3, row.e3, row.x3, row.y3, row.z3])
        # y.append(row[-12:-8])
        # y.append(row[-12:-4])
        # 4 permute 3
        if outputs:
            t = zeros(25)
            t[int(row[-1])] = 1
            # 123 124 132 134 142 143 213 214 231 234 241 243 312 314 321 324 341 342 412 413 421 423 431 432
            y.append(t)
            # Slice the first camera event and third camera event
            # y.append(numpy.array(row[-12:])[[*list(range(-12,-8)),*list(range(-4,0))]])
    y = numpy.array(y)
    x = numpy.array(x)
    # print(x.shape)
    # xr = x.reshape(x.shape[0], 3, 4, 1)
    # print(x[0])
    # y = y.reshape(y.shape[0], 3, 4)
    # print(y.shape)
    # print(y[0])
    return x, y

def saveData(data, outName=None, scaler=None):
    if outName is None:
        outputName = "output.out.npy"
        inputName  = "input.in.npy"
    else:
        outputName = "{}.out.npy".format(outName)
        inputName  = "{}.in.npy".format(outName)
    # catch it before it breaks.
    assert not numpy.any(numpy.isnan(data[0]))
    assert not numpy.any(numpy.isnan(data[1]))
    numpy.save(inputName, data[0])
    numpy.save(outputName, data[1])
    if scaler is not None:
        with open("{}.scaler".format(outName), "wb") as outfile:
            outfile.write(scaler.getbuffer())

def normalizeData(data):
    from sklearn.preprocessing import StandardScaler
    import dill as pickle
    from io import BytesIO
    # https://stackoverflow.com/questions/41165642/how-to-normalize-input-of-neural-network-predicting-stock-market-python
    scale  = StandardScaler(with_mean=0, with_std=1)
    xscale = scale.fit_transform(data[0])
    scale_bytes = BytesIO()
    pickle.dump(scale, scale_bytes)
    # xscale = xscale.reshape(xscale.shape[0], 3, 4, 1)
    print(xscale.shape)
    print(xscale[0])
    return xscale, data[1], scale_bytes

if __name__ == "__main__":
    from numpy import load, save, concatenate
    import dill as pickle
    # rootPath = "/umbc/xfs1/gobbert/common/research/oncology/maggi/data/40Gprocessed_permute_labels-3"
    rootPath = "/umbc/xfs1/gobbert/common/research/oncology/maggi/data/40Gprocessed_permute_labels-2"
    # All files 
    allFiles = ["tTriples_all.csv","tDoubles_all.csv", "fTriples_all.csv",
            "fDoubles_all.csv", "aDtoT_all.csv"]
    allData = list(map(lambda x: pd.read_csv("{}/{}".format(rootPath, x)), allFiles))
    allBits = list(map(lambda x: toOneHotArray(x), allData))
    bigData = concatenate([v[0] for v in allBits])

    # normalize
    _, _, scaler_bytes = normalizeData((bigData, array([])))
    scaler = pickle.loads(scaler_bytes.getbuffer())
    allBits_scaled = list(map(lambda x: (scaler.transform(x[0]), x[1]), allBits))
    # Save all data
    with open("{}/newest.scaler".format(rootPath), "wb") as w:
        w.write(scaler_bytes.getbuffer())
    for i, filename in enumerate(allFiles):
        save("{}/inputs.{}.npy".format(rootPath, filename.split(".")[0]),
                allBits_scaled[i][0])
        save("{}/outputs.{}.npy".format(rootPath, filename.split(".")[0]),
                allBits_scaled[i][1])
