from timeit import timeit

setup = """
import numpy as np
import pandas as pd
x = np.arange(1000)
"""

number = 1000

stmt1 = """
y = np.array(list(map(np.binary_repr, x)))
"""

stmt2 = """
y = pd.Series(x).map(np.binary_repr)
"""

stmt3 = """
a = np.asarray(x)
a = a.astype(f">u{a.itemsize}")

bits = a[np.newaxis].view(np.uint8)
bits = bits.reshape((*a.shape, a.itemsize))
bits = np.unpackbits(bits, axis = -1)
"""

kwargs  = {"setup" : setup, "number" : number}
kwargs1 = {**kwargs, "stmt" : stmt1}
kwargs2 = {**kwargs, "stmt" : stmt2}
kwargs3 = {**kwargs, "stmt" : stmt3}

t1 = timeit(**kwargs1)
t2 = timeit(**kwargs2)
t3 = timeit(**kwargs3)

T1 = t1 / number
T2 = t2 / number
T3 = t3 / number

print("numpy.binary_repr:     {:.8f} seconds".format(T1))
print("pandas.Series mapping: {:.8f} seconds".format(T2))
print("numpy.unpackbits:      {:.8f} seconds".format(T3))
