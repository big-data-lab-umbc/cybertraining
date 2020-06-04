from timeit import timeit

number = 10000

setup = f"""
import numpy as np
"""

stmt1 = """
np.logspace(0, 15, num = 16, base = 2, dtype = np.uint16)
"""

stmt2 = """
2 ** np.arange(16)
"""

kwargs  = {"setup" : setup, "number" : number}
kwargs1 = {**kwargs, "stmt" : stmt1}
kwargs2 = {**kwargs, "stmt" : stmt2}

t1 = timeit(**kwargs1)
t2 = timeit(**kwargs2)

T1 = t1 / number
T2 = t2 / number

print("numpy.binary_repr:     {:.8f} seconds".format(T1))
print("pandas.Series mapping: {:.8f} seconds".format(T2))
