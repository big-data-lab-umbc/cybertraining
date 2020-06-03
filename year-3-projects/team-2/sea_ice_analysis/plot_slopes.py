import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm

df = pd.read_csv("2007_over-water_worldview_anomalies_with_sea_ice_and_slopes.csv")

t = df.timestamp
cth = df.cloud_top_height
d = df.distance
m = df.cloud_top_slope
b = df.cloud_top_intercept

for i in tqdm(range(df.index.size), total = df.index.size):
    
    x1 = np.array(d[i][1:-1].split()).astype(np.float32)
    y1 = np.array(cth[i][1:-1].split()).astype(np.float32)
    
    x2 = np.linspace(x1.min(), x1.max(), 1000)
    y2 = m[i] * x2 + b[i]
    
    plt.figure()
    plt.title(t[i])
    plt.plot(x1, y1, c = "b", zorder = 2)
    plt.scatter(x1, y1, c = "b", edgecolors = "k", zorder = 3)
    plt.plot(x2, y2, ls = "--", c = "k", label = f"y = {m[i]:0.4f}x + {b[i]:0.4f}")
    plt.axvline(x1[5], c = "r", ls = "--")
    plt.xlim(x1.min(), x1.max())
    plt.xlabel("Displacement from Anomaly Along Track (m)")
    plt.ylabel("Cloud Top Height (m)")
    plt.legend()
    plt.savefig(f"test_{i:02d}.png")
    plt.close()

