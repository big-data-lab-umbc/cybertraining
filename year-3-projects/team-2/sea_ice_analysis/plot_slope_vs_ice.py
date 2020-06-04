import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

df = pd.read_csv("2007_over-water_worldview_anomalies_with_sea_ice_and_slopes.csv")


plt.figure()

plt.scatter(df["sea_ice_concentration"], np.abs(df["cloud_top_slope"]), c = df["sza"], cmap = "jet", s = 6, vmin = 45, vmax = 90)
plt.ylim(0, 2.5)
plt.xlim(-1, 101)
plt.xticks(list(range(0, 101, 10)))
plt.xlabel("Sea Ice Concentration (%)")
plt.ylabel("Absolute Cloud Top Height Slope")
cbar = plt.colorbar()
cbar.set_label("Solar Zenith Angle")

nan_mask = (~df["sea_ice_concentration"].isna()) & (~np.abs(df["cloud_top_slope"]).isna())

plt.figure()
plt.hist2d(df["sea_ice_concentration"][nan_mask], np.abs(df["cloud_top_slope"][nan_mask]), norm = matplotlib.colors.LogNorm(), vmin = 1e0, vmax = 1e4, cmap = "jet")
cbar = plt.colorbar()
cbar.set_label("Observations")
plt.ylim(0, 3)
plt.xlim(0, 100)
plt.xlabel("Sea Ice Concentration (%)")
plt.ylabel("Absolute Cloud Top Height Slope")

plt.show()
