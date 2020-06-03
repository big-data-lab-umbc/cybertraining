import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

df = pd.read_csv("2007_over-water_worldview_anomalies_with_sea_ice.csv")

ice = df["sea_ice_concentration"]
ice.loc[ice.isnull()] = -20

bin_edges = np.arange(-20, 105, 5)

H, bins = np.histogram(ice, bins = bin_edges)


# This is achieved by dividing the count by the number of
# observations times the bin width and not dividing by the total number of observations.
plt.figure()
plt.title("Anomalies Over Sea Ice")
p = plt.hist(ice, bins = bin_edges, edgecolor = "black", linewidth = 1.2)
#p = plt.bar(H / 

xlocs, xlabels = plt.xticks()
plt.xticks(xlocs, [xlocs[0], "No Data", *xlocs[2:]])
plt.xlim(-20, 100)
plt.xlabel("Sea Ice Concentration (%)")
plt.ylabel("Anomaly Frequency (%)")
#plt.ylabel("Anomaly Counts")
p[-1][0].set_facecolor("red")
p[-1][0].set_edgecolor("black")

for a in p[-1]:

    a.set_height(a.get_height() / ice.size)
    
plt.ylim(0, 1)

plt.yticks(np.arange(0, 1.01, 0.1))
ylocs, ylabels = plt.yticks()
plt.yticks(ylocs, [f"{int(n * 100)}" for n in ylocs])

#p.patches[0].set_facecolor("red")
#p.patches[0].set_edgecolor("black")
#p[-1][0].set_linewidth(1.2)
plt.savefig("anomalies_over_sea_ice.png")
