# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import matplotlib.pyplot as plt

e5_LR001 = [0.825,0.840,0.846,0.7982, 0.8508,0.7463]
e10_LR001 = [0.8978,0.8680,0.8503,0.8072,0.8362,0.8153]
e15_LR001 = [0.9502,0.8905,0.8592,0.8495,0.8322,0.7922]

e5_LR001_std = [0.0620, 0.0499, 0.0362, 0.1052, 0.0291, 0.0832]
e10_LR001_std = [0.0379, 0.0476, 0.0464, 0.0724, 0.0290, 0.0515]
e15_LR001_std = [0.0220, 0.0553, 0.0493, 0.0293, 0.0492, 0.0652]

e5_LR005 = [0.9040,0.9002, 0.8823, 0.8409, 0.8060, 0.8048]
e10_LR005 = [0.9497, 0.9491, 0.9278, 0.8542, 0.8503, 0.8457]
e15_LR005 = [0.9512, 0.9621, 0.9626, 0.9469, 0.8618, 0.8372]

e5_LR005_std = [0.0372, 0.0365, 0.0285, 0.0508, 0.0872, 0.0676]
e10_LR005_std = [0.0220, 0.0113, 0.0248, 0.1003, 0.0513, 0.0403]
e15_LR005_std = [0.0233, 0.0122, 0.0114, 0.0213, 0.0859, 0.0502]

e5_LR01 = [0.8747, 0.9146, 0.8706, 0.8595, 0.8107, 0.7879]
e10_LR01 = [0.9097, 0.9396, 0.8949, 0.8735, 0.8557, 0.8424]
e15_LR01 = [0.9295, 0.9241, 0.9428, 0.9421, 0.8396, 0.8695]

e5_LR01_std = [0.0335, 0.0324, 0.0528, 0.0455, 0.0756, 0.0752]
e10_LR01_std = [0.0400, 0.0166, 0.0398, 0.0614, 0.0363, 0.0447]
e15_LR01_std = [0.0200, 0.0232, 0.0141, 0.0155, 0.1210, 0.0272]

batch_size = [128,256,512,1024,2048, 4096]
batch_time = [3864.4712, 3319.0052, 2949.0057, 2779.8445, 2678.4618, 2643.0482]

batch_size_2018 = [512, 1024, 2048, 4096]
batch_time_2018 = [1662.0, 1225, 1304, 1674]

plt.figure(figsize=(10,7))
plt.subplot(2,1,1)
plt.title("Epochs = 5", fontsize = 18)
plt.ylabel("Mean Accuracy", fontsize = 16)
plt.semilogx(batch_size, e5_LR001, label = "Learning Rate = 0.0001", color = "black")
plt.semilogx(batch_size, e5_LR005, label = "Learning Rate = 0.0005", color = "blue")
plt.semilogx(batch_size, e5_LR01, label = "Learning Rate = 0.001", color = "red")
plt.xticks(batch_size, labels = batch_size)
plt.legend(prop={'size': 14})

plt.subplot(2,1,2)
plt.semilogx(batch_size, e5_LR001_std, color = "black", marker='o', linewidth=0)
plt.semilogx(batch_size, e5_LR005_std, color = "blue", marker='o', linewidth=0)
plt.semilogx(batch_size, e5_LR01_std, color = "red", marker='o', linewidth=0)
plt.xlabel("Batch Size", fontsize = 16)
plt.ylabel("Standard Deviation", fontsize = 16)
plt.xticks(batch_size, labels = batch_size)
plt.tight_layout()
plt.savefig('/Users/Sarah/Desktop/E5.png')


plt.figure(figsize=(10,7))
plt.subplot(2,1,1)
plt.title("Epochs = 10", fontsize = 18)
plt.ylabel("Mean Accuracy", fontsize = 16)
plt.semilogx(batch_size, e10_LR001, label = "Learning Rate = 0.0001", color = "black")
plt.semilogx(batch_size, e10_LR005, label = "Learning Rate = 0.0005", color = "blue")
plt.semilogx(batch_size, e10_LR01, label = "Learning Rate = 0.001", color = "red")
plt.legend(prop={'size': 14})
plt.xticks(batch_size, labels = batch_size)

plt.subplot(2,1,2)
plt.semilogx(batch_size, e10_LR001_std, color = "black", marker='o', linewidth=0)
plt.semilogx(batch_size, e10_LR005_std, color = "blue", marker='o', linewidth=0)
plt.semilogx(batch_size, e10_LR01_std, color = "red", marker='o', linewidth=0)
plt.xlabel("Batch Size", fontsize = 16)
plt.ylabel("Standard Deviation", fontsize = 16)
plt.xticks(batch_size, labels = batch_size)
plt.tight_layout()
plt.savefig('/Users/Sarah/Desktop/E10.png')

plt.figure(figsize=(10,7))
plt.subplot(2,1,1)
plt.title("Epochs = 15", fontsize = 18)
plt.ylabel("Mean Accuracy", fontsize = 16)
plt.semilogx(batch_size, e15_LR001, label = "Learning Rate = 0.0001", color = "black")
plt.semilogx(batch_size, e15_LR005, label = "Learning Rate = 0.0005", color = "blue")
plt.semilogx(batch_size, e15_LR01, label = "Learning Rate = 0.001", color = "red")
plt.legend(prop={'size': 14})
plt.xticks(batch_size, labels = batch_size)

plt.subplot(2,1,2)
plt.semilogx(batch_size, e15_LR001_std, color = "black", marker='o', linewidth=0)
plt.semilogx(batch_size, e15_LR005_std, color = "blue", marker='o', linewidth=0)
plt.semilogx(batch_size, e15_LR01_std, color = "red", marker='o', linewidth=0)
plt.xlabel("Batch Size", fontsize = 16)
plt.ylabel("Standard Deviation", fontsize = 16)
plt.xticks(batch_size, labels = batch_size)
plt.tight_layout()
plt.savefig('/Users/Sarah/Desktop/E15.png')

plt.figure(figsize=(10,5))
plt.title("Train Time per Batch", fontsize = 18)
plt.ylabel("Train Time (s)", fontsize = 16)
plt.xlabel("Batch Size", fontsize = 16)
plt.semilogx(batch_size, batch_time, color = "red", marker='o', label = "2013")
plt.semilogx(batch_size_2018, batch_time_2018, color = "blue", marker='o', label = "2018")
plt.xticks(batch_size, labels = batch_size)
plt.legend(prop={'size': 14})

plt.tight_layout()
plt.savefig('/Users/Sarah/Desktop/TrainTime.png')

