#!/usr/bin/env python
# coding: utf-8

# In[47]:

import matplotlib
matplotlib.use('Agg')

import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

sns.set()


dataFiles = ["Data_Bcast.txt", "Data_Reduce.txt", "Data_Gather.txt", "Data_Alltoallv.txt"]
plotName = ["plot_Bcast.jpg", "plot_Reduce.jpg", "plot_Gather.jpg", "plot_Alltoallv.jpg"]

for i in range(4):
    try:
        demo_input_format = pd.read_csv(dataFiles[i], sep=", ", header=None)
        demo_input_format.columns = ["D", "P", "ppn", "mode", "Time"]
        demo_input_format["(P, ppn)"] = list(map(lambda x, y: ("(" + x + ", " + y + ")"), map(str, demo_input_format["P"]), map(str, demo_input_format["ppn"])))
        print(demo_input_format)
        sns.barplot(x="(P, ppn)", y="Time", data=demo_input_format,  hue="mode")
        filename = str(plotName[i])
        plt.savefig(filename)
        plt.cla()
    except:
        print("Error Occured")



