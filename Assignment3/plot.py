import matplotlib
matplotlib.use('Agg')

import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

sns.set()


datafile = "Temp_output.txt"

for i in range(1):
    try:
        demo_input_format = pd.read_csv(datafile, sep=" ", header=None)
        demo_input_format.columns = ["P", "ppn", "Time"]
        
        sns.barplot(x="ppn", y="Time", data=demo_input_format,  hue="P")
        filename = str("plot.jpg")
        plt.savefig(filename)
        plt.cla()
    except:
        print("Error Occured")
