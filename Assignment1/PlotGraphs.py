import matplotlib
matplotlib.use('Agg')

import seaborn as sns
import csv 
import pandas as pd 
  
import matplotlib.pyplot as plt

import numpy as np

#sns.set_theme(style="whitegrid")



def readFromFile(fname, arr):
	with open(fname, 'r') as file:
	    reader = csv.reader(file, delimiter = ',')
	    for row in reader:
	        new = []
	        for i in row:
	        	try:
	        		new.append(int(i))
	        	except:
	        		new.append(float(i))
	        try:
	        	arr[new[0]].append([new[1], new[2]])
	        except:
	        	print("Outlier")


def complete(dis, vac, plots):
	for  i in dis[plots]:
		vac[i[0]].append(i[1])




plots = [16, 36, 49,64]
position = [16, 32, 64, 128, 256, 512, 1024]

vec = {} ; pack = {}; norml = {}

for i in plots:
	vec[i] = []; pack[i] = [] ; norml[i] = []
readFromFile("Data_VectorVersion.txt", vec)
readFromFile("Data_NormalVersion.txt", norml)
readFromFile("Data_PackVersion.txt", pack)



VposiTime = {}; NposiTime = {}; PposiTime = {}
for i in position:
	VposiTime[i] = []; NposiTime[i] = []; PposiTime[i] = []
for pl in plots:
	complete(vec, VposiTime, pl)
	complete(pack, PposiTime, pl)
	complete(norml, NposiTime, pl)

	posi = []
	for i in position:
		posi.append(str(i))


	Vdata_to_plot = [ VposiTime[16], VposiTime[32], VposiTime[64], VposiTime[128], VposiTime[256], VposiTime[512], VposiTime[1024]]
	Pdata_to_plot = [ PposiTime[16], PposiTime[32], PposiTime[64], PposiTime[128], PposiTime[256], PposiTime[512], PposiTime[1024]]
	Ndata_to_plot = [ NposiTime[16], NposiTime[32], NposiTime[64], NposiTime[128], NposiTime[256], NposiTime[512], NposiTime[1024]]

	van = []
	for i in VposiTime.keys():
		for j in VposiTime[i]:
			van.append([str(i),j])

	nor = []
	for i in NposiTime.keys():
		for j in NposiTime[i]:
			nor.append([str(i),j])
	pak = []
	for i in PposiTime.keys():
		for j in PposiTime[i]:
			pak.append([str(i),j])

	dfV = pd.DataFrame(van)
	dfN = pd.DataFrame(nor)
	dfP = pd.DataFrame(pak)
	dfV[3] = 'Vector'
	dfN[3] = 'Normal'
	dfP[3] = 'Pack'
	frames = [dfN, dfV, dfP]
	df = pd.concat(frames)
	df.sort_values(by = [0])
	ax = sns.boxplot(x=0, y=1, hue = 3,data=df).set(xlabel='Size of Row/Column', ylabel='Time')
	filename = str("Plot") + str(pl) + str(".png")
	plt.savefig(filename)
	plt.cla()
