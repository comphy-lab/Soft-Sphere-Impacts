import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 
import os
import sys

We = int(sys.argv[1]) #take in We as argument
El = float(sys.argv[2]) #take in El as an arument
filename = "%3.2f_El_pforce.dat" % El
df = pd.read_csv(filename, delimiter=' ') #read file

pforce = np.array(df['pforce'])
t = np.array(df['t'])

#Plot
plt.plot(t, pforce, 'b')
plt.title("Force vs time")
plt.xlabel("time")
plt.ylabel("Force")

plt.savefig('pforce.png')
#plt.show() ##debugging

#Get and write max pforce and We
f_1 = max(pforce)
f = open("../force_max.dat", "a")
f.write(str(We)+" "+str(f_1)+" "+ str(El)+'\n')
f.close()
