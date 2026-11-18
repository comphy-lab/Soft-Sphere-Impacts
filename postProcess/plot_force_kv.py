import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 
import os
import sys

from snapshot_tables import force_summary_path, weber_label

We = weber_label(sys.argv[1])
Ec = sys.argv[2] #take in Ec as an arument
filename = f"{We}_We_pforce.dat"
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
f = force_summary_path("force_max", "bounce-kv").open("a")
f.write(str(We)+" "+str(f_1)+" "+ str(Ec)+'\n')
f.close()
