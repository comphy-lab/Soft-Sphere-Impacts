import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 
import sys


We = int(sys.argv[1]) #take in We as argument
Ec = sys.argv[2] #take in Ec as an arument
#import the data log file
filename = "log"
df = pd.read_csv(filename, delimiter=' ', skiprows=1)

#initialise dVcm_dt and force
dVcm = np.zeros(len(df))
F = np.zeros(len(df))
df_f = df.assign(dVcm_dt = dVcm, force = F)

#calculating force
for i in range(1,len(df_f)-1, 1):
    df_f.loc[i,"dVcm_dt"] = 0.5*((df_f.loc[i,"Vcm"]-df_f.loc[i-1, "Vcm"])/(df_f.loc[i,"t"]-df_f.loc[i-1, "t"]) + (df_f.loc[i+1,"Vcm"]-df_f.loc[i, "Vcm"])/(df_f.loc[i+1,"t"]-df_f.loc[i, "t"]))
    #print(i)
df_f["force"] = df_f["dVcm_dt"].multiply(4/3*np.pi) #F=4/3*pi*r^3*rho*dVcm_dt

#drop rows
df_f=df_f.drop(len(df_f)-1) 

for j in range(0,15,1):
    df_f=df_f.drop(j)

#Simple moving average
df_f['force_sma'] = df_f['force'].rolling(17).mean()
df_f.dropna(inplace=True) #Drop rows with NA
force_sma = np.array(df_f['force_sma'])
#print(df_f)# debugging
print("DONE! Plotting force vs time")

##PLOT
force = np.array(df_f['force'])
t = np.array(df_f['t'])
force_sma = np.array(df_f['force_sma'])

plt.title('Force vs Time')
plt.ylabel("Force")
plt.xlabel("time")
plt.plot(t, force, 'ro', ms = 0.5, label = "Actual data")
plt.plot(t, force, 'k--', linewidth = 0.2)
plt.plot(t[5:], force_sma[5:], 'b-', linewidth = 1, label= "Moving average")
#plt.legend(loc="best")
plt.savefig('force_vcm.png')
plt.savefig('force_vcm.pdf', bbox_inches='tight')

#Get and write max force and We
f_1 = max(force_sma)
f = open("../force_max_Vcm.dat", "a")
f.write(str(We)+" "+str(f_1)+" "+ str(Ec)+'\n')
f.close()