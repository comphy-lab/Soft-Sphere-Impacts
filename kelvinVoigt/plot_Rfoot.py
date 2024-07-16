#import libraries
import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 
import sys

We = float(sys.argv[1])
Ec = float(sys.argv[2])

#import the data log file
filename = "log"
df1 = pd.read_csv(filename, delimiter=' ', skiprows=1)

#initialise dVcm_dt and force
dVcm = np.zeros(len(df1))
F = np.zeros(len(df1))
df_f = df1.assign(dVcm_dt = dVcm, force = F)

#Get tmax
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
df_f.dropna(inplace=True)

#print(df_f) #debug
max_index = df_f['force_sma'].idxmax() #f_max
#print(max_index)
tmax = df_f['t'].iloc[max_index]
print("tmax = ", tmax)

#Read Rfoot
filename = "Rfoot.dat"
df = pd.read_csv(filename, delimiter=' ') #read file

df = df[df['t'] < 0.75*tmax]
t = df['t'].iloc[1:]
Rfoot = df['Rfoot'].iloc[1:]
#Fit the curve
x = np.arange(0.01,0.75*tmax,0.001)

#calculate best fit - F1 vs We
from scipy.optimize import curve_fit
#objective function
def func(t, a):
 return np.sqrt(a*(t-0.01))

popt, pcov = curve_fit(func, t, Rfoot)
print(popt, pcov)
a = popt[0]
Rfoot_fit = func(x,*popt)

#PLOT
plt.rcParams["figure.figsize"] = [6,6]
plt.rcParams['figure.autolayout'] = True
plt.rcParams['axes.linewidth'] = 1
plt.rc('xtick',labelsize=18)
plt.rc('ytick',labelsize=18)
#
plt.xscale("linear")
plt.yscale("linear")
#
h1, = plt.plot(t,Rfoot, 'bo', label = "Simulation data")
#plt.plot(t,Rfoot, 'k--')
#plt.plot(x,y, 'k-')
h2, = plt.plot(x, Rfoot_fit, 'r-', label = '$\sqrt{aVR_0t}$')

#plt.title("Drop spread with time", fontsize = 20)
plt.xlabel("$t/t_c$", fontsize = 20)
plt.ylabel("$R_{foot}$", fontsize = 20)
first_legend = plt.legend(handles=[h1, h2], loc='best', fontsize = 14)

#plt.savefig('Rfootvstime.png')
plt.savefig('Rfootvstime.jpg')
plt.savefig('Rfootvstime.pdf')

#plt.show()

#
#Get and write max force and We
f = open("../Rfoot_coeff.dat", "a")
#f.write(str(We)+" "+str(a)+ '\n')
f.write(str(We)+" "+str(a)+" "+str(Ec)+'\n')
f.close()
