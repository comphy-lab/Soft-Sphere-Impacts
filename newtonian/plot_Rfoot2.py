#import libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
FigureWidth = 3.375
factor = 2
# Set global parameters
params = {
          'lines.linewidth': 3,
          'axes.labelsize': int(11*factor),
          'legend.fontsize': int(7*factor),
          'xtick.labelsize': int(10*factor),
          'ytick.labelsize': int(10*factor),
          'text.usetex': True,
          'font.family': 'serif'}
plt.rcParams.update(params)
plt.rcParams['text.latex.preamble'] = r"\usepackage{amsmath}"

We = float(sys.argv[1])
#Ec = float(sys.argv[2])

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
t = df['t']
Rfoot = df['Rfoot']
#Fit the curve
x = np.arange(0.0,0.75*tmax,0.001)

#calculate best fit - F1 vs We
from scipy.optimize import curve_fit
#objective function
def func(t, a, b):
    return np.sqrt(a*np.maximum((t-b),0))

popt, pcov = curve_fit(func, t, Rfoot, p0=[0,0])
print(popt, pcov)
a = popt[0]
b = popt[1]
Rfoot_fit = func(x,*popt)

#PLOT
plt.rcParams["figure.figsize"] = [6,6]
plt.rcParams['figure.autolayout'] = True
plt.xscale("linear")
plt.yscale("linear")
plt.xlim(0, 0.75*tmax)
plt.ylim(0,1)

#
h1, = plt.plot(t-b,Rfoot, 'bo', label = "Simulation data")
h2, = plt.plot(x-b, Rfoot_fit, 'r-', label = r'$\sqrt{aV_0R_0(t-t_0)}$')


#plt.title("Drop spread with time", fontsize = 20)
plt.xlabel(r"$(t-t_0)/t_c$", fontsize = 20)
plt.ylabel("$R_{foot}$", fontsize = 20)
first_legend = plt.legend(handles=[h1, h2], loc='best', fontsize = 14)

#plt.savefig('Rfootvstime.png')
plt.savefig('Rfootvstime.jpg', bbox_inches='tight')
plt.savefig('Rfootvstime.pdf', bbox_inches='tight')
#
#Get and write coeffs a and b(t_0) and We
f = open("../Rfoot_coeff_2.dat", "a")
f.write(str(We)+" "+str(a)+" "+str(b)+'\n')
#f.write(str(We)+" "+str(a)+" "+str(b)+" " + str(Ec)+'\n')
f.close()