#LSMPSElasticSolver v.1.0
#Created by Devvskiii

import numpy as np
import array as arr
import scipy as sc
import scipy.sparse.linalg as SClinalg
import math
from plotter import plotterrainbow
from plotter import plottergray
from plotter import plotterviridis
from plotter import plotterjet
from neighbourfind import neighbourfind
from neighbourfind import neighbourfindlimited
import matplotlib.pyplot as plt
from LSMPSgeneral import LSMPS
from LSMPSgeneral import LSMPSconst
from LSMPSgeneral import LSMPSHrs
import datetime
plt.close("all")

start = datetime.datetime.now()

#%% Opening the geometry data file

file = open("Hasil Elips 2.txt","r")
readdata = file.readlines()
file.close()

x = arr.array('d',[])
y = arr.array('d',[])
z = arr.array('d',[])
vonmises = arr.array('d',[])
sigmax = arr.array('d',[])
sigmay = arr.array('d',[])
sigmaxy = arr.array('d',[])
Ux = arr.array('d',[])
Uy = arr.array('d',[])
bx = arr.array('d',[])
by = arr.array('d',[])
deleter = arr.array('d',[])
deleter = []
for i in range(1,len(readdata)):
    if readdata[i] != 'nan' + '\n':
        temp = readdata[i].split("\t")
        x.append(float(temp[0]))
        y.append(float(temp[1]))
        vonmises.append(float(temp[2]))
        sigmax.append(float(temp[3]))
        sigmay.append(float(temp[4]))
        sigmaxy.append(float(temp[5]))
        Ux.append(float(temp[6]))
        Uy.append(float(temp[7]))
    else:
        deleter.append(i-1)

del readdata,temp

boundaryplotmethod = 1
# 1 for plotting all the nodes
# 2 for plotting without dirichlet(displacement) boundary nodes
# 3 for plotting without both dirichlet and neumann boundary nodes

dispscale = 2000
# The displacement scale for plotting

contourscheme = 3
# 1 for rainbow color scheme
# 2 for viridis color scheme
# 3 for jet color scheme
# 4 for gray color scheme

datatype = 'd'
# 'f' for float
# 'd' for double

scattersize = 250**2/200**2+1
#%% Calculating the final position
    
Xnew = arr.array(datatype,[])
Ynew = arr.array(datatype,[])
        
if boundaryplotmethod == 1:
    for i in range(len(x)):
        Xnew.append(x[i] + Ux[i]*dispscale)
        Ynew.append(y[i] + Uy[i]*dispscale)

#%% Converting to SI units (m and Pa)
            
calcendtime = datetime.datetime.now()

#%% Plotting The Contour

if contourscheme == 1 and z==0:
    plotterrainbow(Xnew,Ynew,vonmises,scattersize,r'$\sigma_{VonMises},\delta_{Scale}:$'+f'{dispscale}')
    plotterrainbow(Xnew,Ynew,sigmax,scattersize,r'$\sigma_{xx},\delta_{Scale}:$'+f'{dispscale}')
    plotterrainbow(Xnew,Ynew,sigmay,scattersize,r'$\sigma_{yy},\delta_{Scale}:$'+f'{dispscale}')
    plotterrainbow(Xnew,Ynew,sigmaxy,scattersize,r'$\sigma_{xy},\delta_{Scale}:$'+f'{dispscale}')
elif contourscheme == 2:
    plotterviridis(Xnew,Ynew,vonmises,scattersize,r'$\sigma_{VonMises},\delta_{Scale}:$'+f'{dispscale}')
    plotterviridis(Xnew,Ynew,sigmax,scattersize,r'$\sigma_{xx},\delta_{Scale}:$'+f'{dispscale}')
    plotterviridis(Xnew,Ynew,sigmay,scattersize,r'$\sigma_{yy},\delta_{Scale}:$'+f'{dispscale}')
elif contourscheme == 3:
    plotterjet(Xnew,Ynew,vonmises,scattersize,r'$\sigma_{VonMises},\delta_{Scale}:$'+f'{dispscale}')
    plotterjet(Xnew,Ynew,sigmax,scattersize,r'$\sigma_{xx},\delta_{Scale}:$'+f'{dispscale}')
    plotterjet(Xnew,Ynew,sigmay,scattersize,r'$\sigma_{yy},\delta_{Scale}:$'+f'{dispscale}')
    plotterjet(Xnew,Ynew,sigmaxy,scattersize,r'$\sigma_{xy},\delta_{Scale}:$'+f'{dispscale}')
elif contourscheme == 4:
    plottergray(Xnew,Ynew,vonmises,scattersize,r'$\sigma_{VonMises},\delta_{Scale}:$'+f'{dispscale}')
    plottergray(Xnew,Ynew,sigmax,scattersize,r'$\sigma_{xx},\delta_{Scale}:$'+f'{dispscale}')
    plottergray(Xnew,Ynew,sigmay,scattersize,r'$\sigma_{yy},\delta_{Scale}:$'+f'{dispscale}')
    plottergray(Xnew,Ynew,sigmaxy,scattersize,r'$\sigma_{xy},\delta_{Scale}:$'+f'{dispscale}')

print(f"Max Positive X displacement: {max(Ux)}")
print(f"Max Negative X displacement: {min(Ux)}")
print(f"Max Positive Y displacement: {max(Uy)}")
print(f"Max Negative Y displacement: {min(Uy)}")
print(f"Max Positive Sigma Von Mises: {max(vonmises)}")
print(f"Max Negative Sigma Von Mises: {min(vonmises)}")
print(f"Max Positive Sigma XX: {max(sigmax)}")
print(f"Max Negative Sigma XX: {min(sigmax)}")
print(f"Max Positive Sigma XY: {max(sigmaxy)}")
print(f"Max Negative Sigma XY: {min(sigmaxy)}")
print(f"Max Positive Sigma YY: {max(sigmay)}")
print(f"Max Negative Sigma YY: {min(sigmay)}")





