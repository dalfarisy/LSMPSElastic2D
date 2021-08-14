"""
LSMPS2D Viewer v1.0
by Dew
"""
import array as arr
import matplotlib.pyplot as plt
import numpy as np
from plotter import plotterrainbow
from plotter import plottergray
from plotter import plotterviridis
from plotter import plotterjet
from neighbourfind import sorting

plt.close('all')

#%% Opening the result file

file = open("result.txt","r")
readdata = file.readlines()
file.close()

h = float(readdata[0])

x = arr.array('d',[])
y = arr.array('d',[])
vonmises = arr.array('d',[])
sigmax = arr.array('d',[])
sigmay = arr.array('d',[])
sigmaxy = arr.array('d',[])
Ux = arr.array('d',[])
Uy = arr.array('d',[])

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
maxnumx = (max(x) - min(x))/h
maxnumy = (max(y) - min(y))/h
maxnum = max([maxnumx,maxnumy])

#%% Color Scheme Setting

scattersize = 250**(2)/(maxnum)**(2) + 1
# The sizing of each node in the scatter plot

dispscale = 100
# The displacement scale for plotting

contourscheme = 3
# 1 for rainbow color scheme
# 2 for viridis color scheme
# 3 for jet color scheme
# 4 for gray color scheme

convertunits = 2
# 1 to keep current units
# 2 to convert unit from mm to m

#%% Calculating Final Position and Plotting the result

Xnew = arr.array('d',[])
Ynew = arr.array('d',[])

for i in range(len(x)):
    Xnew.append(x[i] + Ux[i]*dispscale)
    Ynew.append(y[i] + Uy[i]*dispscale)
    
if convertunits == 2:
    # Converting to SI units (m and Pa)
    for i in range(len(Xnew)):
        Xnew[i] = Xnew[i]/1000
        Ynew[i] = Ynew[i]/1000
        sigmax[i] = sigmax[i]*1000
        sigmay[i] = sigmay[i]*1000
        sigmaxy[i] = sigmaxy[i]*1000
        vonmises[i] = vonmises[i]*1000
        Ux[i] = Ux[i]/1000
        Uy[i] = Uy[i]/1000
    for i in range(len(x)):
        x[i] = x[i]/1000
        y[i] = y[i]/1000

if contourscheme == 1:
    plotterrainbow(Xnew,Ynew,vonmises,scattersize,r'$\sigma_{VonMises},\delta_{Scale}:$'+f'{dispscale}')
    plotterrainbow(Xnew,Ynew,sigmax,scattersize,r'$\sigma_{xx},\delta_{Scale}:$'+f'{dispscale}')
    plotterrainbow(Xnew,Ynew,sigmay,scattersize,r'$\sigma_{yy},\delta_{Scale}:$'+f'{dispscale}')
    plotterrainbow(Xnew,Ynew,sigmaxy,scattersize,r'$\sigma_{xy},\delta_{Scale}:$'+f'{dispscale}')
elif contourscheme == 2:
    plotterviridis(Xnew,Ynew,vonmises,scattersize,r'$\sigma_{VonMises},\delta_{Scale}:$'+f'{dispscale}')
    plotterviridis(Xnew,Ynew,sigmax,scattersize,r'$\sigma_{xx},\delta_{Scale}:$'+f'{dispscale}')
    plotterviridis(Xnew,Ynew,sigmay,scattersize,r'$\sigma_{yy},\delta_{Scale}:$'+f'{dispscale}')
    plotterviridis(Xnew,Ynew,sigmaxy,scattersize,r'$\sigma_{xy},\delta_{Scale}:$'+f'{dispscale}')
elif contourscheme == 3:
    plotterjet(Xnew,Ynew,vonmises,scattersize,r'$\sigma_{VonMises},\delta_{Scale}:$'+f'{dispscale}')
    plotterjet(Xnew,Ynew,sigmax,scattersize,r'$\sigma_{xx},\delta_{Scale}:$'+f'{dispscale}')
    plotterjet(Xnew,Ynew,sigmay,scattersize,r'$\sigma_{yy},\delta_{Scale}:$'+f'{dispscale}')
    plotterjet(Xnew,Ynew,sigmaxy,scattersize,r'$\sigma_{xy},\delta_{Scale}:$'+f'{dispscale}')
    plotterjet(Xnew,Ynew,Ux,scattersize,r'$Ux,\delta_{Scale}:$'+f'{dispscale}')
    plotterjet(Xnew,Ynew,Uy,scattersize,r'$Uy,\delta_{Scale}:$'+f'{dispscale}')
elif contourscheme == 4:
    plottergray(Xnew,Ynew,vonmises,scattersize,r'$\sigma_{VonMises},\delta_{Scale}:$'+f'{dispscale}')
    plottergray(Xnew,Ynew,sigmax,scattersize,r'$\sigma_{xx},\delta_{Scale}:$'+f'{dispscale}')
    plottergray(Xnew,Ynew,sigmay,scattersize,r'$\sigma_{yy},\delta_{Scale}:$'+f'{dispscale}')
    plottergray(Xnew,Ynew,sigmaxy,scattersize,r'$\sigma_{xy},\delta_{Scale}:$'+f'{dispscale}')

print(f"Max Positive x displacement: {max(Ux)}")
print(f"Max Negative x displacement: {min(Ux)}")
print(f"Max Positive Y displacement: {max(Uy)}")
print(f"Max Negative Y displacement: {min(Uy)}")

#%% Probing and Analytical Comparison

# Opening the Porbe Data File

file = open("geom_probe.txt","r")
readdata = file.readlines()
file.close()

probe = arr.array('I',[])

for i in range(len(readdata)):
    probe.append(int(readdata[i]))
    
for i in range(len(deleter)):
    probe.pop(deleter[i])
    for j in range(len(deleter)):
        deleter[j] = deleter[j] - 1
    
# Probe Plot Setting

xprobe = []
yprobe = []

for i in range(len(x)):
    if probe[i] == 1:
        xprobe.append(x[i])
        yprobe.append(Uy[i])
        
# Rearranging the probe array
xprobesorted = sorting(xprobe,xprobe)
yprobesorted = sorting(xprobe,yprobe)

fig,ax = plt.subplots()
line1, = ax.plot(xprobesorted,yprobesorted)
ax.set_title(f'Displacement Y')
line1.set_label('')

# Probe Plot Comparison

# Beam Cantilever One Side 
xcomp = []
ycomp = []
dy = []
L = 1
h = 0.05
E = 200e9
I = 1/12 *h**(3)

xcomp = xprobesorted
for i in range(len(xprobesorted)):
    #ycomp.append(-200000*h*xcomp[i]**(2)/(6*E*I) *(3*L-xcomp[i]))                                                      # Beam Cantilever on One Side
    ycomp.append(-20000*xcomp[i]**(2)/(24*E*I) *((xcomp[i]**2)-4*L*(xcomp[i])+6*L**2))                                                      # Beam Cantilever on One Side
    #ycomp.append((-200000*xcomp[i])/(24*E*I) *(L**(3) - 2*L*xcomp[i]**(2) + xcomp[i]**(3)))                            # Beam Pin on One Side
    #ycomp.append(200000*(-(L - xcomp[i])**(5) + 2*L**(2)*(L - xcomp[i])**(3) - L**(4)*(L - xcomp[i]))/(120*E*I*L))     # Beam Pin Cantilever
    dy.append(abs(ycomp[i] - yprobesorted[i])**(2))
    
# Analytical Ploting
line2, = ax.plot(xcomp,ycomp)
line2.set_label("")
ax.legend((line1,line2),("LSMPS","Euler's Beam Theory"))
Error = abs(min(ycomp) - min(yprobesorted))/abs(min(ycomp))
RSME = np.sqrt(sum(dy)/len(xprobesorted))

yError = []
for i in range(len(xprobesorted)):
    yError.append(abs((yprobesorted[i]-ycomp[i]))/L)
    
fig1,ax1 = plt.subplots()
ax1.plot(xprobesorted,yError)
ax1.set_title(r'$|\frac{\Delta Y_{LSMPS} - \Delta Y_{Euler}}{L}|$')

print(f"Max Euler's Negative Y displacement: {(min(ycomp))}")
print(f"Max Error: {Error*100} %")
print(f"RSME: {RSME}")
