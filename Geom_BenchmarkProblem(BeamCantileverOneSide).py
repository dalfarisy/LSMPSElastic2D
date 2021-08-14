
"""
LSMPSGeom v1.0
by Dew
"""

from plotter import boundaryplot
from plotter import plotterviridis
import numpy as np

E = 2e8
v = 0.3
rho = 0#70
go = 0#-9.8

a = 1000
Nx = 400
Ny = 20
h = a/(Nx-1)

#%% Outer Rectangular Boundary

xleft = []
yleft = []
nxleft = []
nyleft = []
dispXleft = []
dispYleft = [] 
forceXleft = []
forceYleft = []
bxleft = []
byleft = []
probeleft = []

for i in range(Ny):
    xleft.append(0.0)
    yleft.append(i*h)
    nxleft.append(-1.0)
    nyleft.append(0.0)
    dispXleft.append(0.0)
    dispYleft.append(0.0)
    forceXleft.append('nan')
    forceYleft.append('nan')
    bxleft.append(0.0)
    byleft.append(-rho*go)
    if i == int(Ny/2):
        probeleft.append(1)
    else:
        probeleft.append(0)
    
xright = []
yright = []
nxright = []
nyright = []
dispXright = []
dispYright = [] 
forceXright = []
forceYright = []
bxright = []
byright = []
proberight = []

for i in range(Ny):
    xright.append((Nx-1)*h)
    yright.append(i*h)
    nxright.append(1.0)
    nyright.append(0.0)
    dispXright.append('nan')
    dispYright.append('nan')
    forceXright.append(0.0)
    forceYright.append(0.0)
    bxright.append(0.0)
    byright.append(-rho*go)
    if i == int(Ny/2):
        proberight.append(1)
    else:
        proberight.append(0)
    
xbottom = []
ybottom = []
nxbottom = []
nybottom = []
dispXbottom = []
dispYbottom = []
forceXbottom = []
forceYbottom = []
bxbottom = []
bybottom = []
probebottom = []

for i in range(1,Nx-1):
    xbottom.append(i*h)
    ybottom.append(0.0)
    nxbottom.append(0.0)
    nybottom.append(-1.0)
    dispXbottom.append('nan')
    dispYbottom.append('nan')
    forceXbottom.append(0.0)
    forceYbottom.append(0.0)
    bxbottom.append(0.0)
    bybottom.append(-rho*go)
    probebottom.append(0)
    
xtop = []
ytop = []
nxtop = []
nytop = []
dispXtop = []
dispYtop = []
forceXtop = []
forceYtop = []
bxtop = []
bytop = []
probetop = []

for i in range(1,Nx-1):
    xtop.append(i*h)
    ytop.append((Ny-1)*h)
    nxtop.append(0.0)
    nytop.append(1.0)
    dispXtop.append('nan')
    dispYtop.append('nan')
    forceXtop.append(0.0)
    forceYtop.append(-20.0)
    bxtop.append(0.0)
    bytop.append(-rho*go)
    probetop.append(0)

#%% Inner nodes
    
xinner = []
yinner = []
nxinner = []
nyinner = []
dispXinner = []
dispYinner = []
forceXinner = []
forceYinner = []
bxinner = []
byinner = []
probeinner = []

for i in range(1,Nx-1):
    for j in range(1,Ny-1):
        xtemp = i*h
        ytemp = j*h
        xinner.append(xtemp)
        yinner.append(ytemp)
        nxinner.append('nan')
        nyinner.append('nan')
        dispXinner.append('nan')
        dispYinner.append('nan')
        forceXinner.append('nan')
        forceYinner.append('nan')
        bxinner.append(0.0)
        byinner.append(-rho*go)
        if j == int(Ny/2):
            probeinner.append(1)
        else:
            probeinner.append(0)
        
    
#%% Combining the boundary and inner nodes list
    
x = xtop + xbottom + xright + xleft + xinner
y = ytop + ybottom + yright + yleft + yinner
nx = nxtop + nxbottom + nxright + nxleft + nxinner
ny = nytop + nybottom + nyright + nyleft + nyinner
dispX = dispXtop + dispXbottom + dispXright + dispXleft + dispXinner
dispY = dispYtop + dispYbottom + dispYright + dispYleft + dispYinner
forceX = forceXtop + forceXbottom + forceXright + forceXleft + forceXinner
forceY = forceYtop + forceYbottom + forceYright + forceYleft + forceYinner
bx = bxtop + bxbottom + bxright + bxleft + bxinner
by = bytop + bybottom + byright + byleft + byinner
probe = probetop + probebottom + proberight + probeleft + probeinner
boundaryplot(x,y,dispX,dispY,forceX,forceY,nx,ny,10)
plotterviridis(x,y,probe,10,'Probe Plot')

#print(f"{np.sum(probe)}")

#%% Writing the txt data

file = open("geom.txt","w")

file.write(f"{E}")
file.write("\n")
file.write(f"{v}")
file.write("\n")
file.write(f"{h}")
file.write("\n")

for i in range(len(x)):
    file.write(f"{x[i]}" + "\t" + f"{y[i]}" + "\t" + f"{nx[i]}" + "\t" + f"{ny[i]}" + "\t" + f"{dispX[i]}" + "\t" + f"{dispY[i]}" + "\t" + f"{forceX[i]}" + "\t" + f"{forceY[i]}" + "\t" + f"{bx[i]}" + "\t" + f"{by[i]}")
    file.write("\n")
    
file.close()

file = open("geom_probe.txt","w")

for i in range(len(x)):    
    file.write(f"{probe[i]}")
    file.write("\n")

file.close()




