"""
LSMPSGeom v1.0
by Dew
"""
from plotter import boundaryplot
import numpy as np

E = 200e6
v = 0.3
rho = 0#70
go = 0.0
a = 600
Nx = 450
Ny = 225
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
for i in range(Ny):
    xleft.append(0.0)
    yleft.append(i*h)
    nxleft.append('nan')
    nyleft.append('nan')
    dispXleft.append(0.0)
    dispYleft.append(0.0)
    forceXleft.append('nan')
    forceYleft.append('nan')
    bxleft.append(0.0)
    byleft.append(-rho*go)
    
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

for i in range(Ny):
    xright.append((Nx-1)*h)
    yright.append(i*h)
    nxright.append(1.0)
    nyright.append(0.0)
    dispXright.append('nan')
    dispYright.append('nan')
    forceXright.append(10000)
    forceYright.append(0.0)
    bxright.append(0.0)
    byright.append(-rho*go)
    
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

for i in range(1,Nx-1):
    xtop.append(i*h)
    ytop.append((Ny-1)*h)
    nxtop.append(0.0)
    nytop.append(1.0)
    dispXtop.append('nan')
    dispYtop.append('nan')
    forceXtop.append(0.0)
    forceYtop.append(0.0)
    bxtop.append(0.0)
    bytop.append(-rho*go)

#%% Inner Circle Boundary
    
xuppercircle = []
yuppercircle = []
nxuppercircle = []
nyuppercircle = []
dispXuppercircle = []
dispYuppercircle = []
forceXuppercircle = []
forceYuppercircle = []
bxuppercircle = []
byuppercircle = []

xlowercircle = []
ylowercircle = []
nxlowercircle = []
nylowercircle = []
dispXlowercircle = []
dispYlowercircle = []
forceXlowercircle = []
forceYlowercircle = []
bxlowercircle = []
bylowercircle = []

xmid = Nx/2*h
ymid = Ny/3*h
radius = Ny/4*h

for i in range(int(Nx/2)):
    xuppercircle.append(radius*np.cos((i)*np.pi/(Nx/2-1)) + xmid)
    yuppercircle.append(radius*np.sin((i)*np.pi/(Nx/2-1)) + ymid)
    nxtemp = xmid - xuppercircle[i]
    nytemp = ymid - yuppercircle[i]
    nxtemp1 = float(nxtemp/np.sqrt(nxtemp**(2) + nytemp**(2)))
    nytemp1 = float(nytemp/np.sqrt(nxtemp**(2) + nytemp**(2)))
    nxuppercircle.append(nxtemp1)
    nyuppercircle.append(nytemp1)
    dispXuppercircle.append('nan')
    dispYuppercircle.append('nan')
    forceXuppercircle.append(0.0)
    forceYuppercircle.append(0.0)
    bxuppercircle.append(0.0)
    byuppercircle.append(-rho*go)
    
for i in range(1,int(Nx/2)-1):
    xlowercircle.append(radius*np.cos((i)*np.pi/(Nx/2-1)) + xmid)
    ylowercircle.append(-radius*np.sin((i)*np.pi/(Nx/2-1)) + ymid)
    nxtemp = xmid - xlowercircle[i-1]
    nytemp = ymid - ylowercircle[i-1]
    nxtemp1 = float(nxtemp/np.sqrt(nxtemp**(2) + nytemp**(2)))
    nytemp1 = float(nytemp/np.sqrt(nxtemp**(2) + nytemp**(2)))
    nxlowercircle.append(nxtemp1)
    nylowercircle.append(nytemp1)
    dispXlowercircle.append('nan')
    dispYlowercircle.append('nan')
    forceXlowercircle.append(0.0)
    forceYlowercircle.append(0.0)
    bxlowercircle.append(0.0)
    bylowercircle.append(-rho*go)

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

for i in range(1,Nx-1):
    for j in range(1,Ny-1):
        xtemp = i*h
        ytemp = j*h
        if (xtemp - xmid)**(2) + (ytemp - ymid)**(2) > radius**(2):
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
    
#%% Combining the boundary and inner nodes list
    
x = xtop + xbottom + xright + xleft + xuppercircle + xlowercircle + xinner
y = ytop + ybottom + yright + yleft + yuppercircle + ylowercircle + yinner
nx = nxtop + nxbottom + nxright + nxleft + nxuppercircle + nxlowercircle + nxinner
ny = nytop + nybottom + nyright + nyleft + nyuppercircle + nylowercircle + nyinner
dispX = dispXtop + dispXbottom + dispXright + dispXleft + dispXuppercircle + dispXlowercircle + dispXinner
dispY = dispYtop + dispYbottom + dispYright + dispYleft + dispYuppercircle + dispYlowercircle + dispYinner
forceX = forceXtop + forceXbottom + forceXright + forceXleft + forceXuppercircle + forceXlowercircle + forceXinner
forceY = forceYtop + forceYbottom + forceYright + forceYleft + forceYuppercircle + forceYlowercircle + forceYinner
bx = bxtop + bxbottom + bxright + bxleft + bxuppercircle + bxlowercircle + bxinner
by = bytop + bybottom + byright + byleft + byuppercircle + bylowercircle + byinner 
    
boundaryplot(x,y,dispX,dispY,forceX,forceY,nx,ny,10)

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
