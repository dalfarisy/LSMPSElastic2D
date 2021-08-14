
"""
LSMPSGeom v1.0
by Dew
"""

from plotter import boundaryplot
import numpy as np

E = 200e6
v = 0.3
rho = 70
go = 0.0
a = 300
Nx = 300
Ny = 150
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
    nxleft.append(-1.0)
    nyleft.append(0.0)
    dispXleft.append('nan')
    dispYleft.append('nan')
    forceXleft.append(0.0)
    forceYleft.append(0.0)
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
    forceXright.append(0.0)
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
    dispXbottom.append(0.0)
    dispYbottom.append(0.0)
    forceXbottom.append('nan')
    forceYbottom.append('nan')
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
    dispXtop.append(0.0)
    dispYtop.append(0.02)
    forceXtop.append('nan')
    forceYtop.append('nan')#100000)
    bxtop.append(0.0)
    bytop.append(-rho*go)

#%% Inner Circle Boundary 1
    
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
ymid = Ny/2*h
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

#%% Inner Circle Boundary 2
    
xuppercircle2 = []
yuppercircle2 = []
nxuppercircle2 = []
nyuppercircle2 = []
dispXuppercircle2 = []
dispYuppercircle2 = []
forceXuppercircle2 = []
forceYuppercircle2 = []
bxuppercircle2 = []
byuppercircle2 = []

xlowercircle2 = []
ylowercircle2 = []
nxlowercircle2 = []
nylowercircle2 = []
dispXlowercircle2 = []
dispYlowercircle2 = []
forceXlowercircle2 = []
forceYlowercircle2 = []
bxlowercircle2 = []
bylowercircle2 = []

xmid2 = 0.8*Nx/4*h
ymid2 = 1.7*Ny/5*h
radius2 = Ny/4*h

for i in range(int(Nx/2)):
    xuppercircle2.append(radius2*np.cos((i)*np.pi/(Nx/2-1)) + xmid2)
    yuppercircle2.append(radius2*np.sin((i)*np.pi/(Nx/2-1)) + ymid2)
    nxtemp = xmid2 - xuppercircle2[i]
    nytemp = ymid2 - yuppercircle2[i]
    nxtemp1 = float(nxtemp/np.sqrt(nxtemp**(2) + nytemp**(2)))
    nytemp1 = float(nytemp/np.sqrt(nxtemp**(2) + nytemp**(2)))
    nxuppercircle2.append(nxtemp1)
    nyuppercircle2.append(nytemp1)
    dispXuppercircle2.append('nan')
    dispYuppercircle2.append('nan')
    forceXuppercircle2.append(0.0)
    forceYuppercircle2.append(0.0)
    bxuppercircle2.append(0.0)
    byuppercircle2.append(-rho*go)
    
for i in range(1,int(Nx/2)-1):
    xlowercircle2.append(radius2*np.cos((i)*np.pi/(Nx/2-1)) + xmid2)
    ylowercircle2.append(-radius2*np.sin((i)*np.pi/(Nx/2-1)) + ymid2)
    nxtemp = xmid2 - xlowercircle2[i-1]
    nytemp = ymid2 - ylowercircle2[i-1]
    nxtemp1 = float(nxtemp/np.sqrt(nxtemp**(2) + nytemp**(2)))
    nytemp1 = float(nytemp/np.sqrt(nxtemp**(2) + nytemp**(2)))
    nxlowercircle2.append(nxtemp1)
    nylowercircle2.append(nytemp1)
    dispXlowercircle2.append('nan')
    dispYlowercircle2.append('nan')
    forceXlowercircle2.append(0.0)
    forceYlowercircle2.append(0.0)
    bxlowercircle2.append(0.0)
    bylowercircle2.append(-rho*go)

#%% Inner Circle Boundary 3

xuppercircle3 = []
yuppercircle3 = []
nxuppercircle3 = []
nyuppercircle3 = []
dispXuppercircle3 = []
dispYuppercircle3 = []
forceXuppercircle3 = []
forceYuppercircle3 = []
bxuppercircle3 = []
byuppercircle3 = []

xlowercircle3 = []
ylowercircle3 = []
nxlowercircle3 = []
nylowercircle3 = []
dispXlowercircle3 = []
dispYlowercircle3 = []
forceXlowercircle3 = []
forceYlowercircle3 = []
bxlowercircle3 = []
bylowercircle3 = []

xmid3 = 3.2*Nx/4*h
ymid3 = 3*Ny/5*h
radius3 = Ny/4*h

for i in range(int(Nx/2)):
    xuppercircle3.append(radius3*np.cos((i)*np.pi/(Nx/2-1)) + xmid3)
    yuppercircle3.append(radius3*np.sin((i)*np.pi/(Nx/2-1)) + ymid3)
    nxtemp = xmid3 - xuppercircle3[i]
    nytemp = ymid3 - yuppercircle3[i]
    nxtemp1 = float(nxtemp/np.sqrt(nxtemp**(2) + nytemp**(2)))
    nytemp1 = float(nytemp/np.sqrt(nxtemp**(2) + nytemp**(2)))
    nxuppercircle3.append(nxtemp1)
    nyuppercircle3.append(nytemp1)
    dispXuppercircle3.append('nan')
    dispYuppercircle3.append('nan')
    forceXuppercircle3.append(0.0)
    forceYuppercircle3.append(0.0)
    bxuppercircle3.append(0.0)
    byuppercircle3.append(-rho*go)
    
for i in range(1,int(Nx/2)-1):
    xlowercircle3.append(radius3*np.cos((i)*np.pi/(Nx/2-1)) + xmid3)
    ylowercircle3.append(-radius3*np.sin((i)*np.pi/(Nx/2-1)) + ymid3)
    nxtemp = xmid3 - xlowercircle3[i-1]
    nytemp = ymid3 - ylowercircle3[i-1]
    nxtemp1 = float(nxtemp/np.sqrt(nxtemp**(2) + nytemp**(2)))
    nytemp1 = float(nytemp/np.sqrt(nxtemp**(2) + nytemp**(2)))
    nxlowercircle3.append(nxtemp1)
    nylowercircle3.append(nytemp1)
    dispXlowercircle3.append('nan')
    dispYlowercircle3.append('nan')
    forceXlowercircle3.append(0.0)
    forceYlowercircle3.append(0.0)
    bxlowercircle3.append(0.0)
    bylowercircle3.append(-rho*go)

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
        if (xtemp - xmid)**(2) + (ytemp - ymid)**(2) > radius**(2) and (xtemp - xmid2)**(2) + (ytemp - ymid2)**(2) > radius2**(2) and (xtemp - xmid3)**(2) + (ytemp - ymid3)**(2) > radius3**(2):
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
    
x = xtop + xbottom + xright + xleft + xinner + xuppercircle + xlowercircle + xuppercircle2 + xlowercircle2 + xuppercircle3 + xlowercircle3
y = ytop + ybottom + yright + yleft + yinner + yuppercircle + ylowercircle + yuppercircle2 + ylowercircle2 + yuppercircle3 + ylowercircle3
nx = nxtop + nxbottom + nxright + nxleft + nxinner + nxuppercircle + nxlowercircle + nxuppercircle2 + nxlowercircle2 + nxuppercircle3 + nxlowercircle3
ny = nytop + nybottom + nyright + nyleft + nyinner + nyuppercircle + nylowercircle + nyuppercircle2 + nylowercircle2 + nyuppercircle3 + nylowercircle3
dispX = dispXtop + dispXbottom + dispXright + dispXleft + dispXinner + dispXuppercircle + dispXlowercircle + dispXuppercircle2 + dispXlowercircle2 + dispXuppercircle3 + dispXlowercircle3
dispY = dispYtop + dispYbottom + dispYright + dispYleft + dispYinner + dispYuppercircle + dispYlowercircle + dispYuppercircle2 + dispYlowercircle2 + dispYuppercircle3 + dispYlowercircle3
forceX = forceXtop + forceXbottom + forceXright + forceXleft + forceXinner + forceXuppercircle + forceXlowercircle + forceXuppercircle2 + forceXlowercircle2 + forceXuppercircle3 + forceXlowercircle3
forceY = forceYtop + forceYbottom + forceYright + forceYleft + forceYinner + forceYuppercircle + forceYlowercircle + forceYuppercircle2 + forceYlowercircle2 + forceYuppercircle3 + forceYlowercircle3
bx = bxtop + bxbottom + bxright + bxleft + bxinner + bxuppercircle + bxlowercircle + bxuppercircle2 + bxlowercircle2 + bxuppercircle3 + bxlowercircle3
by = bytop + bybottom + byright + byleft + byinner + byuppercircle + bylowercircle + byuppercircle2 + bylowercircle2 + byuppercircle3 + bylowercircle3
boundaryplot(x,y,dispX,dispY,forceX,forceY,nx,ny,10)

print(f"Particle Numbers: {len(x)}")

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
