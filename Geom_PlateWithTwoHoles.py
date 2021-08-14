
"""
LSMPSGeom v1.0
by Dew
"""

from plotter import boundaryplot
import numpy as np

E = 2e8
v = 0.3
rho = 0.0#70
go = 0.0
a = 600
Nx = 600
Ny = 300
h = a/(Nx-1)

htol = 0.3

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
for i in range(1,Ny-1):
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

for i in range(1,Ny-1):
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

for i in range(Nx):
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

for i in range(Nx):
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

xmid2 = 2*Nx/4*h
ymid2 = 2*Ny/3*h
radius2 = Ny/12*h

for i in range(int(Nx/6)):
    xuppercircle2.append(radius2*np.cos((i)*np.pi/(Nx/6-1)) + xmid2)
    yuppercircle2.append(radius2*np.sin((i)*np.pi/(Nx/6-1)) + ymid2)
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
    
for i in range(1,int(Nx/6)-1):
    xlowercircle2.append(radius2*np.cos((i)*np.pi/(Nx/6-1)) + xmid2)
    ylowercircle2.append(-radius2*np.sin((i)*np.pi/(Nx/6-1)) + ymid2)
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

xmid3 = 2*Nx/4*h
ymid3 = 1*Ny/3*h
radius3 = Ny/12*h

for i in range(int(Nx/6)):
    xuppercircle3.append(radius3*np.cos((i)*np.pi/(Nx/6-1)) + xmid3)
    yuppercircle3.append(radius3*np.sin((i)*np.pi/(Nx/6-1)) + ymid3)
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
    
for i in range(1,int(Nx/6)-1):
    xlowercircle3.append(radius3*np.cos((i)*np.pi/(Nx/6-1)) + xmid3)
    ylowercircle3.append(-radius3*np.sin((i)*np.pi/(Nx/6-1)) + ymid3)
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
        if (xtemp - xmid2)**(2) + (ytemp - ymid2)**(2) > (radius2+htol*h)**(2) and (xtemp - xmid3)**(2) + (ytemp - ymid3)**(2) > (radius3+htol*h)**(2):
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
    
x = xtop + xbottom + xright + xleft + xinner + xuppercircle2 + xlowercircle2 + xuppercircle3 + xlowercircle3
y = ytop + ybottom + yright + yleft + yinner + yuppercircle2 + ylowercircle2 + yuppercircle3 + ylowercircle3
nx = nxtop + nxbottom + nxright + nxleft + nxinner + nxuppercircle2 + nxlowercircle2 + nxuppercircle3 + nxlowercircle3
ny = nytop + nybottom + nyright + nyleft + nyinner + nyuppercircle2 + nylowercircle2 + nyuppercircle3 + nylowercircle3
dispX = dispXtop + dispXbottom + dispXright + dispXleft + dispXinner + dispXuppercircle2 + dispXlowercircle2 + dispXuppercircle3 + dispXlowercircle3
dispY = dispYtop + dispYbottom + dispYright + dispYleft + dispYinner + dispYuppercircle2 + dispYlowercircle2 + dispYuppercircle3 + dispYlowercircle3
forceX = forceXtop + forceXbottom + forceXright + forceXleft + forceXinner + forceXuppercircle2 + forceXlowercircle2 + forceXuppercircle3 + forceXlowercircle3
forceY = forceYtop + forceYbottom + forceYright + forceYleft + forceYinner + forceYuppercircle2 + forceYlowercircle2 + forceYuppercircle3 + forceYlowercircle3
bx = bxtop + bxbottom + bxright + bxleft + bxinner + bxuppercircle2 + bxlowercircle2 + bxuppercircle3 + bxlowercircle3
by = bytop + bybottom + byright + byleft + byinner + byuppercircle2 + bylowercircle2 + byuppercircle3 + bylowercircle3
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
