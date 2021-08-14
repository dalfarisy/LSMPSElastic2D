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

file = open("geom.txt","r")
readdata = file.readlines()
file.close()
E = float(readdata[0])
v = float(readdata[1])
h = float(readdata[2])

x = []
y = []
nx = []
ny = []
dispX = []
dispY = []
forceX = []
forceY = []
bx = []
by = []

for i in range(3,len(readdata)):
    temp = readdata[i].split("\t")
    x.append(float(temp[0]))
    y.append(float(temp[1]))
    if temp[2] != 'nan':
        nx.append(float(temp[2]))
    elif temp[2] == 'nan':
        nx.append('nan')
    if temp[3] != 'nan':
        ny.append(float(temp[3]))
    elif temp[3] == 'nan':
        ny.append('nan')
    if temp[4] != 'nan':
        dispX.append(float(temp[4]))
    elif temp[4] == 'nan':
        dispX.append('nan')
    if temp[5] != 'nan':
        dispY.append(float(temp[5]))
    elif temp[5] == 'nan':
        dispY.append('nan')
    if temp[6] != 'nan':
        forceX.append(float(temp[6]))
    elif temp[6] == 'nan':
        forceX.append('nan')
    if temp[7] != 'nan':
        forceY.append(float(temp[7]))
    elif temp[7] == 'nan':
        forceY.append('nan')
    if temp[8] != 'nan':
        bx.append(float(temp[8]))
    elif temp[8] == 'nan':
        bx.append('nan')
    if temp[9] != 'nan' + '\n':
        by.append(float(temp[9]))
    elif temp[9] == 'nan' + '\n':
        by.append('nan')

maxnumx = (max(x) - min(x))/h
maxnumy = (max(y) - min(y))/h
maxnum = max([maxnumx,maxnumy])

print(f"Number of Particles: {len(x)}")
    
#%% Solver Control

scattersize = 250**(2)/(maxnum)**(2)
# The sizing of each node in the scatter plot

neighboursearchmethod = 2
# 1 for constant cutoff radius
# 2 for constant number of neighbours

scaledcutoffradius = 6.5001
# The cutoff radius for neighbour search

targetneighbournumber = 30
# Number of neighbour (for neighboursearchmethod = 2 only)

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

convertunits = 2
# 1 to keep current units
# 2 to convert unit from mm to m

writeresult = 2
# 1 to not write the result file
# 2 to write the result file

#%% Material data calculation

lamda = v*E/((1+v)*(1-2*v))
Mu = E/(2*(1+v))

#%% Boundary Search

neighbourstarttime = datetime.datetime.now()

if neighboursearchmethod == 1:
    neighbour = neighbourfind(x,y,scaledcutoffradius*h)
elif neighboursearchmethod == 2:
    neighbour = neighbourfindlimited(x,y,scaledcutoffradius*h,targetneighbournumber)
else:
    print("Neighbour search method selection INVALID")

for i in range(len(x)):
    neighbour[i].append(i)
  
neighbourendtime = datetime.datetime.now()

#%% Solving for LSMPS coefficients

calcstarttime = datetime.datetime.now()
e = h
# P = [1,x,y,xy,x^2,y^2,x^2y,xy^2,x^2y^2]
LSMPSconstant = LSMPSconst(x,y,neighbour,e)
LSMPSHRS = LSMPSHrs(x,y,neighbour,e)
 
#%% Rearranging The Coefficients in Matrix Form

print('Rearranging the Matrix ...')

b = np.zeros(len(x)*2)
row = arr.array('I',[])
col = arr.array('I',[])
data = arr.array(datatype,[])

for i in range(len(x)):
    
    if dispX[i] != 'nan':
               
       #A[2*i][2*i] = 1.0
        row.append(2*i)
        col.append(2*i)
        data.append(1.0)
        
        b[2*i] = dispX[i]
        
    if dispY[i] != 'nan':
        
       #A[2*i+1][2*i+1] = 1.0
        row.append(2*i+1)
        col.append(2*i+1)
        data.append(1.0)
        
        b[2*i+1] = dispY[i]
        
    if forceX[i] != 'nan':
        
        temp1 = 0.0
        temp2 = 0.0
        
        for j in range(len(neighbour[i])):
            
            temp1 = temp1 + (lamda + 2*Mu)*nx[i]*LSMPSconstant[i][1][j]*LSMPSHRS[1] + Mu*ny[i]*LSMPSconstant[i][2][j]*LSMPSHRS[2]
            
           #A[2*i][2*i+1] = A[2*i][2*i+1] + lamda*nx[i]*LSMPSconstant[2][j]/e**(1) + Mu*ny[i]*LSMPSconstant[1][j]/e**(1)
            temp2 = temp2 + lamda*nx[i]*LSMPSconstant[i][2][j]*LSMPSHRS[2] + Mu*ny[i]*LSMPSconstant[i][1][j]*LSMPSHRS[1]
            
           #A[2*i][2*neighbour[i][j]] = (lamda + 2*Mu)*nx[i]*LSMPSconstant[1][j]/e**(1) + Mu*ny[i]*LSMPSconstant[2][j]/e**(1)
            row.append(2*i)
            col.append(2*neighbour[i][j])
            data.append((lamda + 2*Mu)*nx[i]*LSMPSconstant[i][1][j]*LSMPSHRS[1] + Mu*ny[i]*LSMPSconstant[i][2][j]*LSMPSHRS[2])
        
           #A[2*i][2*neighbour[i][j]+1] = lamda*nx[i]*LSMPSconstant[2][j]/e**(1) + Mu*ny[i]*LSMPSconstant[1][j]/e**(1)
            row.append(2*i)
            col.append(2*neighbour[i][j]+1)
            data.append(lamda*nx[i]*LSMPSconstant[i][2][j]*LSMPSHRS[2] + Mu*ny[i]*LSMPSconstant[i][1][j]*LSMPSHRS[1])
        
            
            if j == len(neighbour[i])-1:
                
                row.append(2*i)
                col.append(2*i)
                data.append(temp1)
                
                row.append(2*i)
                col.append(2*i+1)
                data.append(temp2)
            
        b[2*i] = forceX[i]
        
    if forceY[i] != 'nan':
        
        temp1 = 0.0
        temp2 = 0.0      
        
        for j in range(len(neighbour[i])):
            
           #A[2*i+1][2*i] =  A[2*i+1][2*i] + lamda*ny[i]*LSMPSconstant[1][j]/e**(1) + Mu*nx[i]*LSMPSconstant[2][j]/e**(1)
            #temp1 = temp1 + lamda*ny[i]*DCPSEEtaDx[i][j]/e**(1) + Mu*nx[i]*DCPSEEtaDy[i][j]/e**(1)
            temp1 = temp1 + lamda*ny[i]*LSMPSconstant[i][1][j]*LSMPSHRS[1] + Mu*nx[i]*LSMPSconstant[i][2][j]*LSMPSHRS[2]
            
           #A[2*i+1][2*i+1] = A[2*i+1][2*i+1] + (lamda + 2*Mu)*ny[i]*LSMPSconstant[2][j]/e**(1) + Mu*nx[i]*LSMPSconstant[1][j]/e**(1)
            #temp2 = temp2 + (lamda + 2*Mu)*ny[i]*DCPSEEtaDy[i][j]/e**(1) + Mu*nx[i]*DCPSEEtaDx[i][j]/e**(1)
            temp2 = temp2 + (lamda + 2*Mu)*ny[i]*LSMPSconstant[i][2][j]*LSMPSHRS[2] + Mu*nx[i]*LSMPSconstant[i][1][j]*LSMPSHRS[1]
            
           #A[2*i+1][2*neighbour[i][j]] = lamda*ny[i]*LSMPSconstant[i][1][j]/e**(1) + Mu*nx[i]*LSMPSconstant[i][2][j]/e**(1)
            row.append(2*i+1)
            col.append(2*neighbour[i][j])
            data.append(lamda*ny[i]*LSMPSconstant[i][1][j]*LSMPSHRS[1] + Mu*nx[i]*LSMPSconstant[i][2][j]*LSMPSHRS[2])
            
           #A[2*i+1][2*neighbour[i][j]+1] = (lamda + 2*Mu)*ny[i]*LSMPSconstant[2][j]/e**(1) + Mu*nx[i]*LSMPSconstant[1][j]/e**(1)
            row.append(2*i+1)
            col.append(2*neighbour[i][j]+1)
            data.append((lamda + 2*Mu)*ny[i]*LSMPSconstant[i][2][j]*LSMPSHRS[2] + Mu*nx[i]*LSMPSconstant[i][1][j]*LSMPSHRS[1])

            if j == len(neighbour[i])-1:
                                
                row.append(2*i+1)
                col.append(2*i)
                data.append(temp1)
                
                row.append(2*i+1)
                col.append(2*i+1)
                data.append(temp2)
                
        b[2*i+1] = forceY[i]
            
    if forceX[i] == 'nan' and forceY[i] == 'nan' and dispX[i] == 'nan' and dispY[i] == 'nan':
                
        temp1 = 0.0
        temp2 = 0.0
        temp3 = 0.0
        temp4 = 0.0
        
        for j in range(len(neighbour[i])):
                        
           #A[2*i][2*i] = A[2*i][2*i] + (lamda + 2*Mu)*LSMPSconstant[4][j]/e**(2) + Mu*DCPSEEtaDy2[i][j]/e**(2)
           #temp1 = temp1 + (lamda + 2*Mu)*DCPSEEtaDx2[i][j]/e**(2) + Mu*DCPSEEtaDy2[i][j]/e**(2)
            temp1 = temp1 + (lamda + 2*Mu)*LSMPSconstant[i][4][j]*LSMPSHRS[4] + Mu*LSMPSconstant[i][5][j]*LSMPSHRS[5]
            
           #A[2*i][2*i+1] = A[2*i][2*i+1] + (lamda + Mu)*LSMPSconstant[3][j]/e**(2)
            #temp2 = temp2 + (lamda + Mu)*DCPSEEtaDxDy[i][j]/e**(2)
            temp2 = temp2 + (lamda + Mu)*LSMPSconstant[i][3][j]*LSMPSHRS[3]
            
           #A[2*i][2*neighbour[i][j]] = (lamda + 2*Mu)*LSMPSconstant[i][4][j]/e**(2) + Mu*LSMPSconstant[i][5][j]/e**(2)
            row.append(2*i)
            col.append(2*neighbour[i][j])
            data.append((lamda + 2*Mu)*LSMPSconstant[i][4][j]*LSMPSHRS[4] + Mu*LSMPSconstant[i][5][j]*LSMPSHRS[5])
            
           #A[2*i][2*neighbour[i][j]+1] = (lamda + Mu)*LSMPSconstant[3][j]/e**(2)
            row.append(2*i)
            col.append(2*neighbour[i][j]+1)
            data.append((lamda + Mu)*LSMPSconstant[i][3][j]*LSMPSHRS[3])
            
            
           #A[2*i+1][2*i] = A[2*i+1][2*i] + (lamda + Mu)*LSMPSconstant[3][j]/e**(2)
            #temp3 = temp3 + (lamda + Mu)*DCPSEEtaDxDy[i][j]/e**(2)
            temp3 = temp3 + (lamda + Mu)*LSMPSconstant[i][3][j]*LSMPSHRS[3]
            
           #A[2*i+1][2*i+1] = A[2*i+1][2*i+1] + (lamda + 2*Mu)*LSMPSconstant[5][j]/e**(2) + Mu*LSMPSconstant[4][j]/e**(2)
            #temp4 = temp4 + (lamda + 2*Mu)*DCPSEEtaDy2[i][j]/e**(2) + Mu*DCPSEEtaDx2[i][j]/e**(2)
            temp4 = temp4 + (lamda + 2*Mu)*LSMPSconstant[i][5][j]*LSMPSHRS[5] + Mu*LSMPSconstant[i][4][j]*LSMPSHRS[4]
            
           #A[2*i+1][2*neighbour[i][j]] = (lamda + Mu)*LSMPSconstant[3][j]/e**(2)
            row.append(2*i+1)
            col.append(2*neighbour[i][j])
            data.append((lamda + Mu)*LSMPSconstant[i][3][j]*LSMPSHRS[3])
            
           #A[2*i+1][2*neighbour[i][j]+1] = (lamda + 2*Mu)*LSMPSconstant[5][j]/e**(2) + Mu*LSMPSconstant[4][j]/e**(2)
            row.append(2*i+1)
            col.append(2*neighbour[i][j]+1)
            data.append((lamda + 2*Mu)*LSMPSconstant[i][5][j]*LSMPSHRS[5] + Mu*LSMPSconstant[i][4][j]*LSMPSHRS[4])
            
            
            if j == len(neighbour[i])-1:
                
                row.append(2*i)
                col.append(2*i)
                data.append(temp1)
                
                row.append(2*i)
                col.append(2*i+1)
                data.append(temp2)
                
                row.append(2*i+1)
                col.append(2*i)
                data.append(temp3)
                
                row.append(2*i+1)
                col.append(2*i+1)
                data.append(temp4)
                
        b[2*i] = bx[i]
        b[2*i+1] = by[i]
            
#%% Calculating Displacement Vector and Stresses
print('Matrix Rearrangement Complete!')
A = sc.sparse.csr_matrix((data, (row, col)), shape=(len(x)*2, len(x)*2))
print('Solving...')
startsolve = datetime.datetime.now()
U = SClinalg.spsolve(A,b,use_umfpack=True)
endsolve = datetime.datetime.now()

print('Solving Complete!')
Ux = arr.array(datatype,[])
Uy = arr.array(datatype,[])
for i in range(len(U)):
    if i%2 == 0:
        Ux.append(U[i])
    if i%2 == 1:
        Uy.append(U[i])
    
dUx = np.zeros((len(x),9))
#dUxdY = np.zeros(len(x))
dUy = np.zeros((len(x),9))
#dUydY = np.zeros(len(x))

epsxx = arr.array(datatype,[])
epsyy = arr.array(datatype,[])
epsxy = arr.array(datatype,[])
sigmax = arr.array(datatype,[])
sigmay = arr.array(datatype,[])
sigmaxy = arr.array(datatype,[])
vonmises = arr.array(datatype,[])

dUx = LSMPS(x,y,Ux,neighbour,e)
#dUxdY = LSMPS(x,y,Ux,neighbour,e)
dUy = LSMPS(x,y,Uy,neighbour,e)
#dUydY = LSMPS(x,y,Uy,neighbour,e)

# =============================================================================
# print(f"Max dUxdX:{max(dUxdX)}")
# print(f"Min dUxdX:{min(dUxdX)}")
# print(f"Max dUxdY:{max(dUxdY)}")
# print(f"Min dUxdY:{min(dUxdY)}")
# print(f"Max dUydX:{max(dUydX)}")
# print(f"Min dUydX:{min(dUydX)}")
# print(f"Max dUydY:{max(dUydY)}")
# print(f"Min dUydY:{min(dUydY)}")
# =============================================================================

j = 0
if boundaryplotmethod == 1:
    for i in range(len(x)):
        epsxx.append(dUx[i][1])
        epsyy.append(dUy[i][2])
        epsxy.append(0.5*(dUx[i][2] + dUy[i][1]))
        
        sigmax.append(lamda*(epsxx[j] + epsyy[j]) + 2*Mu*(epsxx[j]))
        sigmay.append(lamda*(epsxx[j] + epsyy[j]) + 2*Mu*(epsyy[j]))
        sigmaxy.append(2*Mu*epsxy[j])
        vonmises.append(np.sqrt(sigmax[j]**(2) + sigmay[j]**(2) - sigmax[j]*sigmay[j] + 3*sigmaxy[j]**(2)))
        j = j + 1   

elif boundaryplotmethod == 2:
    for i in range(len(x)):
        if dispX[i] == 'nan' and dispY[i] == 'nan':
            epsxx.append(dUx[i][1])
            epsyy.append(dUy[i][2])
            epsxy.append(0.5*(dUx[i][2] + dUy[i][1]))
            
            sigmax.append(lamda*(epsxx[j] + epsyy[j]) + 2*Mu*(epsxx[j]))
            sigmay.append(lamda*(epsxx[j] + epsyy[j]) + 2*Mu*(epsyy[j]))
            sigmaxy.append(2*Mu*epsxy[j])
            vonmises.append(np.sqrt(sigmax[j]**(2) + sigmay[j]**(2) - sigmax[j]*sigmay[j] + 3*sigmaxy[j]**(2)))
            j = j + 1

elif boundaryplotmethod == 3: 
    for i in range(len(x)):
        if nx[i] == 'nan' and ny[i] == 'nan' and dispX[i] == 'nan' and dispY[i] == 'nan':
            epsxx.append(dUx[i][1])
            epsyy.append(dUy[i][2])
            epsxy.append(0.5*(dUx[i][2] + dUy[i][1]))
            
            sigmax.append(lamda*(epsxx[j] + epsyy[j]) + 2*Mu*(epsxx[j]))
            sigmay.append(lamda*(epsxx[j] + epsyy[j]) + 2*Mu*(epsyy[j]))
            sigmaxy.append(2*Mu*epsxy[j])
            vonmises.append(np.sqrt(sigmax[j]**(2) + sigmay[j]**(2) - sigmax[j]*sigmay[j] + 3*sigmaxy[j]**(2)))
            j = j + 1

else:
    print("Boundary plot method selection INVALID")
    
#%% Writing the result file

if writeresult == 2:
    
    file = open("result.txt","w")
    file.write(f"{h}")
    file.write("\n")
    j = 0
    
    if boundaryplotmethod == 1:
        for i in range(len(x)):
            file.write(f"{x[i]}" + "\t" + f"{y[i]}" + "\t" + f"{vonmises[j]}" + "\t" + f"{sigmax[j]}" + "\t" + f"{sigmay[j]}" + "\t" + f"{sigmaxy[j]}" + "\t" + f"{Ux[i]}" + "\t" + f"{Uy[i]}")
            file.write("\n")
            j = j + 1
                
    elif boundaryplotmethod == 2:
        for i in range(len(x)):
            if dispX[i] == 'nan' and dispY[i] == 'nan':
                file.write(f"{x[i]}" + "\t" + f"{y[i]}" + "\t" + f"{vonmises[j]}" + "\t" + f"{sigmax[j]}" + "\t" + f"{sigmay[j]}" + "\t" + f"{sigmaxy[j]}" + "\t" + f"{Ux[i]}" + "\t" + f"{Uy[i]}")
                j = j + 1
            else:
                file.write('nan')
            file.write("\n")
            
    elif boundaryplotmethod == 3:
        for i in range(len(x)):
            if nx[i] == 'nan' and ny[i] == 'nan' and dispX[i] == 'nan' and dispY[i] == 'nan':
                file.write(f"{x[i]}" + "\t" + f"{y[i]}" + "\t" + f"{vonmises[j]}" + "\t" + f"{sigmax[j]}" + "\t" + f"{sigmay[j]}" + "\t" + f"{sigmaxy[j]}" + "\t" + f"{Ux[i]}" + "\t" + f"{Uy[i]}")
                j = j + 1
            else:
                file.write('nan')
            file.write("\n")
    
    file.close()
    
#%% Calculating the final position
    
Xnew = arr.array(datatype,[])
Ynew = arr.array(datatype,[])
        
if boundaryplotmethod == 1:
    for i in range(len(x)):
        Xnew.append(x[i] + Ux[i]*dispscale)
        Ynew.append(y[i] + Uy[i]*dispscale)
elif boundaryplotmethod == 2:
    for i in range(len(x)):
        if dispX[i] == 'nan' and dispY[i] == 'nan':
            Xnew.append(x[i] + Ux[i]*dispscale)
            Ynew.append(y[i] + Uy[i]*dispscale)
elif boundaryplotmethod == 3:
    for i in range(len(x)):
        if nx[i] == 'nan' and ny[i] == 'nan' and dispX[i] == 'nan' and dispY[i] == 'nan':
            Xnew.append(x[i] + Ux[i]*dispscale)
            Ynew.append(y[i] + Uy[i]*dispscale)

#%% Converting to SI units (m and Pa)
            
if convertunits == 2:
    for i in range(len(Xnew)):
        Xnew[i] = Xnew[i]/1000
        Ynew[i] = Ynew[i]/1000
        sigmax[i] = sigmax[i]*1000
        sigmay[i] = sigmay[i]*1000
        sigmaxy[i] = sigmaxy[i]*1000
        vonmises[i] = vonmises[i]*1000
    for i in range(len(x)):
        x[i] = x[i]/1000
        y[i] = y[i]/1000
        Ux[i] = Ux[i]/1000
        Uy[i] = Uy[i]/1000

calcendtime = datetime.datetime.now()

#%% Plotting The Contour

if contourscheme == 1:
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

neighbourtime = neighbourendtime - neighbourstarttime
calculationtime = calcendtime - calcstarttime
solvetime = endsolve - startsolve
totaltime = calcendtime - start
print(f"Neighbour Search Elapsed Time: {int(neighbourtime.total_seconds()*1000)} ms")
print(f"Calculation Elapsed Time: {int(calculationtime.total_seconds()*1000)} ms")
print(f"Solving Elapsed Time: {int(solvetime.total_seconds()*1000)} ms")
print(f"Total Elapsed Time: {int(totaltime.total_seconds()*1000)} ms")
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





