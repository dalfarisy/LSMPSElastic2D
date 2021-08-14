
"""
LSMPSPlotter v1.0
by Dew
"""

import matplotlib.pyplot as plt

def plotter(x,y,w,scattersize,title):
    fig,ax = plt.subplots()
    a = ax.scatter(x,y,c=w,s=scattersize)
    ax.axis('equal')
    plt.colorbar(a)
    ax.set_title(f'{title}')
    
def plotterrainbow(x,y,w,scattersize,title):
    fig,ax = plt.subplots()
    a = ax.scatter(x,y,c=w,s=scattersize,cmap='rainbow')
    ax.axis('equal')
    plt.colorbar(a)
    ax.set_title(f'{title}')
    
def plotterjet(x,y,w,scattersize,title):
    fig,ax = plt.subplots()
    a = ax.scatter(x,y,c=w,s=scattersize,cmap='jet')
    ax.axis('equal')
    plt.colorbar(a)
    ax.set_title(f'{title}')
    
def plottergray(x,y,w,scattersize,title):
    fig,ax = plt.subplots()
    a = ax.scatter(x,y,c=w,s=scattersize,cmap='gray')
    ax.axis('equal')
    plt.colorbar(a)
    ax.set_title(f'{title}')

def plotterviridis(x,y,w,scattersize,title):
    fig,ax = plt.subplots()
    a = ax.scatter(x,y,c=w,s=scattersize,cmap='viridis')
    ax.axis('equal')
    plt.colorbar(a)
    ax.set_title(f'{title}')
    
def plotgraph(x,y,title):
    fig,ax = plt.subplots()
    ax.plot(x,y)
    ax.set_title(f'{title}')

def boundaryplot(x,y,dispX,dispY,Fx,Fy,nx,ny,ScatterSize):
    fig1,ax1 = plt.subplots()
    ax1.scatter(x,y,s = ScatterSize)
    
    xdirichlet = []
    ydirichlet = []
    xneumann = []
    yneumann = []
    nxboundary = []
    nyboundary = []
    for i in range(len(x)):
        if Fx[i] != 'nan' and Fy[i] != 'nan':
            if Fx[i] != 0 or Fy[i] != 0:
                xneumann.append(x[i])
                yneumann.append(y[i])
                nxboundary.append(nx[i])
                nyboundary.append(ny[i])
        elif dispX[i] != 'nan' or dispY[i] != 'nan':
            xdirichlet.append(x[i])
            ydirichlet.append(y[i])
    
    ax1.scatter(xneumann,yneumann,s = ScatterSize, c = 'red')
    ax1.quiver(xneumann,yneumann,nxboundary,nyboundary)
    ax1.scatter(xdirichlet,ydirichlet,s = ScatterSize, c = 'yellow')
    ax1.axis('equal')
    
    xfree = []
    yfree = []
    nxfree = []
    nyfree = []
    for i in range(len(x)):
        if Fx[i] == 0 and Fy[i] == 0:
                xfree.append(x[i])
                yfree.append(y[i])
                nxfree.append(nx[i])
                nyfree.append(ny[i])
    
    ax1.quiver(xfree,yfree,nxfree,nyfree)
