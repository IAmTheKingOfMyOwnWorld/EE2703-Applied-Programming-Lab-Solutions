"""
Course: EE2703-Applied Programming Lab
Name: Nihal Gajjala
Roll Number: EE19B044
Assignment 5
"""

import argparse
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import numpy as np
import pylab

parse=argparse.ArgumentParser()
parse.add_argument('--Nx',default=25,type=int,help='Size along the x-axis')
parse.add_argument('--Ny',default=25,type=int,help='Size along the y-axis')
parse.add_argument('--Radius',default=8,type=float,help='Radius of central lead')
parse.add_argument('--Niter',default=1500,type=int,help='Number of iterations to be performed')
args=parse.parse_args()
Nx,Ny,Radius,Niter=args.Nx,args.Ny,args.Radius,args.Niter

# Function To Update phi
def update(phi,oldpfi):
    phi[1:-1,1:-1]=0.25*(oldphi[0:-2,1:-1]+oldphi[2:,1:-1]+oldphi[1:-1,0:-2]+oldphi[1:-1,2:])
    return phi

# Function To Assert Boundaries
def boundary(phi,index):
    phi[1:-1,0]=phi[1:-1,1]
    phi[0,1:-1]=phi[1,1:-1]
    phi[1:-1,-1]=phi[1:-1,-2]
    phi[-1,1:-1]=0
    phi[index]=1.0
    return phi

# Function To Evaluate The Parameter Of An Exponent
def error(a,b):
    logb=np.log(b)                                          # Finding Logarithm
    vec=np.zeros((len(a),2))                                # Creating An Array Filled With Zeros
    vec[:,0]=a
    vec[:,1]=1
    B,logA=np.linalg.lstsq(vec,np.transpose(logb))[0]       # Finding The Least Square Solution
    return (np.exp(logA),B)                                 # Finding Exponent

# Function To Calculate The Exponent
def fit(x,A,B):
    return A*np.exp(B*x)

phi=np.zeros((Ny,Nx))                       # Creating An Array Filled With Zeros
x=np.linspace(-(Ny-1)/2,(Ny-1)/2,Ny)        # Constructing An Array Of Evenly Spaced Numbers Over A Specified Interval
y=np.linspace(-(Nx-1)/2,(Nx-1)/2,Nx)        # Constructing An Array Of Evenly Spaced Numbers Over A Specified Interval
Y,X=np.meshgrid(y,x)                        # Creating Coordinate Matrices From Coordinate Vectors
index=np.where((X**2+Y**2)<=Radius**2)      # Finding Elements From x Or y Depending On The Condition
phi[index]=1.0
errors=np.zeros(Niter)                      # Creating An Array Filled With Zeros

# Plotting The Contour Of phi
pylab.figure(1)
pylab.contour(Y,X,phi,cmap="winter")                        # Plotting Contours
pylab.xlabel(r'x$\rightarrow$',fontsize=15)                 # Setting The Label For The x-axis
pylab.ylabel(r'y$\rightarrow$',fontsize=15)                 # Setting The Label For The y-axis
pylab.title('Contour Of phi',fontsize=15)                   # Setting Title Of The Graph
pylab.show()                                                # Displaying The Figure

# Calculating The Potential
for i in range(Niter):
    oldphi=phi.copy()                               # Copying phi
    phi=update(phi,oldphi)                          # Updating phi
    phi=boundary(phi,index)                         # Asserting Boundaries
    errors[i]=(abs(phi-oldphi)).max()               # Errors

# Plotting Error vs Number of Iterations
pylab.figure(2)
pylab.plot(range(Niter),errors,'-r')                            # Plotting y vs x As Lines And Markers
pylab.xlabel(r'Niter$\rightarrow$',fontsize=15)                 # Setting The Label For The x-axis
pylab.ylabel(r'Error$\rightarrow$',fontsize=15)                 # Setting The Label For The y-axis
pylab.title('Error vs Number of Iterations',fontsize=15)        # Setting Title Of The Graph
pylab.show()                                                    # Displaying The Figure

pylab.figure(3)
pylab.semilogy(range(Niter)[::50],errors[::50],'or',markersize=2.5) # Making A Plot With Log Scaling On The y-axis
pylab.xlabel(r'Niter$\rightarrow$',fontsize=15)                     # Setting The Label For The x-axis
pylab.ylabel(r'Error$\rightarrow$',fontsize=15)                     # Setting The Label For The y-axis
pylab.title('Error vs Number of Iterations',fontsize=15)            # Setting Title Of The Graph
pylab.show()                                                        # Displaying The Figure

pylab.figure(4)
pylab.loglog(range(Niter)[::50],errors[::50],'or',markersize=2.5)   # Making A Plot With Log Scaling On Both The Axis
pylab.xlabel(r'Niter$\rightarrow$',fontsize=15)                     # Setting The Label For The x-axis
pylab.ylabel(r'Error$\rightarrow$',fontsize=15)                     # Setting The Label For The y-axis
pylab.title('Error vs Number of Iterations',fontsize=15)            # Setting Title Of The Graph
pylab.show()                                                        # Displaying The Figure

# Errors
A,B=error(range(Niter),errors)
A500,B500=error(range(Niter)[500:],errors[500:])
pylab.figure(5)
pylab.semilogy(range(Niter)[::50],errors[::50],'ro',label='errors',markersize=2.5)                      # Making A Plot With Log Scaling On The y-axis
pylab.semilogy(range(Niter)[::50],fit(range(Niter)[::50],A,B),'go',label='fit1',markersize=2.5)         # Making A Plot With Log Scaling On The y-axis
pylab.semilogy(range(Niter)[::50],fit(range(Niter)[::50],A500,B500),'bo',label='fit2',markersize=2.5)   # Making A Plot With Log Scaling On The y-axis
pylab.legend(loc='upper right')                                                                         # Placing A Legend On The Top Right Corner Of The Graph
pylab.xlabel(r'Niter$\rightarrow$',fontsize=15)                                                         # Setting The Label For The x-axis
pylab.ylabel(r'Error$\rightarrow$',fontsize=15)                                                         # Setting The Label For The y-axis
pylab.title('Semilog Plot Of Error vs Number of Iterations',fontsize=15)                                # Setting Title Of The Graph
pylab.show()                                                                                            # Displaying The Figure

# 3-D Plot
fig=pylab.figure(6)
ax=p3.Axes3D(fig)                                                                   # Plotting 3-D data
pylab.title('3-D Plot Of Potential',fontsize=15)                                    # Setting Title Of The Graph
surface=ax.plot_surface(Y,X,phi.T,rstride=1,cstride=1,cmap="winter")            # Creating A Surface Plot
pylab.xlabel(r'x$\rightarrow$',fontsize=15)                                         # Setting The Label For The x-axis
pylab.ylabel(r'y$\rightarrow$',fontsize=15)                                         # Setting The Label For The y-axis
ax.set_zlabel(r'$\phi\rightarrow$',fontsize=15)                                     # Setting The Label For The z-axis
pylab.show()                                                                        # Displaying The Figure

# Plotting The Contour Of Potential
pylab.figure(7)
pylab.contour(Y,X[::-1],phi,cmap="inferno")                 # Plotting Contours
pylab.title('Contour Of Potential',fontsize=15)             # Setting Title Of The Graph
pylab.xlabel(r'x$\rightarrow$',fontsize=15)                 # Setting The Label For The x-axis
pylab.ylabel(r'y$\rightarrow$',fontsize=15)                 # Setting The Label For The y-axis
pylab.plot(index[1]-(Nx-1)/2,index[0]-(Ny-1)/2,'ro')        # Plotting y vs x As Lines And Markers
pylab.show()                                                # Displaying The Figure

# Calculating The Values Of Current Density
Jx=np.zeros((Nx,Ny))                                                # Creating An Array Filled With Zeros
Jy=np.zeros((Nx,Ny))                                                # Creating An Array Filled With Zeros
Jx[:,1:-1]=0.5*(phi[:,0:-2]-phi[:,2:])
Jy[1:-1,:]=0.5*(phi[2:,:]-phi[0:-2,:])
fig,ax=plt.subplots()                                               # Creating A Figure And A Set Of Subplots
fig=ax.quiver(Y[2:-1],X[::-1][2:-1],Jx[2:-1],Jy[2:-1],scale=8)      # Plotting A 2-D Field Of Arrows
pylab.plot(index[1]-(Nx-1)/2,index[0]-(Ny-1)/2,'ro')                # Plotting y vs x As Lines And Markers
pylab.title('Contour Of Current Density',fontsize=15)               # Setting Title Of The Graph
pylab.xlabel(r'x$\rightarrow$',fontsize=15)                         # Setting The Label For The x-axis
pylab.ylabel(r'y$\rightarrow$',fontsize=15)                         # Setting The Label For The y-axis
pylab.show()                                                        # Displaying The Figure