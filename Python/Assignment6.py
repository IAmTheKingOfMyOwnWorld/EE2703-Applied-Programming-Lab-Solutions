"""
Course: EE2703-Applied Programming Lab
Name: Nihal Gajjala
Roll Number: EE19B044
Assignment 6
"""

import argparse
import numpy as np
import pandas
import pylab

"""
n: Length Of Tubelight
M: Average Number Of Electrons Generated Per Turn
nk: Total Number Of Turns To Stimulate
u0: Threshold Velocity For Ionization
p: Probability Of Ionization (Electron Is Travelling With A Velocity Greater Than The Threshold Velocity)
Msig: Stddev Of Number Of Electrons Generated Per Turn
"""
parse=argparse.ArgumentParser()
parse.add_argument('--n',default=100,type=int,help='Spatial Grid Size')
parse.add_argument('--M',default=5,type=int,help='Number Of Electrons Injected Per Turn')
parse.add_argument('--nk',default=500,type=int,help='Number Of Turns To Stimulate')
parse.add_argument('--u0',default=5,type=int,help='Threshold Velocity')
parse.add_argument('--p',default=0.25,type=float,help='Probability Of Ionization')
parse.add_argument('--Msig',default=0.2,type=float,help='Variance of probability distribution')
args=parse.parse_args()
n,M,nk,u0,p,Msig=args.n,args.M,args.nk,args.u0,args.p,args.Msig

# Electron Position
xx=np.zeros(n*M)
# Electron Velocity
u=np.zeros(n*M)
# Displacement In Current Turn
dx=np.zeros(n*M)
# Intensity Of Emitted Light
I=[]
# Electron Position
X=[]
# Electron Velocity
V=[]

for i in range(nk):
    # Electrons Ejected From Cathode
    ii=np.where(xx>0)[0]        # Returns Elements Where xx>0
    dx[ii]=u[ii]+0.5            # Updating Displacement
    xx[ii]=xx[ii]+dx[ii]        # Updating Position
    u[ii]=u[ii]+1               # Updating Velocity

    # Electrons Absorbed At Anode
    jj=np.where(xx>n)[0]        # Returns Elements Where xx>n
    xx[jj]=0                    # Setting Position To 0
    u[jj]=0                     # Setting Velocity To 0
    dx[jj]=0                    # Setting Displacement To 0

    # Electrons With Velocity Greater Than Threshold
    kk=np.where(u>=u0)[0]       # Returns Elements Where u>u0 (Threshold)
    ll=np.where(np.random.rand(len(kk))<=p)[0]
    kl=kk[ll]

    # Electrons After Collision
    xx[kl]=xx[kl]-dx[kl]*np.random.rand()           # Setting Position
    u[kl]=0                                         # Setting Velocity

    I.extend(xx[kl].tolist())           # Excited Electrons
    m=int(np.random.randn()*Msig+M)     # Injected Electrons
    empty=np.where(xx==0)[0]            # Empty Slots
    if len(empty)>=m:
        random_start=np.random.randint(len(empty))
        xx[empty[random_start:m+random_start]]=1
        u[empty[random_start-m:random_start]]=0
    else:
        xx[empty]=1
        u[empty]=0
    occupied=np.where(xx>0)[0]
    X.extend(xx[occupied].tolist())
    V.extend(u[occupied].tolist())

# Light Intensity Histogram
pylab.figure(0)                                                 
a,bins,c=pylab.hist(I,bins=np.arange(1,100),ec='black')         # Plotting Histogram
pylab.title("Light Intensity Histogram")                        # Setting Title Of The Graph
pylab.show()                                                    # Displaying The Figure

# Electron Density Histogram
pylab.figure(1)                                                 
pylab.hist(X,bins=np.arange(1,100),ec='black')                  # Plotting Histogram
pylab.title("Electron Density Histogram")                       # Setting Title Of The Graph
pylab.show()                                                    # Displaying The Figure

xpos=0.5*(bins[0:-1]+bins[1:])
d={'Position':xpos,'Count':a}
p=pandas.DataFrame(data=d)
pandas.set_option("display.max_rows", None, "display.max_columns", None)
print(p)

# Electron Phase Space
pylab.figure(2)                                                 
pylab.plot(xx,u,'or')                                           # Plotting y vs x As Lines And Markers
pylab.title("Electron Phase Space")                             # Setting Title Of The Graph
pylab.xlabel("Position")                                        # Setting The Label For The x-axis
pylab.ylabel("Velocity")                                        # Setting The Label For The y-axis
pylab.show()                                                    # Displaying The Figure