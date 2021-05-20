"""
Course: EE2703-Applied Programming Lab
Name: Nihal Gajjala
Roll Number: EE19B044
Assignment 10
"""

import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import numpy as np
import sympy
import scipy.signal as sp
import pylab
from scipy.integrate import quad

# Function cos(wt+d)
def coswd(t,w,d):
    y=np.cos(w*t+d)         # Cosine Function
    return y

# Function cos(wt+d)*(random number)
def coswdr(t,w,d):
    N=np.size(t)                                    # Returns Number Of Elements In An Array
    y=np.cos(w*t+d)+0.1*np.random.randn(N)          # Cosine Function, Random Number
    return y

# Discrete Time Fourier Transform
def Dft(y):
    N=np.size(y)                                    # Returns Number Of Elements In An Array
    y[0]=0
    y=np.fft.fftshift(y)                            # Shifting The Zero-Frequency Component To The Center Of The Spectrum
    Y=np.fft.fftshift(np.fft.fft(y))/N              # Computing The One-Dimensional Discrete Fourier Transform
    return Y

# Function To Find wd
def findwd(Y,w):
    y=np.abs(Y)                     # Returns Absolute Value
    m=np.max(y)                     # Returns Maximum Value
    mval=np.where(y==m)             # Finding Elements From x Or y Depending On The Condition
    wdi=w[mval]
    di=np.angle(Y[mval])            # Returns The Angle Of The Complex Argument
    wd=np.abs(wdi[0])               # Returns Absolute Value
    d=np.abs(di[0])                 # Returns Absolute Value
    return wd,d

# Function To Find fmax
def fm(t1,t2):
    dt=t1-t2
    fmax=1/dt
    return fmax

# Function To Calculate Data For 3-D Plot
def freqtime(t):
    ncol=len(t)                                 # Returns Length Of Array
    Z=np.zeros((64,ncol),dtype=complex)         # Creating An Array Filled With Zeros
    X=np.zeros((64,ncol),dtype=float)           # Creating An Array Filled With Zeros
    Y=np.zeros((64,ncol),dtype=float)           # Creating An Array Filled With Zeros
    # Creating Time Array For Calculating Frequency
    time=linexc(-4*np.pi/ncol,4*np.pi/ncol,64)
    # Assigning Values Of fmax
    fmax=fm(time[1],time[0])
    # Calculating The Value Of wn
    wn=wnd(64)
    # Calculating The Values Of Frequency
    w=linexc(-fmax*np.pi,fmax*np.pi,64)
    # Calculating Y Array Which Contains The Values Of Frequency
    for j in range(64):
            Y[j,:]=w[j]
    for i in range(ncol):
        # Creating Arrays For Time
        time=linexc(t[i]+(0.8-1)*(4*np.pi/ncol),t[i]+(0.8+1)*(4*np.pi/ncol),64)
        # Calculating The Function Values
        y1=x3(time)*wn
        # Calculating The DFTs
        Y1=Dft(y1)
        Z[:,i]=Y1
        X[:,i]=t[i]
    return Z,X,Y

# Function To Construct An Array Of Evenly Spaced Numbers Over A Specified Interval
def linexc(t1,t2,N):
    t=np.linspace(t1,t2,N+1)            # Constructing An Array Of Evenly Spaced Numbers Over A Specified Interval
    t=t[:-1]
    return t

# Function To Calculate Data For 3-D Plot
def trimx(x,y,z,wmin):
    row=len(x)                                      # Returns Length Of Array
    col=len(x[0])                                   # Returns Length Of Array
    row2=0
    for i in range(row):
        if(abs(x[i,0])<=wmin):
            row2+=1
    z1=np.zeros((row2,col),dtype=complex)           # Creating An Array Filled With Zeros
    x1=np.zeros((row2,col),dtype=float)             # Creating An Array Filled With Zeros
    y1=np.zeros((row2,col),dtype=float)             # Creating An Array Filled With Zeros
    j=0
    for i in range(row):
        if(abs(x[i,0])<=wmin):                      # Returns Absolute Value
            x1[j,:]=x[i,:]
            y1[j,:]=y[i,:]
            z1[j,:]=z[i,:]
            j+=1
    return(x1,y1,z1)

# Function To Find w
def wnd(N):
    a=0.54
    b=0.46
    n=np.arange(N)                                      # Returns Evenly Spaced Values Within A Given Interval
    w=np.fft.fftshift(a+b*np.cos(2*np.pi*n/(N-1)))      # Shifting The Zero-Frequency Component To The Center Of The Spectrum
    return(w)

# Functions To Plot Different Kinds Of Data
def pltdiff(t1,t2,t3,y1,y2,y3,title,xlabel,ylabel,test):
    pylab.figure()                              # Creating Figure
    if test=='l':
        pylab.plot(t1,y1,'b')                   # Plotting y vs x As Lines And Markers
        pylab.plot(t2,y2,'r')                   # Plotting y vs x As Lines And Markers
        pylab.plot(t3,y3,'r')                   # Plotting y vs x As Lines And Markers
    if test=='d':
        pylab.plot(t1,y1,'bo')                  # Plotting y vs x As Lines And Markers
        pylab.plot(t2,y2,'ro')                  # Plotting y vs x As Lines And Markers
        pylab.plot(t3,y3,'ro')                  # Plotting y vs x As Lines And Markers
    pylab.ylabel(ylabel)                        # Setting The Label For The y-axis
    pylab.xlabel(xlabel)                        # Setting The Label For The x-axis
    pylab.title(title)                          # Setting Title For The Graph
    pylab.grid(True)                            # Displaying The Grid
    pylab.show()                                # Displaying The Figure

def pltfunc(t,y,xlabel,ylabel,title):
    pylab.figure()                              # Creating Figure
    pylab.plot(t,y)                             # Plotting y vs x As Lines And Markers
    pylab.ylabel(ylabel)                        # Setting The Label For The y-axis
    pylab.xlabel(xlabel)                        # Setting The Label For The x-axis
    pylab.title(title)                          # Setting Title For The Graph
    pylab.grid(True)                            # Displaying The Grid
    pylab.show()                                # Displaying The Figure

def pltsame(t1,t2,t3,y,title,xlabel,ylabel,test):
    pylab.figure()                              # Creating Figure
    if test=='l':
        pylab.plot(t1,y,'b')                    # Plotting y vs x As Lines And Markers
        pylab.plot(t2,y,'r')                    # Plotting y vs x As Lines And Markers
        pylab.plot(t3,y,'r')                    # Plotting y vs x As Lines And Markers
    if test=='d':
        pylab.plot(t1,y,'bo')                   # Plotting y vs x As Lines And Markers
        pylab.plot(t2,y,'ro')                   # Plotting y vs x As Lines And Markers
        pylab.plot(t3,y,'ro')                   # Plotting y vs x As Lines And Markers
    pylab.xlabel(xlabel)                        # Setting The Label For The x-axis
    pylab.title(title)                          # Setting Title For The Graph
    pylab.grid(True)                            # Displaying The Grid
    pylab.show()                                # Displaying The Figure

def subp(w,Y,xlim1,xlim2,title,xlabel,ylabel1,ylabel2,phasemin):
    pylab.figure()                              # Creating Figure
    pylab.subplot(2,1,1)                        # Adding A Subplot To The Current Figure
    pylab.plot(w,abs(Y))                        # Plotting y vs x As Lines And Markers
    pylab.ylabel(ylabel1)                       # Setting The Label For The y-axis
    pylab.xlim([xlim1,xlim2])                   # Setting The x-limits Of The x-axis
    pylab.title(title)                          # Setting Title For The Graph
    pylab.grid(True)                            # Displaying The Grid
    pylab.subplot(2,1,2)                        # Adding A Subplot To The Current Figure
    ii=np.where(abs(Y)>phasemin)                # Finding Elements From x Or y Depending On The Condition
    pylab.plot(w[ii],np.angle(Y[ii]),'go')      # Plotting y vs x As Lines And Markers
    pylab.xlim([xlim1,xlim2])                   # Setting The x-limits Of The x-axis
    pylab.ylim([-4,4])                          # Setting The y-limits Of The y-axis
    pylab.xlabel(xlabel)                        # Setting The Label For The x-axis
    pylab.ylabel(ylabel2)                       # Setting The Label For The y-axis
    pylab.grid(True)                            # Displaying The Grid
    pylab.show()                                # Displaying The Figure

def calculated_subp(w,wc,dp,dc,xlabel,ylabel1,ylabel2,title):
    pylab.figure()                                      # Creating Figure
    pylab.subplot(2,1,1)                                # Adding A Subplot To The Current Figure
    pylab.plot(w,wc,lw=2)                               # Plotting y vs x As Lines And Markers
    pylab.ylabel(ylabel1)                               # Setting The Label For The y-axis
    pylab.xlim([0.25,1.75])                             # Setting The x-limits Of The x-axis
    pylab.title(title)                                  # Setting Title For The Graph
    pylab.grid(True)                                    # Displaying The Grid
    pylab.subplot(2,1,2)                                # Adding A Subplot To The Current Figure
    pylab.plot(w,dp,'ro',label='True value')            # Plotting y vs x As Lines And Markers
    pylab.plot(w,dc,'go',label='Calculated value')      # Plotting y vs x As Lines And Markers
    pylab.legend(loc='upper right')                     # Placing A Legend On The Top Right Corner Of The Graph
    pylab.xlim([0.25,1.75])                             # Setting The x-limits Of The x-axis
    pylab.ylim([-4,4])                                  # Setting The y-limits Of The y-axis
    pylab.xlabel(xlabel)                                # Setting The Label For The x-axis
    pylab.ylabel(ylabel2)                               # Setting The Label For The y-axis
    pylab.grid(True)                                    # Displaying The Grid
    pylab.show()                                        # Displaying The Figure

def error_subp(w,wc,dp,dc,xlabel,ylabel1,ylabel2,title):
    pylab.figure()                              # Creating Figure
    pylab.subplot(2,1,1)                        # Adding A Subplot To The Current Figure
    pylab.plot(w,abs(wc-w),lw=2)                # Plotting y vs x As Lines And Markers
    pylab.ylabel(ylabel1)                       # Setting The Label For The y-axis
    pylab.xlim([0.25,1.75])                     # Setting The x-limits Of The x-axis
    pylab.title(title)                          # Setting Title For The Graph
    pylab.grid(True)                            # Displaying The Grid
    pylab.subplot(2,1,2)                        # Adding A Subplot To The Current Figure
    pylab.plot(w,abs(dp-dc),'ro')               # Plotting y vs x As Lines And Markers
    pylab.xlim([0.25,1.75])                     # Setting The x-limits Of The x-axis
    pylab.xlabel(xlabel)                        # Setting The Label For The x-axis
    pylab.ylabel(ylabel2)                       # Setting The Label For The y-axis
    pylab.grid(True)                            # Displaying The Grid
    pylab.show()                                # Displaying The Figure

# Declaring Functions
x1 = lambda t : np.sin(np.sqrt(2)*t)
x2 = lambda t : (np.cos(0.86*t))**3
x3 = lambda t : np.cos((24*t)+(8*t*t/np.pi))

#y=sin(sqrt(2)*t)
# Creating Arrays For Time
t11=linexc(-np.pi,np.pi,64)
t12=linexc(-16*np.pi,16*np.pi,1024)
t13=linexc(-3*np.pi,-np.pi,64)
t14=linexc(np.pi,3*np.pi,64)
# Assigning Values Of fmax
fmax11=fm(t11[1],t11[0])
fmax12=fm(t12[1],t12[0])
# Calculating The Value Of wn
wnd11=wnd(64)
wnd12=wnd(1024)
# Calculating The Function Values
y11=x1(t11)
y12=x1(t11)*wnd11
y13=x1(t12)*wnd12
# Calculating The Values Of Frequency
w11=linexc(-fmax11*np.pi,fmax11*np.pi,64)
w12=linexc(-fmax12*np.pi,fmax12*np.pi,1024)
# Calculating The DFTs
Y11=Dft(y11)
Y12=Dft(y12)
Y13=Dft(y13)

#y=(cos(0.86*t))**3
#These values tell the no of divisions
k21=8
N21=64
N22=N21*k21
# Creating Arrays For Time
t21=linexc(-np.pi,np.pi,64)
t22=linexc(-8*np.pi,8*np.pi,512)
# Assigning Values Of fmax
fmax21=fm(t21[1],t21[0])
fmax22=fm(t22[1],t22[0])
# Calculating The Value Of wn
wnd21=wnd(64)
wnd22=wnd(512)
# Calculating The Function Values
y21=x2(t21)
y22=x2(t21)*wnd21
y23=x2(t22)*wnd22
# Calculating The Values Of Frequency
w21=linexc(-fmax21*np.pi,fmax21*np.pi,64)
w22=linexc(-fmax22*np.pi,fmax22*np.pi,512)
# Calculating The DFTs
Y21=Dft(y21)
Y22=Dft(y22)
Y23=Dft(y23)

#y=(cos(w*t+d))
w1=0.86
d1=np.pi/4
w2=np.linspace(0.5,1.5,20)
d2=np.pi/6
w2c=np.zeros(20)
d2c=np.zeros(20)
d2p=np.zeros(20)
w3=np.linspace(0.5,1.5,20)
d3=np.pi/6
w3c=np.zeros(20)
d3c=np.zeros(20)
d3p=np.zeros(20)
# Creating Arrays For Time
t31=linexc(-np.pi,np.pi,128)
t32=linexc(-4*np.pi,4*np.pi,512)
# Assigning Values Of fmax
fmax31=fm(t31[1],t31[0])
fmax32=fm(t32[1],t32[0])
# Calculating The Value Of wn
wnd31=wnd(128)
wnd32=wnd(512)
# Calculating The Function Values
y31=coswd(t31,w1,d1)
y32=coswd(t31,w1,d1)*wnd31
y33=coswd(t32,w1,d1)*wnd32
# Calculating The Values Of Frequency
w31=linexc(-fmax31*np.pi,fmax31*np.pi,128)
w32=linexc(-fmax32*np.pi,fmax32*np.pi,512)
# Calculating The DFTs
Y31=Dft(y31)
Y32=Dft(y32)
Y33=Dft(y33)

print("Exact values of wo and delta(in radians) are:")
print((w1,d1))                  # Printing Exact Values Of wo And delta(in radians)
print("Calculated values of wo and delta(in radians) are:")
print(findwd(Y33,w32))          # Printing Exact Values Of wo And delta(in radians)

for i in range(20):
    y34=coswd(t32,w2[i],d2)*wnd32
    Y34=Dft(y34)
    w2c[i],d2c[i]=findwd(Y34,w32)
    d2p[i]=d2

for i in range(20):
    y35=(coswdr(t32,w3[i],d3))*wnd32
    Y35=Dft(y35)
    w3c[i],d3c[i]=findwd(Y35,w32)
    d3p[i]=d3

#y=(cos((24*t)+(8*t*t/np.pi)))
# Creating Arrays For Time
t41=linexc(-np.pi,np.pi,1024)
t42=linexc(-np.pi,np.pi,16)
# Assigning Values Of fmax
fmax41=fm(t41[1],t41[0])
# Calculating The Value Of wn
wnd41=wnd(1024)
# Calculating The Function Values
y41=x3(t41)
y42=x3(t41)*wnd41
# Calculating The Values Of Frequency
w41=linexc(-fmax41*np.pi,fmax41*np.pi,1024)
# Calculating The DFTs
Y41=Dft(y41)
Y42=Dft(y42)
Y43,XX,YY=freqtime(t42)

#Figure 1
subp(w11,Y11,-10,10,r"Spectrum of $\sin\left(\sqrt{2}t\right)$",r"$\omega$",r"$|Y|$","Phase of Y",1e-2)
#Figure 2
pltdiff(t11,t13,t14,np.sin(np.sqrt(2)*t11),np.sin(np.sqrt(2)*t13),np.sin(np.sqrt(2)*t14),r"$\sin\left(\sqrt{2}t\right)$",r"$t$",r"$y$",'l')
#Figure 3
pltsame(t11,t13,t14,y11,r"$\sin\left(\sqrt{2}t\right)$ with $t$ wrapping every $2\pi$ ",r"$t$",r"$y$",'d')
#Figure 4
pltsame(t11,t13,t14,y12,r"$\sin\left(\sqrt{2}t\right)\times w(t)$ with $t$ wrapping every $2\pi$ ",r"$t$",r"$y$",'d')
#Figure 5
subp(w11,Y12,-8,8,r"Spectrum of $\sin\left(\sqrt{2}t\right)\times w(t)$",r"$\omega$",r"$|Y|$","Phase of Y",1e-2)
#Figure 6
subp(w12,Y13,-4,4,r"Spectrum of $\sin\left(\sqrt{2}t\right)\times w(t)$",r"$\omega$",r"$|Y|$","Phase of Y",1e-2)

#Figure 7
subp(w21,Y21,-10,10,r"Spectrum of $\cos^3\left(0.86t\right)$",r"$\omega$",r"$|Y|$","Phase of Y",1e-2)
#Figure 8
subp(w21,Y22,-8,8,r"Spectrum of $\cos^3\left(0.86t\right)\times w(t)$",r"$\omega$",r"$|Y|$","Phase of Y",1e-2)
#Figure 9
subp(w22,Y23,-4,4,r"Spectrum of $\cos^3\left(0.86t\right)\times w(t)$",r"$\omega$",r"$|Y|$","Phase of Y",1e-2)

#Figure 10
subp(w31,Y31,-10,10,r"Spectrum of $\cos\left(0.86t+\delta\right)$",r"$\omega$",r"$|Y|$","Phase of Y",1e-2)
#Figure 11
subp(w31,Y32,-8,8,r"Spectrum of $\cos\left(0.86t+\delta\right)\times w(t)$",r"$\omega$",r"$|Y|$","Phase of Y",1e-2)
#Figure 12
subp(w32,Y33,-4,4,r"Spectrum of $\cos\left(0.86t+\delta\right)\times w(t)$",r"$\omega$",r"$|Y|$","Phase of Y",1e-2)

#Figure 13
calculated_subp(w2,w2c,d2p,d2c,r"Exact value of w",r"Calculated value of w",r"Calculated value of d",r"Calculating values of w and d")
#Figure 14
error_subp(w2,w2c,d2p,d2c,r"Exact value of w",r"Error in value of w",r"Error in value of d",r"Calculating values of w and d")
#Figure 15
calculated_subp(w3,w3c,d3p,d3c,r"Exact value of w",r"Calculated value of w",r"Calculated value of d",r"Calculating values of w and d")
#Figure 16
error_subp(w3,w3c,d3p,d3c,r"Exact value of w",r"Error in value of w",r"Error in value of d",r"Calculating values of w and d")

#Figure 17
pltfunc(t41,y41,r"t",r"y",r"y vs t")
#Figure 18
subp(w41,Y42,-60,60,r"Spectrum of $\cos\left(16(1.5+t/2\pi)t\right)\times w(t)$",r"$\omega$",r"$|Y|$","Phase of Y",1e-2)
#Plotting Figure 19
Yp,Xp,Zp=trimx(YY,XX,Y43,60)
ax=p3.Axes3D(pylab.figure())
ax.set_xlim(-60,60)
pylab.title(r"3-D surface plot of DFT")
surface=ax.plot_surface(Yp,Xp,abs(Zp),rstride=1,cstride=1,cmap='winter')
pylab.xlabel(r"w")
pylab.ylabel(r"t")
pylab.show()