"""
Course: EE2703-Applied Programming Lab
Name: Nihal Gajjala
Roll Number: EE19B044
Assignment 7
"""

import numpy as np
import scipy.signal as sp
import pylab

# Initializing Constants
decay1=0.5                                                              # Decay 1
decay2=0.05                                                             # Decay 2
freq=1.5                                                                # Frequency

# Laplace Transform Of Function f(t)
def Flti(freq,decay):
    poly1=np.poly1d([1,decay])                                          # 1-dimensional Polynomial
    poly2=np.poly1d([1,(2*decay),((decay*decay)+(freq*freq))])          # 1-dimensional Polynomial
    return poly1,poly2

# Function f(t)
def flti(freq,decay,time):
    x1=np.cos(freq*time)                                                # Cosine Function
    x2=np.multiply(np.exp(-decay*time),np.heaviside(time,0.5))          # Multiplication
    return np.multiply(x1,x2)

# Plotting Time Response Of Spring For Decay=0.5
time=np.linspace(0,50,1000)                                             # Constructing An Array Of Evenly Spaced Numbers Over A Specified Interval
numerator1,denominator1=Flti(freq,decay1)
denominator1=np.polymul([1,0,freq*freq],denominator1)                   # Multiplication Of Two Polynomials
Fspr1=sp.lti(numerator1,denominator1)                                   # Continuous Time LTI System
imp1=sp.impulse(Fspr1,T=time)                                           # Impulse Response Of Continuous Time System
pylab.figure(0)                                                         # Creating A New Figure
pylab.plot(imp1[0],imp1[1],'-b')                                        # Ploting y vs x As Lines And Markers
pylab.title('Time Response Of Spring For Decay=0.5')                    # Setting Title Of The Graph
pylab.xlabel(r'time')                                                   # Setting The Label For The x-axis
pylab.ylabel(r'x')                                                      # Setting The Label For The y-axis
pylab.grid()                                                            # Displaying The Grid
pylab.show()                                                            # Displaying The Figure

# Plotting Time Response Of Spring For Decay=0.05
numerator2,denominator2=Flti(freq,decay2)
denominator2=np.polymul([1,0,freq*freq],denominator2)                   # Multiplication Of Two Polynomials
Fspr2=sp.lti(numerator2,denominator2)                                   # Continuous Time LTI System
imp2=sp.impulse(Fspr2,T=time)                                           # Impulse Response Of Continuous Time System
pylab.figure(1)                                                         # Creating A New Figure
pylab.plot(imp2[0],imp2[1],'-b')                                        # Ploting y vs x As Lines And Markers
pylab.title('Time Response Of Spring For Decay=0.05')                   # Setting Title Of The Graph
pylab.xlabel(r'time')                                                   # Setting The Label For The x-axis
pylab.ylabel(r'x')                                                      # Setting The Label For The y-axis
pylab.grid()                                                            # Displaying The Grid
pylab.show()                                                            # Displaying The Figure

# Plotting Time Response Of Spring For Different Frequencies
j=1
color={1:'b',2:'g',3:'r',4:'y',5:'m',6:'c'}
frequency={1:'Frequency=0.140',2:'Frequency=0.145',3:'Frequency=0.150',4:'Frequency=0.155',5:'Frequency=0.160'}
pylab.figure(2)                                                         # Creating A New Figure
for i in np.arange(1.4,1.6,0.05):                                       # Evenly Spaced Numbers Over A Specified Interval
    Fspr=sp.lti([1],[1,0,freq*freq])                                    # Continuous Time LTI System
    f,x,_=sp.lsim(Fspr,flti(i,0.05,time),time)                          # Simulate Output Of A Continuous Time LTI system
    pylab.plot(time,x,color[j],label=frequency[j])                      # Ploting y vs x As Lines And Markers
    j+=1
pylab.title('Time Response Of Spring')                                  # Setting Title Of The Graph
pylab.xlabel(r'time')                                                   # Setting The Label For The x-axis
pylab.ylabel(r'x')                                                      # Setting The Label For The y-axis
pylab.legend(loc='upper left',fontsize=8)                               # Placing A Legend On The Top Left Corner Of The Graph
pylab.grid()                                                            # Displaying The Grid
pylab.show()                                                            # Displaying The Figure

# Plotting Time Response Of Coupled Spring
time=np.linspace(0,20,1000)                                             # Constructing An Array Of Evenly Spaced Numbers Over A Specified Interval
Fsprx=sp.lti(np.poly1d([1,0,2]),np.poly1d([1,0,3,0]))                   # Continuous Time LTI System Of Two 1-dimensional Polynomial
impx=sp.impulse(Fsprx,T=time)                                           # Impulse Response Of Continuous Time System
Fspry=sp.lti(np.poly1d([2]),np.poly1d([1,0,3,0]))                       # Continuous Time LTI System Of Two 1-dimensional Polynomial
impy=sp.impulse(Fspry,T=time)                                           # Impulse Response Of Continuous Time System
pylab.figure(3)                                                         # Creating A New Figure
pylab.plot(impx[0],impx[1],'-r',label='x')                              # Ploting y vs x As Lines And Markers
pylab.plot(impy[0],impy[1],'-g',label='y')                              # Ploting y vs x As Lines And Markers
pylab.title('Time Response Of Coupled Spring')                          # Setting Title Of The Graph
pylab.xlabel(r'time')                                                   # Setting The Label For The x-axis
pylab.ylabel(r'x')                                                      # Setting The Label For The y-axis
pylab.legend(loc='upper right')                                         # Placing A Legend On The Top Right Corner Of The Graph
pylab.grid()                                                            # Displaying The Grid
pylab.show()                                                            # Displaying The Figure

# Plotting Bode Graph
Hv=sp.lti(np.poly1d([1000000]),np.poly1d([0.000001,100,1000000]))       # Continuous Time LTI System Of Two 1-dimensional Polynomial
w,S,phi=Hv.bode()                                                       # Bode Graph
pylab.figure(4)                                                         # Creating A New Figure
pylab.subplot(2,1,1)                                                    # Adding A Subplot To The Current Figure
pylab.semilogx(w,S)                                                     # Making A Plot With Log Scaling On The x-axis
pylab.xlabel(r'w')                                                      # Setting The Label For The x-axis
pylab.ylabel(r'$|H(s)|$')                                               # Setting The Label For The y-axis
pylab.subplot(2,1,2)                                                    # Adding A Subplot To The Current Figure
pylab.semilogx(w,phi)                                                   # Making A Plot With Log Scaling On The x-axis
pylab.xlabel(r'w')                                                      # Setting The Label For The x-axis
pylab.ylabel(r'$\angle(H(s))$')                                         # Setting The Label For The y-axis
pylab.show()                                                            # Displaying The Figure

# Plotting Output Voltage (e-6)
time=np.linspace(0,0.00003,1000)                                        # Constructing An Array Of Evenly Spaced Numbers Over A Specified Interval
Vi=np.multiply(np.cos(1000*time)-np.cos(1000000*time),np.heaviside(time,0.5))
_,y1,vec1=sp.lsim(Hv,Vi,time)                                           # Simulate Output Of A Continuous Time LTI system
pylab.figure(5)                                                         # Creating A New Figure
pylab.plot(time,y1,'-b')                                                # Ploting y vs x As Lines And Markers
pylab.title('Output Voltage')                                           # Setting Title Of The Graph
pylab.xlabel(r'time')                                                   # Setting The Label For The x-axis
pylab.ylabel(r'$V_{o}(t)$')                                             # Setting The Label For The y-axis
pylab.grid()                                                            # Displaying The Grid
pylab.show()                                                            # Displaying The Figure

# Plotting Output Voltage (e-3)
time=np.linspace(0,0.01,100000)                                         # Constructing An Array Of Evenly Spaced Numbers Over A Specified Interval
Vi=np.multiply(np.cos(1000*time)-np.cos(1000000*time),np.heaviside(time,0.5))
_,y2,vec2=sp.lsim(Hv,Vi,time)                                           # Simulate Output Of A Continuous Time LTI system
pylab.figure(6)                                                         # Creating A New Figure
pylab.plot(time,y2,'-b')                                                # Ploting y vs x As Lines And Markers
pylab.title('Output Voltage')                                           # Setting Title Of The Graph
pylab.xlabel(r'time')                                                   # Setting The Label For The x-axis
pylab.ylabel(r'$V_{o}(t)$')                                             # Setting The Label For The y-axis
pylab.grid()                                                            # Displaying The Grid
pylab.show()                                                            # Displaying The Figure

# Plotting Output Voltage For Low Frequency Component (e-6)
time=np.linspace(0,0.00003,1000)                                        # Constructing An Array Of Evenly Spaced Numbers Over A Specified Interval
Vi=np.multiply(np.cos(1000*time),np.heaviside(time,0.5))
_,yl1,vecl1=sp.lsim(Hv,Vi,time)                                         # Simulate Output Of A Continuous Time LTI system
pylab.figure(7)                                                         # Creating A New Figure
pylab.plot(time,yl1,'-b')                                               # Ploting y vs x As Lines And Markers
pylab.title('Output Voltage For Low Frequency Component')               # Setting Title Of The Graph
pylab.xlabel(r'time')                                                   # Setting The Label For The x-axis
pylab.ylabel(r'$V_{o}(t)$')                                             # Setting The Label For The y-axis
pylab.grid()                                                            # Displaying The Grid
pylab.show()                                                            # Displaying The Figure

# Plotting Output Voltage For Low Frequency Component (e-3)
time=np.linspace(0,0.01,1000)                                           # Constructing An Array Of Evenly Spaced Numbers Over A Specified Interval
Vi=np.multiply(np.cos(1000*time),np.heaviside(time,0.5))
_,yl2,vecl2=sp.lsim(Hv,Vi,time)                                         # Simulate Output Of A Continuous Time LTI system
pylab.figure(8)                                                         # Creating A New Figure
pylab.plot(time,yl2,'-b')                                               # Ploting y vs x As Lines And Markers
pylab.title('Output Voltage For Low Frequency Component')               # Setting Title Of The Graph
pylab.xlabel(r'time')                                                   # Setting The Label For The x-axis
pylab.ylabel(r'$V_{o}(t)$')                                             # Setting The Label For The y-axis
pylab.grid()                                                            # Displaying The Grid
pylab.show()                                                            # Displaying The Figure

# Plotting Output Voltage For High Frequency Component (e-6)
time=np.linspace(0,0.00003,1000)                                        # Constructing An Array Of Evenly Spaced Numbers Over A Specified Interval
Vi=np.multiply(-np.cos(1000000*time),np.heaviside(time,0.5))
_,yh1,vech1=sp.lsim(Hv,Vi,time)                                         # Simulate Output Of A Continuous Time LTI system
pylab.figure(9)                                                         # Creating A New Figure
pylab.plot(time,yh1,'-b')                                               # Ploting y vs x As Lines And Markers
pylab.title('Output Voltage For High Frequency Component')              # Setting Title Of The Graph
pylab.xlabel(r'time')                                                   # Setting The Label For The x-axis
pylab.ylabel(r'$V_{o}(t)$')                                             # Setting The Label For The y-axis
pylab.grid()                                                            # Displaying The Grid
pylab.show()                                                            # Displaying The Figure

# Plotting Output Voltage For High Frequency Component (e-3)
time=np.linspace(0,0.01,1000)                                           # Constructing An Array Of Evenly Spaced Numbers Over A Specified Interval
Vi=np.multiply(-np.cos(1000000*time),np.heaviside(time,0.5))
_,yh2,vech2=sp.lsim(Hv,Vi,time)                                         # Simulate Output Of A Continuous Time LTI system
pylab.figure(10)                                                        # Creating A New Figure
pylab.plot(time,yh2,'-b')                                               # Ploting y vs x As Lines And Markers
pylab.title('Output Voltage For High Frequency Component')              # Setting Title Of The Graph
pylab.xlabel(r'time')                                                   # Setting The Label For The x-axis
pylab.ylabel(r'$V_{o}(t)$')                                             # Setting The Label For The y-axis
pylab.grid()                                                            # Displaying The Grid
pylab.show()                                                            # Displaying The Figure