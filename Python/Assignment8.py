"""
Course: EE2703-Applied Programming Lab
Name: Nihal Gajjala
Roll Number: EE19B044
Assignment 8
"""

import numpy as np
import sympy
import scipy.signal as sp
import pylab

s=sympy.symbols('s')            # Declaring symbols

# Low Pass Filter
def lowpass(R1,R2,C1,C2,G,Vi):
    s=sympy.symbols('s')                                                                                    # Declaring symbols
    A=sympy.Matrix([[0,0,1,-1/G],[-1/(1+s*R2*C2),1,0,0],[0,-G,G,1],[-1/R1-1/R2-s*C1,1/R2,0,s*C1]])          # Declaring Matrices
    b=sympy.Matrix([0,0,0,-Vi/R1])                                                                          # Declaring Matrices
    V=A.inv()*b                                                                                             # Inverse Matrix
    return (A,b,V)

# High Pass Filter
def highpass(R1,R2,C1,C2,G,Vi):
    s=sympy.symbols('s')                                                                                    # Declaring symbols
    A=sympy.Matrix([[0,0,1,-1/G],[-1/(1+(1/s*R2*C2)),1,0,0], [0,-G,G,1],[-1/R1-s*C1-s*C2,s*C2,0,1/R1]])     # Declaring Matrices
    b=sympy.Matrix([0,0,0,-Vi*s*C1])                                                                        # Declaring Matrices
    V=A.inv()*b                                                                                             # Inverse Matrix
    return (A,b,V)

#  Convert Sympy Transfer Function Polynomial To Scipy LTI
def sympy_lti(S):
    S=sympy.simplify(S)                                                                                     # Simplifies The Expression
    numerator,denominator=sympy.fraction(S)                                                                 # Splits Numerator And Denominator
    numerator,denominator=sympy.Poly(numerator,s),sympy.Poly(denominator,s)                                 # Polynomial In 's'
    numerator,denominator=numerator.all_coeffs(),denominator.all_coeffs()                                   # Storing Coefficents
    numerator,denominator=[float(f) for f in numerator],[float(f) for f in denominator]
    return numerator,denominator

# Step Function
def u(t):
    func=[]
    for i in t:
        if i<0:
            func.append(0)              # Adds An Item To The Existing List
        else:
            func.append(1)              # Adds An Item To The Existing List
    return func

# Constants
k=0.9
freq=2e3*np.pi
decay=freq/(2*np.pi)*np.log(1/k)

# Input Voltage 1
vi=lambda t:(np.sin(2e3*np.pi*t)+np.cos(2e6*np.pi*t))
vilf=lambda t:(np.sin(2e3*np.pi*t))
vihf=lambda t:np.cos(2e6*np.pi*t)

# Damped Input (Input 2 & Input 3)
dampsin=lambda t:(np.sin(freq*t)*np.exp(-1*decay*t))
dampcos=lambda t:(np.cos(freq*t)*np.exp(-1*decay*t))

# Values For Low Pass Filter
R1lp=1e4
R2lp=1e4
C1lp=1e-9
C2lp=1e-9
Glp=1.586

# Values For High Pass Filter
R1hp=1e4
R2hp=1e4
C1hp=1e-9
C2hp=1e-9
Ghp=1.586

# Creating Time Response Arrays
time1=np.linspace(0,0.0005,10000)                   # Constructing An Array Of Evenly Spaced Numbers Over A Specified Interval
time2=np.linspace(0,0.01,100000)                    # Constructing An Array Of Evenly Spaced Numbers Over A Specified Interval
time3=np.linspace(0,0.01,10000)                     # Constructing An Array Of Evenly Spaced Numbers Over A Specified Interval
time4=np.linspace(0,0.001,10000)                    # Constructing An Array Of Evenly Spaced Numbers Over A Specified Interval
time5=np.linspace(0,0.00001,10000)                  # Constructing An Array Of Evenly Spaced Numbers Over A Specified Interval

# Creating Frequency Array
w=np.logspace(0,8,801)                              # Constructing An Array Of Evenly Spaced Numbers In Logarithmic Scale
jw=1j*w                                             # Complex Value

# Low Pass Filter
Alp,blp,Vlp=lowpass(R1lp,R2lp,C1lp,C2lp,Glp,1) 
Hlps=Vlp[3]
numerator1,denominator1=sympy_lti(Hlps)                 
Hlp=sp.lti(numerator1,denominator1)                     # Continuous Time LTI System
hlpf=sympy.lambdify(s,Hlps,'numpy')                     # Translating SymPy Expressions Into Python Functions
Hlpjw=hlpf(jw)
t1,voutlp1,svec=sp.lsim(Hlp,u(time1),time1)             # Simulate Output Of A Continuous Time LTI system
t2,voutlp2,svec=sp.lsim(Hlp,vi(time2),time2)            # Simulate Output Of A Continuous Time LTI system
t3,voutlp3,svec=sp.lsim(Hlp,dampsin(time3),time3)       # Simulate Output Of A Continuous Time LTI system
t3,voutlp4,svec=sp.lsim(Hlp,dampcos(time3),time3)       # Simulate Output Of A Continuous Time LTI system

# High Pass Filter
Ahp,bhp,Vhp=highpass(R1hp,R2hp,C1hp,C2hp,Ghp,1)
Hhps=Vhp[3]
numerator2,denominator2=sympy_lti(Hhps)                 
Hhp=sp.lti(numerator2,denominator2)                     # Continuous Time LTI System
hhpf=sympy.lambdify(s,Hhps,'numpy')                     # Translating SymPy Expressions Into Python Functions
Hhpjw=hhpf(jw)
t1,vouthp1,svec=sp.lsim(Hhp,u(time1),time1)             # Simulate Output Of A Continuous Time LTI system
t2,vouthp2,svec=sp.lsim(Hhp,vi(time2),time2)            # Simulate Output Of A Continuous Time LTI system
t3,vouthp3,svec=sp.lsim(Hhp,dampsin(time3),time3)       # Simulate Output Of A Continuous Time LTI system
t3,vouthp4,svec=sp.lsim(Hhp,dampcos(time3),time3)       # Simulate Output Of A Continuous Time LTI system

# Magnitude Response Of Low Pass Filter
pylab.figure(0)                                         
pylab.loglog(w,abs(Hlpjw),lw=2)                         # Plotting y vs x As Lines And Markers
pylab.title('Magnitude Response Of Low Pass Filter')    # Setting Title Of The Graph
pylab.xlabel(r'freq$\rightarrow$')                      # Setting The Label For The x-axis
pylab.ylabel(r'|H|$\rightarrow$')                       # Setting The Label For The y-axis
pylab.grid()                                            # Displaying The Grid
pylab.show()                                            # Displaying The Figure

# Step Response Of Low Pass Filter
pylab.figure(1)                                         
pylab.plot(time1*1e6,voutlp1)                           # Plotting y vs x As Lines And Markers
pylab.title('Step Response Of Low Pass Filter')         # Setting Title Of The Graph
pylab.xlabel(r'time$\rightarrow$')                      # Setting The Label For The x-axis
pylab.ylabel(r'Vout$\rightarrow$')                      # Setting The Label For The y-axis
pylab.grid()                                            # Displaying The Grid
pylab.show()                                            # Displaying The Figure

# Magnitude Response Of High Pass Filter
pylab.figure(2)                                         
pylab.loglog(w,abs(Hhpjw),lw=2)                         # Plotting y vs x As Lines And Markers
pylab.title('Magnitude Response Of High Pass Filter')   # Setting Title Of The Graph
pylab.xlabel(r'freq$\rightarrow$')                      # Setting The Label For The x-axis
pylab.ylabel(r'|H|$\rightarrow$')                       # Setting The Label For The y-axis
pylab.grid()                                            # Displaying The Grid
pylab.show()                                            # Displaying The Figure

# Step Response Of High Pass Filter
pylab.figure(3)                                         
pylab.plot(time1*1e6,vouthp1)                           # Plotting y vs x As Lines And Markers
pylab.title('Step Response Of High Pass Filter')        # Setting Title Of The Graph
pylab.xlabel(r'time$\rightarrow$')                      # Setting The Label For The x-axis
pylab.ylabel(r'Vout$\rightarrow$')                      # Setting The Label For The y-axis
pylab.grid()                                            # Displaying The Grid
pylab.show()                                            # Displaying The Figure

# Input Voltage 1
pylab.figure(4)                                         
pylab.plot(time4*1e3,vi(time4))                         # Plotting y vs x As Lines And Markers
pylab.title('Input Voltage 1')                          # Setting Title Of The Graph
pylab.xlabel(r'time$\rightarrow$')                      # Setting The Label For The x-axis
pylab.ylabel(r'Vin$\rightarrow$')                       # Setting The Label For The y-axis
pylab.grid()                                            # Displaying The Grid
pylab.show()                                            # Displaying The Figure

pylab.figure(5)                                                 
pylab.plot(time5*1e6,vihf(time5))                               # Plotting y vs x As Lines And Markers
pylab.title('High Frequency Component of Input Voltage 1')      # Setting Title Of The Graph
pylab.xlabel(r'time$\rightarrow$')                              # Setting The Label For The x-axis
pylab.ylabel(r'Vin$\rightarrow$')                               # Setting The Label For The y-axis
pylab.grid()                                                    # Displaying The Grid
pylab.show()                                                    # Displaying The Figure

pylab.figure(6)                                                 
pylab.plot(time4*1e3,vilf(time4))                               # Plotting y vs x As Lines And Markers
pylab.title('Low Frequency Component of Input Voltage 1')       # Setting Title Of The Graph
pylab.xlabel(r'time$\rightarrow$')                              # Setting The Label For The x-axis
pylab.ylabel(r'Vin$\rightarrow$')                               # Setting The Label For The y-axis
pylab.grid()                                                    # Displaying The Grid
pylab.show()                                                    # Displaying The Figure

# Input Voltage 2
pylab.figure(7)                                         
pylab.plot(time3*1e3,dampsin(time3))                    # Plotting y vs x As Lines And Markers
pylab.title('Input Voltage 2')                          # Setting Title Of The Graph
pylab.xlabel(r'time$\rightarrow$')                      # Setting The Label For The x-axis
pylab.ylabel(r'Vin$\rightarrow$')                       # Setting The Label For The y-axis
pylab.grid()                                            # Displaying The Grid
pylab.show()                                            # Displaying The Figure

# Input Voltage 3
pylab.figure(8)                                         
pylab.plot(time3*1e3,dampcos(time3))                    # Plotting y vs x As Lines And Markers
pylab.title('Input Voltage 3')                          # Setting Title Of The Graph
pylab.xlabel(r'time$\rightarrow$')                      # Setting The Label For The x-axis
pylab.ylabel(r'Vin$\rightarrow$')                       # Setting The Label For The y-axis
pylab.grid()                                            # Displaying The Grid
pylab.show()                                            # Displaying The Figure

# Low Pass Filter Response
pylab.figure(9)                                         
pylab.plot(time2*1e3,voutlp2)                           # Plotting y vs x As Lines And Markers
pylab.title('Input 1 Response For Low Pass Filter')     # Setting Title Of The Graph
pylab.xlabel(r'time$\rightarrow$')                      # Setting The Label For The x-axis
pylab.ylabel(r'Vout$\rightarrow$')                      # Setting The Label For The y-axis
pylab.grid()                                            # Displaying The Grid
pylab.show()                                            # Displaying The Figure

pylab.figure(10)                                        
pylab.plot(time3*1e3,voutlp3)                           # Plotting y vs x As Lines And Markers
pylab.title('Input 2 Response For Low Pass Filter')     # Setting Title Of The Graph
pylab.xlabel(r'time$\rightarrow$')                      # Setting The Label For The x-axis
pylab.ylabel(r'Vout$\rightarrow$')                      # Setting The Label For The y-axis
pylab.grid()                                            # Displaying The Grid
pylab.show()                                            # Displaying The Figure

pylab.figure(11)                                        
pylab.plot(time3*1e3,voutlp4)                           # Plotting y vs x As Lines And Markers
pylab.title('Input 3 Response For Low Pass Filter')     # Setting Title Of The Graph
pylab.xlabel(r'time$\rightarrow$')                      # Setting The Label For The x-axis
pylab.ylabel(r'Vout$\rightarrow$')                      # Setting The Label For The y-axis
pylab.grid()                                            # Displaying The Grid
pylab.show()                                            # Displaying The Figure

# High Pass Filter Response
pylab.figure(12)                                        
pylab.plot(time2*1e3,vouthp2)                           # Plotting y vs x As Lines And Markers
pylab.title('Input 1 Response For High Pass Filter')    # Setting Title Of The Graph
pylab.xlabel(r'time$\rightarrow$')                      # Setting The Label For The x-axis
pylab.ylabel(r'Vout$\rightarrow$')                      # Setting The Label For The y-axis
pylab.grid()                                            # Displaying The Grid
pylab.show()                                            # Displaying The Figure

pylab.figure(13)                                        
pylab.plot(time3*1e3,vouthp3)                           # Plotting y vs x As Lines And Markers
pylab.title('Input 2 Response For High Pass Filter')    # Setting Title Of The Graph
pylab.xlabel(r'time$\rightarrow$')                      # Setting The Label For The x-axis
pylab.ylabel(r'Vout$\rightarrow$')                      # Setting The Label For The y-axis
pylab.grid()                                            # Displaying The Grid
pylab.show()                                            # Displaying The Figure

pylab.figure(14)                                        
pylab.plot(time3*1e3,vouthp4)                           # Plotting y vs x As Lines And Markers
pylab.title('Input 3 Response For High Pass Filter')    # Setting Title Of The Graph
pylab.xlabel(r'time$\rightarrow$')                      # Setting The Label For The x-axis
pylab.ylabel(r'Vout$\rightarrow$')                      # Setting The Label For The y-axis
pylab.grid()                                            # Displaying The Grid
pylab.show()                                            # Displaying The Figure