"""
Course: EE2703-Applied Programming Lab
Name: Nihal Gajjala
Roll Number: EE19B044
Assignment 9
"""

import numpy as np
import pylab

# Plotting Spectrum
def subp(w,Y,xlim,title,xlabel,ylabel1,ylabel2,phasemin):
    pylab.figure()
    pylab.title(title)                              # Setting Title Of The Graph
    pylab.subplot(2,1,1)                            # Adding A Subplot To The Current Figure
    pylab.plot(w,abs(Y))                            # Plotting y vs x As Lines And Markers
    pylab.ylabel(ylabel1)                           # Setting The Label For The y-axis
    pylab.xlim([-xlim,xlim])                        # Setting The x-limits Of The x-axis
    pylab.grid()                                    # Displaying The Grid
    pylab.subplot(2,1,2)                            # Adding A Subplot To The Current Figure
    ii=np.where(abs(Y)>phasemin)                    # Finding Elements From x Or y Depending On The Condition
    pylab.plot(w[ii],np.angle(Y[ii]),'go')          # Plotting y vs x As Lines And Markers
    pylab.xlim([-xlim,xlim])                        # Setting The x-limits Of The x-axis
    pylab.xlabel(xlabel)                            # Setting The Label For The x-axis
    pylab.ylabel(ylabel2)                           # Setting The Label For The y-axis
    pylab.grid()                                    # Displaying The Grid
    pylab.show()                                    # Displaying The Figure

func1=lambda t:(np.sin(t))**3                       # Function 1
func2=lambda t:(np.cos(t))**3                       # Function 2
func3=lambda t:np.cos(20*t+5*np.cos(t))             # Function 3
func4=lambda t:(1/2)*np.exp((t**2)*-0.5)            # Function 4
func5=lambda t:np.sin(5*t)                          # Function 5
func6=lambda t:(1+0.1*np.cos(t))*np.cos(10*t)       # Function 6

# Spectrum of sin(5t)
N=512
t=np.linspace(-4*np.pi,4*np.pi,N+1);t=t[:-1]        # Constructing An Array Of Evenly Spaced Numbers Over A Specified Interval
w=np.linspace(-64,64,N+1);w=w[:-1]                  # Constructing An Array Of Evenly Spaced Numbers Over A Specified Interval
y=func5(t)
Y=np.fft.fftshift(np.fft.fft(y))/(N)                # Computing The One-Dimensional Discrete Fourier Transform
subp(w,Y,15,"Spectrum of sin(5t)","k","|Y|","Phase of Y",1e-3)

# Spectrum of (1+0.1cos(t))cos(10t)
N=512
t=np.linspace(-4*np.pi,4*np.pi,N+1);t=t[:-1]        # Constructing An Array Of Evenly Spaced Numbers Over A Specified Interval
w=np.linspace(-64,64,N+1);w=w[:-1]                  # Constructing An Array Of Evenly Spaced Numbers Over A Specified Interval
y=func6(t)
Y=np.fft.fftshift(np.fft.fft(y))/(N)                # Computing The One-Dimensional Discrete Fourier Transform
subp(w,Y,15,"Spectrum of (1+0.1cos(t))cos(10t)","k","|Y|","Phase of Y",1e-3)

# Spectrum of sin^3(t)
N=512
t=np.linspace(-4*np.pi,4*np.pi,N+1);t=t[:-1]        # Constructing An Array Of Evenly Spaced Numbers Over A Specified Interval
w=np.linspace(-64,64,N+1);w=w[:-1]                  # Constructing An Array Of Evenly Spaced Numbers Over A Specified Interval
y=func1(t)
Y=np.fft.fftshift(np.fft.fft(y))/(N)                # Computing The One-Dimensional Discrete Fourier Transform
subp(w,Y,15,"Spectrum of sin^3(t)","k","|Y|","Phase of Y",1e-3)

# Spectrum of cos^3(t)
N=512
t=np.linspace(-4*np.pi,4*np.pi,N+1);t=t[:-1]        # Constructing An Array Of Evenly Spaced Numbers Over A Specified Interval
w=np.linspace(-64,64,N+1);w=w[:-1]                  # Constructing An Array Of Evenly Spaced Numbers Over A Specified Interval
y=func2(t)
Y=np.fft.fftshift(np.fft.fft(y))/(N)                # Computing The One-Dimensional Discrete Fourier Transform
subp(w,Y,15,"Spectrum of cos^3(t)","k","|Y|","Phase of Y",1e-3)

# Spectrum of cos(20*t+5*cos(t))
N=512
t=np.linspace(-4*np.pi,4*np.pi,N+1);t=t[:-1]        # Constructing An Array Of Evenly Spaced Numbers Over A Specified Interval
w=np.linspace(-64,64,N+1);w=w[:-1]                  # Constructing An Array Of Evenly Spaced Numbers Over A Specified Interval
y=func3(t)
Y=np.fft.fftshift(np.fft.fft(y))/(N)                # Computing The One-Dimensional Discrete Fourier Transform
subp(w,Y,40,"Spectrum of cos(20*t+5*cos(t))","k","|Y|","Phase of Y",1e-3)

# Gauss Function For N=256
N=256
t=np.linspace(-4*np.pi,4*np.pi,N+1);t=t[:-1]        # Constructing An Array Of Evenly Spaced Numbers Over A Specified Interval
w=np.linspace(-N/2,N/2,N+1);w=w[:-1]                # Constructing An Array Of Evenly Spaced Numbers Over A Specified Interval
y=func4(t)
Y1=np.fft.fftshift(abs(np.fft.fft(y)))/N            # Computing The One-Dimensional Discrete Fourier Transform
Y1=Y1*np.sqrt(2*np.pi)/max(Y1)                      # Gauss Function
subp(w,Y1,15,"Spectrum of Gauss function","k","|Y|","Phase of Y",1e-3)
t=np.linspace(-8*np.pi,8*np.pi,N+1);t=t[:-1]        # Constructing An Array Of Evenly Spaced Numbers Over A Specified Interval
w=np.linspace(-64,64,N+1);w=w[:-1]                  # Constructing An Array Of Evenly Spaced Numbers Over A Specified Interval
y=func4(t)
Y2=np.fft.fftshift(abs(np.fft.fft(y)))/N            # Computing The One-Dimensional Discrete Fourier Transform
Y2=Y2*np.sqrt(2*np.pi)/max(Y2)                      # Gauss Function
subp(w,Y2,15,"Spectrum of Gauss function","k","|Y|","Phase of Y",1e-3)

# Gauss Function For N=512
N=512
t=np.linspace(-4*np.pi,4*np.pi,N+1);t=t[:-1]        # Constructing An Array Of Evenly Spaced Numbers Over A Specified Interval
w=np.linspace(-N/2,N/2,N+1);w=w[:-1]                # Constructing An Array Of Evenly Spaced Numbers Over A Specified Interval
y=func4(t)
Y1=np.fft.fftshift(abs(np.fft.fft(y)))/N            # Computing The One-Dimensional Discrete Fourier Transform
Y1=Y1*np.sqrt(2*np.pi)/max(Y1)                      # Gauss Function
subp(w,Y1,15,"Spectrum of Gauss function","k","|Y|","Phase of Y",1e-3)
t=np.linspace(-8*np.pi,8*np.pi,N+1);t=t[:-1]        # Constructing An Array Of Evenly Spaced Numbers Over A Specified Interval
w=np.linspace(-64,64,N+1);w=w[:-1]                  # Constructing An Array Of Evenly Spaced Numbers Over A Specified Interval
y=func4(t)
Y2=np.fft.fftshift(abs(np.fft.fft(y)))/N            # Computing The One-Dimensional Discrete Fourier Transform
Y2=Y2*np.sqrt(2*np.pi)/max(Y2)                      # Gauss Function
subp(w,Y2,15,"Spectrum of Gauss function","k","|Y|","Phase of Y",1e-3)

# Gauss Function For N=1024
N=1024
t=np.linspace(-4*np.pi,4*np.pi,N+1);t=t[:-1]        # Constructing An Array Of Evenly Spaced Numbers Over A Specified Interval
w=np.linspace(-N/2,N/2,N+1);w=w[:-1]                # Constructing An Array Of Evenly Spaced Numbers Over A Specified Interval
y=func4(t)
Y1=np.fft.fftshift(abs(np.fft.fft(y)))/N            # Computing The One-Dimensional Discrete Fourier Transform
Y1=Y1*np.sqrt(2*np.pi)/max(Y1)                      # Gauss Function
subp(w,Y1,15,"Spectrum of Gauss function","k","|Y|","Phase of Y",1e-3)
t=np.linspace(-8*np.pi,8*np.pi,N+1);t=t[:-1]        # Constructing An Array Of Evenly Spaced Numbers Over A Specified Interval
w=np.linspace(-64,64,N+1);w=w[:-1]                  # Constructing An Array Of Evenly Spaced Numbers Over A Specified Interval
y=func4(t)
Y2=np.fft.fftshift(abs(np.fft.fft(y)))/N            # Computing The One-Dimensional Discrete Fourier Transform
Y2=Y2*np.sqrt(2*np.pi)/max(Y2)                      # Gauss Function
subp(w,Y2,15,"Spectrum of Gauss function","k","|Y|","Phase of Y",1e-3)