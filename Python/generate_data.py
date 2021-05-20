# Script To Generate Data Files For The Least Squares Assignment
from pylab import c_, diag, dot, grid, linspace, logspace, meshgrid, ones, plot, randn, savetxt, show, title, xlabel, ylabel
import scipy.special as sp
N=101                           # Number Of Data Points
k=9                             # Number Of Sets Of Data With Varying Noise
# Generating The Data Points And Adding Noise
t=linspace(0,10,N)                          # Vector t
f=1.05*sp.jv(2,t)-0.105*t                   # Vector f(t)
func=meshgrid(f,ones(k),indexing='ij')[0]   # Generating Signal
scl=logspace(-1,-3,k)                       # Noise
noise=dot(randn(N,k),diag(scl))             # Generating k Vectors
y=func+noise                                # Adding Noise To The Signal
# Shadow Plot
plot(t,y)                                   # Ploting Time And Signal With Noise
xlabel(r'$t$',size=20)                      # X-axis -> Time
ylabel(r'$f(t)+n$',size=20)                 # Y-axis -> Signal With Noise
title(r'Plot Of The Data To Be Fitted')     # Title
grid(True)                                  # Show Grid
# Write Out The Matrix Into A File (fitting.dat)
savetxt("fitting.dat",c_[t,y]) 
show()