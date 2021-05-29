"""
Course: EE2703-Applied Programming Lab
Name: Nihal Gajjala
Roll Number: EE19B044
End Semester Exam
"""

# Importing Libraries
import numpy as np
import pylab

# Function To Calculate Rijkl
def calc(l):
    return np.linalg.norm(np.tile(rijk,(100,1,1,1)).reshape((100,3,3,3,1000))-np.hstack((r_vector,np.zeros((100,1)))).reshape((100,3,1,1,1)),axis=1)
''' 
Rijkl = | rijk - rl |
.reshape - Gives A New Shape To An Array Without Changing Its Data
np.hstack - Stacks Arrays In Sequence Horizontally (Column Wise)
np.linalg.norm - Returns The Norm Of rijk - rl
np.zeros - Returns An Array Of Given Shape And Type, Filled With Zeros
'''

# Function To Calculate The Current Along x-direction And y-direction
def current(x,y):
    return np.array([-np.sin(np.arctan(y/x)),np.cos(np.arctan(y/x))])
'''
np.array - Constructs An Array
np.arctan - Tangent Inverse Function
np.cos - Cosine Function
'''

# Creating Meshgrid
# Size Along x Direction = 3
# Size Along y Direction = 3
# Size Along z Direction = 1000
x=np.linspace(0,2,num=3)            # Constructing An Array Of Evenly Spaced Numbers Over A Specified Interval
y=np.linspace(0,2,num=3)            # Constructing An Array Of Evenly Spaced Numbers Over A Specified Interval
z=np.linspace(0,999,num=1000)       # Constructing An Array Of Evenly Spaced Numbers Over A Specified Interval
X,Y,Z=np.meshgrid(x,y,z)            # Creating Coordinate Matrices From Coordinate Vectors

# Declaring The Radius Of Wire And Number Of Sections It Is Divided Into
radius=10                   # Radius Of Wire
sections=100                # Number Of Sections

# Radius Vector -> (x - cos(phi), y - sin(phi))
r_vector=np.vstack((radius*np.cos(np.linspace(0,2*np.pi,sections)).T,radius*np.sin(np.linspace(0,2*np.pi,sections)).T)).T
# dl Vector -> (x - Acos(phi), y - Asin(phi))
# A = np.pi/5
dl_vector=2*np.pi*radius/sections*np.vstack((np.cos(np.linspace(0,2*np.pi,sections)).T,np.sin(np.linspace(0,2*np.pi,sections)).T)).T
# Phi
phi=np.linspace(0,2*np.pi,sections)
# Position Vector
rijk=np.array((X,Y,Z))
'''
.T - Transposes The Array
np.array - Constructs An Array
np.cos - Cosine Function
np.sin - Sine Function
np.linspace - Constructs An Array Of Evenly Spaced Numbers Over A Specified Interval
np.vstack - Stacks Arrays In Sequence Vertically (Row Wise)
'''

# Plotting Current Elements In x-y Plane
pylab.figure(1)
pylab.plot(r_vector[:,0],r_vector[:,1],'ro',markersize=2.5)             # Plotting y vs x As Lines And Markers
pylab.xlabel(r'x-axis$\rightarrow$',fontsize=15)                        # Setting The Label For The x-axis
pylab.ylabel(r'y-axis$\rightarrow$',fontsize=15)                        # Setting The Label For The y-axis
pylab.title('Current Elements In x-y Plane',fontsize=15)                 # Setting Title Of The Graph
pylab.grid()                                                            # Displaying The Grid
pylab.show()                                                            # Displaying The Figure

# Currents Along x-direction And y-direction
ix,iy=current(r_vector[:,0],r_vector[:,1])          # Calculating x And y Component Of Current Using 'current' Function

# Plotting Current Direction In x-y Plane
pylab.figure(2)
pylab.quiver(r_vector[:,0],r_vector[:,1],ix,iy,scale=50,headwidth=10,headlength=10,width=0.001)
'''
Plotting A 2-D Field Of Arrows
r_vector[:,0] - The x Coordinates Of The Arrow Locations
r_vector[:,1] - The y Coordinates Of The Arrow Locations
ix - The x-direction Components Of The Arrow Vectors
iy - The y-direction Components Of The Arrow Vectors
scale - The Arrow Length Unit
headwidth - Head Width As Multiple Of Shaft Width
headlength - Head Length As Multiple Of Shaft Width
width - Shaft Width In Arrow Units
'''
pylab.xlabel(r'x-axis$\rightarrow$',fontsize=15)                        # Setting The Label For The x-axis
pylab.ylabel(r'y-axis$\rightarrow$',fontsize=15)                        # Setting The Label For The y-axis
pylab.title('Current In Wire In x-y Plane',fontsize=15)                 # Setting Title Of The Graph
pylab.grid()                                                            # Displaying The Grid
pylab.show()                                                            # Displaying The Figure

Rijkl=calc(0)                                                   # Calculating Rijkl Using 'calc' Function
cosine=np.cos(phi).reshape((100,1,1,1))                         # Construction Cosine Vector
dl=dl_vector[:,0].reshape((100,1,1,1))                          # dl Vector Component Along x-direction
Ax=np.sum(cosine*dl*np.exp(1j*Rijkl/10)*dl/Rijkl,axis=0)        # Potential Along x-direction
dl=dl_vector[:,1].reshape((100,1,1,1))                          # dl Vector Component Along y-direction
Ay=np.sum(cosine*dl*np.exp(1j*Rijkl/10)*dl/Rijkl,axis=0)        # Potential Along y-direction
'''
j - Imaginary Unit
np.cos - Cosine Function
np.exp - Exponential Function
np.sum - Summation Function
'''

Bz=(Ay[1,0,:]-Ax[0,1,:]-Ay[-1,0,:]+Ax[0,-1,:])/(4)              # Calculating Magnetic Field
magneticField=(np.zeros(1000).T,np.zeros(1000).T,Bz.T)          # Vectorizing Magnetic Field
'''
Magnetic Field Is Zero Along x-direction
Magnetic Field Is Zero Along y-direction
Magnetic Field Is Bz Along z-direction
'''

# Plotting Magnetic Field Along z-axis
pylab.figure(3)
pylab.loglog(z,np.abs(Bz))                                              # Making A Plot With Log Scaling On Both The Axis
pylab.xlabel(r'z-axis$\rightarrow$',fontsize=15)                        # Setting The Label For The x-axis
pylab.ylabel(r'B(Magnetic Field)$\rightarrow$',fontsize=15)             # Setting The Label For The y-axis
pylab.title('Magnetic Field Along z Axis',fontsize=15)                  # Setting Title Of The Graph
pylab.grid()                                                            # Displaying The Grid
pylab.show()                                                            # Displaying The Figure

# Solving For B Using Least Squares
A=np.hstack([np.ones(len(Bz[300:]))[:,np.newaxis],np.log(z[300:])[:,np.newaxis]])
log_c,b=np.linalg.lstsq(A,np.log(np.abs(Bz[300:])),rcond=None) [0]      # Returns log(c) and b
c=np.exp(log_c)                 # Exponential Function
print("The Value Of b is:")
print(b)                        # Printing The Value Of 'b' In The Terminal
print("The Value Of c is:")
print(c)                        # Printing The Value Of 'c' In The Terminal