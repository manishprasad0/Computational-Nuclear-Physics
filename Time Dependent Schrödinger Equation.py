# Solving the TDSE using Cayley's operator. The code contains the solutions both the infinite box and the harmonic oscillator. 
# To get the required plot, add and remove the commented parts accordingly. For example, remove the comments from the harmonic oscillator cell
# and add comments to the infinite box cell to get the animated plot for the harmonic oscillator.

import numpy as np
from math import pi
import matplotlib.pyplot as plt
from scipy.sparse import spdiags
from scipy.sparse.linalg import splu

# 2nd order derivative (3pt central difference)
#-----------------------------------------------------------------------
def diff2(f,h,J):
    d2f=(np.roll(f,-1)+np.roll(f,+1)-2*f)/h**2
    d2f[0]=d2f[J-1]=0
    return d2f

# LU matrix on the left: tridiagonal case
#-----------------------------------------------------------------------
def lhs_lumatrix(J,dx,dt,V,hb,hb2m):
    a = 1 + 1j*dt*(2*hb2m/dx**2 + V)/(2*hb)
    b = (-1j*hb2m*dt/(2*hb*dx**2))*np.ones((J),float)
    return splu(spdiags(np.array([b,a,b]),np.array([-1,0,+1]),J,J).tocsc())

# zeta vector on the right: tridiagonal case
#-----------------------------------------------------------------------
def zeta(J,psi,dx,dt,V,hb,hb2m):
    return psi - 1j*(dt/hb)*(-hb2m*diff2(psi,dx,J)+V*psi)/2

hb = 1                                     #Planck's constant, hbar
m =	1                                      #mass of particle, m

""" 
# defining the system for harmonic oscillator
#-----------------------------------------------------------------------

xmin,xmax,dx = -25,+25,0.01                #x-limits of simulation box
x = np.arange(xmin,xmax+dx,dx)             #defining the position grid
J = len(x)                                 #dimension of position grid

w = 0.1                                    #freq of harmonic oscillator
V = 1/2*m*w**2*x**2                        #harmonic oscillator potential

x0,p0,sig = -10,0,1                        #initial position,momentum, position spread

psi = np.exp( -((x-x0)/(2*sig))**2 + 1j*p0*(x-x0) )/np.sqrt( sig*np.sqrt(2*pi) )

tmax,dt,plot_steps =  4*(2*pi/w),0.01,10   #time limit, time step, and interval b/w 2 successive plots
 """


# defining the system for particle in a box
#-----------------------------------------------------------------------

xmin,xmax,dx = -10,+10,0.01                #x-limits of simulation box
x = np.arange(xmin,xmax+dx,dx)             #defining the position grid
J = len(x)                                 #dimension of position grid

V = np.zeros(J,float)                      #zero potential
V[0],V[J-1] = np.tan(pi/2),np.tan(pi/2)    #infinite potential at the boundaries

x0,p0,sig = 0,1,1                          #initial position,momentum, position spread

psi = np.exp( -((x-x0)/(2*sig))**2 + 1j*p0*(x-x0) )/np.sqrt( sig*np.sqrt(2*pi) )

tmax,dt,plot_steps = 100,0.01,10           #time limit, time step, and interval b/w 2 successive plots



# solving the TDSE
#-----------------------------------------------------------------------

hb2m = hb**2/(2*m)                         #value of hbar^2/2m
lhs_lu = lhs_lumatrix(J,dx,dt,V,hb,hb2m)   #LU decomposition for the LHS matrix
t_range=np.linspace(0,tmax+dt,len(x)) 
t = 0
ymax = 1.01*np.max(np.abs(psi)**2)

while t < tmax:

    for j in range(plot_steps):            #evolve plot_steps times
        psi = lhs_lu.solve(zeta(J,psi,dx,dt,V,hb,hb2m))   
    
    t = t + plot_steps*dt
    
    plt.cla()
    plt.xlim(0,tmax)
    plt.ylim(0,1)
    plt.xlabel('$x$',fontsize=14)
    plt.ylabel('$|\psi|^2$',fontsize=14)
    plt.title(f'$t$ = {t:.1f}')
    plt.grid()
    plt.plot(t_range,np.abs(psi)**2,c='red')
    plt.pause(0.001)