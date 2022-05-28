import numpy as np
from math import atan
import matplotlib.pyplot as plt

def Spherical2Cartesian(r, theta, phi):
    x = r* np.sin(theta)* np.cos(phi)
    y = r* np.sin(theta)* np.sin(phi)
    z = r* np.cos(theta)
    
    return x, y, z

def factorial(n):
	fct=[1.0]
	for i in range(n+1):
		fct.append((i+1)*fct[i])
	return fct[n]

# Double Factorial: Product of all odd integers less than or equal to n
def dblfact(n):
    if n==0 or n==1:
        return 1
    if n%2==1:
        return n*dblfact(n-2)
    if n%2==0:
        return (n-1)*dblfact(n-3)

#Associated Legendre Polynomials
def legendre_polynomial(l, m, x):
    pmm = 1.0
    if m > 0:
        sign = 1.0 if m % 2 == 0 else -1.0
        pmm = sign * dblfact(2 * m - 1)*pow((1.0 - x * x), ((m / 2)) )

    if l == m:
        return pmm

    pmm1 = x * (2 * m + 1) * pmm
    if l == m + 1:
        return pmm1

    for n in range(m + 2, l + 1):
        pmn = (x * (2 * n - 1) * pmm1 - (n + m - 1) * pmm) / (n - m)
        pmm = pmm1
        pmm1 = pmn

    return pmm1

#Azimuthal and Magnetic quantum nubers
l = int(input('Azimuthal quantum number "l": '))
if l < 0:
    print('"l" can not be a negative number')

m = int(input('Magnetic quantum number "m": '))
if m > l or m < (-l):
    print('quntum number "m" belong to "[-l,l]"')

pi=4*atan(1.0)

#Normalization constant
K = np.sqrt( ((2*l + 1)* factorial(l - abs(m)))/ (4* pi* factorial(l + abs(m))) )

#Value of Phi and Theta
phi = np.linspace(0, 2* pi, 181)
tht = np.linspace(0, pi, 91)
Phi, Tht = np.meshgrid(phi, tht)

#Value of Y
p, q = np.shape(Phi)
Y    = np.zeros([p, q],dtype=np.complex_)

for i in range(0, p):
    for j in range(0, q):
        
        if m > 0:      
            #Y[i, j] = np.sqrt(2) * K * np.exp(1j*m * Phi[i, j]) * legendre_polynomial(l, m, np.cos(Tht[i, j]))
            Y[i, j] = np.sqrt(2) * K * np.cos(m * Phi[i, j]) * legendre_polynomial(l, m, np.cos(Tht[i, j]))
        elif m < 0:
            Y[i, j] = np.sqrt(2) * K * np.sin(abs(m) * Phi[i, j]) * legendre_polynomial(l, abs(m), np.cos(Tht[i, j]))
        
        else:
            Y[i, j] = K * legendre_polynomial(l, 0, np.cos(Tht[i, j]))

#Convert data to cartesian form
x, y, z = Spherical2Cartesian(abs((Y)), Tht, Phi)

#Plotting
fig = plt.figure('Harmonics')
ax = fig.add_subplot( 111 , projection='3d')

surf=ax.plot_surface(x, y, z, cmap = 'jet', edgecolor = 'k')
fig.colorbar(surf)
#plt.axis('auto')
plt.title('Spherical Harmonics', fontsize = 14, fontweight = 'bold')
#ax.view_init(azim=0, elev=0)
ax.set_xlabel('X', fontsize = 12, fontweight = 'bold')
ax.set_ylabel('Y', fontsize = 12, fontweight = 'bold')
ax.set_zlabel('Z', fontsize = 12, fontweight = 'bold')
plt.show()
