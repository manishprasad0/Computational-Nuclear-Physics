from math import atan
import numpy as np

def factorial(n):
	fct=[1.0]
	for i in range(n+1):
		fct.append((i+1)*fct[i])
	return fct[n]

def dblfact(n):
    if n==0 or n==1:
        return 1
    if n%2==1:
        return n*dblfact(n-2)
    if n%2==0:
        return (n-1)*dblfact(n-3)
    

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

def simpInt(h, fnc):
	I = (fnc[0] + fnc[len(fnc)-1])
	for i in range(1,len(fnc)-2,2):
		I = I + 4*fnc[i] + 2*fnc[i+1]
	I = I + 4*fnc[i+2]
	I = I*h/3
	return I

pi= 4*atan(1.0) 
L=2
M=0
R=1
theta = np.linspace(0, pi, 1000)
phi = np.linspace(0, 2* pi, 181)
r=np.linspace(0,R,1000)

#The Integral is:
#    Qlm = Integral(r^l Ylm(theta,phi) Rho(r) dTau)
#        = Integral((r^l Rho(r) 4*pi*r^2 dr) (Ylm_theta sin(theta) d_theta)  (Ylm_phi d_phi))
#        = R(r)  *   THETA(theta)    *   PHI(phi)    Using Separation of Variables

K = np.sqrt( ((2*L + 1)* factorial(L - abs(M)))/ (4* pi* factorial(L + abs(M))) )
Ylm_theta = np.sin(theta)*np.sin(theta)*legendre_polynomial(L,M,np.cos(theta))
Qlm_theta = simpInt(theta[1]-theta[0],Ylm_theta)
Ylm_r = r*r*r*r
Qlm_r = simpInt(r[1]-r[0],Ylm_r)
Ylm_phi=np.cos(M*phi)
Qlm_phi = simpInt(phi[1]-phi[0],Ylm_phi)
print(Qlm_r*Qlm_theta*Qlm_phi*K)

#Hence we can evaluate the integrals separately. And we can see below that the THETA(theta) Integral evaluates to zero. Hence Q20 = 0