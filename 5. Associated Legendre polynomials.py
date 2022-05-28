from matplotlib import pyplot as plt
import numpy as np

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
        pmm = sign * dblfact(2*m - 1)*pow((1.0 - x * x), ((m / 2)) )

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

x=np.linspace(-1,1,1000)
L=5
M=3
Plm = legendre_polynomial(L,M,x)
plt.plot(x,Plm)
plt.show()