import matplotlib.pyplot as plt
import numpy as np

N=8
x = np.linspace(-0.8,1, 100)
for n in range(N):
    for nz in range(n+1):
        plt.plot(x,(n+1.5-(x*(3*nz-n)/3)))

plt.ylim((0,7.5))
plt.axvline(x=0, color='k')
plt.xlabel(r'$delta$')
plt.ylabel(r'E/$\hbar\omega$')
plt.title('Single-particle level spectrum')
plt.show()