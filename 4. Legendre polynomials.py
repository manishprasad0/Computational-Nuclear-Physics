import matplotlib.pyplot as plt
import numpy as np
x=np.linspace(-1,1,1000)
p0 = np.ones(np.size(x))
p1 = x
plt.plot(x,p0, label='$P_{0}(x)$')
plt.plot(x,p1, label='$P_{1}(x)$')
n=5
for l in range(2,n+1):
    p2 = ((2*l+1)*x*p1 - l*p0)/(l+1)
    plt.plot(x,p2, label='$P_{{{}}}(x)$'.format(l))
    p0=p1
    p1=p2
plt.legend(loc='lower right')
plt.xlabel('x')
plt.ylabel('$P_{n}(x)$')
plt.title('Legendre polynomials')
plt.show()