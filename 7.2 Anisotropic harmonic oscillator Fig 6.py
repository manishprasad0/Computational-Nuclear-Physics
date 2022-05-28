import matplotlib.pyplot as plt
from fractions import Fraction
Ntotal=6
x = [-4,0,2,8]
lname=['s','p','d','f','g','h','i','j']
for N in range(Ntotal):
    E1=N+1.5
    for l in range(N+1):
        n=((N-l)/2)+1
        if float(n).is_integer():
            j=[l-0.5,l+0.5]
            if j[0]>0:
                E2=N+1.5+0.05*(l+1)-0.0225*l*(l+1)
                plt.plot(x,[E1,E1,E2,E2])
                plt.text(8.7, E2,"{}{}{}".format(int(n),lname[l],str(Fraction(j[0]))),fontsize = 6)
            E2=N+1.5-0.05*(l)-0.0225*l*(l+1)
            plt.plot(x,[E1,E1,E2,E2])
            plt.text(8.7, E2,"{}{}{}".format(int(n),lname[l],str(Fraction(j[1]))),fontsize = 6)
    plt.text(-3,E1+0.1,"N={}".format(int(N)),fontsize = 8)
plt.title('The Single-Particle Energy Level Positions for the Spherical Shell Model')
plt.show()