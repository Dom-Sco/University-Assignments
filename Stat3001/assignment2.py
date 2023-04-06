import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as ss

data = np.array([8.23, 7.58, 7.39, 9.02, 6.69, 8.05, 8.38, 8.03, 9.54, 12.10])

x = np.sort(data)
y = np.arange(len(x))/float(len(x))

x1 = []
y1 = []

for i in range(len(x)):
    x1.append(x[i])
    x1.append(x[i])
    y1.append(y[i])
    y1.append(y[i])
    
y1 = np.delete(y1, -1)
x1 = np.delete(x1, 0)

x = np.linspace(6, 12, 1000)
y = ss.norm.cdf(x, loc=9, scale=1)


plt.plot(x, y, label="CDF")
plt.plot(x1, y1, label="ECDF")
plt.xlabel("x")
plt.ylabel("CDF")
plt.legend()

xnorm = ss.norm.cdf(x1, loc=9, scale=1)

y2 = ss.norm.cdf(x1, loc=9, scale=1)

distances = np.abs(xnorm-y1)
Dn = np.max(distances)


#pvalue
p = (-1)**(0)*np.exp(-2*(0*Dn)**2)


for i in range(1,100000):
    p+=(-1)**(i)*np.exp(-2*(i*Dn)**2)
    p+=(-1)**(-i)*np.exp(-2*(-i*Dn)**2)

p = p/np.sqrt(len(data))

print(Dn, p)