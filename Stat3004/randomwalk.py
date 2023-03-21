import numpy as np

j = 0.3
q = 1-j

x=np.random.binomial(size=1000, n=1, p= j)

u=-100
s = 0

for i in range(len(x)):
    s+=2*x[i]-1
    if s==u:
        print("GOTTEM:",i)

print(u/(j-q))