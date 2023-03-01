import numpy as np

def bp(r,s):
    o = 1-r-s
    Sn = [1]
    e = (r+2*s)**100
    
    i = 0
    while Sn[i]!=0:
        S = np.random.multinomial(Sn[i], [o,r,s], size=1)
        Sn.append(0*S[0][0]+1*S[0][1]+2*S[0][2])
        i+=1
    print(Sn)
    print(e)