import numpy as np

w = np.array([1,0,0,1,0,1,0,1,1,0,0,1,1,1,1,0,0,0,1,0,0,0,0])

def C24IMLD(w):
    B = np.array([[1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1],
                  [1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1],
                  [0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1],
                  [1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1],
                  [1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1],
                  [1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1],
                  [0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1],
                  [0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1],
                  [0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1],
                  [1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1],
                  [0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1],
                  [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0]])

    H = np.vstack((np.eye(12),B))
    t = np.array([0,0,0,0,0,0,0,0,0,0,0,0])
    
    syndrome = (w@H)%2
    if sum(syndrome)<=3:
        e = np.hstack((syndrome,t))
        return (e+w)%2
    
    for i in range(B.shape[0]):
        if sum((syndrome+B[i])%2)<=2:
            t = np.array([0,0,0,0,0,0,0,0,0,0,0,0])
            t[i] += 1
            e = np.hstack(((syndrome+B[i])%2,t))
            return (e+w)%2
    
    s2 = syndrome@B
    if sum(s2)<=3:
        e = np.hstack((t,s2))
        return (e+w)%2
    
    for i in range(B.shape[0]):
        if sum((s2+B[i])%2)<=2:
            t = np.array([0,0,0,0,0,0,0,0,0,0,0,0])
            t[i] += 1
            e = np.hstack((t,(s2+B[i])%2))
            return (e+w)%2
    print("Request retransmission")
    

def C23IMLD(w):
    if sum(w)%2==1:
        ws = np.hstack((w,np.array([0])))
    else:
        ws = np.hstack((w,np.array([1])))
    
    cs = C24IMLD(ws)[:-1]
    return cs
    

print(C23IMLD(w))

w = np.array([1,0,1,1,0,1,0,0,0,0,1,0,1,1,1])

            

def cyclicBurst(w):
    gx = np.poly1d([1,0,0,1,1,1,1]) #the polynomial g(x)
    #example 7.6.3
    H = np.array([[1, 1, 1, 1, 0, 0,],
                  [0, 1, 1, 1, 1, 0],
                  [0, 0, 1, 1, 1, 1],
                  [1, 1, 1, 0, 1, 1],
                  [1, 0, 0, 0, 0, 1],
                  [1, 0, 1, 1, 0, 0],
                  [0, 1, 0, 1, 1, 0],
                  [0, 0, 1, 0, 1, 1],
                  [1, 1, 1, 0, 0, 1]])
    H = np.vstack((H,np.eye(6)))
    s = (w@H)%2
    sp = np.poly1d(np.flip(s))
    
    sol = []
    
    for i in range(9):
        t = np.zeros(i+2)
        t[0]+=1
        t = np.poly1d(t)
        a = np.polymul(sp,t)
        
        #find positions of ones
        pos = np.where(np.flip(np.array(a))==1)[0]
        
        polynomials = []
        
        poly = np.zeros(6)
        
        for j in range(len(pos)):
            if pos[j]>=6:
                polynomials.append(np.poly1d(np.flip(H[pos[j]-6])))
            else:
                poly[pos[j]]+=1
        
        polynomials.append(np.poly1d(np.flip(poly)))
        
        si = np.flip(np.array(sum(polynomials))%2)
        
        
        b = max(np.where(si==1)[0])-min(np.where(si==1)[0])+1
        
        sol.append([si,b])
        
        if b<=3:
            si = [si, i] #si, i, dimension
            break
    
    kt = np.zeros(9-i)
    kt[0]+=1
    kt = np.poly1d(kt)
    e = np.polymul(kt,np.flip(si[0]))
    extend = np.zeros(len(w)-len(e)-1)
    e = np.hstack((e,extend))
    code = (w+e)%2
    syndrome = (code@H)%2
    return code, sol

a = cyclicBurst(w)

print(a[0])
print(a[1])
    
        
        
        
        
    