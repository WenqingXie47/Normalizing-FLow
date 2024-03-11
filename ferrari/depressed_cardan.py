import numpy as np

def power(x, index):
    y = np.sign(x)*np.exp(index*np.log(np.abs(x)))
    return y

def cardan_depressed(p,q):

    D= (q**2)/4 + (p**3)/27
    r=np.sqrt(D)        
    u=(-q/2+r)
    u = power(u,1/3)
    v=(-q/2-r)
    v = power(v,1/3)
    return u+v