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


# Algorithm from
# https://quarticequations.com/
# x^3 + Px + Q = 0
def cardan_depressed(P,Q):
    q = P/3
    r = Q/(-2)
    A = np.sqrt(r**2 + q**3)
    u = power(A+r, 1/3)
    v = power(A-r, 1/3)
    x0 = u-v
    return x0


# Algorithm from
# https://quarticequations.com/
# x^3 + Px + Q = 0
def stable_cubic(P,Q):
    q = P/3
    r = Q/(-2)
    A = power((np.abs(r) + np.sqrt(r**2 + q**3)), 1/3)
    x0 = (A-q/A)*np.sign(r)
    return x0

    
    