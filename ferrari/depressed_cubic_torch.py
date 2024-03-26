import torch

def power(x, index):
    y = torch.sign(x)*torch.exp(index*torch.log(torch.abs(x)))
    return y

# def cardan_depressed(p,q):

#     D= (q**2)/4 + (p**3)/27
#     r=torch.sqrt(D)        
#     u=(-q/2+r)
#     u = power(u,1/3)
#     v=(-q/2-r)
#     v = power(v,1/3)
#     return u+v


# Algorithm from
# https://quarticequations.com/
# x^3 + Px + Q = 0
def cardan_depressed(P,Q):
    q = P/3
    r = Q/(-2)
    A = torch.sqrt(r**2 + q**3)
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
    A = power((torch.abs(r) + torch.sqrt(r**2 + q**3)), 1/3)
    x0 = (A-q/A)*torch.sign(r)
    return x0

    