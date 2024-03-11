import numpy as np


J=np.exp(2j*np.pi/3)
Jc=1/J

def roots2(a,b,c):
    bp=b/2    
    delta=bp*bp-a*c
    u1=(-bp-delta**.5)/a
    u2=-u1-b/a
    return u1,u2  


def cardan(a,b,c,d):
    z0=b/3/a
    a2,b2 = a*a,b*b    
    p=-b2/3/a2 +c/a
    q=(b/27*(2*b2/a2-9*c/a)+d)/a
    D=-4*p*p*p-27*q*q
    r=np.sqrt(-D/27+0j)        
    u=((-q-r)/2)**0.33333333333333333333333
    v=((-q+r)/2)**0.33333333333333333333333
    w=u*v
    w0=np.abs(w+p/3)
    w1=np.abs(w*J+p/3)
    w2=np.abs(w*Jc+p/3)
    v = np.where(
        np.logical_and(np.lt(w0,w1), np.lt(w2,w0)),
        v*Jc, v
    )
    v = np.where(
        np.logical_and(np.lt(w1,w0), np.lt(w2,w1)),
        v*Jc, v
    )
    v = np.where(
        np.logical_and(np.lt(w1,w0), np.lt(w1,w2)),
        v*J, v
    )    
    return u+v-z0, u*J+v*Jc-z0,u*Jc+v*J-z0

def ferrari(a,b,c,d,e):
    "resolution of P=ax^4+bx^3+cx^2+dx+e=0"
    "CN all coeffs real."
    "First shift : x= z-b/4/a  =>  P=z^4+pz^2+qz+r"
    z0=b/4/a
    a2,b2,c2,d2 = a*a,b*b,c*c,d*d 
    p = -3*b2/(8*a2)+c/a
    q = b*b2/8/a/a2 - 1/2*b*c/a2 + d/a
    r = -3/256*b2*b2/a2/a2 +c*b2/a2/a/16-b*d/a2/4+e/a
    "Second find X so P2=AX^3+BX^2+C^X+D=0"
    A=8
    B=-4*p
    C=-8*r
    D=4*r*p-q*q
    y0,y1,y2=cardan(A,B,C,D)
    y0 = np.where(
        np.lt(np.abs(y1.imag),np.abs(y0.imag)),
        y1, y0
    )
    y0 = np.where(
        np.lt(np.abs(y2.imag),np.abs(y0.imag)),
        y2, y0
    )
    a0=(-p+2*y0.real)**.5
    b0 = np.where(
        np.eq(a0,0),
        y0**2-r, -q/2/a0
    )
    r0,r1=roots2(1,a0,y0+b0)
    r2,r3=roots2(1,-a0,y0-b0)
    answer = np.stack((r0-z0,r1-z0,r2-z0,r3-z0))
    return answer




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