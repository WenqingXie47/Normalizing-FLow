import torch



def power(x, index):
    y = torch.sign(x)*torch.exp(index*torch.log(torch.abs(x)))
    return y


def cardan_depressed(p,q):

    D= (q**2)/4 + (p**3)/27
    r=torch.sqrt(D)        
    u=(-q/2+r)
    u = power(u,1/3)
    v=(-q/2-r)
    v = power(v,1/3)
    return (u+v)


def depressed_cubic(x,g):
    y = (1-g)*(x**3) + g*x
    return y


def depressed_cubic_gradient(x,g):
    grad = (1-g)*3*(x**2) + g
    return grad


def inverse_depressed_cubic(y,g):
    x = torch.where(
        torch.isclose(g,1),
        y/g,  cardan_depressed(p=g/(1-g),q=(-y)/(1-g))
    )
    return x

def inverse_depressed_cubic_gradient(y,g):
    e = 0.001
    x = inverse_depressed_cubic(y,g)
    grad_x = depressed_cubic_gradient(x,g)
    grad_y = 1/(grad_x+e)
    return grad_y


def forward_layer(x,a,b,g):
    x = a*(x+b)
    x = depressed_cubic(x,g)
    return x


def get_cd(a,b,g):
    max = forward_layer(1,a,b,g)
    min = forward_layer(-1,a,b,g)
    c = 2/(max-min)
    d = -min*c-1
    return c,d


def normalized_forward_layer(x,a,b,g):
    x = forward_layer(x,a,b,g)
    c,d = get_cd(a,b,g)
    x = x*c+d
    return x


def backward_layer(x,a,b,g):
    c,d = get_cd(a,b,g)
    x = (x-d)/c
    y = inverse_depressed_cubic(y)
    y = (y-b[n_layers-i-1])/a[n_layers-i-1]

def forward_function(x,a,b,c,d,n_layers):
    for i in range(n_layers):   
        x = normalized_forward_layer(x,a[i],b[i],c[i],d[i])
    return x



def backward_function(y,a,b,c,d,n_layers):
    for i in range(n_layers):
        y = (y-d[n_layers-i-1])/c[n_layers-i-1]
        y = inverse_non_linear(y)
        y = (y-b[n_layers-i-1])/a[n_layers-i-1]
    return y

    