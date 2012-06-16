from numpy             import *
from cmath             import exp, cos, sin
from scipy.special     import *
#from scipy.optimize    import curve_fit


from pylab             import *
from matplotlib.pyplot import *
from matplotlib        import rc

h = 1e-8

def J_v(v,x,k,d):
    J_v = x**(d/2.0)*jv(v,k*x)
#    J_v  = x
    return J_v

def Y_v(v,x,k,d):
    Y_v = x**(d/2.0)*yv(v,k*x)
#    Y_v  = x**3
    return Y_v

def dJ_v(v,x,k,d):
    dJ_v = (J_v(v,x+h,k,d) - J_v(v,x-h,k,d)) / (2.0*h) 
    return dJ_v

def dY_v(v,x,k,d):
    dY_v = (Y_v(v,x+h,k,d) - Y_v(v,x-h,k,d)) / (2.0*h) 
    return dY_v

def func(x, v, k, d, a, b):
    return a*J_v(v,x,k,d) + b*Y_v(v,x,k,d)

def fun(x,a,b):
    return a*x**3 + b*x

#x = linspace(0.1,1,50)
#	y = fun(x, 2.5, 1.3)
#yn = y + 0.2*random(size=len(x))

#popt, pcov = curve_fit(func, x, yn)
#print popt, pcov

def Gauss_Newton_CF(x,y,v,k,d,a,b):
    Jr   = zeros([len(x),2])
    Jrt  = zeros([2,len(x)])
    for i in range(0,len(x)):
       Jr[i,0] = dJ_v(v,x[i],k,d)
       Jr[i,1] = dY_v(v,x[i],k,d)
       Jrt[0,i]= Jr[i,0]
       Jrt[1,i]= Jr[i,1]
    J    = zeros([2,2])
    for i in range(0,len(x)):
       J[0,0]  = J[0,0] + Jrt[0,i]*Jr[i,0]
       J[0,1]  = J[0,1] + Jrt[0,i]*Jr[i,1]
       J[1,0]  = J[1,0] + Jrt[1,i]*Jr[i,0]
       J[1,1]  = J[1,1] + Jrt[1,i]*Jr[i,1]
    detJ = J[0,0]*J[1,1] - J[0,1]*J[1,0]
    invJ = zeros([2,2])
    invJ[0,0]  = J[1,1]/detJ
    invJ[0,1]  = -J[0,1]/detJ
    invJ[1,0]  = -J[1,0]/detJ
    invJ[1,1]  = J[0,0]/detJ
    Jy   = zeros(2)
    for i in range(0,len(x)):
       Jy[0]   = Jy[0] + Jr[i,0]*(func(x[i], v, k, d, a, b) - y[i])
       Jy[1]   = Jy[1] + Jr[i,1]*(func(x[i], v, k, d, a, b) - y[i])
    da   = invJ[0,0]*Jy[0] + invJ[0,1]*Jy[1]
    db   = invJ[1,0]*Jy[0] + invJ[1,1]*Jy[1]
    anew = - da + a
    bnew = - db + b
    return [anew , bnew]


def Gauss_Newton_Solver(x,y,v,k,d,a,b):
    n = 100
    delta = 1e-4
    i = 0
    deltap = 1
    while (i < n and deltap > delta ):
        [anew , bnew] = Gauss_Newton_CF(x,y,v,k,d,a,b)
        deltap        = sqrt((anew-a)**2 + (bnew-b)**2) / sqrt(a**2+b**2)
        a  = anew
        b  = bnew
        i += 1
    print deltap
    return [a,b]

'''
a = 1
b = 1
kt  = 1.0j
k   = 1.0
Lm2 = -4.0
A   = 1.0
B   = 0.0
d   = 4
v   = sqrt(Lm2 + d**2/4.0)

[a,b] = Gauss_Newton_Solver(x,y,v,kt,d,a,b)

yp = func(x, v, kt, d, a, b)
figure(1)
plot(x,y,'b')
plot(x,yp,'r')
show()
'''
#def polyFit(xList,yList,v,k,d):
#   '''fit the data using a least squares and polynomial'''
#    fList = [(lambda x,n=n: x**n) for n in range(order,-1,-4)]
#    fList = [(lambda x,n=n: x**n) for n in (1,-1)]
'''   fList = [(lambda x,n=n: x**(d/2.0)*J_v(n*v,k*x)) for n in (1,-1)]
    # build row for each element in y
    bList = []
    A_List = []
    for (thisX,thisY) in zip(xList,yList):
        bList.append(thisY)
        A_Row = [f(thisX) for f in fList]
        A_List.append(A_Row)
    b = matrix(bList).T
    A = matrix(A_List)
    w = inv(A.T*A)*A.T*b
    return w.T.tolist()[0]

def f1(x,d,v,k):
    y = x**(d/2.0)*J_v(v,kx)

def f2(x,d,v,k):
    y = x**(d/2.0)*Y_v(v,kx)
'''
