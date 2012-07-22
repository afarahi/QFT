from numpy             import *
from cmath             import *
from scipy.special     import *
#from scipy.optimize    import curve_fit

from pylab             import *
from matplotlib.pyplot import *
from matplotlib        import rc

h = 1e-8

def J_v(v,x,k,d):
    J_v = x**(d/2.0)*jv(v,k*x)
    return J_v

def Y_v(v,x,k,d):
    Y_v = x**(d/2.0)*yv(v,k*x)
    return Y_v

def dJ_v(v,x,k,d):
    dJ_v = (J_v(v,x+h,k,d) - J_v(v,x-h,k,d)) / (2.0*h) 
    return dJ_v

def dY_v(v,x,k,d):
    dY_v = (Y_v(v,x+h,k,d) - Y_v(v,x-h,k,d)) / (2.0*h) 
    return dY_v

def df1(x,n):
    df_1 = n*x**(n-1) 
    return df_1

def df2(x,n):
    df_2 = n*x**(n-1) 
    return df_2

def f1(x,n,R2,kx,kt,Rm2,muq):
    f_1 = x**n*(1.0 - R2**2*((kt+muq)**2 - kx**2)/(x**2*(n**2-4-Rm2)))
    return f_1

def f2(x,n,R2,kx,kt,Rm2,muq):
    f_2 = x**n*(1.0 - R2**2*((kt+muq)**2 - kx**2)/(x**2*(n**2-4-Rm2)))
    return f_2

def f1n(x,n,R2,kx,kt,Rm2,muq):
    f_1 = x**n
    return f_1

def f2n(x,n,R2,kx,kt,Rm2,muq):
    f_2 = x**n*log(x)
    return f_2

def f12(x,n,R2,kx,kt,Rm2,muq):
    f_1 = x**n
    return f_1

def f22(x,n,R2,kx,kt,Rm2,muq):
    f_2 = x**(n-2)*(1 + R2**2*(kx**2-(kt+muq)**2)*log(x)/(2.0*n) )
    return f_2

def func(x, n, a, b):
    return a*f1(x,n[0]) + b*f2(x,n[1])

def Gauss_Newton_CF(x,y,pow_r,a,b):
    Jr   = zeros([len(x),2],complex)
    Jrt  = zeros([2,len(x)],complex)
    for i in range(0,len(x)):
       Jr[i,0] = df1(x[i],pow_r[0])
       Jr[i,1] = df2(x[i],pow_r[1])
       Jrt[0,i]= Jr[i,0]
       Jrt[1,i]= Jr[i,1]
    J    = zeros([2,2],complex)
    for i in range(0,len(x)):
       J[0,0]  = J[0,0] + Jrt[0,i]*Jr[i,0]
       J[0,1]  = J[0,1] + Jrt[0,i]*Jr[i,1]
       J[1,0]  = J[1,0] + Jrt[1,i]*Jr[i,0]
       J[1,1]  = J[1,1] + Jrt[1,i]*Jr[i,1]
    detJ = J[0,0]*J[1,1] - J[0,1]*J[1,0]
    invJ = zeros([2,2],complex)
    invJ[0,0]  = J[1,1]/detJ
    invJ[0,1]  = -J[0,1]/detJ
    invJ[1,0]  = -J[1,0]/detJ
    invJ[1,1]  = J[0,0]/detJ
    Jy   = zeros(2,complex)
    for i in range(0,len(x)):
       Jy[0]   = Jy[0] + Jr[i,0]*(func(x[i], pow_r, a, b) - y[i])
       Jy[1]   = Jy[1] + Jr[i,1]*(func(x[i], pow_r, a, b) - y[i])
    da   = invJ[0,0]*Jy[0] + invJ[0,1]*Jy[1]
    db   = invJ[1,0]*Jy[0] + invJ[1,1]*Jy[1]
    anew = - da + a
    bnew = - db + b
    return [anew , bnew]

def Gauss_Newton_Solver(x,y,pow_r,a,b):
    n = 100
    delta = 1e-6
    i = 0
    deltap = 1
    while (i < n and deltap > delta ):
        [anew , bnew] = Gauss_Newton_CF(x,y,pow_r,a,b)
        deltap        = sqrt((anew-a)**2 + (bnew-b)**2) / sqrt(a**2+b**2)
        a  = anew
        b  = bnew
        i += 1
    print deltap
    return [a,b]

def two_point_curve_fit(x,y,n,R2,kx,kt,Rm2,rmax,Q,M,muq):
    if ( n[1]!=n[0] and abs(n[0]-n[1])!= 2):
       det = 1.0 / (f1(x[0],n[0],R2,kx,kt,Rm2,muq)*f2(x[1],n[1],R2,kx,kt,Rm2,muq)-f2(x[0],n[1],R2,kx,kt,Rm2,muq)*f1(x[1],n[0],R2,kx,kt,Rm2,muq))
       a = det*(f2(x[1],n[1],R2,kx,kt,Rm2,muq)*y[0] - f2(x[0],n[1],R2,kx,kt,Rm2,muq)*y[1])
       b = det*(-f1(x[1],n[0],R2,kx,kt,Rm2,muq)*y[0]+ f1(x[0],n[0],R2,kx,kt,Rm2,muq)*y[1])
    elif ( n[1] == n[0]):
       det = 1.0 / (f1n(x[0],n[0],R2,kx,kt,Rm2,muq)*f2n(x[1],n[1],R2,kx,kt,Rm2,muq)-f2n(x[0],n[1],R2,kx,kt,Rm2,muq)*f1n(x[1],n[0],R2,kx,kt,Rm2,muq))
       a = det*(f2n(x[1],n[1],R2,kx,kt,Rm2,muq)*y[0] - f2n(x[0],n[1],R2,kx,kt,Rm2,muq)*y[1])
       b = det*(-f1n(x[1],n[0],R2,kx,kt,Rm2,muq)*y[0]+ f1n(x[0],n[0],R2,kx,kt,Rm2,muq)*y[1])
    else:
       np = max(abs(n[0]),abs(n[1]))
       det = 1.0 / (f12(x[0],np,R2,kx,kt,Rm2,muq)*f22(x[1],np,R2,kx,kt,Rm2,muq)-f22(x[0],np,R2,kx,kt,Rm2,muq)*f12(x[1],np,R2,kx,kt,Rm2,muq))
       a = det*(f22(x[1],np,R2,kx,kt,Rm2,muq)*y[0] - f22(x[0],np,R2,kx,kt,Rm2,muq)*y[1])
       b = det*(-f12(x[1],np,R2,kx,kt,Rm2,muq)*y[0]+ f12(x[0],np,R2,kx,kt,Rm2,muq)*y[1])
    return [a,b]

def polyFit(xList,yList,A):
#   '''fit the data using a least squares and polynomial'''
#    fList = [(lambda x,n=n: x**n) for n in range(order,-1,-4)]
    fList = [(lambda x,n=n: x**n) for n in A]
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

