from numpy             import *
from cmath             import *
from scipy.special     import *
#from scipy.optimize    import curve_fit

from pylab             import *
from matplotlib.pyplot import *
from matplotlib        import rc

h = 1e-8

def two_point_curve_fit(r,phip,phin,L,kx,kt,ms,muq):
#     B = ( L*kt*(1.0+muq) - L*kx)/(2.0*ms+1.0)*C
#     D = ( L*kt*(1.0+muq) + L*kx)/(2.0*ms-1.0)*A
#     phip = A*r**(m*L)  + B*r**(-m*L-1.0)
#     phin = C*r**(-m*L) + D*r**(m*L-1.0)

#     phip = A*r**(m*L)  + C*( L*kt*(1.0+muq) - L*kx)/(2.0*ms+1.0)*r**(-m*L-1.0)
#     phin = A*( L*kt*(1.0+muq) + L*kx)/(2.0*ms-1.0)*r**(m*L-1.0) + C*r**(-m*L)
    
    a11 = r**(ms*L)
    a12 = ( L*kt*(1.0+muq) - L*kx) / (2.0*ms+1.0) * r**(-ms*L-1.0)
    a21 = ( L*kt*(1.0+muq) + L*kx) / (2.0*ms-1.0) * r**(ms*L-1.0)
    a22 = r**(-ms*L)

    det = 1.0 / (a11*a22-a12*a21)
    A = det*( a22*phip - a12*phin)
    C = det*(-a21*phip + a11*phin)
     
#    if ( n[1]!=n[0] and abs(n[0]-n[1])!= 2):
#       det = 1.0 / (f1(x[0],n[0],R2,kx,kt,Rm2,muq)*f2(x[1],n[1],R2,kx,kt,Rm2,muq)-f2(x[0],n[1],R2,kx,kt,Rm2,muq)*f1(x[1],n[0],R2,kx,kt,Rm2,muq))
#       a = det*(f2(x[1],n[1],R2,kx,kt,Rm2,muq)*y[0] - f2(x[0],n[1],R2,kx,kt,Rm2,muq)*y[1])
#       b = det*(-f1(x[1],n[0],R2,kx,kt,Rm2,muq)*y[0]+ f1(x[0],n[0],R2,kx,kt,Rm2,muq)*y[1])
#    elif ( n[1] == n[0]):
#       det = 1.0 / (f1n(x[0],n[0],R2,kx,kt,Rm2,muq)*f2n(x[1],n[1],R2,kx,kt,Rm2,muq)-f2n(x[0],n[1],R2,kx,kt,Rm2,muq)*f1n(x[1],n[0],R2,kx,kt,Rm2,muq))
#       a = det*(f2n(x[1],n[1],R2,kx,kt,Rm2,muq)*y[0] - f2n(x[0],n[1],R2,kx,kt,Rm2,muq)*y[1])
#       b = det*(-f1n(x[1],n[0],R2,kx,kt,Rm2,muq)*y[0]+ f1n(x[0],n[0],R2,kx,kt,Rm2,muq)*y[1])
#    else:
#       np = max(abs(n[0]),abs(n[1]))
#       det = 1.0 / (f12(x[0],np,R2,kx,kt,Rm2,muq)*f22(x[1],np,R2,kx,kt,Rm2,muq)-f22(x[0],np,R2,kx,kt,Rm2,muq)*f12(x[1],np,R2,kx,kt,Rm2,muq))
#       a = det*(f22(x[1],np,R2,kx,kt,Rm2,muq)*y[0] - f22(x[0],np,R2,kx,kt,Rm2,muq)*y[1])
#       b = det*(-f12(x[1],np,R2,kx,kt,Rm2,muq)*y[0]+ f12(x[0],np,R2,kx,kt,Rm2,muq)*y[1])
    
    return [A,C]

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

