from numpy             import *
from cmath             import *
from scipy.special     import *
#from scipy.optimize    import curve_fit

from pylab             import *
from matplotlib.pyplot import *
from matplotlib        import rc

h = 1e-8

def two_point_curve_fit_G22(r,phip,phin,L,kx,kt,ms,muq):
# For Finding G_{22}

# For Field with non-Zero Mass
#     B = L^2*( kt+muq - kx)/(2.0*ms+1.0)*C
#     D = L^2*( kt+muq + kx)/(2.0*ms-1.0)*A
#     phip = A*r**(m*L)  + B*r**(-m*L-1.0)
#     phin = C*r**(-m*L) + D*r**(m*L-1.0)

# For Field with Zero Mass
#     A1 = L^2*( kt+muq - kx)*C
#     C1 =-L^2*( kt+muq + kx)*A
#     A2 = L^2*( kt+muq)*C1
#     C2 =-L^2*( kt+muq)*A2
#     phip = A + A1*r^(-1) + A2*r^(-2)
#     phin = C + C1*r^(-1) + C2*r^(-2)

    if (ms != 0.0):
       a11 = r**(ms*L)
       a12 = L**2*( kt+muq - kx) / (2.0*ms+1.0) * r**(-ms*L-1.0)
       a21 = L**2*( kt+muq + kx) / (2.0*ms-1.0) * r**(ms*L-1.0)
       a22 = r**(-ms*L)

       det = 1.0 / (a11*a22-a12*a21)
       A = det*( a22*phip - a12*phin)
       C = det*(-a21*phip + a11*phin)

    else:
       A1 = L**2*( kt+muq - kx)
       C1 =-L**2*( kt+muq + kx)
       A2 = L**2*( kt+muq)*C1
       C2 =-L**2*( kt+muq)*A2

       a11 = 1.0
       a12 = A1/r + A2/(r**2)
       a21 = C1/r + C2/(r**2)
       a22 = 1.0

       det = 1.0 / (a11*a22-a12*a21)
       A = det*( a22*phip - a12*phin)
       C = det*(-a21*phip + a11*phin)
     
    return [A,C]


def two_point_curve_fit_G11(r,phip,phin,L,kx,kt,ms,muq):
# For Finding G_{11}

# For Field with non-Zero Mass
#     B = L^2*( kt+muq + kx)/(2.0*ms+1.0)*C
#     D = L^2*( kt+muq - kx)/(2.0*ms-1.0)*A
#     phip = A*r**(m*L)  + B*r**(-m*L-1.0)
#     phin = C*r**(-m*L) + D*r**(m*L-1.0)

# For Field with Zero Mass
#     A1 = L^2*( kt+muq + kx)*C
#     C1 =-L^2*( kt+muq - kx)*A
#     A2 = L^2*( kt+muq)*C1
#     C2 =-L^2*( kt+muq)*A2
#     phip = A + A1*r^(-1) + A2*r^(-2)
#     phin = C + C1*r^(-1) + C2*r^(-2)

    if (ms != 0.0):
       a11 = r**(ms*L)
       a12 = L**2*( kt+muq + kx) / (2.0*ms+1.0) * r**(-ms*L-1.0)
       a21 = L**2*( kt+muq - kx) / (2.0*ms-1.0) * r**(ms*L-1.0)
       a22 = r**(-ms*L)

       det = 1.0 / (a11*a22-a12*a21)
       A = det*( a22*phip - a12*phin)
       C = det*(-a21*phip + a11*phin)

    else:
       A1 = L**2*( kt+muq + kx)
       C1 =-L**2*( kt+muq - kx)
       A2 = L**2*( kt+muq)*C1
       C2 =-L**2*( kt+muq)*A2

       a11 = 1.0
       a12 = A1/r + A2/(r**2)
       a21 = C1/r + C2/(r**2)
       a22 = 1.0

       det = 1.0 / (a11*a22-a12*a21)
       A = det*( a22*phip - a12*phin)
       C = det*(-a21*phip + a11*phin)
     
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

