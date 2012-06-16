from numpy             import *
from pylab             import *
from cmath             import exp, cos, sin

def polyFit(xList,yList,order=1):
    '''fit the data using a least squares and polynomial'''
#    fList = [(lambda x,n=n: x**n) for n in range(order,-1,-4)]
    fList = [(lambda x,n=n: x**n) for n in (4,0)]
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


def GN_cons(X,u,v):
    f = zeros([len(u),1])
    for i in range(0,len(u)):
       f[i,0] = X[0]*u[i]**3 + X[1]*u[i] - v[i]
    J = zeros([len(u),2])
    for i in range(0,len(u)):
       J[i,0] = 3.0*X[0]*u[i]**2
       J[i,1] = X[1]
    return (f,J)

def GN_solver(u,v,X):
    err = 1e-4
    maxn= 100
    n   = 0
    Aerr=10
    Xnew  = zeros(2)
    while (Aerr > err and n < maxn):
       (fList,JList) = GN_cons(X,u,v)
       f     = matrix(fList)
       J     = matrix(JList)
       delta = -inv(J.T*J)*(J.T*f)
       Xnew[0] = X[0] + delta[0]
       Xnew[1] = X[1] + delta[1]
       Aerr  = sqrt( delta[0]**2+delta[1]**2 )/sqrt( X[0]**2+X[1]**2 )
       X[0]  = Xnew[0]
       X[1]  = Xnew[1]
       n     = n+1
    return X 

