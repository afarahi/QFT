#! /usr/bin/python

from numpy             import *
from pylab             import *
from matplotlib.pyplot import *
from matplotlib        import rc
from cmath             import exp, cos, sin

#
#rc('text', usetex = True)
#rc('font', family = 'serif')

#Constants
kx  = 0.0
kt  = 1.0j
Lm2 = -3.0
A   = 1.0
B   = 0.0
d   = 4

# write hand constraction
def RH_cons(w,f,fp):
    fpp = -( 1.0-d-w**d)*fp/(w*(1.0-w**d) ) - ( kt**2/((1.0-w**d)**2) - kx**2/(1.0-w**d) - Lm2/(w**2*(1.0-w**d)) )*f
#    fpp = -( 1.0-d-w**d)*fp/(w*(1.0-w**d) ) + (kt**2/((1.0-w**d)**2) + kx**2/(1.0-w**d) + Lm2/(w**2*(1.0-w**d)) )*f
    return fpp

#Boundary condition 
def appl_bondary(h):
    w = 1.0+h
    f = A*exp(-1.0j*kt*log(1.0-w)/d) + B*exp(1.0j*kt*log(1.0-w)/d)
    fp= A*exp(-1.0j*kt*log(1.0-w)/d)*1.0j*kt/(d*(1.0-w)) - B*exp(1.0j*kt*log(1.0-w)/d)*1.0j*kt/(d*(1.0-w))
#    f = A*exp(kt*log(1-w)/d) + B*exp(-kt*log(1-w)/d)
#    fp= - A*exp(kt*log(1.0-w)/d)*kt/(d*(1.0-w)) + B*exp(kt*log(1.0-w)/d)*kt/(d*(1.0-w))
    return (f,fp)

#Runge-Kutta coefficient
def RK_cons(w,f,fp,h):
#    k1 = 0.5*h**2*RH_cons( w , f , fp )
#    k2 = 0.5*h**2*RH_cons( w+0.5*h , f+0.5*h*fp+0.25*k1 , fp+k1/h)
#    k3 = 0.5*h**2*RH_cons( w+0.5*h , f+0.5*h*fp+0.25*k1 , fp+k2/h)
#    k4 = 0.5*h**2*RH_cons( w+h , f+h*fp+k3 , fp+2.0*k3/h)
    l1 = h*RH_cons( w , f , fp )
    k1 = h*fp
    l2 = h*RH_cons( w+0.5*h , f+0.5*k1 , fp+0.5*l1)
    k2 = h*(fp+0.5*l1)
    l3 = h*RH_cons( w+0.5*h , f+0.5*k2 , fp+0.5*l2)
    k3 = h*(fp+0.5*l2)
    l4 = h*RH_cons( w+h , f+k3 , fp+l3) 
    k4 = h*(fp+0.5*l3)
    return (k1,k2,k3,k4,l1,l2,l3,l4)

#Left hand constraction
def LH_cons(w,f,fp,h):
#    (k1,k2,k3,k4) = RK_cons(w,f,fp,h)
#    fnew  = f + h*fp + (k1+k2+k3)/3.0
#    fpnew = fp + (k1+2.0*k2+2.0*k3+k4)/(2.0*h)
    (k1,k2,k3,k4,l1,l2,l3,l4) = RK_cons(w,f,fp,h)
    fnew  = f  + (k1+2.0*k2+2.0*k3+k4)/6.0
    fpnew = fp + (l1+2.0*l2+2.0*l3+l4)/6.0
    return (fnew,fpnew)

# w is changing from zero to one w = [0 1[
# w is not defined in w = 1 so we need to use somewhere close to w = 1
# so we would start from point w = 1 - h

#Apply Boundary condition
h      = -1e-7
(f,fp) = appl_bondary(h)
fw     = [f]
w      = [1.0+h]

#Mesh Grid Size
h = -1e-7
wp     = w[0]+h

#Solver
while (wp>0):
   w.append(wp)
   (fnew,fpnew) = LH_cons(wp,f,fp,h)
   fw.append(fnew)
   wp = wp+h
   f  = fnew
   fp = fpnew
w.pop()
fw.pop()

ff = zeros(len(fw))
for i in range(0, len(fw)):
   ff[i] = fw[i].real
print ff

figure(1)
plot(w,ff,'b',label='1e-3')

def polyFit(xList,yList,order=1):
    '''fit the data using a least squares and polynomial'''
#    fList = [(lambda x,n=n: x**n) for n in range(order,-1,-1)]
    fList = [(lambda x,n=n: x**n) for n in range(order,-1,-1)]
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

xList = w[len(w)-400:len(w)]
yList = ff[len(w)-400:len(w)]     
a = polyFit(xList,yList,order=3)
aHat = a[0]
bHat = a[1]
cHat = a[2]
dHat = a[3]

x = linspace(0,0.5,50)
y = aHat*x**3 + bHat*x**2 + cHat*x + dHat

plot(x,y,'r')

print 'Fit equation is :%4.5f*x^3 + %4.5f*x^2 + %4.5f*x + %4.5f' % (aHat,bHat,cHat,dHat)

#xlabel(r'\textbf{$w$}')
#ylabel(r'\textbf{$f(w)$}')
#title(r'\textbf{$L^2m^2 = 0$}')
#legend(loc=2)
show()
#savefig('plot_00.pdf')

