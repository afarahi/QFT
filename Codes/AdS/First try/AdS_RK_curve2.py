#! /usr/bin/python

from numpy             import *
from pylab             import *
from matplotlib.pyplot import *
from matplotlib        import rc
from cmath             import exp, cos, sin
from curvfit           import *

#
#rc('text', usetex = True)
#rc('font', family = 'serif')

#Constants
kx  = 0.0
kt  = 1.0j
k   = 1.0
Lm2 = -4.0
A   = 1.0
B   = 0.0
d   = 4
v   = sqrt(Lm2 + d**2/4.0)
h1  = -1e-5
h2  = -1e-5
mm  = 10000

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
    (k1,k2,k3,k4,l1,l2,l3,l4) = RK_cons(w,f,fp,h)
    fnew  = f  + (k1+2.0*k2+2.0*k3+k4)/6.0
    fpnew = fp + (l1+2.0*l2+2.0*l3+l4)/6.0
    return (fnew,fpnew)

# w is changing from zero to one w = [0 1[
# w is not defined in w = 1 so we need to use somewhere close to w = 1
# so we would start from point w = 1 - h

#Apply Boundary condition
h      = h1
(f,fp) = appl_bondary(h)
fw     = [f]
w      = [1.0+h]

#Mesh Grid Size
h = h2
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
#print ff

figure(1)
plot(w,ff,'b',label='1e-3')


#a = GN_solver(w[len(w)-mm:len(w)],ff[len(w)-mm:len(w)],[-1,1])
#x = linspace(0,0.5,50)
#y = a[0]*x**3 + a[1]*x
#plot(x,y,'m')
#print 'Fit equation is :%4.5f*x^3 + %4.5f*x' % (a[0],a[1])


xList = w[len(w)-mm:len(w)]
yList = ff[len(w)-mm:len(w)]     
#a = polyFit(xList,yList,order=4)
#a = polyFit(xList,yList,v,k,d)
a = 1
b = -1
[a,b] = Gauss_Newton_Solver(xList,yList,v,kt,d,a,b)

aHat = a
bHat = b
#cHat = a[2]
#dHat = a[3]
#eHat = a[4]

x = linspace(0,0.2,50)
#y = aHat*x**3 + bHat*x**2 + cHat*x + dHat
#y = aHat*x**4 + bHat#*x
#y = aHat*x**(d/2.0)*J_v(n*v,k*x) + bHat*x**(d/2.0)*Y_v(n*v,kt*x)
y = aHat*J_v(v,x,kt,d) + bHat*Y_v(v,x,kt,d)

plot(x,y,'r')
show()
#print 'Fit equation is :%4.5f*x^3 + %4.5f*x^2 + %4.5f*x + %4.5f' % (aHat,bHat,cHat,dHat)
#print 'Fit equation is :%4.5f*x^4 + %4.5f*x^3 + %4.5f*x^2 +  %4.5f*x + %4.5f' % (aHat,bHat,cHat,dHat,eHat)
print 'Fit equation is :%4.5f*J_v + %4.5f*Y_v' % (aHat,bHat)


#xlabel(r'\textbf{$w$}')
#ylabel(r'\textbf{$f(w)$}')
#title(r'\textbf{$L^2m^2 = 0$}')
#legend(loc=2)
#show()
#savefig('plot_00.pdf')

#A    = [[ w[len(w)-mm]**3   , w[len(w)-mm] ]
#       ,[ w[len(w)-mm-5]**3 , w[len(w)-mm-5] ]]

#Y    = [[ ff[len(w)-mm] ]
#       ,[ ff[len(w)-mm-5] ]]

#detA = w[len(w)-mm]**3*w[len(w)-mm+10] - w[len(w)-mm]*w[len(w)-mm+10]**3
#x1   = 1.0/detA * (w[len(w)-mm+10]*ff[len(w)-mm] - w[len(w)-mm]*ff[len(w)-mm+10] )
#x2   = 1.0/detA * ( - w[len(w)-mm+10]**3*ff[len(w)-mm] + w[len(w)-mm]**3*ff[len(w)-mm+10] )
#print 'Fit equation is :%4.5f*x^4 + %4.5f ' % (x1,x2)
