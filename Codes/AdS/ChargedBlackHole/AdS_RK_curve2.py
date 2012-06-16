#! /usr/bin/python

from numpy             import *
from pylab             import *
from matplotlib.pyplot import *
from matplotlib        import rc
from cmath             import exp, cos, sin
from curvfit           import *

#Constants
kx  = 0.0
kt  = 1.0j
k   = 1.0
Rm2 = -4.0
R2  = 1.0
d   = 4
v   = sqrt(Rm2 + d**2/4.0)

rs  = 0.25
r0  = 1.0
rmax= 10*r0
h1  = 1.0*r0
h2  = 0.01*r0
h3  = 1e-3*r0
mm  = 10

Q   = sqrt(d/(d-2))*rs**(d-1)
M   = 	r0**d + Q**2/(r0**(d-2))

# write hand constraction
def redshift_f(r):
    f = 1.0 + Q**2/(r**(2*d-2)) - M/r**d
    return f

def RH_cons(r,phi,phip):
#    fpp = -( 1.0-d-w**d)*fp/(w*(1.0-w**d) ) - ( kt**2/((1.0-w**d)**2) - kx**2/(1.0-w**d) - Lm2/(w**2*(1.0-w**d)) )*f
    phipp = ( \
             - ( (d+1)*r - (d-3)*Q**2/r**(2*d-3) - M/r**(d-1) ) * fp \
             - ( - R2*kt**2/(r**2*redshift_f(r)) + kx**2*R2/(r**2) - Rm2 )*f \
            ) / (r**2*redshift_f(r))
    return phipp

#Boundary condition 
def appl_bondary(rstart,h):
    r   = rstart+h
    beta= R2*kt/( d*(1.0 - (rs/r0)**(2*d-2)) )
    phi = exp(-1.0j*beta*log(r-r0))
    phip= -exp(-1.0j*beta*log(r-r0))*1.0j*beta/(r-r0)
    return (phi,phip)

#Runge-Kutta coefficient
def RK_cons(r,phi,phip,h):
    l1 = h*RH_cons( r , phi , phip )
    k1 = h*phip
    l2 = h*RH_cons( r+0.5*h , phi+0.5*k1 , phip+0.5*l1 )
    k2 = h*(phip+0.5*l1)
    l3 = h*RH_cons( r+0.5*h , phi+0.5*k2 , phip+0.5*l2 )
    k3 = h*(phip+0.5*l2)
    l4 = h*RH_cons( r+h , phi+k3 , phip+l3) 
    k4 = h*(phip+0.5*l3)
    return (k1,k2,k3,k4,l1,l2,l3,l4)

#Left hand constraction
def LH_cons(r,phi,phip,h):
    (k1,k2,k3,k4,l1,l2,l3,l4) = RK_cons(r,phi,phip,h)
    phinew  = phi  + (k1+2.0*k2+2.0*k3+k4)/6.0
    phipnew = phip + (l1+2.0*l2+2.0*l3+l4)/6.0
    return (phinew,phipnew)

# r is changing from zero to one w = [r0 inf[
# w is not defined in w = 1 so we need to use somewhere close to w = 1
# so we would start from point w = 1 - h

#Apply Boundary condition
h      = h3
(f,fp) = appl_bondary(r0,h)
fr     = [f]
r      = [r0+h]

#Mesh Grid Size
h = h3
rp     = r[0]+h

#Solver
while (rp<1.5*r0):
   r.append(rp)
   (fnew,fpnew) = LH_cons(rp,f,fp,h)
   fr.append(fnew)
   rp = rp+h
   f  = fnew
   fp = fpnew

h = h2
while (rp<10.0*r0):
   r.append(rp)
   (fnew,fpnew) = LH_cons(rp,f,fp,h)
   fr.append(fnew)
   rp = rp+h
   f  = fnew
   fp = fpnew

h = h1
while (rp<rmax):
   r.append(rp)
   (fnew,fpnew) = LH_cons(rp,f,fp,h)
   fr.append(fnew)
   rp = rp+h
   f  = fnew
   fp = fpnew

r.pop()
fr.pop()


ff = zeros(len(fr))
for i in range(0, len(fr)):
   ff[i] = fr[i].real
print ff

figure(1)
semilogx(r,ff,'b',label='Scale 1')


###########################################################################
###########################################################################
###########################################################################

del fr,r,ff

#Apply Boundary condition
h      = h3
(f,fp) = appl_bondary(r0,h)
fr     = [f]
r      = [r0+h]

#Mesh Grid Size
h = h3/10
rp     = r[0]+h

#Solver
while (rp<1.5*r0):
   r.append(rp)
   (fnew,fpnew) = LH_cons(rp,f,fp,h)
   fr.append(fnew)
   rp = rp+h
   f  = fnew
   fp = fpnew

h = h2/10
while (rp<10.0*r0):
   r.append(rp)
   (fnew,fpnew) = LH_cons(rp,f,fp,h)
   fr.append(fnew)
   rp = rp+h
   f  = fnew
   fp = fpnew

h = h1/10
while (rp<rmax):
   r.append(rp)
   (fnew,fpnew) = LH_cons(rp,f,fp,h)
   fr.append(fnew)
   rp = rp+h
   f  = fnew
   fp = fpnew

r.pop()
fr.pop()


ff = zeros(len(fr))
for i in range(0, len(fr)):
   ff[i] = fr[i].real
print ff

figure(1)
semilogx(r,ff,'r',label='Scale 2')

###########################################################################
###########################################################################
###########################################################################

del fr,r,ff

#Apply Boundary condition
h      = h3
(f,fp) = appl_bondary(r0,h)
fr     = [f]
r      = [r0+h]

#Mesh Grid Size
h = h3/100
rp     = r[0]+h

#Solver
while (rp<1.5*r0):
   r.append(rp)
   (fnew,fpnew) = LH_cons(rp,f,fp,h)
   fr.append(fnew)
   rp = rp+h
   f  = fnew
   fp = fpnew

h = h2/10
while (rp<10.0*r0):
   r.append(rp)
   (fnew,fpnew) = LH_cons(rp,f,fp,h)
   fr.append(fnew)
   rp = rp+h
   f  = fnew
   fp = fpnew

h = h1/100
while (rp<rmax):
   r.append(rp)
   (fnew,fpnew) = LH_cons(rp,f,fp,h)
   fr.append(fnew)
   rp = rp+h
   f  = fnew
   fp = fpnew

r.pop()
fr.pop()


ff = zeros(len(fr))
for i in range(0, len(fr)):
   ff[i] = fr[i].real
print ff

figure(1)
semilogx(r,ff,'g',label='Scale 3')
legend(loc=1)
savefig('plot4.pdf')

#show()

#a = GN_solver(w[len(w)-mm:len(w)],ff[len(w)-mm:len(w)],[-1,1])
#x = linspace(0,0.5,50)
#y = a[0]*x**3 + a[1]*x
#plot(x,y,'m')
#print 'Fit equation is :%4.5f*x^3 + %4.5f*x' % (a[0],a[1])


xList = r[len(r)-mm:len(r)]
yList = ff[len(r)-mm:len(r)]     
#a = polyFit(xList,yList,order=4)
#a = polyFit(xList,yList,v,k,d)
a = 1
b = -1
#[a,b] = Gauss_Newton_Solver(xList,yList,v,kt,d,a,b)

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

#plot(x,y,'r')
#show()
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
