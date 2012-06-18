#! /usr/bin/python

from numpy             import *
from pylab             import *
from matplotlib.pyplot import *
from matplotlib        import rc
from cmath             import exp, cos, sin
from curvfit           import *
from readdata         import read_data
from RK_solver        import RK_solver

#Constants
(rs,r0,R2,kx,kt,Rm2,d,h3,h2,h1,rmax,Q,M,v,muq) = read_data('input.xml')
#k   = 1.0
mm    = 10

[r,fr] = RK_solver()

f_real = zeros(len(fr))
for i in range(0, len(fr)):
   f_real[i] = fr[i].real
print f_real

f_imag = zeros(len(fr))
for i in range(0, len(fr)):
   f_imag[i] = fr[i].imag
print f_imag

figure(1)
semilogx(r,f_real,'b',label='Real Part')
semilogx(r,f_imag,'r',label='Imaginary Part')


###########################################################################
###########################################################################
###########################################################################

#del fr,r,f_real,f_imag

#legend(loc=1)
#savefig('plot4.pdf')
show()

#a = GN_solver(w[len(w)-mm:len(w)],ff[len(w)-mm:len(w)],[-1,1])
#x = linspace(0,0.5,50)
#y = a[0]*x**3 + a[1]*x
#plot(x,y,'m')
#print 'Fit equation is :%4.5f*x^3 + %4.5f*x' % (a[0],a[1])


xList = r[len(r)-mm:len(r)]
yList = fr[len(r)-mm:len(r)]     
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

#x = linspace(0,0.2,50)
#y = aHat*x**3 + bHat*x**2 + cHat*x + dHat
#y = aHat*x**4 + bHat#*x
#y = aHat*x**(d/2.0)*J_v(n*v,k*x) + bHat*x**(d/2.0)*Y_v(n*v,kt*x)
#y = aHat*J_v(v,x,kt,d) + bHat*Y_v(v,x,kt,d)

#plot(x,y,'r')
#show()
#print 'Fit equation is :%4.5f*x^3 + %4.5f*x^2 + %4.5f*x + %4.5f' % (aHat,bHat,cHat,dHat)
#print 'Fit equation is :%4.5f*x^4 + %4.5f*x^3 + %4.5f*x^2 +  %4.5f*x + %4.5f' % (aHat,bHat,cHat,dHat,eHat)
#print 'Fit equation is :%4.5f*J_v + %4.5f*Y_v' % (aHat,bHat)


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
