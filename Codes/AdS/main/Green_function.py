from numpy             import *
from curvfit           import *
from cmath             import *
from matplotlib.pyplot import *
from RK_solver         import RK4_solver, Adaptive_RK4_solver

def Green_func(rs,r0,R2,kx,kt,Rm2,d,h3,h2,h1,rmax,Q,M,v,muq,Scale):
   #Constants
   Delta = d/2 + sqrt(d**2/4 + Rm2)
   pow_r = [-Delta,Delta-d]

   err = 1e-10
   [r,fr,tot_err]= Adaptive_RK4_solver(rs,r0,R2,kx,kt,Rm2,d,rmax,Q,M,v,muq,err,Scale)

   rlen   = len(r)
   for i in range(1,len(r)):
       if (r[rlen-i] < (rmax-2)):
          mm = i
          break

#   xList = r[len(r)-mm:len(r)]
#   yList = fr[len(r)-mm:len(r)]     
#   a     = polyFit(xList,yList,pow_r)
#   a    = Gauss_Newton_Solver(xList,yList,pow_r,1.0,1.0)
   x   = [r[len(r)-mm],r[len(r)-1]]
   y   = [fr[len(r)-mm],fr[len(r)-1]]
   a   = two_point_curve_fit(x,y,pow_r,R2,kx,kt,Rm2,rmax,Q,M,muq) 

   aHat = a[0] # Coefficient of r^{-\Delta}
   bHat = a[1] # Coefficient of r^{\Delta-d}

   #Test Part
   '''
   x =linspace(rmax-25,rmax,1000)
   y = aHat*x**pow_r[0] + bHat*x**pow_r[1]

   plot(r,fr,'r')
   plot(x,y,'b')
   show()

   plot(r,fr,'r')
   plot(x,y,'b')
   show()
   '''

   Green_f = aHat/bHat

   return Green_f

