from numpy             import *
from curvfit           import *
from cmath             import *
from matplotlib.pyplot import *
from RK_solver         import RK4_solver, Adaptive_RK4_solver
from RK_solver_T0      import RK4_solver_T0, Adaptive_RK4_solver

def Green_func(rs,r0,L,kx,kt,ms,d,h3,h2,h1,rmax,Q,M,muq,Scale):

#   err = 1e-10
#   [r,fr,tot_err]= Adaptive_RK4_solver(rs,r0,R2,kx,kt,Rm2,d,rmax,Q,M,v,muq,err,Scale)
   if (r0 != rs):
      # Finite Temprature
      [r,fpr,fnr] = RK4_solver(rs,r0,L,kx,kt,ms,d,h3,h2,h1,rmax,Q,M,muq)
   else:
      # Zero   Temprature
#      [r,fpr,fnr] = RK4_solver_T0(r0,L,kx,kt,ms,d,h3,h2,h1,rmax,muq)
      [r,fpr,fnr] = RK4_solver(rs,r0,L,kx,kt,ms,d,h3,h2,h1,rmax,Q,M,muq)


   rs    = r[len(r)-2]
   phi_p = fpr[len(r)-2]
   phi_n = fnr[len(r)-2]
   [A,C] = two_point_curve_fit(rs,phi_p,phi_n,L,kx,kt,ms,muq) 

   aHat = A # Coefficient of r^{mL}
   bHat = C # Coefficient of r^{-mL}

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

   Green_f = bHat/aHat

   return Green_f

