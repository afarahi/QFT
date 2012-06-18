from numpy             import *
from pylab             import *
from matplotlib.pyplot import *
from matplotlib        import rc
from cmath             import exp, cos, sin
from readdata          import read_data

# write hand constraction
def redshift_f(r,rs,r0,R2,kx,kt,Rm2,d,h3,h2,h1,rmax,Q,M,v):
    f = 1.0 + Q**2/(r**(2*d-2)) - M/r**d
    return f

#Right hand constraction
def RH_cons(r,phi,phip,rs,r0,R2,kx,kt,Rm2,d,h3,h2,h1,rmax,Q,M,v,muq):
    phipp = ( \
             - ( (d+1)*r - (d-3)*Q**2/r**(2*d-3) - M/r**(d-1) ) * phip \
             - ( - R2*(kt+muq*(1-(r0/r)**(d-2)))**2/(r**2*redshift_f(r,rs,r0,R2,kx,kt,Rm2,d,h3,h2,h1,rmax,Q,M,v)) \
                 + kx**2*R2/(r**2) \
                 - Rm2             \
               ) * phi               \
            ) / (r**2*redshift_f(r,rs,r0,R2,kx,kt,Rm2,d,h3,h2,h1,rmax,Q,M,v))
    return phipp

#Boundary condition 
def appl_bondary(rstart,h,rs,r0,R2,kx,kt,Rm2,d,h3,h2,h1,rmax,Q,M,v,muq):
    r   = rstart+h
    beta= R2*kt/( d*(1.0 - (rs/r0)**(2*d-2)) )
    phi = exp(-1.0j*beta*log(r-r0))
    phip= -exp(-1.0j*beta*log(r-r0))*1.0j*beta/(r-r0)
    return (phi,phip)

#Runge-Kutta coefficient
def RK_cons(r,phi,phip,h,rs,r0,R2,kx,kt,Rm2,d,h3,h2,h1,rmax,Q,M,v,muq):
    l1 = h*RH_cons( r , phi , phip ,rs,r0,R2,kx,kt,Rm2,d,h3,h2,h1,rmax,Q,M,v,muq)
    k1 = h*phip
    l2 = h*RH_cons( r+0.5*h , phi+0.5*k1 , phip+0.5*l1 ,rs,r0,R2,kx,kt,Rm2,d,h3,h2,h1,rmax,Q,M,v,muq)
    k2 = h*(phip+0.5*l1)
    l3 = h*RH_cons( r+0.5*h , phi+0.5*k2 , phip+0.5*l2 ,rs,r0,R2,kx,kt,Rm2,d,h3,h2,h1,rmax,Q,M,v,muq)
    k3 = h*(phip+0.5*l2)
    l4 = h*RH_cons( r+h , phi+k3 , phip+l3 ,rs,r0,R2,kx,kt,Rm2,d,h3,h2,h1,rmax,Q,M,v,muq) 
    k4 = h*(phip+0.5*l3)
    return (k1,k2,k3,k4,l1,l2,l3,l4)

#Left hand constraction
def LH_cons(r,phi,phip,h,rs,r0,R2,kx,kt,Rm2,d,h3,h2,h1,rmax,Q,M,v,muq):
    (k1,k2,k3,k4,l1,l2,l3,l4) = RK_cons(r,phi,phip,h,rs,r0,R2,kx,kt,Rm2,d,h3,h2,h1,rmax,Q,M,v,muq)
    phinew  = phi  + (k1+2.0*k2+2.0*k3+k4)/6.0
    phipnew = phip + (l1+2.0*l2+2.0*l3+l4)/6.0
    return (phinew,phipnew)

#Solver
def RK_solver(scale=1):
    #Constants
    (rs,r0,R2,kx,kt,Rm2,d,h3,h2,h1,rmax,Q,M,v,muq) = read_data('input.xml')
    #k   = 1.0
    mm  = 10

    # r is changing from r0 to infty w = [r0 inf[
    # r is not defined in r = r0 so we need to use somewhere close to r = r0
    # so we would start from point r = r0 + h

    #Apply Boundary condition
    h      = 5e-5
    (f,fp) = appl_bondary(r0,h,rs,r0,R2,kx,kt,Rm2,d,h3,h2,h1,rmax,Q,M,v,muq)
    fr     = [f]
    r      = [r0+h]

    #Mesh Grid Size
    h = h3/scale
    rp     = r[0]+h

    #Solver
    while (rp<1.5*r0):
      r.append(rp)
      (fnew,fpnew) = LH_cons(rp,f,fp,h,rs,r0,R2,kx,kt,Rm2,d,h3,h2,h1,rmax,Q,M,v,muq)
      fr.append(fnew)
      rp = rp+h
      f  = fnew
      fp = fpnew

    h = h2/scale
    while (rp<10.0*r0):
      r.append(rp)
      (fnew,fpnew) = LH_cons(rp,f,fp,h,rs,r0,R2,kx,kt,Rm2,d,h3,h2,h1,rmax,Q,M,v,muq)
      fr.append(fnew)
      rp = rp+h
      f  = fnew
      fp = fpnew

    h = h1/scale
    while (rp<rmax):
      r.append(rp)
      (fnew,fpnew) = LH_cons(rp,f,fp,h,rs,r0,R2,kx,kt,Rm2,d,h3,h2,h1,rmax,Q,M,v,muq)
      fr.append(fnew)
      rp = rp+h
      f  = fnew
      fp = fpnew

    r.pop()
    fr.pop()
    return(r,fr)
