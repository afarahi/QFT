from numpy             import *
from pylab             import *
from matplotlib.pyplot import *
from matplotlib        import rc
from cmath             import exp, cos, sin, sqrt
from readdata          import read_data

# write hand constraction
def redshift_f(r,d,Q,M):
    f = 1.0 + Q**2/(r**(2*d-2)) - M/r**d
    return f

#Right hand constraction
def RH_cons(r,phi,phip,r0,R2,kx,kt,Rm2,d,Q,M,muq):
    phipp = ( \
             - ( (d+1)/r - (d-3)*Q**2/r**(2*d-1) - M/r**(d+1) ) * phip \
             - ( R2*R2*(kt+muq*(1-(r0/r)**(d-2)))**2/(r**4*redshift_f(r,d,Q,M)) \
                 - kx*kx*R2*R2/(r**4) \
                 - Rm2/r**2         \
               ) * phi               \
            ) / redshift_f(r,d,Q,M)

    return phipp

#Boundary condition 
def appl_bondary(h,r0,rs,R2,kx,kt,Rm2,d,muq):
    r     = r0+h
    #rr = (r*/r0)^6
#    rr    = (rs/r0)**6
    beta  = R2*kt/(d*(1.0 - (rs/r0)**(2*d-2))*r0)
#    lambd1= (kt*muq*R2*R2)/((1.0-rr)*(7.0-rr)) \
#            - (R2*R2*kx*kx/(r0*r0) + Rm2)/(14.0 - 2.0*rr)
#    alpha1= 1.0/8.0 - sqrt(1.0/64.0 - lambd1/4.0)
#    lambd2= 
#    alpha2= -1.0 - sqrt(1.0/4.0 - lambd2)
    phi   = exp(-1.0j*beta*log(r-r0))
    phip  = -exp(-1.0j*beta*log(r-r0))*1.0j*beta/(r-r0)
    return (phi,phip)

#Runge-Kutta coefficient
def RK_cons(r,phi,phip,h,r0,R2,kx,kt,Rm2,d,Q,M,muq):
    l1 = h*RH_cons( r , phi , phip ,r0,R2,kx,kt,Rm2,d,Q,M,muq)
    k1 = h*phip
    l2 = h*RH_cons( r+0.5*h , phi+0.5*k1 , phip+0.5*l1 ,r0,R2,kx,kt,Rm2,d,Q,M,muq)
    k2 = h*(phip+0.5*l1)
    l3 = h*RH_cons( r+0.5*h , phi+0.5*k2 , phip+0.5*l2 ,r0,R2,kx,kt,Rm2,d,Q,M,muq)
    k3 = h*(phip+0.5*l2)
    l4 = h*RH_cons( r+h , phi+k3 , phip+l3 ,r0,R2,kx,kt,Rm2,d,Q,M,muq) 
    k4 = h*(phip+0.5*l3)
    return (k1,k2,k3,k4,l1,l2,l3,l4)

#Left hand constraction
def LH_cons(r,phi,phip,h,r0,R2,kx,kt,Rm2,d,Q,M,muq):
    (k1,k2,k3,k4,l1,l2,l3,l4) = RK_cons(r,phi,phip,h,r0,R2,kx,kt,Rm2,d,Q,M,muq)
    phinew  = phi  + (k1+2.0*k2+2.0*k3+k4)/6.0
    phipnew = phip + (l1+2.0*l2+2.0*l3+l4)/6.0
    return (phinew,phipnew)

#Solver
def RK4_solver(rs,r0,R2,kx,kt,Rm2,d,h3,h2,h1,rmax,Q,M,v,muq,scale=1):
    # r is changing from r0 to infty w = [r0 inf[
    # r is not defined in r = r0 so we need to use somewhere close to r = r0
    # so we would start from point r = r0 + h

    #Apply Boundary condition
    h      = 1e-4
    (f,fp) = appl_bondary(h,r0,rs,R2,kx,kt,Rm2,d,muq)
    fr     = [f]
    r      = [r0+h]
    rp     = r[0]

    #Solver
    #Mesh Grid Size
    h = h3/scale/1000
    while (rp<1.001*r0):
      (fnew,fpnew) = LH_cons(rp,f,fp,h,r0,R2,kx,kt,Rm2,d,Q,M,muq)
      rp = rp+h
      fr.append(fnew)
      r.append(rp)
      f  = fnew
      fp = fpnew

    #Mesh Grid Size
    h = h3/scale/100
    while (rp<1.01*r0):
      (fnew,fpnew) = LH_cons(rp,f,fp,h,r0,R2,kx,kt,Rm2,d,Q,M,muq)
      rp = rp+h
      r.append(rp)
      fr.append(fnew)
      f  = fnew
      fp = fpnew

    h = h3/scale/10
    while (rp<1.1*r0):
      (fnew,fpnew) = LH_cons(rp,f,fp,h,r0,R2,kx,kt,Rm2,d,Q,M,muq)
      rp = rp+h
      r.append(rp)
      fr.append(fnew)
      f  = fnew
      fp = fpnew

    #Mesh Grid Size
    h = h3/scale
    while (rp<2*r0):
      (fnew,fpnew) = LH_cons(rp,f,fp,h,r0,R2,kx,kt,Rm2,d,Q,M,muq)
      rp = rp+h
      r.append(rp)
      fr.append(fnew)
      f  = fnew
      fp = fpnew

    #Mesh Grid Size
    h = h2/scale
    while (rp<10.0*r0):
      (fnew,fpnew) = LH_cons(rp,f,fp,h,r0,R2,kx,kt,Rm2,d,Q,M,muq)
      fr.append(fnew)
      rp = rp+h
      r.append(rp)
      f  = fnew
      fp = fpnew

    #Mesh Grid Size
    h = max(0.1,h1/scale)
    while (rp<100*r0):
      (fnew,fpnew) = LH_cons(rp,f,fp,h,r0,R2,kx,kt,Rm2,d,Q,M,muq)
      rp = rp+h
      r.append(rp)
      fr.append(fnew)
      f  = fnew
      fp = fpnew

    #Mesh Grid Size
    h = max(0.5,h1/scale)
    while (rp<rmax):
      (fnew,fpnew) = LH_cons(rp,f,fp,h,r0,R2,kx,kt,Rm2,d,Q,M,muq)
      rp = rp+h
      r.append(rp)
      fr.append(fnew)
      f  = fnew
      fp = fpnew

    r.pop()
    fr.pop()
    return(r,fr)

#h Calculator
def h_cal(err,hp,h1):
    hs  = h1*hp
    hp  = max(hp,0.5)
    h   = min(hs,hp)/2.0
    return h

#Solver
def Adaptive_RK4_solver(rs,r0,R2,kx,kt,Rm2,d,rmax,Q,M,v,muq,err,scale=1.0):
    # r is changing from r0 to infty w = [r0 inf[
    # r is not defined in r = r0 so we need to use somewhere close to r = r0
    # so we would start from point r = r0 + h
    beta  = R2*kt/(d*(1.0 - (rs/r0)**(2*d-2))*r0)
    hp    = abs( (120.0*err/(sqrt((10.0*beta**2+5.0*beta**4)**2+(24.0*beta+beta**5)**2)))**(0.2) )

    #Apply Boundary condition
    h1     = 1e-4
    (f,fp) = appl_bondary(h1,r0,rs,R2,kx,kt,Rm2,d,muq)
    fr     = [f]
    r      = [r0+h1]
    rold   = r[0]

    #Solver
    while (rold<rmax):
      #Mesh Grid Size
      h = h_cal(err,hp,r[len(r)-1]-r0)/scale
      (fnew,fpnew) = LH_cons(rold,f,fp,h,r0,R2,kx,kt,Rm2,d,Q,M,muq)
      rold = rold+h
      r.append(rold)
      fr.append(fnew)
      f    = fnew
      fp   = fpnew

    r.pop()
    fr.pop()

    tot_err = len(r)*err

    return(r,fr,tot_err)
