from numpy             import *
from pylab             import *
from matplotlib.pyplot import *
from matplotlib        import rc
from cmath             import exp, cos, sin, sqrt
from readdata          import read_data

#Green Function_{11} Solver

#phip = \phi_{+}
#phin = \phi_{-}
#phipp= \phi'_{+}
#phinp= \phi'_{-}

# write hand constraction
def redshift_f(r,d,Q,M):
    f = 1.0 + Q**2/(r**(2*d-2)) - M/r**d
    return f

#Right hand constraction for Phi'_{+}
def RHp_cons(r,phip,phin,r0,L,kx,kt,ms,d,Q,M,muq):
    phipp = L*( ms*phip +  (-L*kx/r - L*kt*(1.0+muq*(1.0-(r0/r)**2))/(r*sqrt(redshift_f(r,d,Q,M))))*phin ) / ( r*sqrt(redshift_f(r,d,Q,M)) )
    return phipp

#Right hand constraction for Phi'_{-}
def RHn_cons(r,phip,phin,r0,L,kx,kt,ms,d,Q,M,muq):
    phinp = L*( -ms*phin + (-L*kx/r + L*kt*(1.0+muq*(1.0-(r0/r)**2))/(r*sqrt(redshift_f(r,d,Q,M))))*phip ) / ( r*sqrt(redshift_f(r,d,Q,M)) )
    return phinp

#Boundary condition 
def appl_bondary(h,r0,rs,L,kx,kt,ms,d,muq):
    r     = r0+h
    if (rs != r0):
       # Finite Temprature
       kappa = 4.0*(r0**6-rs**6)/(r0**7)
       beta  = L*L*kt/(kappa*r0)
       phip  = exp(-1.0j*beta*log(r-r0))
       phin  = 1.0j*exp(-1.0j*beta*log(r-r0))
    else:
       # Zero   Temprature
       beta  = L*L*kt/12.0
       a1p   = ( - 12.0*beta/L * ( -28.0*beta/(r0*L) - 2.0*L*muq + 2.0*L*kx*sqrt(3)/r0 ) \
                 - L*kt * ( -28.0*beta/(r0*L) - 2.0*L*muq - 2.0*L*kx*sqrt(3)/r0        ) \
               ) / (144.0*beta**2/L**2 - L**2*kt**2 )
       a1n   = ( - L*kt * ( -28.0*beta/(r0*L) - 2.0*L*muq + 2.0*L*kx*sqrt(3)/r0        ) \
                 - 12.0*beta/L * ( -28.0*beta/(r0*L) - 2.0*L*muq - 2.0*L*kx*sqrt(3)/r0 ) \
               ) / (144.0*beta**2/L**2 - L**2*kt**2 )
       phip  = exp(1.0j*beta/(r-r0)) * (1.0 + a1p*(r-r0))
       phin  = 1.0j*exp(1.0j*beta/(r-r0)) * (1.0 + a1n*(r-r0))
    return (phip,phin)

#Runge-Kutta coefficient
def RK_cons(r,phip,phin,h,r0,L,kx,kt,ms,d,Q,M,muq):
    kp1 = h*RHp_cons( r       , phip         , phin         ,r0,L,kx,kt,ms,d,Q,M,muq)
    kn1 = h*RHn_cons( r       , phip         , phin         ,r0,L,kx,kt,ms,d,Q,M,muq)
    kp2 = h*RHp_cons( r+0.5*h , phip+0.5*kp1 , phin+0.5*kn1 ,r0,L,kx,kt,ms,d,Q,M,muq)
    kn2 = h*RHn_cons( r+0.5*h , phip+0.5*kp1 , phin+0.5*kn1 ,r0,L,kx,kt,ms,d,Q,M,muq)
    kp3 = h*RHp_cons( r+0.5*h , phip+0.5*kp2 , phin+0.5*kn2 ,r0,L,kx,kt,ms,d,Q,M,muq)
    kn3 = h*RHn_cons( r+0.5*h , phip+0.5*kp2 , phin+0.5*kn2 ,r0,L,kx,kt,ms,d,Q,M,muq)
    kp4 = h*RHp_cons( r+h     , phip+kp3     , phin+kn3     ,r0,L,kx,kt,ms,d,Q,M,muq) 
    kn4 = h*RHn_cons( r+h     , phip+kp3     , phin+kn3     ,r0,L,kx,kt,ms,d,Q,M,muq)
    return (kp1,kp2,kp3,kp4,kn1,kn2,kn3,kn4)

#Left hand constraction
def LH_cons(r,phip,phin,h,r0,L,kx,kt,ms,d,Q,M,muq):
    (kp1,kp2,kp3,kp4,kn1,kn2,kn3,kn4) = RK_cons(r,phip,phin,h,r0,L,kx,kt,ms,d,Q,M,muq)
    phipnew  = phip + (kp1+2.0*kp2+2.0*kp3+kp4)/6.0
    phinnew  = phin + (kn1+2.0*kn2+2.0*kn3+kn4)/6.0
    return (phipnew,phinnew)

#Solver
def RK4_solver_G11(rs,r0,L,kx,kt,ms,d,h3,h2,h1,rmax,Q,M,muq,scale=1):
    # r is changing from r0 to infty w = [r0 inf[
    # r is not defined in r = r0 so we need to use somewhere close to r = r0
    # so we would start from point r = r0 + h

    #Apply Boundary condition
    h       = 1e-2
    (fp,fn) = appl_bondary(h,r0,rs,L,kx,kt,ms,d,muq)
    fpr     = [fp]
    fnr     = [fn]
    r       = [r0+h]
    rp      = r[0]

    #Solver
    #Mesh Grid Size
    h = h3/scale/10
    while (rp<1.1*r0):
      (fpnew,fnnew) = LH_cons(rp,fp,fn,h,r0,L,kx,kt,ms,d,Q,M,muq)
      rp = rp+h
      r.append(rp)
      fpr.append(fpnew)
      fnr.append(fnnew)
      fp  = fpnew
      fn  = fnnew

    #Mesh Grid Size
    h = h3/scale
    while (rp<2*r0):
      (fpnew,fnnew) = LH_cons(rp,fp,fn,h,r0,L,kx,kt,ms,d,Q,M,muq)
      rp = rp+h
      r.append(rp)
      fpr.append(fpnew)
      fnr.append(fnnew)
      fp  = fpnew
      fn  = fnnew

    #Mesh Grid Size
    h = h2/scale
    while (rp<10.0*r0):
      (fpnew,fnnew) = LH_cons(rp,fp,fn,h,r0,L,kx,kt,ms,d,Q,M,muq)
      rp = rp+h
      r.append(rp)
      fpr.append(fpnew)
      fnr.append(fnnew)
      fp  = fpnew
      fn  = fnnew

    #Mesh Grid Size
    h = max(0.1,h1/scale)
    while (rp<100*r0):
      (fpnew,fnnew) = LH_cons(rp,fp,fn,h,r0,L,kx,kt,ms,d,Q,M,muq)
      rp = rp+h
      r.append(rp)
      fpr.append(fpnew)
      fnr.append(fnnew)
      fp  = fpnew
      fn  = fnnew

    #Mesh Grid Size
    h = max(0.5,h1/scale)
    while (rp<rmax):
      (fpnew,fnnew) = LH_cons(rp,fp,fn,h,r0,L,kx,kt,ms,d,Q,M,muq)
      rp = rp+h
      r.append(rp)
      fpr.append(fpnew)
      fnr.append(fnnew)
      fp  = fpnew
      fn  = fnnew

    r.pop()
    fpr.pop()
    fnr.pop()
    return(r,fpr,fnr)

#THIS SECTION NEEDS WORK
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
