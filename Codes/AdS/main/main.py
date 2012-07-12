from matplotlib.pyplot import *
from numpy             import *
from Green_function    import Green_func
from readdata          import read_data


(rs,r0,R2,kx,kt,Rm2,d,h3,h2,h1,rmax,Q,M,v,muq) = read_data()

#kx = [0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0]
#kx = [0.0,0.5,1.5,2.5,5.0,10.0]
kx = 0.0
kt = [1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0]
Scale = 1
Gr_f = zeros(len(kt),complex)
for i in range(0,len(kt)):
   Gr_f[i] = Green_func(rs,r0,R2,kx,kt[i],Rm2,d,h3,h2,h1,rmax,Q,M,v,muq,Scale)
   print Gr_f[i]
print Scale
figure(1)
plot(kt,Gr_f.imag,label='Scale = %f'%(kx))
######################################################################################
kx = 5.0
for i in range(0,len(kt)):
   Gr_f[i] = Green_func(rs,r0,R2,kx,kt[i],Rm2,d,h3,h2,h1,rmax,Q,M,v,muq,Scale)
   print Gr_f[i]
print Scale
plot(kt,Gr_f.imag,label='Scale = %f'%(kx))
######################################################################################
kx = 10.0
for i in range(0,len(kt)):
   Gr_f[i] = Green_func(rs,r0,R2,kx,kt[i],Rm2,d,h3,h2,h1,rmax,Q,M,v,muq,Scale)
   print Gr_f[i]
print Scale
plot(kt,Gr_f.imag,label='Scale = %f'%(kx))
'''
######################################################################################
Scale = 10
for i in range(0,len(kx)):
   Gr_f[i] = Green_func(rs,r0,R2,kx[i],kt,Rm2,d,h3,h2,h1,rmax,Q,M,v,muq,Scale)
   print Gr_f[i]
print Scale
plot(kx,Gr_f.imag,label='Scale = %i'%(Scale))
######################################################################################
Scale = 20
for i in range(0,len(kx)):
   Gr_f[i] = Green_func(rs,r0,R2,kx[i],kt,Rm2,d,h3,h2,h1,rmax,Q,M,v,muq,Scale)
   print Gr_f[i]
print Scale
plot(kx,Gr_f.imag,label='Scale = %i'%(Scale))
######################################################################################
Scale = 50
for i in range(0,len(kx)):
   Gr_f[i] = Green_func(rs,r0,R2,kx[i],kt,Rm2,d,h3,h2,h1,rmax,Q,M,v,muq,Scale)
   print Gr_f[i]
print Scale
plot(kx,Gr_f.imag,label='Scale = %i'%(Scale))
'''
ylabel('Green Function (Imaginary Part)')
xlabel('Omega')
legend(loc=2)
show()
