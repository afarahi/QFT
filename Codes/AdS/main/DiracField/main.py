from matplotlib.pyplot import *
from numpy             import *
from Green_function    import Green_func
from readdata          import read_data

(rs,r0,L,kx,kt,ms,d,h3,h2,h1,rmax,Q,M,muq) = read_data()

#kx = [0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0]
#kx = [0.0,0.5,1.5,2.5,5.0,10.0]
n  = 1.0
kx = 1.2
kt = linspace(-2.0,2.0,40)#[ -n*2.0 , -n*1.95 , -n*1.9 , -n*1.8 , -n*1.85 , -n*1.7 , -n*1.6 , -n*1.5 , -n*1.4 , -n*1.3 , -n*1.2 , -n*1.1 , -n*1.0 , -n*0.9 , -n*0.8 ,  -n*0.7 , -n*0.6 , -n*0.5 , -n*0.4 , -n*0.3 , -n*0.2 , -n*0.1 , n*0.0 , n*0.25 , n*0.5 , n*0.75 , n*1.0 , n*1.25 , n*1.5 , n*1.75 , n*2.0 ]
Scale = 0.001
Gr_f = zeros(len(kt),complex)
for i in range(0,len(kt)):
   Gr_f[i] = Green_func(rs,r0,L,kx,kt[i],ms,d,h3,h2,h1,rmax,Q,M,muq,Scale)
   print Gr_f[i]
print Scale
figure(1)
plot(kt,Gr_f.imag,label='Scale = %f'%Scale)
#######################################################################################
#Scale = 10
#Gr_f = zeros(len(kt),complex)
#for i in range(0,len(kt)):
#   Gr_f[i] = Green_func(rs,r0,L,kx,kt[i],ms,d,h3,h2,h1,rmax,Q,M,v,muq,Scale)
#   print Gr_f[i]
#print Scale
#figure(1)
#plot(kt,Gr_f.imag,label='Scale = %f'%Scale)

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
#legend(loc=2)
show()
