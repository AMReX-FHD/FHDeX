import math

# cgs units

Runiv = 8.31446261815324e7  # universal gas constant
kB = 1.38064852e-16         # Boltzmann constant
Navo = 6.0221409e23         # Avogadro constant

##########

# molecular weights
M1 = 28.01    # CO
M2 = 39.95    # Ar
print "- molecular masses: [%.4f, %.4f]" % (M1,M2)

# mass fractions
Y1 = 0.05
Y2 = 0.95
print "- mass fractions: [%.3f, %.3f]" % (Y1,Y2)

# mole fractions
tmp = Y1/M1+Y2/M2
X1 = Y1/M1/tmp
X2 = Y2/M2/tmp
print "- mole fractions: [%.3f, %.3f]" % (X1,X2)

# average molecular weight
Mavg = 1/(Y1/M1+Y2/M2)
mavg = Mavg/Navo
print "- average molecular weight = %.2f" % Mavg
print "- average mass of a gas molecule = %e\n" % mavg

##########

temp = 1000.                # room temperature
pres = 1.01325e6            # atmospheric pressure
rho = pres*Mavg/Runiv/temp  # total mass density

print "- temperature = %f" % temp
print "- pressure = %e" % pres
print "- total mass density = %e\n" % rho

##########

ntot = pres/kB/temp
n1 = X1*ntot
n2 = X2*ntot

m1 = M1/Navo
m2 = M2/Navo

print "- total number density = %e" % ntot
print "- number density of spec1 = %e" % n1
print "- number density of spec2 = %e\n" % n2

##########

lat_const = 2.8e-8
n1x = 300
n1y = 150

dx = n1x*lat_const
dy = n1y*lat_const*math.sqrt(3)
dz = dx
dv = dx*dy*dz

n2x = 8
n2y = 8
n2z = 8
ncell = n2x*n2y*n2z

Lx = n2x*dx
Ly = n2y*dy
Lz = n2z*dz

print "- lattice constant for Pt(111): %e" % lat_const
print "- dx (FHD) = %e" % dx
print "- dy (FHD) = %e" % dy
print "- dz (FHD) = %e" % dz
print "- dv (FHD) = %e" % dv
print "- Lx = %e" % Lx
print "- Ly = %e" % Ly
print "- Lz = %e\n" % Lz

##########

n1s = 2*n1x*n1y
print "- total number of sites per dxFHD*dyFHD = %d" % n1s
print "- surf_site_num_dens = %e" % (n1s/dx/dy)

Ntot = ntot*dv
N1 = n1*dv
N2 = n2*dv

print "- total number of gas molecules in dv = %.3e" % Ntot
print "- number of gas molecules of spec1 in dv = %.3e" % N1
print "- number of gas molecules of spec2 in dv = %.3e\n" % N2

##########

rcol = math.sqrt(kB*temp/2/math.pi/m1)*n1*(lat_const**2*math.sqrt(3)/2)
k1ads = rcol
k1des = 1.25e15*math.exp(-1.514/(8.617e-5*temp))
theta1 = k1ads/(k1ads+k1des)

print "- rcol = %e" % rcol
print "- k1ads = %e (rate)" % k1ads
print "- k1ads/n1 (rate const) = %e" % (k1ads/n1)
print "- k1des = %e" % k1des
print "- eq coverage of spec1: %e\n" % theta1

##########

dt = 1.e-12

print "- dt (FHD) = %e" % dt
print "- max mean number of ads events per dxFHD*dyFHD per dt = %e" % (n1s*k1ads*dt)
print "- eq mean number of ads events per dxFHD*dyFHD per dt = %e" % (n1s*k1ads*dt*(1-theta1))
print "- max mean number of des events per dxFHD*dyFHD per dt = %e" % (n1s*k1des*dt)
print "- eq mean number of des events per dxFHD*dyFHD per dt = %e\n" % (n1s*k1des*dt*theta1)

##########

rho1 = n1*m1
rho2 = n2*m2
rhotot = rho1+rho2

drho1sq = rho1**2/N1
drho2sq = rho2**2/N2

drhosq = drho1sq+drho2sq
djasq = rho*kB*temp/dv

dof1 = 5
dof2 = 3
e01 = 0.
e02 = 0.

cv1 = 0.5*dof1*Runiv/M1
cv2 = 0.5*dof2*Runiv/M2
cvmix = Y1*cv1+Y2*cv2

e1 = e01+cv1*temp
e2 = e02+cv2*temp

dTsq = kB*temp**2/rhotot/cvmix/dv
dEsq = e1**2*drho1sq+e2**2*drho2sq+rhotot**2*cvmix**2*dTsq

print "drhosq = %e" % drhosq
print "drho1sq = %e" % drho1sq
print "drho2sq = %e" % drho2sq
print "djasq = %e" % djasq
print "dTsq = %e" % dTsq
print "dEsq = %e\n" % dEsq

print "drhosq = %e (size correction)" % (drhosq*(1-1./ncell))
print "drho1sq = %e (size correction)" % (drho1sq*(1-1./ncell))
print "drho2sq = %e (size correction)" % (drho2sq*(1-1./ncell))
print "djasq = %e (size correction - slip jx, jy)" % (djasq*(1-1./ncell))
print "dTsq = %e (size correction)" % (dTsq*(1+1./ncell))
print "dEsq = %e (size correction)\n" % (e1**2*drho1sq*(1-1./ncell)+e2**2*drho2sq*(1-1./ncell)+rhotot**2*cvmix**2*dTsq*(1+1./ncell))

print "drhosq  (+/-50%%) = %e\t%e\t%e\t%e" % (drhosq,2*drhosq,0.5*drhosq,1.5*drhosq)
print "drho1sq (+/-50%%) = %e\t%e\t%e\t%e" % (drho1sq,2*drho1sq,0.5*drho1sq,1.5*drho1sq)
print "drho2sq (+/-50%%) = %e\t%e\t%e\t%e" % (drho2sq,2*drho2sq,0.5*drho2sq,1.5*drho2sq)
print "djasq   (+/- 5%%) = %e\t%e\t%e\t%e" % (djasq,2*djasq,0.95*djasq,1.05*djasq)
print "dEsq    (+/-25%%) = %e\t%e\t%e\t%e" % (dEsq,2*dEsq,0.75*dEsq,1.25*dEsq)
print "dTsq    (+/- 5%%) = %e\t%e\t%e\t%e\n" % (dTsq,2*dTsq,0.95*dTsq,1.05*dTsq)

print "drhosq*dv  (+/-50%%) = %e\t%e\t%e\t%e" % (drhosq*dv,2*drhosq*dv,0.5*drhosq*dv,1.5*drhosq*dv)
print "drho1sq*dv (+/-50%%) = %e\t%e\t%e\t%e" % (drho1sq*dv,2*drho1sq*dv,0.5*drho1sq*dv,1.5*drho1sq*dv)
print "drho2sq*dv (+/-50%%) = %e\t%e\t%e\t%e" % (drho2sq*dv,2*drho2sq*dv,0.5*drho2sq*dv,1.5*drho2sq*dv)
print "djasq*dv   (+/-25%%) = %e\t%e\t%e\t%e" % (djasq*dv,2*djasq*dv,0.75*djasq*dv,1.25*djasq*dv)
print "dEsq*dv    (+/-50%%) = %e\t%e\t%e\t%e" % (dEsq*dv,2*dEsq*dv,0.5*dEsq*dv,1.5*dEsq*dv)
print "dTsq*dv    (+/-25%%) = %e\t%e\t%e\t%e\n" % (dTsq*dv,2*dTsq*dv,0.75*dTsq*dv,1.25*dTsq*dv)

d1 = 3.76e-8
d2 = 3.63e-8
d12 = (d1+d2)/2
D12 = 3./16*math.sqrt(2*math.pi*kB**3*(m1+m2)/m1/m2)/math.pi/d12**2*temp*math.sqrt(temp)/pres

u = 3*math.sqrt(kB*temp/mavg/Ntot)

print "dx = %e " % dx
print "dt = %e " % dt
print "D12 = %e" % D12
print "D12*dt/dx**2 = %e" % (D12*dt/dx**2)
print "u = %e" % u
print "u*dt/dx = %e" % (u*dt/dx)
