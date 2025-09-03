import math

# cgs units

Runiv = 8.31446261815324e7  # universal gas constant
kB = 1.38064852e-16         # Boltzmann constant
Navo = 6.0221409e23         # Avogadro constant

##########

# molecular weights
M1 = 46.0055    # NO2
M2 = 92.0110    # N2O4
print "- molecular masses: [%.4f, %.4f]" % (M1,M2)

# mass fractions
Y1 = 7.267597e-01
Y2 = 2.732403e-01
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

temp = 350.                 # room temperature
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

dx = 8e-6
dy = 8e-6
dz = 8e-6
dv = dx*dy*dz

Nx = 8
Ny = 8
Nz = 8

Lx = Nx*dx
Ly = Ny*dy
Lz = Nz*dz

print "- dx = %e" % dx
print "- dy = %e" % dy
print "- dz = %e" % dz
print "- dv = %e" % dv
print "- Lx = %e" % Lx
print "- Ly = %e" % Ly
print "- Lz = %e\n" % Lz

##########

Ntot = ntot*dv
N1 = n1*dv
N2 = n2*dv

print "- total number of gas molecules in dv = %.3e" % Ntot
print "- number of gas molecules of spec1 in dv = %.3e" % N1
print "- number of gas molecules of spec2 in dv = %.3e\n" % N2

##########

rho1 = n1*m1
rho2 = n2*m2
rhotot = rho1+rho2

drho1sq = rho1**2/N1
drho2sq = rho2**2/N2

drhosq = drho1sq+drho2sq
djasq = rho*kB*temp/dv

dof1 = 7.279582
dof2 = 18.049065
e01 = 4.684542e+09
e02 = -1.730892e+09

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
print "djasq = %e\n" % djasq

print "drhosq  (+/-50%%) = %e\t%e\t%e\t%e" % (drhosq,2*drhosq,0.5*drhosq,1.5*drhosq)
print "drho1sq (+/-50%%) = %e\t%e\t%e\t%e" % (drho1sq,2*drho1sq,0.5*drho1sq,1.5*drho1sq)
print "drho2sq (+/-50%%) = %e\t%e\t%e\t%e" % (drho2sq,2*drho2sq,0.5*drho2sq,1.5*drho2sq)
print "djasq   (+/- 3%%) = %e\t%e\t%e\t%e" % (djasq,2*djasq,0.97*djasq,1.03*djasq)
print "dEsq    (+/-50%%) = %e\t%e\t%e\t%e" % (dEsq,2*dEsq,0.5*dEsq,1.5*dEsq)
print "dTsq    (+/- 3%%) = %e\t%e\t%e\t%e\n" % (dTsq,2*dTsq,0.97*dTsq,1.03*dTsq)

print "drhosq*dv  (+/-50%%) = %e\t%e\t%e\t%e" % (drhosq*dv,2*drhosq*dv,0.5*drhosq*dv,1.5*drhosq*dv)
print "drho1sq*dv (+/-50%%) = %e\t%e\t%e\t%e" % (drho1sq*dv,2*drho1sq*dv,0.5*drho1sq*dv,1.5*drho1sq*dv)
print "drho2sq*dv (+/-50%%) = %e\t%e\t%e\t%e" % (drho2sq*dv,2*drho2sq*dv,0.5*drho2sq*dv,1.5*drho2sq*dv)
print "djasq*dv   (+/-10%%) = %e\t%e\t%e\t%e" % (djasq*dv,2*djasq*dv,0.9*djasq*dv,1.1*djasq*dv)
print "dEsq*dv    (+/-50%%) = %e\t%e\t%e\t%e" % (dEsq*dv,2*dEsq*dv,0.5*dEsq*dv,1.5*dEsq*dv)
print "dTsq*dv    (+/-10%%) = %e\t%e\t%e\t%e\n" % (dTsq*dv,2*dTsq*dv,0.9*dTsq*dv,1.1*dTsq*dv)

d1 = 3.765e-08
d2 = 4.621e-08
d12 = (d1+d2)/2
D12 = 3./16*math.sqrt(2*math.pi*kB**3*(m1+m2)/m1/m2)/math.pi/d12**2*temp*math.sqrt(temp)/pres

u = 3*math.sqrt(kB*temp/mavg/Ntot)

dt = 1.e-12
print "dx = %e " % dx
print "dt = %e " % dt
print "D12 = %e" % D12
print "D12*dt/dx**2 = %e" % (D12*dt/dx**2)
print "u = %e" % u
print "u*dt/dx = %e" % (u*dt/dx)
