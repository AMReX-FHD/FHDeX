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

drho1sq = rho1**2/N1
drho2sq = rho2**2/N2

drhosq = drho1sq+drho2sq
djasq = rho*kB*temp/dv

dof1 = 6. 
dof2 = 6.

dEsq = (kB*temp)**2/dv*(0.5*dof1*0.5*(dof1+2)*rho1/m1+0.5*dof2*0.5*(dof2+2)*rho2/m2)
dTsq = temp**2/dv/(0.5*dof1*n1+0.5*dof2*n2)

print "drhosq = %e" % drhosq
print "drho1sq = %e" % drho1sq
print "drho2sq = %e" % drho2sq
print "djasq = %e\n" % djasq

print "drhosq*dv  = %e\t%e\t%e\t%e" % (drhosq*dv,2*drhosq*dv,0.5*drhosq*dv,1.5*drhosq*dv)
print "drho1sq*dv = %e\t%e\t%e\t%e" % (drho1sq*dv,2*drho1sq*dv,0.5*drho1sq*dv,1.5*drho1sq*dv)
print "drho2sq*dv = %e\t%e\t%e\t%e" % (drho2sq*dv,2*drho2sq*dv,0.5*drho2sq*dv,1.5*drho2sq*dv)
print "djasq*dv   = %e\t%e\t%e\t%e" % (djasq*dv,2*djasq*dv,0.5*djasq*dv,1.5*djasq*dv)

print "dEsq*dv    = %e\t%e\t%e\t%e" % (dEsq*dv,2*dEsq*dv,0.5*dEsq*dv,1.5*dEsq*dv)
print "dTsq*dv    = %e\t%e\t%e\t%e" % (dTsq*dv,2*dTsq*dv,0.5*dTsq*dv,1.5*dTsq*dv)
