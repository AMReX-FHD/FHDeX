import math
import numpy as np

# cgs units

Runiv = 8.31446261815324e7  # universal gas constant
kB = 1.38064852e-16         # Boltzmann constant
Navo = 6.0221409e23         # Avogadro constant

##########

# molecular weights
M1 =  4.0026    # Helium
M2 = 20.1797    # Neon
M3 = 39.9480    # Argon
M4 = 83.8000    # Krypton
print "- molecular masses: [%.2f, %.2f, %.2f, %.2f]" % (M1,M2,M3,M4)

# mass fractions
Y1 = 0.25
Y2 = 0.25
Y3 = 0.25
Y4 = 0.25
print "- mass fractions: [%.3f, %.3f, %.3f, %.3f]" % (Y1,Y2,Y3,Y4)

# mole fractions
tmp = Y1/M1+Y2/M2+Y3/M3+Y4/M4
X1 = Y1/M1/tmp
X2 = Y2/M2/tmp
X3 = Y3/M3/tmp
X4 = Y4/M4/tmp
print "- mole fractions: [%.3f, %.3f, %.3f, %.3f]" % (X1,X2,X3,X4)

# average molecular weight
Mavg = 1/(Y1/M1+Y2/M2+Y3/M3+Y4/M4)
mavg = Mavg/Navo
print "- average molecular weight = %.2f" % Mavg
print "- average mass of a gas molecule = %e\n" % mavg

##########

temp = 300.                 # room temperature
pres = 1.01325e6            # atmospheric pressure
rho = pres*Mavg/Runiv/temp  # total mass density

print "- temperature = %f" % temp
print "- pressure = %e" % pres
print "- total mass density = %e\n" % rho

##########

ntot = pres/kB/temp
n1 = X1*ntot
n2 = X2*ntot
n3 = X3*ntot
n4 = X4*ntot

m1 = M1/Navo
m2 = M2/Navo
m3 = M3/Navo
m4 = M4/Navo

print "- total number density = %e" % ntot
print "- number density of spec1 = %e" % n1
print "- number density of spec2 = %e" % n2
print "- number density of spec3 = %e" % n3
print "- number density of spec4 = %e\n" % n4

##########

dx = 8.e-6
dt = 1.e-12
dV = dx**3
print "- dx (FHD) = %e" % dx
print "- dV (FHD) = %e" % dV
print "- dt (FHD) = %e" % dt

lat_const = dx/10
print "- lattice constant = %e\n" % lat_const

##########

Ntot = ntot*dV
N1 = n1*dV
N2 = n2*dV
N3 = n3*dV
N4 = n4*dV

print "- total number of gas molecules in dV = %.3e" % Ntot
print "- number of gas molecules of spec1 in dV = %.3e" % N1
print "- number of gas molecules of spec2 in dV = %.3e" % N2
print "- number of gas molecules of spec3 in dV = %.3e" % N3
print "- number of gas molecules of spec4 in dV = %.3e\n" % N4


##########

rho1 = n1*m1
rho2 = n2*m2
rho3 = n3*m3
rho4 = n4*m4

drho1sq = rho1**2/N1
drho2sq = rho2**2/N2
drho3sq = rho3**2/N3
drho4sq = rho4**2/N4

drhosq = drho1sq+drho2sq+drho3sq+drho4sq

djasq = rho*kB*temp/dV

print "drhosq = %e" % drhosq
print "drho1sq = %e" % drho1sq
print "drho2sq = %e" % drho2sq
print "drho3sq = %e" % drho3sq
print "drho4sq = %e" % drho4sq
print "djasq = %e" % djasq
print

print "drhosq*dV = %e" % (drhosq*dV)
print "drho1sq*dV = %e" % (drho1sq*dV)
print "drho2sq*dV = %e" % (drho2sq*dV)
print "drho3sq*dV = %e" % (drho3sq*dV)
print "drho4sq*dV = %e" % (drho4sq*dV)
print "djasq*dV = %e" % (djasq*dV)
print
