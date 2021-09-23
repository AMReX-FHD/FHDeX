import math

# cgs units

Runiv = 8.31446261815324e7  # universal gas constant
kB = 1.38064852e-16         # Boltzmann constant
Navo = 6.0221409e23         # Avogadro constant

##########

# molecular weights
M1 = 28.01      # CO 
M2 =  4.0026    # He
print "- molecular masses: [%.2f, %.2f]" % (M1,M2)

# mass fractions
Y1 = 0.2
Y2 = 1-Y1 
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

dx = 300*2.8e-8 
dy = 150*2.8e-8*math.sqrt(3)
dz = dx
dV = dx*dy*dz

nx = 8
ny = 8
nz = 8

Lx = nx*dx
Ly = ny*dy
Lz = nz*dz

print "- dx = %e" % dx
print "- dy = %e" % dy
print "- dz = %e" % dz
print "- dV = %e" % dV
print "- Lx = %e" % Lx
print "- Ly = %e" % Ly
print "- Lz = %e\n" % Lz

##########

Ntot = ntot*dV
N1 = n1*dV
N2 = n2*dV

print "- total number of gas molecules in dV = %.3e" % Ntot
print "- number of gas molecules of spec1 in dV = %.3e" % N1
print "- number of gas molecules of spec2 in dV = %.3e\n" % N2

##########

rho1 = n1*m1
rho2 = n2*m2

drho1sq = rho1**2/N1
drho2sq = rho2**2/N2

drhosq = drho1sq+drho2sq
djasq = rho*kB*temp/dV

dof1 = 5
dof2 = 3 

dEsq = (kB*temp)**2/dV*((0.5*dof1)*(0.5*dof1+1)*rho1/m1+(0.5*dof2)*(0.5*dof2+1)*rho2/m2)
dTsq = temp**2/dV/((0.5*dof1)*n1+(0.5*dof2)*n2)

print "drhosq  (+/-50%%) = %e\t%e\t%e\t%e" % (drhosq,2*drhosq,0.5*drhosq,1.5*drhosq)
print "drho1sq (+/-50%%) = %e\t%e\t%e\t%e" % (drho1sq,2*drho1sq,0.5*drho1sq,1.5*drho1sq)
print "drho2sq (+/-50%%) = %e\t%e\t%e\t%e" % (drho2sq,2*drho2sq,0.5*drho2sq,1.5*drho2sq)
print "djasq   (+/-10%%) = %e\t%e\t%e\t%e" % (djasq,2*djasq,0.9*djasq,1.1*djasq)
print "dEsq    (+/-50%%) = %e\t%e\t%e\t%e" % (dEsq,2*dEsq,0.5*dEsq,1.5*dEsq)
print "dTsq    (+/-10%%) = %e\t%e\t%e\t%e\n" % (dTsq,2*dTsq,0.9*dTsq,1.1*dTsq)

print "drhosq*dV  (+/-50%%) = %e\t%e\t%e\t%e" % (drhosq*dV,2*drhosq*dV,0.5*drhosq*dV,1.5*drhosq*dV)
print "drho1sq*dV (+/-50%%) = %e\t%e\t%e\t%e" % (drho1sq*dV,2*drho1sq*dV,0.5*drho1sq*dV,1.5*drho1sq*dV)
print "drho2sq*dV (+/-50%%) = %e\t%e\t%e\t%e" % (drho2sq*dV,2*drho2sq*dV,0.5*drho2sq*dV,1.5*drho2sq*dV)
print "djasq*dV   (+/-10%%) = %e\t%e\t%e\t%e" % (djasq*dV,2*djasq*dV,0.9*djasq*dV,1.1*djasq*dV)
print "dEsq*dV    (+/-50%%) = %e\t%e\t%e\t%e" % (dEsq*dV,2*dEsq*dV,0.5*dEsq*dV,1.5*dEsq*dV)
print "dTsq*dV    (+/-10%%) = %e\t%e\t%e\t%e" % (dTsq*dV,2*dTsq*dV,0.9*dTsq*dV,1.1*dTsq*dV)
