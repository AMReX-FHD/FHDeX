import math

# constants in cgs units
Runiv = 8.31446261815324e7  # universal gas constant
kB = 1.38064852e-16         # Boltzmann constant
Navo = 6.0221409e23         # Avogadro constant

# molecular weights
MNO2 = 46.0055
MN2O4 = 2*MNO2
# mass of a molecule
mNO2 = MNO2/Navo
mN2O4 = MN2O4/Navo

# N2O4 <-> 2NO2
#
# some reference values for Kp = P_NO2^2 / P_N2O4
# Temp      Kp (atm)
# 298.15    0.137
# 308.15    0.294
# 318.15    0.642
#
# J. Atmos. Chem. 16, 257 (1993), https://doi.org/10.1007/BF00696899

Temp = 298.15
ptot_atm = 1.

# compute Kp
log10Kp = -4.7040-2691.09/Temp+6.9628e-4*Temp+9.981163*math.log10(Temp)-6.6794e-8*Temp**2-14.284*(1+85./Temp)*math.log10(1+Temp/85.)
Kp = math.pow(10.,log10Kp)

# compute equilibirum pNO2 and pN2O4
# need to solve x^2/(ptot-x)=Kp
pNO2 = 0.5*(-Kp+math.sqrt(Kp*(Kp+4*ptot_atm)))
pN2O4 = ptot_atm-pNO2

print("Temp=      \t%f" % Temp)
print("Kp(atm)=   \t%e" % Kp)
print("pNO2(atm)= \t%e" % pNO2)
print("pN2O4(atm)=\t%e\n" % pN2O4)

# number densities
nNO2 = (1.013250e6*pNO2)/kB/Temp
nN2O4 = (1.013250e6*pN2O4)/kB/Temp
ntot = nNO2+nN2O4
print("nNO2= \t%e" % nNO2)
print("nN2O4=\t%e" % nN2O4)
print("ntot= \t%e" % ntot)
print("mole fractions= %e %e\n" % (nNO2/ntot,nN2O4/ntot))

# mass densities
rhoNO2 = mNO2*nNO2
rhoN2O4 = mN2O4*nN2O4
rhotot = rhoNO2+rhoN2O4
print("rhoNO2= \t%e" % rhoNO2)
print("rhoN2O4=\t%e" % rhoN2O4)
print("rhotot= \t%e" % rhotot)
print("mass fractions= %e %e\n" % (rhoNO2/rhotot,rhoN2O4/rhotot))

# Study of a Reversible Gas Phase Reaction: An Integrated Physical Chemistry Project
# Chem. Educator 9, 32 (2004)
# http://web.ist.utl.pt/ist12219/data/85.pdf
#
# reference values at 298K
# kac = 1.1e9 M^-1 s^-1
# kdc = 5.6e6 s^-1

# association (i.e. dimerization) rate constant conversion
kac = 2.2e3*math.pow(Temp,2.3)  # units: M^-1 s^-1
kan = kac*1e3/Navo              # units: cm^3 s^-1
print("kac=\t%e" % kac)
print("kan=\t%e" % kan)

# dissociation rate constant computed from equilibrium state and forward rate
kdn = kan*nNO2**2/nN2O4         # units: s^-1
kdc = kdn                       # units: s^-1
print("kdc=\t%e" % kdc)
print("kdn=\t%e" % kdn)
