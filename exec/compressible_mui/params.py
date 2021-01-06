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

lat_const = 2.8e-8
n1x = 300
n1y = 150

dx = n1x*lat_const
dy = n1y*lat_const*math.sqrt(3)
dz = dx
dV = dx*dy*dz

n2x = 4
n2y = 4
n2z = 4

Lx = n2x*dx
Ly = n2y*dy
Lz = n2z*dz

print "- lattice constant for Pt(111): %e" % lat_const
print "- dx (FHD) = %e" % dx
print "- dy (FHD) = %e" % dy
print "- dz (FHD) = %e" % dz
print "- dV (FHD) = %e" % dV
print "- Lx = %e" % Lx
print "- Ly = %e" % Ly
print "- Lz = %e\n" % Lz

##########

n1s = 2*n1x*n1y
print "- total number of sites per dxFHD*dyFHD = %d" % n1s 

Ntot = ntot*dV
N1 = n1*dV
N2 = n2*dV

print "- total number of gas molecules in dV = %.3e" % Ntot
print "- number of gas molecules of spec1 in dV = %.3e" % N1
print "- number of gas molecules of spec2 in dV = %.3e\n" % N2

##########

p1torr = pres*X1/1.01325e6*760
k1ads = 2e5*p1torr

print "- p1 in torr = %e" % p1torr
print "- k1ads = %e" % k1ads

k1des = 1.25e15*math.exp(-1.514/(8.617e-5*temp))
print "- temp in K = %e" % temp
print "- k1des = %e" % k1des

theta1 = k1ads/(k1ads+k1des)
print "- eq coverage of spec1: %e\n" % theta1

##########

dt = 1e-12

print "- dt (FHD) = %e" % dt
print "- max mean number of ads events per dxFHD*dyFHD per dt = %e" % (n1s*k1ads*dt)
print "- max mean number of des events per dxFHD*dyFHD per dt = %e" % (n1s*k1des*dt)




