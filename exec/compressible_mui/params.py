import math

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
print "- dx (FHD) = %e" % dx
print "- dt (FHD) = %e" % dt

lat_const = dx/10
print "- lattice constant = %e\n" % lat_const

##########

print "*** spec1 ***"
rcol = math.sqrt(0.5*kB*temp/math.pi/m1)*n1*dx**2
print "- collision rate (onto dx^2) = %e" % rcol
print "- number of collisions during dt onto dx^2 = %e\n" % (rcol*dt)

rcol_site = math.sqrt(0.5*kB*temp/math.pi/m1)*n1*lat_const**2
print "- collision rate (onto a site) = %e" % rcol_site

ads_prob = 0.1
eq_cov = 0.4
rads = ads_prob*rcol_site
rdes = (1-eq_cov)/eq_cov*rads
print "- conversion probability (collision->adsorption) = %f" % ads_prob 
print "- equilibrium coverage (for single component) = %f" % eq_cov 
print "- r_ads = %e" % rads 
print "- r_des = %e" % rdes

cads = rads/n1
print "- c_ads = %e\n" % cads

##########

print "*** spec2 ***"
rcol = math.sqrt(0.5*kB*temp/math.pi/m2)*n2*dx**2
print "- collision rate (onto dx^2) = %e" % rcol
print "- number of collisions during dt onto dx^2 = %e\n" % (rcol*dt)

rcol_site = math.sqrt(0.5*kB*temp/math.pi/m2)*n2*lat_const**2
print "- collision rate (onto a site) = %e" % rcol_site

ads_prob = 0.1
eq_cov = 0.5
rads = ads_prob*rcol_site
rdes = (1-eq_cov)/eq_cov*rads
print "- conversion probability (collision->adsorption) = %f" % ads_prob 
print "- equilibrium coverage (for single component) = %f" % eq_cov 
print "- r_ads = %e" % rads 
print "- r_des = %e" % rdes

cads = rads/n2
print "- c_ads = %e" % cads
