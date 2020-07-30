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

print "- total number density = %.3e" % ntot
print "- number density of spec1 = %.3e" % n1
print "- number density of spec2 = %.3e" % n2
print "- number density of spec3 = %.3e" % n3
print "- number density of spec4 = %.3e\n" % n4

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

Ncol1 = math.sqrt(0.5*kB*temp/math.pi/m1)*n1*dx**2*dt
Ncol2 = math.sqrt(0.5*kB*temp/math.pi/m2)*n2*dx**2*dt
Ncol3 = math.sqrt(0.5*kB*temp/math.pi/m3)*n3*dx**2*dt
Ncol4 = math.sqrt(0.5*kB*temp/math.pi/m4)*n4*dx**2*dt

print "- number of collisions by gas molecules of spec1 per dx^2 per dt = %.3e" % Ncol1
print "- number of collisions by gas molecules of spec2 per dx^2 per dt = %.3e" % Ncol2
print "- number of collisions by gas molecules of spec3 per dx^2 per dt = %.3e" % Ncol3
print "- number of collisions by gas molecules of spec4 per dx^2 per dt = %.3e\n" % Ncol4

##########

lat_const = dx/10
print "- lattice constant = %e\n" % lat_const

##########

ads_prob = 0.1
eq_cov = 0.1
rcol_site = math.sqrt(0.5*kB*temp/math.pi/m1)*n1*lat_const**2
rads = ads_prob*rcol_site
rdes = (1-eq_cov)/eq_cov*rads
cads = rads/n1

print "*** spec1 ***"
print "- collision rate onto a site = %e" % rcol_site
print "- conversion probability (collision->adsorption) = %f" % ads_prob 
print "- equilibrium coverage (for single component) = %f" % eq_cov 
print "- r_ads = %e" % rads 
print "- r_des = %e" % rdes
print "- c_ads = %e\n" % cads

##########

ads_prob = 0.1
eq_cov = 0.1
rcol_site = math.sqrt(0.5*kB*temp/math.pi/m2)*n2*lat_const**2
rads = ads_prob*rcol_site
rdes = (1-eq_cov)/eq_cov*rads
cads = rads/n2

print "*** spec2 ***"
print "- collision rate onto a site = %e" % rcol_site
print "- conversion probability (collision->adsorption) = %f" % ads_prob 
print "- equilibrium coverage (for single component) = %f" % eq_cov 
print "- r_ads = %e" % rads 
print "- r_des = %e" % rdes
print "- c_ads = %e\n" % cads

##########

ads_prob = 0.1
eq_cov = 0.05
rcol_site = math.sqrt(0.5*kB*temp/math.pi/m3)*n3*lat_const**2
rads = ads_prob*rcol_site
rdes = (1-eq_cov)/eq_cov*rads
cads = rads/n3

print "*** spec3 ***"
print "- collision rate (onto a site) = %e" % rcol_site
print "- conversion probability (collision->adsorption) = %f" % ads_prob 
print "- equilibrium coverage (for single component) = %f" % eq_cov 
print "- r_ads = %e" % rads 
print "- r_des = %e" % rdes
print "- c_ads = %e\n" % cads

##########

ads_prob = 0.1
eq_cov = 0.05
rcol_site = math.sqrt(0.5*kB*temp/math.pi/m4)*n4*lat_const**2
rads = ads_prob*rcol_site
rdes = (1-eq_cov)/eq_cov*rads
cads = rads/n4

print "*** spec4 ***"
print "- collision rate (onto a site) = %e" % rcol_site
print "- conversion probability (collision->adsorption) = %f" % ads_prob 
print "- equilibrium coverage (for single component) = %f" % eq_cov 
print "- r_ads = %e" % rads 
print "- r_des = %e" % rdes
print "- c_ads = %e\n" % cads
