import sys
import math
import numpy as np
# this script uses cgs units

Runiv = 8.31446261815324e7  # universal gas constant
kB = 1.381e-16#64852e-16         # Boltzmann constant
AVONUM = 6.02214076e23      # Avogadro constant
h_cgs = 6.6261e-27#76e-27
eV_cgs = 1.602e-12
temp = 600.
pres = 10*1.01325e6            # 1atm

##########

# molecular weights
M1 = 28.01    # CO
M2 = 32.00    # O2
M3 = 44.01    # CO2
M4 = 39.95    # Ar
# molecular masses
m1 = M1/AVONUM
m2 = M2/AVONUM
m3 = M3/AVONUM
m4 = M4/AVONUM

# mass fractions
Y1 = 0.20
Y2 = 0.10
Y3 = 0.05
Y4 = 0.65

tmp = Y1/M1+Y2/M2+Y3/M3+Y4/M4
X1 = Y1/M1/tmp
X2 = Y2/M2/tmp
X3 = Y3/M3/tmp
X4 = Y4/M4/tmp

# average molecular weight
Mavg = 1/(Y1/M1+Y2/M2+Y3/M3+Y4/M4)
mavg = Mavg/AVONUM

print("** composition and molecular masses **")
print("- mass fraction Y1 = %f" % Y1)
print("- mass fraction Y2 = %f" % Y2)
print("- mass fraction Y3 = %f" % Y3)
print("- mass fraction Y4 = %f" % Y4)
print("- mole fraction X1 = %f" % X1)
print("- mole fraction X2 = %f" % X2)
print("- mole fraction X3 = %f" % X3)
print("- mole fraction X4 = %f" % X4)

print("- molecular weight M1 = %f" % M1)
print("- molecular weight M2 = %f" % M2)
print("- molecular weight M3 = %f" % M3)
print("- molecular weight M4 = %f" % M4)
print("- average molecular weight = %f\n" % Mavg)

print("- molecular mass m1 = %e" % m1)
print("- molecular mass m2 = %e" % m2)
print("- molecular mass m3 = %e" % m3)
print("- molecular mass m4 = %e" % m4)
print("- average molecular mass = %e\n" % mavg)

##########
rho = pres*Mavg/Runiv/temp  # total mass density
rho1 = rho*Y1
rho2 = rho*Y2
rho3 = rho*Y3
rho4 = rho*Y4

dens = pres/kB/temp         # total number density
n1 = dens*X1
n2 = dens*X2
n3 = dens*X3
n4 = dens*X4

# partial pressure
p1 = pres*X1
p2 = pres*X2
p3 = pres*X3
p4 = pres*X4

##### constants #####
cm_Hz_factor = 2.997924e10
cm_J_factor = 1.986e-23
J_eV_factor = 6.242e18
eV_J_factor = 1.602e-19

kB_si   = 1.380649e-23 # J/K
h_si    = 6.62607015e-34 # J*Hz^-1
hbar_si = 1.054572e-34 # Planck in J*s

m1_si = 4.65e-26 # 28.0101*10e-3 / 6.022e23 # CO
m2_si = 5.31e-26 # 31.9988*10e-3 / 6.022e23 # O2

p1_si = p1/10
p2_si = p2/10

# O2 parameters
sigma2 = 2
B2 = 1.445622*cm_J_factor # J
omega2 = 1580.1932*cm_Hz_factor # Hz
I2 = 3

# CO
sigma1 = 1 # unitless
B1 = 1.9302*cm_J_factor # J
omega1 = 2169.52*cm_Hz_factor # Hz
I1 = 1 # unitless


delta_mu1 = -1*kB_si*temp*(np.log(np.power(2*math.pi*m1_si/(h_si**2),1.5)*np.power(kB_si*temp,2.5)/p1_si) + np.log(kB_si*temp/(sigma1 * B1)) - np.log(1 - np.exp(-h_si*omega1/(kB_si*temp))) + np.log(I1))
delta_mu2 = -1*kB_si*temp*(np.log(np.power(2*math.pi*m2_si/(h_si**2),1.5)*np.power(kB_si*temp,2.5)/p2_si) + np.log(kB_si*temp/(sigma2 * B2)) - np.log(1 - np.exp(-h_si*omega2/(kB_si*temp))) + np.log(I2))
delta_mu1 = delta_mu1 * J_eV_factor
delta_mu2 = delta_mu2 * J_eV_factor

print("** delta_mu_CO = {}, delta_mu_O2 = {}\n".format(delta_mu1,delta_mu2))

print("** gas state **")

print("- temperature = %f" % temp)
print("- pressure = %e" % pres)

print("- partial pressure of spec1 = %e" % p1)
print("- partial pressure of spec2 = %e" % p2)
print("- partial pressure of spec3 = %e" % p3)
print("- partial pressure of spec4 = %e" % p4)

print("- total mass density = %e" % rho)
print("- mass density of spec1 = %e" % rho1)
print("- mass density of spec2 = %e" % rho2)
print("- mass density of spec3 = %e" % rho3)
print("- mass density of spec4 = %e\n" % rho4)

print("- total number density = %e" % dens)
print("- number density of spec1 = %e" % n1)
print("- number density of spec2 = %e" % n2)
print("- number density of spec3 = %e" % n3)
print("- number density of spec4 = %e\n" % n4)

##########

lat_const_x = 6.43e-8  # lattice constant of Ru(110)
lat_const_y = 3.12e-8  # lattice constant of Ru(110)

Nucx = 150          # number of unit cells in dx
Nucy = 300          # number of unit cells in dy
Ntot = Nucx*Nucy  # total number of sites per dx*dy (2 sites per unit cell)

dx = Nucx*lat_const_x
dy = Nucy*lat_const_y
dz = dx
dv = dx*dy*dz

Ncellx = 4
Ncelly = 4
Ncellz = 4
Ncell = Ncellx*Ncelly*Ncellz

Lx = Ncellx*dx
Ly = Ncelly*dy
Lz = Ncellz*dz

N1 = n1*dv
N2 = n2*dv
N3 = n3*dv
N4 = n4*dv

print("** FHD cells **")

print("- dx dy dz = %e %e %e" % (dx,dy,dz))
print("- dv = %e" % dv)

print("- Ncellx Ncelly Ncellz = %d %d %d" %(Ncellx,Ncelly,Ncellz))
print("- Lx Ly Lz = %e %e %e" % (Lx,Ly,Lz))

print("- total number of gas molecules in dv = %e" % (N1+N2+N3+N4))
print("- number of gas molecules of spec1 in dv = %e" % N1)
print("- number of gas molecules of spec2 in dv = %e" % N2)
print("- number of gas molecules of spec3 in dv = %e" % N3)
print("- number of gas molecules of spec4 in dv = %e\n" % N4)

print("** surface **")

print("- lattice constant for Ru(110): %e, %e\n" % (lat_const_x, lat_const_y))
print("- total number of sites per dx*dy = %d" % Ntot)
print("- surf_site_num_dens = %e\n" % (Ntot/dx/dy))

##########

print("** FHD input **")
print("- Ncellx Ncelly Ncellz = %d %d %d" %(Ncellx,Ncelly,Ncellz))
print("- Lx Ly Lz = %e %e %e" % (Lx,Ly,Lz))
print("- total mass density = %e" % rho)
print("- mass fraction Y1 = %f" % Y1)
print("- mass fraction Y2 = %f" % Y2)
print("- mass fraction Y3 = %f" % Y3)
print("- mass fraction Y4 = %f\n" % Y4)

##########

Atot = lat_const_x*lat_const_y
A1 = Atot/2
A2 = Atot/2

sprob_ads = 1.0
sprob_des = 0.5
sprob_dif = 1.0

#print("** adsorption/desorption **")

kads1 = sprob_ads*p1*A1/(math.sqrt(2*math.pi*m1*kB*temp))
kdes1 = kads1*math.exp((-delta_mu1-1.60)*eV_cgs/(kB*temp))

kads2 = sprob_ads*p1*A1/(math.sqrt(2*math.pi*m1*kB*temp))
kdes2 = kads2*math.exp((-delta_mu1-1.30)*eV_cgs/(kB*temp))

kads3 = sprob_ads*p2*A2/(math.sqrt(2*math.pi*m2*kB*temp))
kdes3 = kads3*math.exp((-delta_mu2-4.60)*eV_cgs/(kB*temp))

kads4 = sprob_ads*p2*A2/(math.sqrt(2*math.pi*m2*kB*temp))
kdes4 = kads4*math.exp((-delta_mu2-3.30)*eV_cgs/(kB*temp))

kads5 = sprob_ads*p2*A2/(math.sqrt(2*math.pi*m2*kB*temp))
kdes5 = kads5*math.exp((-delta_mu2-2.00)*eV_cgs/(kB*temp))

krxn1 = sprob_des*(kB*temp/h_cgs)*math.exp(-1.50*eV_cgs/(kB*temp))
krxn2 = sprob_des*(kB*temp/h_cgs)*math.exp(-0.80*eV_cgs/(kB*temp))
krxn3 = sprob_des*(kB*temp/h_cgs)*math.exp(-1.20*eV_cgs/(kB*temp))
krxn4 = sprob_des*(kB*temp/h_cgs)*math.exp(-0.90*eV_cgs/(kB*temp))

kdif1 = sprob_dif*(kB*temp/h_cgs)*math.exp(-0.60*eV_cgs/(kB*temp))
kdif2 = sprob_dif*(kB*temp/h_cgs)*math.exp(-1.60*eV_cgs/(kB*temp))
kdif3 = sprob_dif*(kB*temp/h_cgs)*math.exp(-1.30*eV_cgs/(kB*temp))
kdif4 = sprob_dif*(kB*temp/h_cgs)*math.exp(-1.70*eV_cgs/(kB*temp))
kdif5 = sprob_dif*(kB*temp/h_cgs)*math.exp(-0.70*eV_cgs/(kB*temp))
kdif6 = sprob_dif*(kB*temp/h_cgs)*math.exp(-2.30*eV_cgs/(kB*temp))
kdif7 = sprob_dif*(kB*temp/h_cgs)*math.exp(-1.00*eV_cgs/(kB*temp))
kdif8 = sprob_dif*(kB*temp/h_cgs)*math.exp(-1.60*eV_cgs/(kB*temp))

print("** Event & Rate command **")
print("# siteA = br / siteB = cus / spec1 = CO / spec2 = O2")
print("# Adsorption")
if len(sys.argv) == 2 and sys.argv[1] == "rate":
    print("event       4 siteA vac   %e spec1 rate" % kads1)
    print("event       4 siteB vac   %e spec1 rate" % kads2)
    print("event       6 siteA siteA vac   vac   %e spec2 spec2 spec2 rate #half" % (kads3/2))
    print("event       6 siteA siteB vac   vac   %e spec2 spec2 spec2 rate" % kads4)
    print("event       6 siteB siteA vac   vac   %e spec2 spec2 spec2 rate" % kads4)
    print("event       6 siteB siteB vac   vac   %e spec2 spec2 spec2 rate #half" % (kads5/2))
else:
    print("event       4 siteA vac   %e spec1 beta 0.5" % (kads1/n1))
    print("event       4 siteB vac   %e spec1 beta 0.5" % (kads2/n1))
    print("event       6 siteA siteA vac   vac   %e spec2 spec2 spec2 beta 0.5 #half" % ((kads3/2)/n2))
    print("event       6 siteA siteB vac   vac   %e spec2 spec2 spec2 beta 0.5" % (kads4/n2))
    print("event       6 siteB siteA vac   vac   %e spec2 spec2 spec2 beta 0.5" % (kads4/n2))
    print("event       6 siteB siteB vac   vac   %e spec2 spec2 spec2 beta 0.5 #half" % ((kads5/2)/n2))
print("\n# Desorption")
print("event       5 siteA spec1 %e vac" % kdes1)
print("event       5 siteB spec1 %e vac" % kdes2)
print("event       7 siteA siteA spec2 spec2 %e vac   vac   spec2 #half" % (kdes3/2))
print("event       7 siteA siteB spec2 spec2 %e vac   vac   spec2" % kdes4)
print("event       7 siteB siteA spec2 spec2 %e vac   vac   spec2" % kdes4)
print("event       7 siteB siteB spec2 spec2 %e vac   vac   spec2 #half" % (kdes5/2))
print("\n# Reaction")
print("event       8 siteA siteA spec1 spec2 %e vac   vac   spec3" % krxn1)
print("event       8 siteA siteB spec1 spec2 %e vac   vac   spec3" % krxn2)
print("event       8 siteB siteA spec1 spec2 %e vac   vac   spec3" % krxn3)
print("event       8 siteB siteB spec1 spec2 %e vac   vac   spec3" % krxn4)
print("\n# Diffusion")
print("event       2 siteA siteA spec1 vac   %e vac   spec1" % kdif1)
print("event       2 siteA siteB spec1 vac   %e vac   spec1" % kdif2)
print("event       2 siteB siteA spec1 vac   %e vac   spec1" % kdif3)
print("event       2 siteB siteB spec1 vac   %e vac   spec1" % kdif4)
print("event       2 siteA siteA spec2 vac   %e vac   spec2" % kdif5)
print("event       2 siteA siteB spec2 vac   %e vac   spec2" % kdif6)
print("event       2 siteB siteA spec2 vac   %e vac   spec2" % kdif7)
print("event       2 siteB siteB spec2 vac   %e vac   spec2" % kdif8)
