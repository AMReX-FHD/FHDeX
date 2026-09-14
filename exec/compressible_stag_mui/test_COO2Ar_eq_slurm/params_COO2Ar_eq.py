import math

# this script uses cgs units

Runiv = 8.31446261815324e7  # universal gas constant
kB = 1.38064852e-16         # Boltzmann constant
AVONUM = 6.02214076e23      # Avogadro constant

##########

# molecular weights
M1 = 28.01    # CO
M2 = 32.00    # O2
M3 = 39.95    # Ar

# molecular masses
m1 = M1/AVONUM
m2 = M2/AVONUM
m3 = M3/AVONUM

# mass fractions
Y1 = 0.05
Y2 = 0.05
Y3 = 0.90

# mole fractions
tmp = Y1/M1+Y2/M2+Y3/M3
X1 = Y1/M1/tmp
X2 = Y2/M2/tmp
X3 = Y3/M3/tmp

# average molecular weight
Mavg = 1/(Y1/M1+Y2/M2+Y3/M3)
mavg = Mavg/AVONUM

print("** composition and molecular masses **")

print("- mass fraction Y1 = %f" % Y1)
print("- mass fraction Y2 = %f" % Y2)
print("- mass fraction Y3 = %f" % Y3)

print("- mole fraction X1 = %f" % X1)
print("- mole fraction X2 = %f" % X2)
print("- mole fraction X3 = %f" % X3)

print("- molecular weight M1 = %f" % M1)
print("- molecular weight M2 = %f" % M2)
print("- molecular weight M3 = %f" % M3)
print("- average molecular weight = %f" % Mavg)

print("- molecular mass m1 = %e" % m1)
print("- molecular mass m2 = %e" % m2)
print("- molecular mass m3 = %e" % m3)
print("- average molecular mass = %e\n" % mavg)

##########

temp = 1000.
pres = 1.01325e6            # 1atm

rho = pres*Mavg/Runiv/temp  # total mass density
rho1 = rho*Y1
rho2 = rho*Y2
rho3 = rho*Y3

dens = pres/kB/temp         # total number density
n1 = dens*X1
n2 = dens*X2
n3 = dens*X3

print("** gas state **")

print("- temperature = %f" % temp)
print("- pressure = %e" % pres)

print("- total mass density = %e" % rho)
print("- mass density of spec1 = %e" % rho1)
print("- mass density of spec2 = %e" % rho2)
print("- mass density of spec3 = %e" % rho3)

print("- total number density = %e" % dens)
print("- number density of spec1 = %e" % n1)
print("- number density of spec2 = %e" % n2)
print("- number density of spec3 = %e\n" % n3)

##########

lat_const = 3.0e-8  # lattice constant of Pt(111)
Nucx = 150          # number of unit cells in dx
Nucy = 300          # number of unit cells in dy
Ntot = Nucx*Nucy  # total number of sites per dx*dy (2 sites per unit cell)

dx = Nucx*lat_const*2
dy = Nucy*lat_const
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

print("** FHD cells **")

print("- dx dy dz = %e %e %e" % (dx,dy,dz))
print("- dv = %e" % dv)

print("- Ncellx Ncelly Ncellz = %d %d %d" %(Ncellx,Ncelly,Ncellz))
print("- Lx Ly Lz = %e %e %e" % (Lx,Ly,Lz))

print("- total number of gas molecules in dv = %e" % (N1+N2+N3))
print("- number of gas molecules of spec1 in dv = %e" % N1)
print("- number of gas molecules of spec2 in dv = %e" % N2)
print("- number of gas molecules of spec3 in dv = %e\n" % N3)

print("** surface **")

print("- lattice constant for Pt(111): %e" % lat_const)
print("- total number of sites per dx*dy = %d" % Ntot)
print("- surf_site_num_dens = %e\n" % (Ntot/dx/dy))

##########

rcol1 = math.sqrt(kB*temp/2/math.pi/m1)*n1*(lat_const**2*2)
rcol2 = math.sqrt(kB*temp/2/math.pi/m1)*n2*(lat_const**2*2)

sprob = 1.

rads1 = sprob*rcol1
kads1 = rads1/n1
kdes1 = 1.25e15*math.exp(-1.514/(8.617e-5*temp))

rads2 = sprob*rcol2
kads2 = rads2/n2
kdes2 = 1.25e15*math.exp(-1.514/(8.617e-5*temp))

print("** adsorption/desorption **")

print("- collision rate rcol1 = %e" % rcol1)
print("- collision rate rcol2 = %e" % rcol2)

print("- sticking prob = %f" % sprob)

print("- rads1 = %e (rate)" % rads1)
print("- kads1 = rads1/n1 (rate const) = %e" % kads1)
print("- kdes1 = %e" % kdes1)

print("- rads2 = %e (rate)" % rads2)
print("- kads2 = rads2/n2 (rate const) = %e" % kads2)
print("- kdes2 = %e" % kdes2)

theta1_eq = rads1/(rads1+kdes1)
var_theta1_eq = theta1_eq*(1-theta1_eq)/Ntot
std_theta1_eq_cell = math.sqrt(var_theta1_eq)
std_theta1_eq_sys = math.sqrt(var_theta1_eq/Ncellx/Ncelly)

theta2_eq = rads2/(rads2+math.sqrt(kdes2*rads2))
var_theta2_eq = theta2_eq*(1-theta2_eq)/Ntot
std_theta2_eq_cell = math.sqrt(var_theta2_eq)
std_theta2_eq_sys = math.sqrt(var_theta2_eq/Ncellx/Ncelly)

print("- E[theta1]: %e" % theta1_eq)
print("- Var[theta1]: %e" % var_theta1_eq)
print("- E[theta2]: %e" % theta2_eq)
print("- Var[theta2]: %e" % var_theta2_eq)

print("- std of eq fluct in theta1 (cell): %e (%f%%)" % (std_theta1_eq_cell,std_theta1_eq_cell/theta1_eq*100))
print("- std of eq fluct in theta1 (system): %e (%f%%)\n" % (std_theta1_eq_sys,std_theta1_eq_sys/theta1_eq*100))

print("- std of eq fluct in theta2 (cell): %e (%f%%)" % (std_theta2_eq_cell,std_theta2_eq_cell/theta2_eq*100))
print("- std of eq fluct in theta2 (system): %e (%f%%)\n" % (std_theta2_eq_sys,std_theta2_eq_sys/theta2_eq*100))

##########

dof1 = 5
dof2 = 3
e01 = 0.
e02 = 0.

cv1 = 0.5*dof1*Runiv/M1
cv2 = 0.5*dof2*Runiv/M2
cvmix = Y1*cv1+Y2*cv2

e1 = e01+cv1*temp
e2 = e02+cv2*temp

print("** thermo **")

print("- dof1 dof2 = %f %f" % (dof1,dof2))
print("- e01 e02 = %e %e" % (e01,e02))

print("- cv1 cv2 = %e %e" % (cv1,cv2))
print("- e1 e2 = %e %e\n" % (e1,e2))

##########

drho1sq = rho1**2/N1
drho2sq = rho2**2/N2

drhosq = drho1sq+drho2sq
djasq = rho*kB*temp/dv

dTsq = kB*temp**2/rho/cvmix/dv
drhoEsq = e1**2*drho1sq+e2**2*drho2sq+rho**2*cvmix**2*dTsq

print("** eq fluct **")

print("drhosq = %e" % drhosq)
print("drho1sq = %e" % drho1sq)
print("drho2sq = %e" % drho2sq)
print("drhoEsq = %e" % drhoEsq)
print("dTsq = %e" % dTsq)
print("djasq = %e\n" % djasq)

print("** for amrvis colorbar ranges (cell variance)")

print("drhosq  (+/-50%%) = %e\t%e\t%e\t%e" % (drhosq,2*drhosq,0.5*drhosq,1.5*drhosq))
print("drho1sq (+/-50%%) = %e\t%e\t%e\t%e" % (drho1sq,2*drho1sq,0.5*drho1sq,1.5*drho1sq))
print("drho2sq (+/-50%%) = %e\t%e\t%e\t%e" % (drho2sq,2*drho2sq,0.5*drho2sq,1.5*drho2sq))
print("drhoEsq (+/-25%%) = %e\t%e\t%e\t%e" % (drhoEsq,2*drhoEsq,0.75*drhoEsq,1.25*drhoEsq))
print("dTsq    (+/- 5%%) = %e\t%e\t%e\t%e" % (dTsq,2*dTsq,0.95*dTsq,1.05*dTsq))
print("djasq   (+/- 5%%) = %e\t%e\t%e\t%e\n" % (djasq,2*djasq,0.95*djasq,1.05*djasq))

print("** for amrvis colorbar ranges (structure factor)")

print("drhosq*dv  (+/-50%%) = %e\t%e\t%e\t%e" % (drhosq*dv,2*drhosq*dv,0.5*drhosq*dv,1.5*drhosq*dv))
print("drho1sq*dv (+/-50%%) = %e\t%e\t%e\t%e" % (drho1sq*dv,2*drho1sq*dv,0.5*drho1sq*dv,1.5*drho1sq*dv))
print("drho2sq*dv (+/-50%%) = %e\t%e\t%e\t%e" % (drho2sq*dv,2*drho2sq*dv,0.5*drho2sq*dv,1.5*drho2sq*dv))
print("djasq*dv   (+/-25%%) = %e\t%e\t%e\t%e" % (djasq*dv,2*djasq*dv,0.75*djasq*dv,1.25*djasq*dv))
print("drhoEsq*dv (+/-50%%) = %e\t%e\t%e\t%e" % (drhoEsq*dv,2*drhoEsq*dv,0.5*drhoEsq*dv,1.5*drhoEsq*dv))
print("dTsq*dv    (+/-25%%) = %e\t%e\t%e\t%e\n" % (dTsq*dv,2*dTsq*dv,0.75*dTsq*dv,1.25*dTsq*dv))

##########

# FHD time step size
dt = 1.e-12

# diameters
d1 = 3.46e-8
d2 = 3.63e-8
d12 = (d1+d2)/2

# mass diffusion coefficient
D12 = 3./16*math.sqrt(2*math.pi*kB**3*(m1+m2)/m1/m2)/math.pi/d12**2*temp*math.sqrt(temp)/pres

# magnitude of random advection
u = 3*math.sqrt(kB*temp/mavg/(N1+N2))

print("** time step size and condition numbers **")

print("dt = %e " % dt)

print("D12 = %e" % D12)
print("D12*dt/dx**2 = %e" % (D12*dt/dx**2))

print("u = %e" % u)
print("u*dt/dx = %e\n" % (u*dt/dx))

print("** reaction time scale **")

print("- max mean number of ads events per dxFHD*dyFHD per dt = %e" % (Ntot*rads1*dt))
print("- eq mean number of ads events per dxFHD*dyFHD per dt = %e" % ((1-theta1_eq)*Ntot*rads1*dt))

print("- max mean number of des events per dxFHD*dyFHD per dt = %e" % (Ntot*kdes1*dt))
print("- eq mean number of des events per dxFHD*dyFHD per dt = %e\n" % (theta1_eq*Ntot*kdes1*dt))

print("- 1/(kdes1*theta1_eq) = %e" % (1/(kdes1*theta1_eq)))
print("- 1/(kdes1*theta1_eq) = %e dt" % (1/(kdes1*theta1_eq)/dt))
