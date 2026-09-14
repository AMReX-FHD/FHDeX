import math

# NIST Chemistry WebBook, SRD 69
# NO2 (nitrogen dioxide):
# https://webbook.nist.gov/cgi/cbook.cgi?Formula=NO2&NoIon=on&Units=SI&cTG=on#Thermo-Gas
# N2O4 (dinitrogen tetroxide):
# https://webbook.nist.gov/cgi/cbook.cgi?Formula=N2O4&NoIon=on&Units=SI&cTG=on#Thermo-Gas
# Cp = heat capacity (J/mol*K)
# H = standard enthalpy (kJ/mol)
# S = standard entropy (J/mol*K)

def compute_thermochem_NO2(temp):

    if (temp<298 or temp>6000):
        print("Warning: temp is not in the range of NO2 thermochem data")

    t = temp/1000.

    if t<1.2:
        A = 16.10857
        B = 75.89525
        C = -54.38740
        D = 14.30777
        E = 0.239423
        F = 26.17464
        G = 240.5386
        H = 33.09502
    else:
        A = 56.82541
        B = 0.738053
        C = -0.144721
        D = 0.009777
        E = -5.459911
        F = 2.846456
        G = 290.5056
        H = 33.09502

    Cp = A + B*t + C*t**2 + D*t**3 + E/t**2
    H = A*t + B*t**2/2 + C*t**3/3 + D*t**4/4 - E/t + F
    S = A*math.log(t) + B*t + C*t**2/2 + D*t**3/3 - E/2/t**2 + G

    return [Cp,H,S]

def compute_thermochem_N2O4(temp):

    if (temp<500 or temp>6000):
        print("Warning: temp is not in the range of N2O4 thermochem data")

    t = temp/1000.

    if t<1:
        A = 34.05274
        B = 191.9845
        C = -151.0575
        D = 44.39350
        E = -0.158949
        F = -8.893428
        G = 293.7724
        H = 9.078988
    else:
        A = 128.6220
        B = 2.524345
        C = -0.520883
        D = 0.036630
        E = -11.55704
        F = -59.22619
        G = 417.0444
        H = 9.078988

    Cp = A + B*t + C*t**2 + D*t**3 + E/t**2
    H = A*t + B*t**2/2 + C*t**3/3 + D*t**4/4 - E/t + F
    S = A*math.log(t) + B*t + C*t**2/2 + D*t**3/3 - E/2/t**2 + G

    return [Cp,H,S]

# NIST Kinetics Database Resources
# N2O4 -> 2NO2
# https://kinetics.nist.gov/kinetics/ReactionSearch?r0=10544726&r1=0&r2=0&r3=0&r4=0&p0=10102440&p1=10102440&p2=0&p3=0&p4=0&expandResults=true&
# Squib: 2004ATK/BAU1461-1738
# Author: Atkinson, R.; Baulch, D.L.; Cox, R.A.; Crowley, J.N.; Hampson, R.F.; Hynes, R.G.; Jenkin, M.E.; Rossi, M.J.; Troe, J.
# Title: Evaluated kinetic and photochemical data for atmospheric chemistry: Volume I - gas phase reactions of Ox, HOx, NOx and SOx species
# Journal: Atmos. Chem. Phys. 4, 1461-1738 (2004)
# Temp   rate (s^-1)
# 250    6.89e4
# 275    7.22e5
# 300    5.12e6

def compute_rate_N2O4_2NO2(temp):

    A = 1.15e16             # s^-1
    Ea = 53.71              # kJ/mol
    R = 8.31446261815324e-3 # kJ/mol*K

    rate = A*math.exp(-Ea/R/temp)
    return rate

# equilibrium constant for N2O4 <-> 2NO2
# Kp = P_NO2^2 / P_N2O4
# some reference values
# Temp      Kp (atm)
# 298.15    0.137
# 308.15    0.294
# 318.15    0.642
# J. Atmos. Chem. 16, 257 (1993), https://doi.org/10.1007/BF00696899
# ----> for comparison (i.e. not used for actual parameter evalution)

def compute_Kp_ref(temp):
    log10Kp = -4.7040-2691.09/temp+6.9628e-4*temp+9.981163*math.log10(temp)-6.6794e-8*temp**2-14.284*(1+85./temp)*math.log10(1+temp/85.)
    Kp = math.pow(10.,log10Kp)
    return Kp

# Study of a Reversible Gas Phase Reaction: An Integrated Physical Chemistry Project
# Chem. Educator 9, 32 (2004)
# http://web.ist.utl.pt/ist12219/data/85.pdf
# kac = 2.2e3*temp^2.3                  # units: M^-1*s^-1
# kd = 2.8e13*temp^1.3)*exp(-6790/temp) # units: s^-1
# reference values at 298K
# k(assoc) = 1.1e9 (M^-1*s^-1)
# k(dissoc) = 5.6e6 (s^-1)
# ----> for comparison (i.e. not used for actual parameter evalution)

def compute_rate_N2O4_2NO2_ref(temp):

    rate = 2.8e13*math.pow(temp,1.3)*math.exp(-6790./temp)
    return rate

# constants in cgs units
Runiv = 8.31446261815324e7  # universal gas constant
kB = 1.38064852e-16         # Boltzmann constant
Navo = 6.0221409e23         # Avogadro constant

# state
Temp = 350.
p_atm = 1.

print("Temp=  \t%f" % Temp)
print("p(atm)=\t%f\n" % p_atm)

# molecular weights
MNO2 = 46.0055
MN2O4 = 2*MNO2
# mass of a molecule
mNO2 = MNO2/Navo
mN2O4 = MN2O4/Navo

print("M_NO2= \t%f" % MNO2)
print("M_N2O4=\t%f" % MN2O4)
print("m_NO2= \t%e" % mNO2)
print("m_N2O4=\t%e\n" % mN2O4)

#####

[CpNO2,HNO2,SNO2] = compute_thermochem_NO2(Temp)
[CpN2O4,HN2O4,SN2O4] = compute_thermochem_N2O4(Temp)

print("Cp_NO2(J/mol*K)= \t%f" % CpNO2)
print("Cp_N2O4(J/mol*K)=\t%f" % CpN2O4)
print("H_NO2(kJ/mol)=   \t%f" % HNO2)
print("H_N2O4(kJ/mol)=  \t%f" % HN2O4)
print("S_NO2(J/mol*K)=  \t%f" % SNO2)
print("S_N2O4(J/mol*K)= \t%f\n" % SN2O4)

fNO2 = 2*CpNO2/(Runiv/1e7)-2
fN2O4 = 2*CpN2O4/(Runiv/1e7)-2
eT0_NO2 = (HNO2-(Runiv/1e10)*Temp)*1e10/MNO2
eT0_N2O4 = (HN2O4-(Runiv/1e10)*Temp)*1e10/MN2O4

print("f_NO2=   \t%f" % fNO2)
print("f_N2O4=  \t%f" % fN2O4)
print("eT0_NO2= \t%e" % eT0_NO2)
print("eT0_N2O4=\t%e\n" % eT0_N2O4)

e0_NO2 = eT0_NO2-0.5*fNO2*Runiv*Temp/MNO2
e0_N2O4 = eT0_N2O4-0.5*fN2O4*Runiv*Temp/MN2O4
cvNO2 = 0.5*fNO2*Runiv/MNO2
cvN2O4 = 0.5*fN2O4*Runiv/MN2O4

print("e0_NO2= \t%e" % e0_NO2)
print("e0_N2O4=\t%e" % e0_N2O4)
print("cv_NO2= \t%e" % cvNO2)
print("cv_N2O4=\t%e\n" % cvN2O4)

#####

dH = 2*HNO2-HN2O4
dS = 2*SNO2-SN2O4
dG = dH-Temp*dS/1000.

print("N2O4 -> 2NO2")
print("dH(kJ/mol)= \t%f" % dH)
print("dS(J/mol*K)=\t%f" % dS)
print("dG(kJ/mol=  \t%f\n" % dG)

Kp = math.exp(-dG/(Runiv/1e10)/Temp)
Kp_ref = compute_Kp_ref(Temp)
Kc = Kp/0.0820574614/Temp

print("Kp(atm)=     \t%f" % Kp)
print("ref: Kp(atm)=\t%f" % Kp_ref)
print("Kc(M)=       \t%f\n" % Kc)

# compute equilibirum pNO2 and pN2O4
# need to solve x^2/(ptot-x)=Kp
pNO2 = 0.5*(-Kp+math.sqrt(Kp*(Kp+4*p_atm)))
pN2O4 = p_atm-pNO2

print("eq pNO2(atm)= \t%e" % pNO2)
print("eq pN2O4(atm)=\t%e\n" % pN2O4)

# number densities
nNO2 = (1.013250e6*pNO2)/kB/Temp
nN2O4 = (1.013250e6*pN2O4)/kB/Temp
ntot = nNO2+nN2O4
xNO2 = nNO2/ntot
xN2O4 = nN2O4/ntot
print("eq nNO2= \t%e" % nNO2)
print("eq nN2O4=\t%e" % nN2O4)
print("eq ntot= \t%e" % ntot)
print("eq mole fractions= %e %e\n" % (xNO2,xN2O4))

# mass densities
rhoNO2 = mNO2*nNO2
rhoN2O4 = mN2O4*nN2O4
rhotot = rhoNO2+rhoN2O4
print("eq rhoNO2= \t%e" % rhoNO2)
print("eq rhoN2O4=\t%e" % rhoN2O4)
print("eq rhotot= \t%e" % rhotot)
print("eq mass fractions= %e %e\n" % (rhoNO2/rhotot,rhoN2O4/rhotot))

k_diss = compute_rate_N2O4_2NO2(Temp)
k_diss_ref = compute_rate_N2O4_2NO2_ref(Temp)
k_asso_c = k_diss/Kc
k_asso_n = k_asso_c*1e3/Navo

print("k_diss=     \t%e" % k_diss)
print("ref: k_diss=\t%e\n" % k_diss_ref)
print("k_asso_c(M^-1*s^-1)=\t%e" % k_asso_c)
print("k_asso_n(cm^3*s^-1)=\t%e\n" % k_asso_n)

k_asso_x = k_asso_n*nNO2**2/xNO2**2
k_diss_x = k_diss*nN2O4/xN2O4

print("k_asso_x=\t%e" % k_asso_x)
print("k_diss_x=\t%e\n" % k_diss_x)

#####

# ** Technical note: Determination of binary gas-phase diffusion
# coefficients of unstable and adsorbing atmospheric trace gases at
# low temperature - arrested flow and twin tube method
# Atmos. Chem. Phys. 20, 3669-3682 (2020)
# https://doi.org/10.5194/acp-20-3669-2020
# ** Viscosity and Thermal Conductivity of the N2O4<->2NO2 System
# J. Chem. Phys. 44, 4643 (1966)
# https://doi.org/10.1063/1.1726692

diamNO2 = 3.765e-8
diamN2O4 = 4.621e-8

print("diam_NO2= \t%e" % diamNO2)
print("diam_N2O4=\t%e\n" % diamN2O4)

#####

dx = 8e-6
dy = 8e-6
dz = 8e-6

dt = 1e-12

Nx = 8
Ny = 8
Nz = 8

dv = dx*dy*dz

print("dx=\t%e\tLx=\t%e" % (dx,Nx*dx))
print("dy=\t%e\tLy=\t%e" % (dy,Ny*dy))
print("dz=\t%e\tLz=\t%e" % (dz,Nz*dz))
print("dv=\t%e" % dv)
print("dt=\t%e\n" % dt)

print("n_NO2*dv= \t%e" % (nNO2*dv))
print("n_N2O4*dv=\t%e\n" % (nN2O4*dv))

rate_diss = k_diss*nN2O4
rate_asso = k_asso_n*nNO2**2

print("eq diss rate (cm^-3*s^-1)=\t%e" % rate_diss)
print("eq asso rate (cm^-3*s^-1)=\t%e\n" % rate_asso)

print("No diss rxn during dt=\t%e" % (rate_diss*dv*dt))
print("No asso rxn during dt=\t%e" % (rate_asso*dv*dt))
