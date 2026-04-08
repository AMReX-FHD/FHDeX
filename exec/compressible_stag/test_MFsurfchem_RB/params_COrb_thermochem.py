import math

# NIST Chemistry WebBook, SRD 69
# CO (carbon monoxide):
# https://webbook.nist.gov/cgi/cbook.cgi?ID=C630080&Mask=1
# Cp0 = heat capacity (J/mol*K)
# H0 = standard enthalpy (kJ/mol)
# S0 = standard entropy (J/mol*K)

def compute_thermochem_CO(temp):

    if (temp<298 or temp>6000):
        print("Warning: temp is not in the range of CO thermochem data")

    t = temp/1000.

    if t<1.3:
        A = 25.56759
        B = 6.096130
        C = 4.054656
        D = -2.671301
        E = 0.131021
        F = -118.0089
        G = 227.3665
        H = -110.5271
    else:
        A = 35.15070
        B = 1.300095
        C = -0.205921
        D = 0.013550
        E = -3.282780
        F = -127.8375
        G = 231.7120
        H = -110.5271

    Cp = A + B*t + C*t**2 + D*t**3 + E/t**2
    H = A*t + B*t**2/2 + C*t**3/3 + D*t**4/4 - E/t + F
    S = A*math.log(t) + B*t + C*t**2/2 + D*t**3/3 - E/2/t**2 + G

    return [Cp,H,S]

# constants in cgs units
Runiv = 8.31446261815324e7  # universal gas constant
kB = 1.38064852e-16         # Boltzmann constant
Navo = 6.0221409e23         # Avogadro constant

# state
Temp = 1000.
p_atm = 1.

print("Temp=  \t%f" % Temp)
print("p(atm)=\t%f\n" % p_atm)

# molecular weight
MCO = 28.0101

# mass of a molecule
mCO = MCO/Navo

print("M_CO = \t%f" % MCO)
print("m_CO = \t%e" % mCO)

#####

[CpCO,HCO,SCO] = compute_thermochem_CO(Temp)

print("Cp_CO(J/mol*K)= \t%f" % CpCO)
print("H_CO(kJ/mol)=   \t%f" % HCO)
print("S_CO(J/mol*K)=  \t%f" % SCO)

fCO = 2*CpCO/(Runiv/1e7)-2
eT0_CO = (HCO-(Runiv/1e10)*Temp)*1e10/MCO

print("f_CO =   \t%f" % fCO)
print("eT0_CO = \t%e" % eT0_CO)

e0_CO = eT0_CO-0.5*fCO*Runiv*Temp/MCO
cvCO = 0.5*fCO*Runiv/MCO

print("e0_CO = \t%e" % e0_CO)
print("cv_CO= \t%e" % cvCO)
