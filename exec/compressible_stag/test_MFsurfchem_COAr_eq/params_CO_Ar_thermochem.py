import math

# NIST Chemistry WebBook, SRD 69
# CO (Carbon monoxide)
# https://webbook.nist.gov/cgi/cbook.cgi?Formula=CO&NoIon=on&Units=SI&cTG=on#Thermo-Gas
# Ar (Argon)
# https://webbook.nist.gov/cgi/cbook.cgi?ID=C7440371&Mask=1#Thermo-Gas
# Cp = heat capacity (J/mol*K)
# H = standard enthalpy (kJ/mol)
# S = standard entropy (J/mol*K)

def compute_thermochem(formula,temp):

    if (temp<298 or temp>6000):
        print("Warning: temp is not in the range of thermochem data")
    t = temp/1000.

    if (formula == "CO"):
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
        H_T = A*t + B*t**2/2 + C*t**3/3 + D*t**4/4 - E/t + F #- H
        S_T = A*math.log(t) + B*t + C*t**2/2 + D*t**3/3 - E/2/t**2 + G
        return [Cp,H_T,S_T]

    elif (formula == "Ar"):
        A = 20.78600
        B = 2.825911e-07
        C = -1.464191e-07
        D = 1.092131e-08
        E = -3.661371e-08
        F = -3.661371
        G = 179.9990
        H = 0.000000

        Cp = A + B*t + C*t**2 + D*t**3 + E/t**2
        H_T = A*t + B*t**2/2 + C*t**3/3 + D*t**4/4 - E/t + F #- H
        S_T = A*math.log(t) + B*t + C*t**2/2 + D*t**3/3 - E/2/t**2 + G
        return [Cp,H_T,S_T]

    else:
        print("Warning: not in the range of thermochem data formula")

# constants in cgs units
Runiv = 8.31446261815324e7  # universal gas constant
kB = 1.38064852e-16         # Boltzmann constant
Navo = 6.0221409e23         # Avogadro constant

# state
Temp = 700.
print("Temp   =\t%f" % Temp)

# molecular weights
MCO = 28.0101
MAr = 39.948

# mass of a molecule
mCO = MCO/Navo
mAr = MAr/Navo

print("M_CO  = \t%f" % MCO)
print("M_Ar  = \t%f" % MAr)
print("m_CO  = \t%e" % mCO)
print("m_Ar  = \t%e\n" % mAr)

#####

[CpCO,HCO,SCO] = compute_thermochem("CO",Temp)
[CpAr,HAr,SAr] = compute_thermochem("Ar",Temp)

print("Cp_CO(J/mol*K)  = \t%f" % CpCO)
print("Cp_Ar(J/mol*K)  = \t%f" % CpAr)
print("H_CO(kJ/mol)    = \t%f" % HCO)
print("H_Ar(kJ/mol)    = \t%f" % HAr)
print("S_CO(J/mol*K)   = \t%f" % SCO)
print("S_Ar(J/mol*K)   = \t%f\n" % SAr)

fCO = 2*CpCO/(Runiv/1e7)-2
fAr = 2*CpAr/(Runiv/1e7)-2
eT0_CO = (HCO-(Runiv/1e10)*Temp)*1e10/MCO
eT0_Ar = (HAr-(Runiv/1e10)*Temp)*1e10/MAr

print("f_CO    = \t%f" % fCO)
print("f_Ar    = \t%f" % fAr)
print("eT0_CO  = \t%e" % eT0_CO)
print("eT0_Ar  = \t%e\n" % eT0_Ar)

e0_CO = eT0_CO-0.5*fCO*Runiv*Temp/MCO
e0_Ar = eT0_Ar-0.5*fAr*Runiv*Temp/MAr
cvCO = 0.5*fCO*Runiv/MCO
cvAr = 0.5*fAr*Runiv/MAr

print("e0_CO  = \t%e" % e0_CO)
print("e0_Ar  = \t%e" % e0_Ar)
print("cv_CO  = \t%e" % cvCO)
print("cv_Ar  = \t%e\n" % cvAr)

#####

# https://en.wikipedia.org/wiki/Kinetic_diameter

diamCO = 3.76e-8
diamAr = 3.40e-8

print("diam_CO  = \t%e" % diamCO)
print("diam_Ar  = \t%e\n" % diamAr)

#####
