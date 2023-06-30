import math
import numpy as np

kB = 1.38064852e-16 

# spec1 = NO2
# spec2 = N2O4
# fwd reaction: 2 NO2 -> N2O4
# bwd reaction: N2O4 -> 2 NO2

m1 = 7.639393e-23
m2 = 2*m1 

T_eq = 350.

e0_NO2  =  4.684542e+09
e0_N2O4 = -1.730892e+09
cv_NO2  =  6.578106e+06
cv_N2O4 =  8.154909e+06

e0hat_NO2 = m1/kB*e0_NO2
e0hat_N2O4 = m2/kB*e0_N2O4
cphat_NO2 = m1/kB*cv_NO2+1
cphat_N2O4 = m2/kB*cv_N2O4+1

print("e0hat_NO2 = %e" % e0hat_NO2)
print("e0hat_N2O4 = %e" % e0hat_N2O4)
print("cphat_NO2 = %e" % cphat_NO2)
print("cphat_N2O4 = %e" % cphat_N2O4)

kf_eq = 5.194896e+26
kb_eq = 2.326169e+27

rho_eq = 1.855335e-03
Y1_eq = 7.267597e-01
Y2_eq = 1-Y1_eq
rho1_eq = rho_eq*Y1_eq
rho2_eq = rho_eq*Y2_eq
p1_eq = rho1_eq/m1*kB*T_eq
p2_eq = rho2_eq/m2*kB*T_eq
p_std = 1.013250e6

print("fwd rate at T_eq = %e" % (kf_eq*(p1_eq/p_std)**2))
print("bwd rate at T_eq = %e" % (kb_eq*(p2_eq/p_std)))

def kf(T):
    return kf_eq*math.exp(2*e0hat_NO2*(1/T-1/T_eq))*math.pow(T/T_eq,-2*cphat_NO2)

def kb(T):
    return kb_eq*math.exp(e0hat_N2O4*(1/T-1/T_eq))*math.pow(T/T_eq,-cphat_N2O4)

for T in np.linspace(0.95*T_eq,1.05*T_eq,num=1000):
    print("%e\t%e\t%e" % (T,kf(T)/kf_eq,kb(T)/kb_eq))
