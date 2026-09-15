import sys
import numpy as np
import math

file = sys.argv[1]

Ncellx = 16
Ncelly = 16
Ncellz = 16

dx = 9.360000e-06
dy = 9.360000e-06
dz = 9.360000e-06
dv = dx*dy*dz

Ntot = 90000

meandens = 6.163095e+18
meantheta = 7.466307e-01

data = np.loadtxt(file)
surfcov0 = data[0][1]

dtheta0 = surfcov0-meantheta

r = Ntot/Ncellz/dv
a = r*meantheta
b = -(meandens+r*meantheta*(1-meantheta)+r*meantheta*dtheta0)
c = r*meantheta*(1-meantheta)*dtheta0

# sol1 is the correct one if dtheta0 is small
sol1 = (-b-math.sqrt(b**2-4*a*c))/2/a
sol2 = (-b+math.sqrt(b**2-4*a*c))/2/a

dtheta_eq = Ntot/(Ntot+Ncellz*meandens*dv/meantheta/(1-meantheta))*dtheta0

print("theta0    = %f" % surfcov0)
print("theta_bar = %f" % meantheta)

print("dtheta0   = %f (%.2f%%)" % (dtheta0,dtheta0/meantheta*100))
print("dtheta_eq = %f (%.2f%%)" % (dtheta_eq,dtheta_eq/meantheta*100))
print("sol1 sol2 = %f %f" % (sol1,sol2))

print("theta_eq  = %f" % (meantheta+dtheta_eq))

neq = meandens-Ntot/(Ncellz*dv)*(dtheta_eq-dtheta0)
m1 =  4.651170e-23
print("n_eq = %e" % neq)
print("rhoY1_eq = %e" % (m1*neq))

print("total CO (initial) = %e" % (meandens*dv*Ncellx*Ncelly*Ncellz+Ntot*surfcov0*Ncellx*Ncelly))
print("total CO (eq)      = %e" % (neq*dv*Ncellx*Ncelly*Ncellz+Ntot*(meantheta+dtheta_eq)*Ncellx*Ncelly))
