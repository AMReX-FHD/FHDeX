import sys
import numpy as np
import math

Ncellx = 4
Ncelly = 4
Ncellz = 512

dx = 8.400000e-06
dy = 7.274613e-06
dz = 8.400000e-06
dv = dx*dy*dz

Ntot = 90000

mean_rho1 = 2.383472e-05
mean_rho2 = 4.528596e-04
mean_rho = mean_rho1 + mean_rho2
print("mean_rho1 = %e" % mean_rho1)
print("mean_rho2 = %e" % mean_rho2)
print("mean_rho  = %e" % mean_rho)

mean_theta1 = 2.053333e-01
m1 = 4.651170e-23
mean_dens1 = mean_rho1/m1
init_dens1 = mean_dens1+Ntot*mean_theta1/Ncellz/dv
print("mean_theta1 = %e" % mean_theta1)
print("mean_dens1 = %e" % mean_dens1)
print("init_dens1 = %e" % init_dens1)

init_N1 = init_dens1*dv*Ncellx*Ncelly*Ncellz
eq_N1 = mean_dens1*dv*Ncellx*Ncelly*Ncellz+Ntot*mean_theta1*Ncellx*Ncelly
print("init_N1 = %e" % init_N1)
print("eq_N1   = %e" % eq_N1)

init_rho1 = init_dens1*m1
init_rho = init_rho1 + mean_rho2
print("init_rho1 = %e" % init_rho1)
print("init_rho  = %e" % init_rho)

init_Y1 = init_rho1/init_rho
init_Y2 = 1 - init_Y1
print("Y1 Y2 (init) = %e %e" % (init_Y1,init_Y2))

