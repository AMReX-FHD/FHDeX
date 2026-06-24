import sys
import numpy as np
import math

file = sys.argv[1]

Ncellx = 16
Ncelly = 16
Ncellz = 16

dx = 9.36e-6
dy = dx
dz = dx
dv = dx*dy*dz

Ntot = 90000

mass = 28.01/6.02214076e23

[surfcov,rhoY] = np.loadtxt(file,unpack=True,usecols=(1,3))

# compute the total number of CO molecules in the system (gas+surface)
print(Ncellx*Ncelly*Ntot*surfcov[0]+sum(rhoY)/mass*dv*Ncellx*Ncelly)
