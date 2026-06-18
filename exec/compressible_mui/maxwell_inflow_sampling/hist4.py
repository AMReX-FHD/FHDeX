import numpy as np
import math
import sys

datafile = sys.argv[1]
histfile = sys.argv[2]
Vzratio = float(sys.argv[3])

kB = 1.38064852e-16;
Navo = 6.0221409e+23;
Temp = 300;
M = 83.8000;
m = M/Navo;
vT = math.sqrt(2*kB*Temp/m);

[vec1,vec2,vec3] = np.loadtxt(datafile,unpack=True)

Nbin = 100

[hist1,bins1] = np.histogram(vec1,bins=Nbin,range=(-3*vT,3*vT),density=True)
[hist2,bins1] = np.histogram(vec2,bins=Nbin,range=(-3*vT,3*vT),density=True)
[hist3,bins1] = np.histogram(vec3,bins=Nbin,range=(-3*vT,3*vT),density=True)

bins1 = 0.5*(bins1[:-1]+bins1[1:])

nbin = len(hist1)
hist4 = np.zeros(nbin)
V = Vzratio*vT

for i in range(nbin):
    v = bins1[i]
    if v<0:
        hist4[i] =  -v*math.exp(-((v-V)/vT)**2)

hist4 /= sum(hist4)*(bins1[1]-bins1[0])

np.savetxt(histfile,np.transpose([bins1,hist1,hist2,hist3,hist4]))
