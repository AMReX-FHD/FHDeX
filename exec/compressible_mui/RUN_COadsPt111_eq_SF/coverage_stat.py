import sys
import numpy as np

assert len(sys.argv)==2
infile = sys.argv[1]

[tt,nvac,nCO] = np.loadtxt(infile,usecols=(0,1,2),unpack=True)

assert np.min(nvac+nCO) == np.max(nvac+nCO)
Nsite = nvac[0]+nCO[0]
#t0 = int(len(tt)*0.1)
t0 = 0
theta = nCO[t0:]/float(Nsite)
Ntpt = len(theta)
dt = tt[1]-tt[0];
print "No sites = %d" % Nsite
print "No time points: %d (skipping initial %d points)" % (Ntpt,t0)
print "dt = %e\n" % dt

ra = 7.562496e+06
rd = 2.926784e+07
meancov = ra/(ra+rd)

T = 1e-12*1e6
a = (ra+rd)*T
ratio = 1-2*(a-1+np.exp(-a))/a**2

print "mean(theta) = %f" % (np.average(theta))
print "expected    = %f\n" % meancov

print "var(theta)  = %e" % (np.var(theta,ddof=1))
print "expected    = %e" % (ratio*meancov*(1-meancov)/Nsite)
print "** ratio    = %f\n" % ratio

print "var2(theta) = %e" % (np.average(np.square(theta-meancov)))
print "expected    = %e" % (meancov*(1-meancov)/Nsite)
