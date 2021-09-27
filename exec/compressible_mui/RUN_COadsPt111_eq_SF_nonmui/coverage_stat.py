import sys
import numpy as np

assert len(sys.argv)==2
infile = sys.argv[1] 

[tt,nvac,nCO] = np.loadtxt(infile,usecols=(0,1,2),unpack=True)

assert np.min(nvac+nCO) == np.max(nvac+nCO)
Nsite = nvac[0]+nCO[0]

t0 = int(len(tt)*0.1)
nCO_ss = nCO[t0:]
cov =  np.average(nCO_ss/float(Nsite))

print "Nsite = %d" % Nsite
print "mean cov = %e" % cov
print "var of nCO = %e" % np.var(nCO_ss,ddof=1)
print "Nsite*cov = %e" % (Nsite*cov)
print "Nsite*cov*(1-cov) = %e" % (Nsite*cov*(1-cov))
