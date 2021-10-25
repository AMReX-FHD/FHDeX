import numpy as np

vv = np.loadtxt("res.frac",usecols=(0),unpack=True)

meancov = 2.053333e-01
Ntot = 300*300

print "mean     = %f" % (np.mean(vv))
print "expected = %f" % (Ntot*meancov)
print "std      = %f" % (np.std(vv,ddof=1))
print "expected = %f" % (np.sqrt(meancov*(1-meancov)*Ntot))
