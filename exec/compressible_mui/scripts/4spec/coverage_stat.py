import numpy as np

nrun = 64
cc1 = []
cc2 = []
cc3 = []
cc4 = []

for i in range(nrun):
  filename = "RUN%d/res.coverage" % (i+1)
  [tt,c1,c2,c3,c4] = np.loadtxt(filename,usecols=(0,1,2,3,4),unpack=True)
  cc1.append(c1)
  cc2.append(c2)
  cc3.append(c3)
  cc4.append(c4)

mc1 = np.mean(cc1,axis=0)
mc2 = np.mean(cc2,axis=0)
mc3 = np.mean(cc3,axis=0)
mc4 = np.mean(cc4,axis=0)

vc1 = np.std(cc1,axis=0,ddof=1)
vc2 = np.std(cc2,axis=0,ddof=1)
vc3 = np.std(cc3,axis=0,ddof=1)
vc4 = np.std(cc4,axis=0,ddof=1)

np.savetxt("res.coverage_stat",np.transpose([tt,mc1,vc1,mc2,vc2,mc3,vc3,mc4,vc4]))
