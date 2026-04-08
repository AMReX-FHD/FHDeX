import math
import numpy as np

outfile = "res.coverage_det"

ka = 7.562496e+06
kd = 2.926784e+07

cov0 = 0.

dt = 1e-12
Nstep = 100000

tt = np.zeros(Nstep+1)
cc = np.zeros(Nstep+1)

for i in range(Nstep+1):
  t = i*dt
  tt[i] = t
  cc[i] = math.exp(-(ka+kd)*t)*cov0 + ka/(ka+kd)*(1-math.exp(-(ka+kd)*t))

np.savetxt(outfile,np.transpose([tt,cc]))
