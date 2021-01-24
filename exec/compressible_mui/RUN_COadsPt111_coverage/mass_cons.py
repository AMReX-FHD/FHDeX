import numpy as np

plot_int = 200

mCO = 28.01/6.0221409e23
dV = 5.132967e-16
nx = 4
ny = 4

print("** plot_int=%d" % plot_int)
print("** nx=%d" % nx)
print("** ny=%d" % ny)

Nsurf = np.loadtxt("res.coverage",unpack=True,usecols=(2,))
npts = len(Nsurf)/plot_int

M1 = np.zeros(npts)
M2 = np.zeros(npts)

for i in range(npts):
    a = "%d" % (plot_int*(i+1))
    filename = "havg"+a.zfill(9)
    rhoY1z = np.loadtxt(filename,unpack=True,usecols=(6,))
    M1[i] = int(Nsurf[plot_int*(i+1)-1])
    M2[i] = int(nx*ny*dV*sum(rhoY1z)/mCO)

np.savetxt("res.mass_cons",np.transpose([M1,M2,M1+M2]))
