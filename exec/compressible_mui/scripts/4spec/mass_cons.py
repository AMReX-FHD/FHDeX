import numpy as np

infile0 = "res.coverage"
infile1 = "res.havg_spec1"
infile2 = "res.havg_spec2"
infile3 = "res.havg_spec3"
infile4 = "res.havg_spec4"

outfile1 = "res.mass_spec1"
outfile2 = "res.mass_spec2"
outfile3 = "res.mass_spec3"
outfile4 = "res.mass_spec4"

navo = 6.02e23
m1 = 4.00/navo
m2 = 20.18/navo
m3 = 39.95/navo
m4 = 83.80/navo

stepincr = 10
stepmax  = 10000

dx = 8.e-6
dV = dx**3
nc = 8*8

# surface

data0 = np.loadtxt(infile0)
data0 = np.transpose(data0)

N1s = data0[1][stepincr-1:stepmax:stepincr]
N2s = data0[2][stepincr-1:stepmax:stepincr]
N3s = data0[3][stepincr-1:stepmax:stepincr]
N4s = data0[4][stepincr-1:stepmax:stepincr]

m1s = N1s*m1
m2s = N2s*m2
m3s = N3s*m3
m4s = N4s*m4

# gas

data1 = np.loadtxt(infile1)
m1g = data1*dV*nc

data2 = np.loadtxt(infile2)
m2g = data2*dV*nc

data3 = np.loadtxt(infile3)
m3g = data3*dV*nc

data4 = np.loadtxt(infile4)
m4g = data4*dV*nc

# output

tmp = np.column_stack((m1s,m1g))
np.savetxt(outfile1+"_final_vert",np.transpose(tmp[-1]))
m1sum = np.sum(tmp,axis=1)
tmp = np.column_stack((tmp,m1sum))
np.savetxt(outfile1,tmp)

tmp = np.column_stack((m2s,m2g))
np.savetxt(outfile2+"_final_vert",np.transpose(tmp[-1]))
m2sum = np.sum(tmp,axis=1)
tmp = np.column_stack((tmp,m2sum))
np.savetxt(outfile2,tmp)

tmp = np.column_stack((m3s,m3g))
np.savetxt(outfile3+"_final_vert",np.transpose(tmp[-1]))
m3sum = np.sum(tmp,axis=1)
tmp = np.column_stack((tmp,m3sum))
np.savetxt(outfile3,tmp)

tmp = np.column_stack((m4s,m4g))
np.savetxt(outfile4+"_final_vert",np.transpose(tmp[-1]))
m4sum = np.sum(tmp,axis=1)
tmp = np.column_stack((tmp,m4sum))
np.savetxt(outfile4,tmp)
