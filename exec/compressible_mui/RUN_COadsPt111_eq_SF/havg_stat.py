import numpy as np

a0=2000
a1=2000
a2=200000

nz=8
nf=7

tmp1 = np.zeros((nz,nf+1))
tmp2 = np.zeros((nz,nf+1))
cnt = 0

for i in range(a0+a1,a2+a1,a1):
    datafile = "havg%09d" %i
    print datafile
    data = np.loadtxt(datafile)

    tmp1 = tmp1 + data
    tmp2 = tmp2 + np.power(data,2)
    cnt = cnt + 1

tmp1 = tmp1/cnt
tmp2 = tmp2/cnt
tmp2 = tmp2 - np.power(tmp1,2)
tmp2[tmp2<0] = 0
tmp2 = np.sqrt(tmp2/cnt)

tmp1 = tmp1.T
tmp2 = tmp2.T

tmp = []
tmp.append(tmp1[0])
for k in range(nf):
    tmp.append(tmp1[k+1])
    tmp.append(tmp2[k+1])

tmp = np.array(tmp)
tmp = tmp.T
np.savetxt("havg_stat",tmp)

