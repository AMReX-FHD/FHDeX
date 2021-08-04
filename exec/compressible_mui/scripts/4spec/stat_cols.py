import numpy as np

Nsample = 64
filefmt = "RUN%d/res.havg_spec4"
outfile = "res.havg_spec4_stat"
ntavg = 10

# check ncol and npt
infile = filefmt % 1
data = np.transpose(np.loadtxt(infile))
ncol = len(data)
npt = len(data[0])

print "ncol = %d" % ncol
print "npt  = %d" % npt

data2 = []
for j in range(ncol):
    data2.append([])

for i in range(Nsample):
    infile = filefmt % (i+1)
    print "reading %s" % infile
    data = np.transpose(np.loadtxt(infile))

    for j in range(ncol):

        # time average
        vec = np.zeros(npt/ntavg)
        for k in range(npt/ntavg):
            vec[k] = np.mean(data[j][k*ntavg:(k+1)*ntavg])
        data2[j].append(vec)

data3 = []
for j in range(ncol):
    data3.append(np.mean(data2[j],axis=0))
    data3.append(np.std(data2[j],axis=0,ddof=1))

np.savetxt(outfile,np.transpose(data3))
