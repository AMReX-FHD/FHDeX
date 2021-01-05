import numpy as np

stepmax  = 5000
stepincr = 10

# figure out sizes
infile = "havg"+str(stepincr).zfill(9)
data = np.loadtxt(infile)
data = np.transpose(data)

zvec = data[0]
nzpt = len(zvec)
ntpt = stepmax/stepincr

# 0: z coordinate
# 1: j_x^2
# 2: j_y^2
# 3: j_z^2

# init arrays
res = np.zeros((3,ntpt,nzpt))

# copy data
for i in range(ntpt):
    step = (i+1)*stepincr
    infile = "havg"+str(step).zfill(9)
    data = np.loadtxt(infile)
    data = np.transpose(data)

    for n in range(3):
        res[n][i] = data[n+1]

# file output
outfile = "res.havg_jxsq"
np.savetxt(outfile,res[0])

outfile = "res.havg_jysq"
np.savetxt(outfile,res[1])

outfile = "res.havg_jzsq"
np.savetxt(outfile,res[2])
