import numpy as np

stepmax  = 10000
stepincr = 10

# figure out sizes
infile = "havg"+str(stepincr).zfill(9)
data = np.loadtxt(infile)
data = np.transpose(data)

nspec = len(data)-6
zvec = data[0]
nzpt = len(zvec)
ntpt = stepmax/stepincr

# 0: z coordinate
# 1: rho
# 2: j_x
# 3: j_y
# 4: j_z
# 5: E
# 6...(5+nspec): rho_k

# init arrays
res = np.zeros((5+nspec,ntpt,nzpt))

# copy data
for i in range(ntpt):
    step = (i+1)*stepincr
    infile = "havg"+str(step).zfill(9)
    data = np.loadtxt(infile)
    data = np.transpose(data)

    for n in range(5+nspec):
        res[n][i] = data[n+1]

# file output
outfile = "res.havg_rho"
np.savetxt(outfile,res[0])

outfile = "res.havg_jx"
np.savetxt(outfile,res[1])

outfile = "res.havg_jy"
np.savetxt(outfile,res[2])

outfile = "res.havg_jz"
np.savetxt(outfile,res[3])

outfile = "res.havg_E"
np.savetxt(outfile,res[4])

for n in range(nspec):
    outfile = "res.havg_spec%d" % (n+1)
    np.savetxt(outfile,res[5+n])
