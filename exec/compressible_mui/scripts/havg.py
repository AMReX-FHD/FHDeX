import numpy as np

stepmax  = 10000
stepincr = 10

# figure out sizes
infile = "havg"+str(stepincr).zfill(9)
data = np.loadtxt(infile)
data = np.transpose(data)

nspec = len(data)-1
zvec = data[0]
nzpt = len(zvec)
ntpt = stepmax/stepincr

# init arrays
res = np.zeros((nspec,ntpt,nzpt))

# copy data
for i in range(ntpt):
    step = (i+1)*stepincr
    infile = "havg"+str(step).zfill(9)
    data = np.loadtxt(infile)
    data = np.transpose(data)

    for n in range(nspec):
        res[n][i] = data[n+1]

# file output
for n in range(nspec):
    outfile = "res.havg_spec%d" % (n+1)
    np.savetxt(outfile,res[n])
