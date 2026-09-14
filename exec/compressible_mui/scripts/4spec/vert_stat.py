import numpy as np

NSAMPLE=64
NSPEC=4

for spec in range(NSPEC):

    list = []

    for i in range(NSAMPLE):
        infile = "RUN%d/res.mass_spec%d_final_vert" % (i+1,spec+1)
        a = np.loadtxt(infile)
        list.append(a)

    list = np.transpose(np.array(list))

    mean = np.mean(list,axis=1)
    std = np.std(list,axis=1,ddof=1)

    res_stat = np.transpose([mean,std])

    outfile = "res.mass_spec%d_final_vert_stat" % (spec+1)
    np.savetxt(outfile,res_stat)
