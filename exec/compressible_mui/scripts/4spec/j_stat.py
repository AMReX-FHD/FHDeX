import numpy as np

NSAMPLE=64

for dim in ["x","y","z"]:

    list = []

    for i in range(NSAMPLE):
        infile = "RUN%d/res.havg_j%s" % (i+1,dim)
        a = np.loadtxt(infile)
        list.append(np.transpose(a))

    nz = len(list[0])

    res_stat = []

    for n in range(nz):
        list2 = []
        for i in range(NSAMPLE):
            list2.append(list[i][n])

        mean = np.mean(list2,axis=0)
        std = np.std(list2,axis=0,ddof=1)

        res_stat.append(mean)
        res_stat.append(std)

    res_stat = np.transpose(res_stat)

    outfile = "res.havg_j%s_stat" % (dim)
    np.savetxt(outfile,res_stat)
