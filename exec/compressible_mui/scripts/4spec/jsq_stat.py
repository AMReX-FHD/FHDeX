import numpy as np

NSAMPLE=4

for dim in ["x","y","z"]:

    list = []

    for i in range(NSAMPLE):
        infile = "m01_cutemp_16x16x8_5000_RUN%d/res.havg_j%ssq" % (i+1,dim)
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

    outfile = "res.havg_j%ssq_stat" % (dim)
    np.savetxt(outfile,res_stat)
