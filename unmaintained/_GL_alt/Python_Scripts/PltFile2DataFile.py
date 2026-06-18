import sys, os, os.path
import glob

import scipy as sp
import numpy as np

import matplotlib
import matplotlib.pyplot as pp
from mpl_toolkits.mplot3d import Axes3D

import yt
from yt.frontends.boxlib.data_structures import AMReXDataset

from tempfile import TemporaryFile


def substitute_header(plt_file, source="phi", target="con"):

    # load header file
    header_file = os.path.join(plt_file, "Header")
    with open(header_file, "r") as f:
        header_orig = f.readlines()

    # select variable lables
    n_lables   = int(header_orig[1])
    l_offset   = 2

    # make a backup copy(iff the source was found in original)
    if source+"\n" in header_orig:
        header_cpy  = os.path.join(plt_file, "Header.backup")
        with open(header_cpy, "w") as f:
            for line in header_orig:
                f.write(line)

    # replace source with target
    for i in range(l_offset, n_lables+l_offset):
        if header_orig[i] == source+"\n":
            header_orig[i] = target+"\n"

    # save substituted file in place of original
    with open(header_file, "w") as f:
        for line in header_orig:
            f.write(line)


data_root = "."
data_dir  = "."



data_path = os.path.join(data_root, data_dir)
print(data_path)

n_fill   = 5
prefix   = "plt"
file_fmt = prefix + "{:0" + str(n_fill) + "d}"


data_glob  = os.path.join(data_path, prefix + "*")
data_files = glob.glob(data_glob)
data_files.sort()
n = 100 #sample number per window
N= 55 #number of umbrellas
phi_init0=0.05
phi_n_inc=0.01
for y in range (0,N):
    values = [0]*n
    for x in range (0,n):
        substitute_header(data_files[x+y*n])
        ds = yt.load(data_files[x+y*n])
        square = ds.covering_grid(level=0, fields=["con"],left_edge=ds.domain_left_edge, dims=ds.domain_dimensions)
        data = square["con"]
        #matshow(data[:,:,0])
        new_data = data[:,:,0]
        values[x] = np.sum(new_data)*(1./128.)*(1./128.)
        print("here",y,x,np.sum(new_data)*(1./128.)*(1./128.)) # note that np.sum sums all the data in the arrays if you don't specify an axis
    print(values)
    phi_0_loc =str(format(((y)*phi_increment+phi_initial), '.3f'))
    np.savetxt('data'+phi_0_loc, values)
