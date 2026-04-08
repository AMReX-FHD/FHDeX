import os
import numpy as np
from matplotlib import pyplot as plt
import yt

# This script checks if there is conservation along x-direction as
# y-direction is an ensemble direction , i.e., no flux exchange.

parent_folder = "../"

fld_path = parent_folder
# Get plt files
onlyplotfiles = [fl for fl in os.listdir(fld_path) if "plt" in fl]

# Ensure to remove files that contain "old" string
onlyplotfiles = [fl for fl in onlyplotfiles if "old" not in fl]
onlyplotfiles.sort()
n_files = len(onlyplotfiles)

# Load the first plotfiles
plotfile_nm = os.path.join(fld_path, onlyplotfiles[0])
ds = yt.load(plotfile_nm)

print(plotfile_nm)
print(float(ds.current_time))
print(ds.field_list)

all_data = ds.covering_grid(
    level=0, left_edge=ds.domain_left_edge, dims=ds.domain_dimensions)
phi_start = all_data['boxlib', 'phi0'].to_ndarray()[:,:]

# Load the last file
plotfile_nm = os.path.join(fld_path, onlyplotfiles[-2])
ds = yt.load(plotfile_nm)
print(plotfile_nm)
print(float(ds.current_time))
print(ds.field_list)

all_data = ds.covering_grid(
    level=0, left_edge=ds.domain_left_edge, dims=ds.domain_dimensions)

phi_end = all_data['boxlib', 'phi0'].to_ndarray()[:,:]

# Sum along x-dimension as y-dimension is ensemble direction
aa_start = np.sum(phi_start, axis=0)
aa_end   = np.sum(phi_start, axis=0)

print(aa_start[0], aa_end[0])
print(aa_start[-1], aa_end[-1])

print(np.max(np.abs(aa_start - aa_end)))
