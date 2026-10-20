import numpy as np
import sys

if len(sys.argv)!=3:
  print("Usage: python %s xhi yhi" % sys.argv[0])
  sys.exit()

xhi = int(sys.argv[1])
yhi = int(sys.argv[2])

a1 = 3.12e-08
a2 = 3.12e-08

nsites = xhi*yhi

file_name = "data.strips"
site_file = open(file_name,"w")


site_file.write("Site file written by create_site_file.py\n\n")
site_file.write("{} sites\n".format(nsites))
site_file.write("{} max neighbors\n".format(4))
site_file.write("id site i1 i2 values\n\n")

site_file.write("0 {:.6e} xlo xhi\n".format(a1*xhi))
site_file.write("0 {:.6e} ylo yhi\n".format(a2*yhi))
site_file.write("-{:.6e} {:.6e} zlo zhi\n".format(a1, a1))

site_file.write("\n")

site_file.write("Sites\n\n")

for j in range(0,xhi):
    for i in range(1,yhi+1):
        site_file.write("{} {:.6e} {:.6e} 0.0\n".format(yhi*j+i,a1*j,a2*(i-1)))

site_file.write("\n")
site_file.write("Neighbors\n\n")

for i in range(1,nsites+1):
    neigh1 = i-1
    if neigh1 % yhi == 0:
        neigh1 += yhi

    neigh2 = i+1
    if (neigh2-1) % yhi == 0:
        neigh2 -= yhi

    neigh3 = i+yhi
    if neigh3 > nsites:
        neigh3 -= nsites

    neigh4 = i-yhi
    if neigh4 <= 0:
        neigh4 += nsites

    site_file.write("{} {} {} {} {}\n".format(i,neigh1,neigh2,neigh3,neigh4))

site_file.write("\n")

"""
site_file.write("Values\n\n")

for i in range(1,nsites+1):
    site_file.write("{} 1 1 0\n".format(i))
"""
site_file.close()

print("%s generated" % file_name)
