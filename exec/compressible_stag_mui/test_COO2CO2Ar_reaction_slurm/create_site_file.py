import numpy as np
import sys

if len(sys.argv)!=3:
  print("Usage: python %s xhi yhi" % sys.argv[0])
  sys.exit()

xhi = int(sys.argv[1])
yhi = int(sys.argv[2])
a1 = 6.43
a2 = 3.12

nsites = xhi*yhi

file_name = "data.strips"
site_file = open(file_name,"w")


site_file.write("Site file written by create_site_file.py\n\n")
# site_file.write("{} dimension\n".format(2))
site_file.write("{} sites\n".format(2*nsites))
site_file.write("{} max neighbors\n".format(4))
site_file.write("id site i1 i2 values\n\n")

site_file.write("0 {:.3f}e-08 xlo xhi\n".format(a1*xhi))
site_file.write("0 {:.3f}e-08 ylo yhi\n".format(a2*yhi))
site_file.write("-{:.3f}e-08 {:.3f}e-08 zlo zhi\n".format(a1, a1))

site_file.write("\n")

site_file.write("Sites\n\n")

for i in range(0,yhi):
    for j in range(0,xhi):
        site_file.write("{} {:.3f}e-08 {:.3f}e-08 0.0\n".format(2*xhi*i+2*j+1,a1*j,i*a2))
        site_file.write("{} {:.3f}e-08 {:.3f}e-08 0.0\n".format(2*xhi*i+2*j+2,a1*(j+1/2),i*a2))

site_file.write("\n")
site_file.write("Neighbors\n\n")

for i in range(1,2*nsites+1):
    if i + 2*xhi > 2*nsites:
        a = i + 2*xhi - 2*nsites
    else:
        a = i + 2*xhi

    if i <= 2*xhi:
        b = i + 2*nsites - 2*xhi
    else:
        b = i - 2*xhi

    if (i -1) % (2*xhi) == 0:
        c = i + 2*xhi - 1
    else:
        c = i - 1

    if i % (2*xhi) == 0:
        d = i - 2*xhi + 1
    else:
        d = i + 1

    site_file.write("{} {} {} {} {}\n".format(i,a,b,c,d,))

site_file.write("\n")
site_file.write("Values\n\n")

for i in range(0,nsites):
    site_file.write("{} 1 1 0\n".format(2*i+1))
    site_file.write("{} 2 2 0\n".format(2*i+2))

site_file.close()

print("%s generated" % file_name)
