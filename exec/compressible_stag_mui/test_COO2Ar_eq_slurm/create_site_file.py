import numpy as np
import sys

if len(sys.argv)!=4:
  print("Usage: python %s xhi yhi scale" % sys.argv[0])
  sys.exit()

xhi = int(sys.argv[1])
yhi = int(sys.argv[2])
scale = str(sys.argv[3])
a1 = 6.0
a2 = 3.0

nsites = xhi*yhi

file_name = "data.strips"
site_file = open(file_name,"w")


site_file.write("Site file written by create_site_file.py\n\n")
# site_file.write("{} dimension\n".format(2))
site_file.write("{} sites\n".format(2*nsites))
site_file.write("{} max neighbors\n".format(2))
site_file.write("id site values\n\n")

site_file.write("0 3.60e-05 xlo xhi\n")
site_file.write("0 3.60e-05 ylo yhi\n")
site_file.write("-1.80e-08 1.80e-08 zlo zhi\n")

site_file.write("\n")

site_file.write("Sites\n\n")

for i in range(0,xhi):
    for j in range(0,yhi):
        site_file.write("{} {} {} 0.0\n".format(str(i*yhi+j+1),str(a1*i)+scale,str(a2*j)+scale))

for i in range(0,xhi):
    for j in range(0,yhi):
        site_file.write("{} {} {} 0.0\n".format(str(nsites+i*yhi+j+1),str(a1*(i+1/2))+scale,str(a2*(j+1/2))+scale))

m = xhi

site_file.write("\n")
site_file.write("Neighbors\n\n")

for i in range(0,xhi):
    site_file.write("{} {} {}\n".format(1+i*yhi,2+i*yhi,yhi*(i+1)))
    for j in range(2,yhi):
        site_file.write("{} {} {}\n".format(i*yhi+j,i*yhi+j-1,i*yhi+j+1))
    site_file.write("{} {} {}\n".format(yhi*(i+1),yhi*(i+1)-1,1+yhi*i))

for i in range(0,xhi):
    site_file.write("{} {} {}\n".format(nsites+1+i*yhi,nsites+2+i*yhi,nsites+yhi*(i+1)))
    for j in range(2,yhi):
        site_file.write("{} {} {}\n".format(nsites+i*yhi+j,nsites+i*yhi+j-1,nsites+i*yhi+j+1))
    site_file.write("{} {} {}\n".format(nsites+yhi*(i+1),nsites+yhi*(i+1)-1,nsites+1+yhi*i))

site_file.write("\n")
site_file.write("Values\n\n")

for i in range(1,nsites+1):
    site_file.write("{} 1\n".format(i))

for i in range(1,nsites+1):
    site_file.write("{} 2\n".format(nsites+i))

site_file.close()

print("%s generated" % file_name)
