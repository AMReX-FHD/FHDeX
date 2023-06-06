set terminal png enhanced size 1200, 960
set output "zavg.png"

val1=3.093067e-11
val2=2.384732e+09
val3=1.061837e+02
val4=9.181427e-02

dz = 8.4e-6
ncell = 8

set multiplot
set size 0.5,0.5
set origin 0,0.5
plot [0:dz*ncell][0.80*val1:1.20*val1] "res.zavg" u 1:2:($3/ncell*3) w e, val1, val1*(1-1./ncell**3)
set origin 0.5,0.5
plot [0:dz*ncell][0.80*val2:1.20*val2] "res.zavg" u 1:4:($5/ncell*3) w e, val2, val2*(1-1./ncell**3)
set origin 0,0
plot [0:dz*ncell][0.95*val3:1.05*val3] "res.zavg" u 1:6:($7/ncell*3) w e, val3, val3*(1+1./ncell**3)
set origin 0.5,0
plot [0:dz*ncell][0.95*val4:1.05*val4] "res.zavg" u 1:8:($9/ncell*3) w e,'' u 1:10:($11/8*3) w e,'' u ($1-dz/2):12:($13/8*3) w e, val4, val4*(1-1./ncell**3)
