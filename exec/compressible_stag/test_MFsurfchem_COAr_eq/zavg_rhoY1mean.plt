set terminal png enhanced size 600, 480
set output "zavg_rhoY1mean.png"

val=2.866559e-04

ncellx = 16
ncelly = 16
ncellz = 16
dz = 9.36e-6
ncell = ncellx*ncelly*ncellz

set bars 0.3
plot [-dz:dz*(ncellz+1)][0.99*val:1.01*val] "res.zavg_rhoY1mean" u 1:2:($3/sqrt(ncellx*ncelly)*3) w e, val
