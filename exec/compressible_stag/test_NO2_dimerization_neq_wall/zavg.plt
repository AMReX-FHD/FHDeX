set terminal png enhanced size 1200, 960
set output "zavg.png"

ncellx = 16
ncelly = 16
ncellz = 16
dz = 8e-6

set multiplot
set size 0.5,0.5
set origin 0,0.5
set bars 0.3
plot [-dz:dz*(ncellz+1)][320:380] "res.zavg" u 1:2:($3/sqrt(ncellx*ncelly)*3) w e t "tInstant"
set origin 0.5,0.5
plot [-dz:dz*(ncellz+1)][:] "res.zavg" u 1:4:($5/sqrt(ncellx*ncelly)*3) w e t "rhoInstant"
set origin 0,0
plot [-dz:dz*(ncellz+1)][:] "res.zavg" u 1:6:($7/sqrt(ncellx*ncelly)*3) w e t "rhoYkInstant\\\_0"
set origin 0.5,0
plot [-dz:dz*(ncellz+1)][:] "res.zavg" u 1:8:($9/sqrt(ncellx*ncelly)*3) w e t "rhoYkInstant\\\_1"
