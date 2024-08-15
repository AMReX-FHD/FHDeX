set terminal png enhanced size 1200, 960
set output "zavg.png"

#val1=6.068741e-11  # for rho (total mass density)
val1=1.625906e-11   # for rho1
val2=2.285787e+10   # for rhoE
val3=2.550839e+01   # for temp
val4=6.756854e-02   # for jx, jy, jz

ncellx = 16
ncelly = 16
ncellz = 16
dz = 8.4e-6
ncell = ncellx*ncelly*ncellz

set multiplot
set size 0.5,0.5
set origin 0,0.5
set bars 0.3
plot [-dz:dz*(ncellz+1)][0.95*val1:1.05*val1] "res.zavg" u 1:2:($3/sqrt(ncellx*ncelly)*3) w e, val1, val1*(1-1./ncell)
set origin 0.5,0.5
plot [-dz:dz*(ncellz+1)][0.95*val2:1.05*val2] "res.zavg" u 1:4:($5/sqrt(ncellx*ncelly)*3) w e, val2, val2*(1-1./ncell)
set origin 0,0
plot [-dz:dz*(ncellz+1)][0.95*val3:1.05*val3] "res.zavg" u 1:6:($7/sqrt(ncellx*ncelly)*3) w e, val3, val3*(1+1./ncell)
set origin 0.5,0
plot [-dz:dz*(ncellz+1)][0.95*val4:1.05*val4] "res.zavg" u 1:8:($9/sqrt(ncellx*ncelly)*3) w e,'' u 1:10:($11/sqrt(ncellx*ncelly)*3) w e,'' u ($1-dz/2):12:($13/sqrt(ncellx*ncelly)*3) w e, val4, val4*(1-1./ncell)
