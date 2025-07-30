set terminal png enhanced size 1200, 960
set output "surfcoVar.png"

val1 = 0
val2 = 0
val3 = 0
val4 = 0

ncellx = 16
ncelly = 16

set multiplot
set size 0.5, 0.5
set bars 0.3

set origin 0,0.5
plot "res.surfcoVar_rhoYk-T" u 0:3:($5/sqrt(ncellx*ncelly)*3) w e, val1
set origin 0.5,0.5
plot "res.surfcoVar_theta-vx" u 0:3:($5/sqrt(ncellx*ncelly)*3) w e, "res.surfcoVar_theta-vy" u 0:3:($5/sqrt(ncellx*ncelly)*3) w e, "res.surfcoVar_theta-vz" u 0:3:($5/sqrt(ncellx*ncelly)*3) w e, val2
set origin 0,0
plot "res.surfcoVar_theta-rhoYk" u 0:3:($5/sqrt(ncellx*ncelly)*3) w e, val3
set origin 0.5,0
plot "res.surfcoVar_theta-T" u 0:3:($5/sqrt(ncellx*ncelly)*3) w e, val4
