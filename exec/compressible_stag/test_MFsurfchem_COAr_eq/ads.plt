set terminal png enhanced size 1200, 480

val1 = 7.466307e-01
val2 = 2.101926e-06

set term png
set output "ads.png"

set multiplot
set origin 0,0
set size 0.5,1
set key bottom
plot [:][0.99*val1:1.01*val1] "res.ads_mean" u 0:3:($5*3) w e, val1
set origin 0.5,0
set size 0.5,1
set key bottom
plot [:][0.95*val2:1.05*val2] "res.ads_var" u 0:3:($5*3) w e, val2
