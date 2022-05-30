val = 1.813017e-06

set term png
set output "ads_var.png"
set key bottom

plot "res.ads_var" u 0:3:($5*3) w e,val
