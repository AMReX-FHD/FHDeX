set terminal png enhanced size 1200, 960
set output "zavg.png"

# 300x150
val1=6.068741e-11
val2=1.117179e+09
val2c=1.116669e+09
val3=1.691011e+02
val4=1.282197e-01
# 300x173
#val1 = 5.261914e-11
#val2 = 9.686527e+08
#val3 = 1.466195e+02
#val4 = 0.1111731

set multiplot
set size 0.5,0.5
set origin 0,0.5
plot [:][0.9*val1:1.1*val1] "res.zavg" u 1:2:($3/8*3) w e, val1, val1*511/512 
set origin 0.5,0.5
plot [:][0.9*val2:1.1*val2] "res.zavg" u 1:4:($5/8*3) w e, val2, val2c, val2*511/512
set origin 0,0
plot [:][0.95*val3:1.05*val3] "res.zavg" u 1:6:($7/8*3) w e, val3, val3*513/512
set origin 0.5,0
plot [:][0.975*val4:1.025*val4] "res.zavg" u 1:8:($9/8*3) w e,'' u 1:10:($11/8*3) w e,'' u 1:12:($13/8*3) w e, val4, val4*511/512
#    EOF
