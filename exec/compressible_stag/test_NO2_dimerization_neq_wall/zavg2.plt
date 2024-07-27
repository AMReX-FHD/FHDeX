set terminal png enhanced size 1200, 960
set output "zavg.png"

ncellx = 16
ncelly = 16
ncellz = 16
dx = 8e-6
dy = 8e-6 
dz = 8e-6
dv = dx*dy*dz

M1 = 46.0055
M2 = 92.0110
cvhat1 = 3.64 
cvhat2 = 9.02 
Navo = 6.02214e23

rhoeq = 1.858338e-03
rho1eq = rhoeq*0.723969 
rho2eq = rhoeq*0.276031
Teq = 350

dTsq = Teq**2/Navo/(cvhat1/M1*rho1eq+cvhat2/M2*rho2eq)/dv

set multiplot
set size 0.5,0.5
set origin 0,0.5
set bars 0.3
plot [0:dz*ncellz][0.5*dTsq:1.5*dTsq] "res.zavg" u 1:10:($11/sqrt(ncellx*ncelly)*3) w e t "TVarDirect",'' u 1:($2**2/Navo/(cvhat1/M1*$6+cvhat2/M2*$8)/dv) w lp t "eq fluct", dTsq t "eq sys"
set origin 0.5,0.5
plot [0:dz*ncellz] "res.zavg" u 1:($10/($2**2/Navo/(cvhat1/M1*$6+cvhat2/M2*$8)/dv)) w lp t "ratio to eq fluct"
set origin 0,0
plot [0:dz*ncellz][320:380] "res.zavg" u 1:2:($3/sqrt(ncellx*ncelly)*3) w e t "tMean", Teq
set origin 0.5,0
plot [0:dz*ncellz][0:1.5*rhoeq] "res.zavg" u 1:6:($7/sqrt(ncellx*ncelly)*3) w e t "rhoYkMean\\\_0", rho1eq, '' u 1:8:($9/sqrt(ncellx*ncelly)*3) w e t "rhoYkMean\\\_1", rho2eq, '' u 1:4:($5/sqrt(ncellx*ncelly)*3) w e t "rhoMean\\\_1", rhoeq
