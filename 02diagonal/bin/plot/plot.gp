set term png fontscale 1.2
set output "lrmesh.png"
set xlabel "xi"
set ylabel "eta"
set title "LR mesh"
set size ratio -1
p "../LRmesh2.res"w l, "../LRmesh.res"w l

set output "surface.png"
set xlabel "xi"
set ylabel "eta"
set title "LR-spline surface"
set size ratio -1
sp "../surface2.res", "../surface.res"
