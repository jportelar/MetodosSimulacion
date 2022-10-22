reset
clear
set term pdf
set output "PresionTemperaturaprueba0.pdf"
set xrange [0:20]
set xlabel "Energía térmica kT" font "cambria,14"
set yrange [200:800]
set ylabel "Presión" font "cambria,14"
set grid
show grid
unset key

R(x)=a*x

fit R(x) 'punto5i.dat' u 3:2 via a

plot 'punto5i.dat' u 3:2 with point pointtype 7 ps 0.5 lc rgb "red" lw 1, R(x) w l lw 2 lc "blue" 