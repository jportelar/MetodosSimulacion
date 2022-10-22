reset
clear
set output 'AlturaPromedio.pdf'
set term pdf
unset key
set xrange [0:200]
set yrange [25:70]
set xlabel 'Tiempo' font "cambria,14"
set ylabel 'Altura Promedio' font "cambria,14"
set grid
show grid

f(x)=55

plot 'punto5e.dat' w l lw 2 lc rgb "red",f(x) dashtype "--" lc rgb "blue" 