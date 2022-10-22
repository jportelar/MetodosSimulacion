clear
reset
set output "Oscilaciones.pdf"
set term pdf
set xrange [0:100]
set xlabel "Tiempo" font "cambria,14"
set ylabel "Posición en x" font "cambria,14"
set grid
unset key
show grid

plot 'punto5a.dat' w l lw 2 lc rgb "blue"

