reset
clear
set terminal png
set output "OscilacionAmortiguada.png"
set xlabel "Tiempo t" font "cambria,14"
set ylabel "Posición de la masa 1" font "cambria,14"
set yrange [8:12]
set grid
show grid
unset key


plot 'OscilacionAmortiguada.dat' u 1:2 w l lw 1 lc rgb "red"