reset
clear
set terminal png
set output "OscilacionForzada.png"
set xlabel "Tiempo t" font "cambria,14"
set ylabel "Posición de la masa 1" font "cambria,14"
set yrange [9.92:10.08]
set grid
show grid
unset key


plot 'OscilacionForzada.dat' u 1:2 w l lw 1 lc rgb "green"