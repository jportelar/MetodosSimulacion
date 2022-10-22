reset
clear
set terminal png
set output "AmplitudFrecuencia.png"
set xlabel "Frecuencia" font "cambria,14"
set ylabel "Posición de la masa 1" font "cambria,14"
set yrange [9.2:10.8]
set grid
show grid
unset key


plot 'AmplitudFrecuencia.dat' u 2:3 w l lw 1 lc rgb "red"