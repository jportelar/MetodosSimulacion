reset
clear
set terminal png
set output "Oscilaciones3.png"
set xlabel "Tiempo t" font "cambria,14"
set ylabel "Posición de la masa 1" font "cambria,14"
set yrange [7.5:12.5]
set grid
show grid
unset key

G(x) = a*cos(b*x+c)+d
a=2
b=0.75
c=2
d=10

fit G(x) 'OscilacionLibre3.dat' u 1:2 via a,b,c,d

plot 'OscilacionLibre3.dat' u 1:2 with point pointtype 3 ps 1 lc rgb "green" lw 1, G(x) lt 12 lw 2 lc "blue"
