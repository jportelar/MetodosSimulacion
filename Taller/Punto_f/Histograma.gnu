set xlabel "Velocidades"
set xrange [-11:11]

set ylabel "Frecuencia"

set grid
set title "Histograma de Velocidades"

set term pdf
set output "Histograma.pdf"
Min = -10.2348
Max = 10.0187
prom = 0.0099
desv = 2.9789
n = 50
width = (Max-Min)/n

hist(x,width) = width/2.0 + width*floor(x/width)

plot 'punto5f.dat' u (hist($1,width)):(1.0) smooth freq w boxes lw 2 lc rgb "red" notitle, 210000*exp(-((x-prom)**2)/(2*desv**2)) title "Gaussiana" lw 3 lc rgb "blue"


show grid
