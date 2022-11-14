reset 
clear 
set xrange [100:180]
set yrange [-9:9]
set terminal jpeg enhanced 
set output "LongitudDeOnda.jpg"

R(x)= a*sin(b*x+c)
a=8
b=6.366
c=-101

fit R(x) "LongitudDeOnda.dat" via a,b,c

plot "LongitudDeOnda.dat" w l, R(x) w l
