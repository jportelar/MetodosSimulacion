reset 
clear 
set xrange [100:180]
set yrange [0:200]
set terminal jpeg enhanced 
set output "LongitudDeOnda20.jpg"

R(x)= a*(x+b)
a=-5.2
b=-174.5

T(x)= a*(x+c)
c=-170

V(x)= a*(x+d)
d=-165

#fit R(x) "LongitudDeOnda20.dat" via a,b
plot "LongitudDeOnda20.dat", R(x) w l, T(x) w l,V(x) w l
