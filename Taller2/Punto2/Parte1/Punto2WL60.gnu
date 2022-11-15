reset 
clear 
set xrange [100:180]
set yrange [150:200]
set terminal jpeg enhanced 
set output "LongitudDeOnda60.jpg"

R(x)= a*(x+b)
a=-1.4
b=-288

T(x)= a*(x+c)
c=-282

V(x)= a*(x+d)
d=-276

#fit R(x) "LongitudDeOnda60.dat" via a,b
plot "LongitudDeOnda60.dat", R(x) w l, T(x) w l,V(x) w l
