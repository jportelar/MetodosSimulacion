reset 
clear 
set xrange [100:180]
set yrange [75:200]
set terminal jpeg enhanced 
set output "LongitudDeOnda30.jpg"

R(x)= a*(x+b)
a=-3.5
b=-195.25

T(x)= a*(x+c)
c=-190.25

V(x)= a*(x+d)
d=-185.25

S(x)= a*(x+e)
e=-180.25

U(x)= a*(x+f)
f=-175.25


#fit R(x) "LongitudDeOnda30.dat" via a,b
plot "LongitudDeOnda30.dat", R(x) w l, T(x) w l,V(x) w l,S(x) w l, U(x) w l
