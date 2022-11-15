reset 
clear 
set xrange [100:180]
set yrange [120:200]
set terminal jpeg enhanced 
set output "LongitudDeOnda45.jpg"

R(x)= a*(x+b)
a=-2.1
b=-235.25

T(x)= a*(x+c)
c=-230

V(x)= a*(x+d)
d=-225

S(x)= a*(x+e)
e=-220

U(x)= a*(x+f)
f=-215


#fit R(x) "LongitudDeOnda45.dat" via a,b
plot "LongitudDeOnda45.dat", R(x) w l, T(x) w l,V(x) w l,S(x) w l, U(x) w l
