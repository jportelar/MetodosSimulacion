#include <iostream>
#include <cmath>
using namespace std;

double f(double t, double x){
  return x;
}

void UnPasoDeEuler(double & t, double & x, double dt){
  double dx;
  dx=dt*f(t,x);
  t+=dt;
  x+=dx;
}

int main(){
  double t,x; double dt=0.1;

  for (t=0;x=1;t<2+dt/2;){
    cout << t <<"\t" << x << "\t" << t*t/2 <<endl;
    UnPasoDeEuler(t,x,dt);
  }

  return 0;
}
