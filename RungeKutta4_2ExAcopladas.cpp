#include <iostream>
#include <cmath>
using namespace std;

const double omega=3.0;

double f1(double t, double x1, double x2){
  return -omega*omega*x2;
}

double f2(double t, double x1, double x2){
  return x1;
}

void UnPasoDeRungeKutta4(double & t0, double & x10, double & x20, double dt){
  double dx11,dx21,dx31,dx41;
  double dx12,dx22,dx32,dx42;

  dx11=dt*f1(t0,x10,x20);  dx12=dt*f2(t0,x10,x20);
  dx21=dt*f1(t0+dt/2,x10+dx11/2,x20+dx12/2);  dx22=dt*f2(t0+dt/2,x10+dx11/2,x20+dx12/2);
  dx31=dt*f1(t0+dt/2,x10+dx21/2,x20+dx22/2);  dx32=dt*f2(t0+dt/2,x10+dx21/2,x20+dx22/2);
  dx41=dt*f1(t0+dt,x10+dx31,x20+dx32);  dx42=dt*f2(t0+dt,x10+dx31,x20+dx32);
  t0+=dt; x10+=(dx11+2*(dx21+dx31)+dx41)/6;  x20+=(dx12+2*(dx22+dx32)+dx42)/6; 
}

int main(){
  double t,x1,x2; double dt=0.1;

  for (t=0,x1=1, x2=0; t<7;){
    cout<<t<<" "<<x1<<" "<<x2<<endl; 
    UnPasoDeRungeKutta4(t,x1,x2,dt);
  }

  return 0;
}
