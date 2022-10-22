#include <iostream>
#include <cmath>
#include "Random64.h"
using namespace std;

const int Lx=1024;
const double p=0.5;
const int Q=2;//Número de estados

//-----Clase LatticeGas----//
class LatticeGas{
private:
  int V[Q]; //V[i] i=0 (derecha) i=1 (izquierda)
  int n[Lx][Q], nnew[Lx][Q]; //n[ix][i]
public:
  LatticeGas(void);
  void Borrese(void);
  void Inicie(int N, double mu, double sigma, Crandom & ran64);
  void Show(void);
  void Colisione(Crandom & ran64);
  void Adveccione(void);
  double rho(int ix);
  double Varianza(void);
};
LatticeGas::LatticeGas(void){
  V[0]=1; V[1]=-1;
}
void LatticeGas::Borrese(void){
  for(int ix=0;ix<Lx;ix++){
    for(int i=0;i<Q;i++){
      n[ix][i]=nnew[ix][i]=0;
    }
  }
}
void LatticeGas::Inicie(int N, double mu, double sigma, Crandom & ran64){
  int ix,i;
  while(N>0){//Mientras (N>0)
    ix=(int) ran64.gauss(mu,sigma);//se convierte en entereo (int)// Escoger una celda al azar
    if(ix<0) ix=0; if(ix>Lx-1) ix=Lx-1; //Corregir en los bordes si es necesario
    i=(int) Q*ran64.r();// Escoger una dirección al azar
    if(n[ix][i]==0){// si (ese sitio está vacio)
      n[ix][i]=1;N--;//   ponga una bolita ahí y decrezca
    }
  }
}
void LatticeGas::Show(void){
  for(int i=0;i<Q;i++){
    for(int ix=0;ix<Lx;ix++){
      cout<<n[ix][i];
    }
    cout<<endl;
  }
  cout<<endl;
}
void LatticeGas::Colisione(Crandom & ran64){
  for(int ix=0;ix<Lx;ix++){ //para cada celda
    if(ran64.r()>p){//Genero un número al azar, y si es menor que p
      nnew[ix][0]=n[ix][1]; nnew[ix][1]=n[ix][0]; //intercambio el contenido
    }
    else{
      nnew[ix][0]=n[ix][0]; nnew[ix][1]=n[ix][1];
    }
  }
}
void LatticeGas::Adveccione(void){
  for(int ix=0;ix<Lx;ix++){
    for(int i=0;i<Q;i++){
      n[(ix+V[i]+Lx)%Lx][i]=nnew[ix][i];
    }
  }
}
double LatticeGas::rho(int ix){
  return n[ix][0]+n[ix][1];
}
double LatticeGas::Varianza(void){
  int ix; double N, Xprom, sigma2;
  for(N=0, ix=0;ix<Lx;ix++){
    N+=rho(ix);
  }
  for(Xprom=0,ix=0;ix<Lx;ix++){
    Xprom+=ix*rho(ix);
  }
  Xprom/=N;
  for(sigma2=0,ix=0;ix<Lx;ix++){
    sigma2+=pow(ix-Xprom,2.0)*rho(ix);
  }
  sigma2/=(N-1);

  return sigma2;
}

//------Programa Principal------//
int main(void){
  LatticeGas Difusion;
  Crandom ran64(1);
  int N=400; double mu=Lx/2, sigma=Lx/8;
  int t, tmax=400;
  
  Difusion.Borrese();
  Difusion.Inicie(N,mu,sigma,ran64);
  for(t=0;t<tmax;t++){
    cout<<t<<" "<<Difusion.Varianza()<<endl;
    Difusion.Colisione(ran64);
    Difusion.Adveccione();
  }
  return 0;
}
