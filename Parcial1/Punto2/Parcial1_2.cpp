#include <iostream>
#include <cmath>
#include "vector.h"
using namespace std;

//--Declaración de las variables globales del sistema--//
const double k=1;
const int N=3;
const double Lx=40;
const double x10=10, x20=20, x30=30; //Posiciones de equilibrio//
const double gam=0.05; //Constante de Fricción//
const double amp=0.1; //Amplitud del forzamiento //

//--Declaración de las constantes del método PEFRL --//
const double epsilon=0.1786178958448091e00;
const double lambda=-0.2123418310626054e00;
const double chi=-0.6626458266981849e-1;
const double lambda2=(1.0-2.0*lambda)/2.0;
const double chiepsilon=1.0-2.0*(chi+epsilon);

//--Decalaración de las clases---//
class Cuerpo;
class Colisionador;

//--Interfase e implementación de las clases ----//
//--Clase Cuerpo donde se establecen las magnitudes asociadas a las tres masas--//
class Cuerpo{
private:
  vector3D r, V, F; double m,R,id;
public:
  void Inicie(double x0, double y0, double Vx0, double Vy0, double m0, double R0);
  void BorreFuerza(){F.load(0,0,0);};
  void AdicioneFuerza(vector3D F0){F+=F0;};
  void Mueva_r(double dt, double Coeficiente);
  void Mueva_V(double dt, double Coeficiente);
  void Dibujese(void);
  double Getx(void){return r.x();}; //inline
  double Gety(void){return r.y();}; //inline
  double GetVx(void){return V.x();};
  double Getm(void){return m;};
  friend class Colisionador;
};
void Cuerpo::Inicie(double x0,double y0,double Vx0,double Vy0,double m0,double R0){
  r.load(x0,y0,0); V.load(Vx0,Vy0,0);  m=m0;  R=R0;
} 
void Cuerpo::Mueva_r(double dt, double Coeficiente){
  r+=V*(Coeficiente*dt);  
}
void Cuerpo::Mueva_V(double dt, double Coeficiente){
  V+=F*(Coeficiente*dt/m);
}
void Cuerpo::Dibujese(void){
  cout<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t)";
}
//--Clase colisionador donde se calculan las interacciones entre las masas--//
class Colisionador{
private:
public:
  void CalculeFuerzas(Cuerpo * Grano, double t, double w);
};

void Colisionador::CalculeFuerzas(Cuerpo * Grano, double t, double w){
  int i;
  double m=Grano[0].Getm();
  double x1=Grano[0].Getx()-x10, x2=Grano[1].Getx()-x20, x3=Grano[2].Getx()-x30;
  double x4=amp*sin(w*t); //Se agrega el movimiento de forzamiento//
  double v1=Grano[0].GetVx(),v2=Grano[1].GetVx(),v3=Grano[2].GetVx();
  vector3D F1,F2,F3;
  double f1 = -k*x1+k*(x2-x1)-gam*m*v1;
  double f2 = -k*(x2-x1)+k*(x3-x2)-gam*m*v2;
  double f3 = -k*(x3-x2)+k*(x4-x3)-gam*m*v3;

  for(i=0;i<N;i++){Grano[i].BorreFuerza();}
  
  F1.load(f1,0,0);
  Grano[0].AdicioneFuerza(F1);
  F2.load(f2,0,0);
  Grano[1].AdicioneFuerza(F2);
  F3.load(f3,0,0);
  Grano[2].AdicioneFuerza(F3);
}


//---Funciones de animación y dibujo---//
void InicieAnimacion(void){
  cout<<"set terminal gif animate"<<endl; 
  cout<<"set output 'OscilaciónAmortiguada.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange[-5:"<<Lx+5<<"]"<<endl;
  cout<<"set yrange[-5:5]"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set trange [0:7]"<<endl;
  cout<<"set isosamples 12"<<endl;  
}
void InicieCuadro(void){
    cout<<"plot 0,0 ";
    cout<<" , 0,"<<-10/7<<"*t";        //pared de la izquierda//
    cout<<" , "<<Lx<<","<<-10/7<<"*t"; //pared de la derecha//
    cout<<" , 0,"<<10/7<<"*t";        //pared de la izquierda//
    cout<<" , "<<Lx<<","<<10/7<<"*t"; //pared de la derecha//
}
void TermineCuadro(void){
    cout<<endl;
}

//-------Programa principal--------//
int main(void){
  int i;
  Cuerpo Grano[N];
  Colisionador Hooke;
  double m=1, R=1;
  double t, tdibujo, dt=0.01;
  double tmax=300;
  double tcuadro=tmax/100;
  double tequilibrio=200;
  double omega=0.5;
  
  //InicieAnimacion();
  
  //--Para inicializar las moléculas--//

  Grano[0].Inicie(x10+2,0,0,0,m,R); //Las posiciones iniciales se miden desde el punto de equilibrio//
  Grano[1].Inicie(x20-2*sqrt(2),0,0,0,m,R);
  Grano[2].Inicie(x30+2,0,0,0,m,R);

  //--Se inicia el ciclo principal--//
  for(t=0, tdibujo=0 ; t<tmax ; t+=dt, tdibujo +=dt){

    //---Dibujo para cada cuadro de la animación--//
    if(tdibujo>tcuadro){
      
      //InicieCuadro();
      //for(i=0;i<N;i++) Grano[i].Dibujese();
      //TermineCuadro();
      
      tdibujo=0;
    }    
    //---Muévase por PEFRL----//
    for(i=0;i<N;i++)Grano[i].Mueva_r(dt,epsilon);
    Hooke.CalculeFuerzas(Grano,t,omega);
    for(i=0;i<N;i++)Grano[i].Mueva_V(dt,lambda2);
    for(i=0;i<N;i++)Grano[i].Mueva_r(dt,chi);
    Hooke.CalculeFuerzas(Grano,t,omega);
    for(i=0;i<N;i++)Grano[i].Mueva_V(dt,lambda);
    for(i=0;i<N;i++)Grano[i].Mueva_r(dt,chiepsilon);
    Hooke.CalculeFuerzas(Grano,t,omega);
    for(i=0;i<N;i++)Grano[i].Mueva_V(dt,lambda);
    for(i=0;i<N;i++)Grano[i].Mueva_r(dt,chi);
    Hooke.CalculeFuerzas(Grano,t,omega);
    for(i=0;i<N;i++)Grano[i].Mueva_V(dt,lambda2);
    for(i=0;i<N;i++)Grano[i].Mueva_r(dt,epsilon);

    if(t>tequilibrio){
      cout<<t<<" "<<Grano[0].Getx()<<endl;
    }
  }
  return 0;
}
