
#include <iostream>
#include <cmath>
#include "vector.h"
#include "Random64.h"
using namespace std;

//--------- DECLARAR CONSTANTES GLOBALES ------------//
const double Lx=60, Ly=120;
const int Nx=5, Ny=5, N=Nx*Ny;
const double re=10, epsi=1.0; //valores dentro de la fuerza de Lennard-Jones//

//------------------ CONSTANTES DEL PERFL ------------//
const double epsilon=0.1786178958448091e00;
const double lambda=-0.2123418310626054e00;
const double chi=-0.6626458266981849e-1;
const double lambda2=(1.0-2.0*lambda)/2.0;
const double chiepsilon=1.0-2.0*(chi+epsilon);

//----------------- DECLARAR CLASES ----------------//
class Cuerpo;
class Colisionador;

//---- INTERFASE E IMPLEMENTACION DE LAS CLASES ----//
//---------------- CLASE CUERPO --------------------//
class Cuerpo{
private:
  vector3D r,V,F; double m,R,intensidad; //Spq
public:
  void Inicie(double x0,double y0,double Vx0,double Vy0,double m0,double R0);
  void BorreFuerza(){F.load(0,0,0);};
  void AdicioneFuerza(vector3D F0){F+=F0;};
  void BorreIntensidad(double I0){intensidad=I0;}; //Valor de la intensidad de las colisiones de una particula contra una pared
  void AdicioneIntensidad(double P){intensidad+=P;}; //
  void Mueva_r(double dt, double Coeficiente);
  void Mueva_V(double dt, double Coeficiente);
  void Dibujese(void);
  double Getx(void){return r.x();}; //inline
  double Gety(void){return r.y();}; //inline
  double GetVx(void){return V.x();}; //inline
  double GetIntensidad(void){return intensidad;};
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
  cout<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t)"; //Solo dibujar las paredes y quitar los granos que las representan 
}

//--- clase Colisionador ----
class Colisionador{
private:
public:
  void CalculeFuerzas(Cuerpo * Grano);
  void CalculeFuerzaPared(Cuerpo & Grano);
  void CalculeFuerzaEntre(Cuerpo & Grano1, Cuerpo & Grano2);
};

void Colisionador::CalculeFuerzaEntre(Cuerpo & Grano1, Cuerpo & Grano2){
  vector3D r21, n, F; double d;
  r21 = Grano1.r-Grano2.r; d=r21.norm(); n= r21/d;
  F = n*(12*epsi/d)*(pow((re/d),12)-pow((re/d),6));
  Grano1.AdicioneFuerza(F); Grano2.AdicioneFuerza(F*(-1));
}

void Colisionador::CalculeFuerzaPared(Cuerpo & Grano){
  vector3D Fp; double R=Grano.R,x=Grano.Getx(),y=Grano.Gety();
  int i; double K=1.0e4; 
  Fp.load(0,0,0);
  double dP[4]={R-x, R+x-Lx, R-y, R+y-Ly};
  for(i=0;i<4;i++){
    if(dP[i]>0){
      if(i<2){
	Fp.load(pow(-1,i)*K*pow(dP[i],1.5),0,0);
	Grano.AdicioneFuerza(Fp);
	Grano.AdicioneIntensidad(pow(-1,i)*K*pow(dP[i],1.5)*1e-3); //Spq
      }
      else{
	Fp.load(0,pow(-1,i)*K*pow(dP[i],1.5),0);
	Grano.AdicioneFuerza(Fp);
	Grano.AdicioneIntensidad(pow(-1,i)*K*pow(dP[i],1.5)*1e-3); //Spq
      }
    }
  }
}

void Colisionador::CalculeFuerzas(Cuerpo * Grano){
  int i,j; vector3D Fg;
  //--- Borrar todas las fuerzas ---
  for(i=0;i<N;i++){
    Grano[i].BorreFuerza();
  }
  //--- Calcular Fuerzas entre pares de granos ---
  for(i=0;i<N;i++){
    CalculeFuerzaPared(Grano[i]);
    for(j=i+1;j<N;j++){
      CalculeFuerzaEntre(Grano[i], Grano[j]);
    }
  }
}

//----------------- Funciones de Animacion y Dibujo ----------
void InicieAnimacion(void){
  cout<<"set terminal gif animate"<<endl; 
  cout<<"set output 'GasLennardJones2DLiquido.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange[-10:"<<Lx+10<<"]"<<endl;
  cout<<"set yrange[-10:"<<Ly+10<<"]"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set trange [0:7]"<<endl;
  cout<<"set isosamples 12"<<endl;  
}
void InicieCuadro(void){
    cout<<"plot 0,0 ";
    cout<<" , "<<Lx/7<<"*t,0";        //pared de abajo
    cout<<" , "<<Lx/7<<"*t,"<<Ly;     //pared de arriba
    cout<<" , 0,"<<Ly/7<<"*t";        //pared de la izquierda
    cout<<" , "<<Lx<<","<<Ly/7<<"*t"; //pared de la derecha
}
void TermineCuadro(void){
    cout<<endl;
}

void GuardeFrame(void){
  cout<<"reset"<<endl;
  cout<<"set term png"<<endl;
  cout<<"set output 'EstructuraLiquida.png'"<<endl;
  cout<<"set xrange[-10:"<<Lx+10<<"]"<<endl;
  cout<<"set yrange[-10:"<<Ly+10<<"]"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set trange [0:7]"<<endl;
  cout<<"set isosamples 12"<<endl;
}

//-----------  Programa Principal --------------  
int main(void){
  Cuerpo Grano[N];
  Colisionador Lennard;
  Crandom ran64(1);
  double m0=1, R0=2.5, kT=10.0, V0=sqrt(2*kT/m0);
  int i,ix,iy;
  double t,tdibujo,dt=1e-3;
  double tmax=200; //double tmax=10*(Lx/V0);//
  double tcuadro=tmax/1000;
  double dx=Lx/(Nx+1), dy=10.0;
  double Theta, Vx;
  double sumacuadrada=0;
  double suma=0, count=0;
  double sumay=0, promy;
  double tequilibrio=50.0;
  double destandar;
  double intot=0, presion;

  //InicieAnimacion(); //Dibujar
  
  //-------------------Inicializar las moléculas-------------------//
  for(ix=0;ix<Nx;ix++)
    for(iy=0;iy<Ny;iy++){
      Theta=2*M_PI*ran64.r();
      //--------------------(   x0,   y0,          Vx0,                  Vy0, m0,R0)
      Grano[Nx*iy+ix].Inicie((ix+1)*dx,(iy+1)*dy,V0*cos(Theta),V0*sin(Theta), m0,R0);//OJO
    }
  
  for(i=0;i<N;i++){ Grano[i].BorreIntensidad(0);}
  
  //--------------------Inicia el ciclo principal---------------//
  for(t=0,tdibujo=0 ; t<tmax ; t+=dt,tdibujo+=dt){

    //------------------------Dibujar--------------------------//
    if(tdibujo>tcuadro){
      
      //InicieCuadro();
      //for(i=0;i<N;i++) Grano[i].Dibujese();
      //TermineCuadro();
      
      tdibujo=0;
    }    

    //---------------------- Muevase por PEFRL -----------------//
    for(i=0;i<N;i++)Grano[i].Mueva_r(dt,epsilon);
    Lennard.CalculeFuerzas(Grano);
    for(i=0;i<N;i++)Grano[i].Mueva_V(dt,lambda2);
    for(i=0;i<N;i++)Grano[i].Mueva_r(dt,chi);
    Lennard.CalculeFuerzas(Grano);
    for(i=0;i<N;i++)Grano[i].Mueva_V(dt,lambda);
    for(i=0;i<N;i++)Grano[i].Mueva_r(dt,chiepsilon);
    Lennard.CalculeFuerzas(Grano);
    for(i=0;i<N;i++)Grano[i].Mueva_V(dt,lambda);
    for(i=0;i<N;i++)Grano[i].Mueva_r(dt,chi);
    Lennard.CalculeFuerzas(Grano);
    for(i=0;i<N;i++)Grano[i].Mueva_V(dt,lambda2);
    for(i=0;i<N;i++)Grano[i].Mueva_r(dt,epsilon);  

    //-------Cálculo de la altura prodmedio---//
    sumay=0;
    for(i=0;i<N;i++){
      sumay+=Grano[i].Gety();
    };
    promy=sumay/N;
    //cout<<t<<" "<<promy<<endl;
    
    //---Desviación estandar de la velocidad y datos para el histograma----//
    if(t>tequilibrio){
      for(i=0;i<N;i++){
	Vx=Grano[i].GetVx();
	//cout<<Vx<<endl;
	suma+=Vx;
	sumacuadrada+=Vx*Vx;
	count+=1;
	//----------Cálculo de la presión global---------//
	intot+=Grano[i].GetIntensidad();
	presion=intot/(tmax*Lx);
      }
    }
  }
  destandar=sqrt((1/count)*(sumacuadrada-(1/count)*suma*suma));
  //cout<<"la desviacion estandar de la velocidad es "<<destandar<<endl;
  //cout<<"la temperatura es "<<pow(destandar,2)*m0<<endl;

  //------------------Guardar el último cuadro----------------//
  //GuardeFrame();
  //InicieCuadro();
  //for(i=0;i<N;i++) Grano[i].Dibujese();
  //TermineCuadro();
  //------------------Valores de temperatura y presión----------//
  //cout<<"la presion para kT="<<kT<<" es "<<presion<<endl;
  return 0;
}
