#include <cmath>
#include <iostream>
#include "../vector.h"
using namespace std;

//Constantes globales
const double g=9.8; //Gravedad
const double G=1; //Const gravitacional
const int N = 2; //Numero de planetas 

//Constantes PERFL
const double Zeta = 0.1786178958448091;
const double Lambda = -0.2123418310626054;
const double Chi = -0.6626458266981849e-1;

const double Lambda1 = 0.5 - Lambda;
const double Zetachi = 1 - 2*(Chi+Zeta);

//Declaracion de las clases
class Cuerpo;
class Interaccionador;

//-----------------------CLASE CUERPO------------------------------
class Cuerpo{//Interfaz de la clase
private:
  vector3D r,V,F;
  double m,R;

public:
  void Inicie(vector3D r0, vector3D V0, double m0,double R0); //Declaración de las funciones, solo se declara para que no se puedan urgar
  void BorreFuerza(void);
  void SumeFuerza(vector3D F0);
  void MoveR(double dt,double a);
  void MoveV(double dt, double a);
  double Positionx(void){return r.x();};
  double Positiony(void){return r.y();};
  double Positionz(void){return r.z();};
  void Draw(void);
  friend class Interaccionador;
};
void Cuerpo::Inicie(vector3D r0, vector3D V0, double m0,double R0){
  r=r0 ; V=V0;m=m0;  R=R0;
}
void Cuerpo::BorreFuerza(void){
  F.load(0,0,0);
}
void Cuerpo::SumeFuerza(vector3D F0){
  F+=F0;
}
void Cuerpo::MoveR(double dt, double a){
  r += V*a*dt;
}
void Cuerpo::MoveV(double dt, double a) {
  V += a*F*dt/m;
}
void Cuerpo::Draw(void){
  cout << " , " << r.x() << "+" << R << "*cos(t)"  << ","<< r.y() << "+"
  << R << "*sin(t)";
}
void StartFrame(void){cout << "plot 0,0 ";}
void EndFrame(void){cout << "\n";}


//----------------CLASE INTERACCIONADOR-----------------
class Interaccionador{
  private:
  public:
    void CalcFuerza(Cuerpo * Planet);
    void CalcFuerzaEntre(Cuerpo & Planet1, Cuerpo & Planet2);
};

void Interaccionador::CalcFuerzaEntre(Cuerpo & Planet1, Cuerpo & Planet2){
  vector3D r21,n; double d,F;
  r21 = Planet2.r-Planet1.r; d = r21.norm(); n=r21/d;
  F = G*Planet1.m*Planet2.m/pow(d,2);
  Planet1.SumeFuerza(F*n); Planet2.SumeFuerza(-F*n);
}
void Interaccionador::CalcFuerza(Cuerpo * Planet){
  int i,j;
  //Set all forces for all planets to cero
  for (i=0;i<N;i++){Planet[i].BorreFuerza();}
  //Calculate total force for the planet
  for (i=0;i<N;i++){
    for (j=i+1;j<N;j++){
      CalcFuerzaEntre(Planet[i],Planet[j]);
    }
  }
}



//---------------FUNCIONES GLOBALES--------------------
void Animation(void){
  cout << "set terminal gif animate" << "\n";
  cout << "set output '2Planetas.gif'"<< "\n";
  cout << "unset key" << "\n";
  cout << "set xrange [-12:12]"<<"\n";
  cout << "set yrange [-12:12]"<<"\n";
  cout << "set size ratio -1" << "\n"; //x-y ratio = 1
  cout << "set parametric" <<  "\n"; //Is used to set a parameter to draw curves
  cout << "set trange[0:7]" << "\n"; //The range of the parameter
  cout << "set isosamples 12" << "\n"; //Take 12 examples, in that way we're drawing a dodecahedron
}


int main (){
  int i;
  Cuerpo Planet[N]; //Instance(ejemplares) of Cuerpo
  Interaccionador Newton;
  vector3D r_0,r_1,V0,V1;
  double m0 = 10000, m1 = 1000, r = 11;
  double M = m0 + m1, x0 = -m1*r/M, x1=m0*r/M; 
  double omega = sqrt(G*M/pow(r,3)),v0 = omega*x0, v1 = omega*x1, T;
  double t, dt=0.01;
  double tdraw,tframe = T/100;

  T=2*M_PI/omega; 
  r_0.load(x0,0,0);V0.load(0,v0,0);
  r_1.load(x1,0,0);V1.load(0,v1,0);

  Planet[0].Inicie(r_0,0.8*V0,m0,1.0);
  Planet[1].Inicie(r_1,0.8*V1,m1,0.5);
  Newton.CalcFuerza(Planet);


  Animation();

  for (t=0,tdraw=0;t<1.6*T+dt;t+=dt,tdraw+=dt){
    // cout << Planet[1].Positionx() << "\t" << Planet[1].Positiony() << "\t"<<Planet[1].Positionz()<< "\n";
    //Move using PERFL
    for (int i = 0; i<N; i++){Planet[i].MoveR(dt,Zeta);} //1 paso
    Newton.CalcFuerza(Planet);
    for (int i = 0; i<N; i++){Planet[i].MoveV(dt,Lambda1);} //2 Paso
    for (int i = 0; i<N; i++){Planet[i].MoveR(dt,Chi);} //3 Paso
    Newton.CalcFuerza(Planet);
    for (int i = 0; i<N; i++){Planet[i].MoveV(dt,Lambda);} //4 Paso
    for (int i = 0; i<N; i++){Planet[i].MoveR(dt,Zetachi);} //5 Paso
    Newton.CalcFuerza(Planet);
    for (int i = 0; i<N; i++){Planet[i].MoveV(dt,Lambda);} //6 Paso
    for (int i = 0; i<N; i++){Planet[i].MoveR(dt,Chi);} //7 Paso
    Newton.CalcFuerza(Planet);
    for (int i = 0; i<N; i++){Planet[i].MoveV(dt,Lambda1);} //8 Paso
    for (int i = 0; i<N; i++){Planet[i].MoveR(dt,Zeta);} // 9 Paso
    //Draw it!
    if(tdraw > tframe){
      StartFrame();
      for(i=0;i<N;i++){
        Planet[i].Draw();
      }
      EndFrame();
      tdraw=0;
    }
  }
  return 0;
}

