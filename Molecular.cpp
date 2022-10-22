#include <iostream>
#include <cmath>

using namespace std;
//CONSTANTES GLOBALES
const double g=9.8;
//DECLARACIÓN DE CLASES
class Cuerpo; // se declara la interfaz, luego se declara el interior
//--------------Clase Cuerpo--------------------------
class Cuerpo{ 
private: //datos que no son accesibles, solo son importantes las ordenes  
  double x,y,Vx,Vy,Fx,Fy,m,R; 
public:  //datos que son accesibles 
  void Inicie(double x0,double y0,double Vx0,double Vy0,double m0,double R0);
  void CalculeFuerza(void); //para indicar que está vacío
  void Muevase(double dt);
  double Getx(void){return x;}; //Función Inline
  double Gety(void){return y;}; 
};
//aquí acabaría un eventual .h

//----------------------------------------------------
void Cuerpo::Inicie(double x0,double y0,double Vx0,double Vy0,double m0,double R0){// tengo acceso a todas las variables privadas porque estoy con Cuerpo::
  x=x0; y=y0; Vx=Vx0; Vy=Vy0; m=m0; R=R0;
}

void Cuerpo::CalculeFuerza(void){
  Fx=0; Fy=-m*g;
}

void Cuerpo::Muevase(double dt){
   x+=Vx*dt;    y+=Vy*dt;
  Vx+=Fx/m*dt; Vy+=Fy/m*dt;
}
//--------------Funciones Globales--------------------
int main(){
  Cuerpo Balon; //Ejemplares de la clase Cuerpo (Instance)
  double t, dt=0.1;
  //----------(x0,y0,Vx0,Vy0,m0,R0)
  Balon.Inicie( 0, 0, 20, 15, 0.453, 0.15);
  for(t=0;t<3.5;t+=dt){
    cout<<Balon.Getx()<<" "<<Balon.Gety()<<endl;
    Balon.CalculeFuerza();
    Balon.Muevase(dt);
  }
  return 0;
}
