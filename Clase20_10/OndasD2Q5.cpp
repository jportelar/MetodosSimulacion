#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

const int Lx=128;
const int Ly=128;

const int Q=5;
const double W0=1.0/3;

const double C=0.5; // C<0.707 cells/click
const double C2=C*C;
const double AUX0=1-3*C2*(1-W0);

const double tau=0.5;
const double Utau=1.0/tau;
const double UmUtau=1-Utau;

//-----Clase LatticeBoltzmann----//
class LatticeBoltzmann{
private:
  double w[Q]; //pesos
  int Vx[Q],Vy[Q]; //vectores de vel
  double *f,*fnew; //func. de distribución
public:
  LatticeBoltzmann(void);
  ~LatticeBoltzmann(void);
  int n(int ix,int iy,int i){return (ix*Ly+iy)*Q+i;};
  double rho(int ix,int iy,bool UseNew);
};

LatticeBoltzmann::LatticeBoltzmann(void){
  //Set the weights
  w[0]=W0; w[1]=w[2]=w[3]=w[4]=(1.0-W0)/4;
  //Set the velocity vectors
  Vx[0]=0;  Vx[1]=1;  Vx[2]=0;  Vx[3]=-1; Vx[4]=0;
  Vy[0]=0;  Vy[1]=0;  Vy[2]=1;  Vy[3]=0;  Vy[4]=-1;
  //Create the dynamic arrays
  int ArraySize=Lx*Ly*Q;
  f=new double [ArraySize];  fnew=new double [ArraySize];
}
LatticeBoltzmann::~LatticeBoltzmann(void){
    delete[] f;  delete[] fnew;
}

double LatticeBoltzmann::rho(int ix,int iy,bool UseNew){
  double sum; int i,n0;
  for(sum=0,i=0;i<Q;i++){
    n0=n(ix,iy,i);
    if(UseNew) sum+=fnew[n0]; else sum+=f[n0];
  }
  return sum;
}



class LatticeGas{
private:
  int V[Q]; //V[i] i=0 (derecha) i=1 (izquierda)
  double f[Lx][Q], fnew[Lx][Q]; //n[ix][i]
public:
  LatticeGas(void);
  void Inicie(int N, double mu, double sigma);
  void Colisione();
  void GrafiqueRho(void);
  void Adveccione(void);
  double rho(int ix);
  double Varianza(void);
};
LatticeGas::LatticeGas(void){
  V[0]=1; V[1]=-1;
}
void LatticeGas::Inicie(int N, double mu, double sigma){
  for(int ix=0;ix<Lx;ix++){
    double rho0 = (N*1.0/(sigma*sqrt(2*M_PI)))*exp(-0.5*pow((ix-mu)/sigma,2.0));
    for(int i=0; i<Q;i++)
      f[ix][i]=rho0/Q;
  }
}
double LatticeGas::rho(int ix){
  double suma; int i;
  for(suma=0, i=0;i<Q;i++)
    suma+=f[ix][i];
  return suma;      
}
void LatticeGas::GrafiqueRho(void){
  for(int ix=0;ix<Lx;ix++)
    cout<<ix<<" "<<rho(ix)<<endl; 
}


void LatticeGas::Colisione(){
  int ix, i, j;
  for(ix=0;ix<Lx;ix++){ //para cada celda
    for(i=0;i<Q;i++){
      j=(i+1)%Q;
      fnew[ix][i]=fnew[ix][i]+(1-p)*(f[ix][j]-f[ix][i]);
    }
  }
}
void LatticeGas::Adveccione(void){
  for(int ix=0;ix<Lx;ix++){
    for(int i=0;i<Q;i++){
      f[(ix+V[i]+Lx)%Lx][i]=fnew[ix][i];
    }
  }
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
  sigma2/=(1.0*N);

  return sigma2;
}

//------Programa Principal------//
int main(void){
  LatticeGas Difusion;
  int N=6; double mu=Lx/2, sigma=Lx/8;
  int t, tmax=400;

 
  Difusion.Inicie(N,mu,sigma);
  for(t=0;t<tmax;t++){
    cout<<t<<" "<<Difusion.Varianza()<<endl;
    Difusion.Colisione();
    Difusion.Adveccione();
  }
  //Difusion.GrafiqueRho();
  return 0;
}
