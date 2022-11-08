#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

#define Lx 16
#define Nx 8
const int Mx=(Lx+Nx-1)/Nx;

//-----------Programa del device----------
//--------------Kernels-----------------
__global__ void AddTwoVectors(float* d_a, float* d_b,float* d_c){
  //Qué tarea me toca?
  int ix; ix=blockIdx.x*blockDim.x+threadIdx.x;
  d_c[ix]=d_a[ix]+d_b[ix];
}
//-----------Codigo del Host---------------
int main(void){
  int ix;
  //Declarar todas las variables por
  //------en el host--------
  float h_a[Lx],h_b[Lx],h_c[Lx];
  //------en el device-------
  float *d_a; cudaMalloc((void**) &d_a, Lx*sizeof(float));
  float *d_b; cudaMalloc((void**) &d_b, Lx*sizeof(float));
  float *d_c; cudaMalloc((void**) &d_c, Lx*sizeof(float));

  //llenar los datos que vamos a procesar
  for(ix=0;ix<Lx;ix++)
    {h_a[ix]=ix; h_b[ix]=2*ix;}
  
  //Enviar al device
  cudaMemcpy(d_a,h_a,Lx*sizeof(float),cudaMemcpyHostToDevice);
  cudaMemcpy(d_b,h_b,Lx*sizeof(float),cudaMemcpyHostToDevice);

  //Correr el Device
  dim3 ThreadsPerBlock(Nx,0,0);
  dim3 BlocksPerGrid(Mx,0,0);
  AddTwoVectors<<<BlocksPerGrid,ThreadsPerBlock>>>(d_a,d_b,d_c);

  cudaMemcpy(h_c,d_c,Lx*sizeof(float),cudaMemcpyDeviceToHost);

  for(ix=0;ix<Lx;ix++)
    cout<<h_c[ix]<<endl;

  cudaFree(d_a);cudaFree(d_b);cudaFree(d_c);
   
  return 0;
}

