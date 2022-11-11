#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

#define Lx 16
#define Nx 8
const int Mx=(Lx+Nx-1)/Nx;

//--------------- KERNELS ----------------
__global__ void SumarTodos(float *d_a,float *d_Resultado){
  //cargar los datos a la memoria compartida
  int ix; ix=blockIdx.x*blockDim.x+threadIdx.x;
  //d_Resultado[0]+=d_a[ix];
  atomicAdd(&d_Resultado[0],d_a[ix]);
}

int main(void){
  int ix;
  //DECLARAR LAS MATRICES
  float h_a[Lx];         float *d_a;  cudaMalloc((void**) &d_a,Lx*sizeof(float));
  float h_Resultado[1];  float *d_Resultado;  cudaMalloc((void**) &d_Resultado,sizeof(float));

  //INICIALIZAR LOS DATOS
  //Cargarlos en el Host
  for(ix=0;ix<Lx;ix++){
    h_a[ix]=ix;
  }
  h_Resultado[0]=0;

  //IMPRIMIRLOS
  for(ix=0;ix<Lx;ix++)
    cout<<h_a[ix]<<" ";
  cout<<endl;

  //Enviarlos al Device
  cudaMemcpy(d_a,h_a,Lx*sizeof(float),cudaMemcpyHostToDevice);
  cudaMemcpy(d_Resultado,h_Resultado,sizeof(float),cudaMemcpyHostToDevice);
  //PROCESAR EN LA TARJETA GRAFICA
  dim3 ThreadsPerBlock(Nx,1,1);
  dim3 BlocksPerGrid(Mx,1,1);
  SumarTodos<<<BlocksPerGrid,ThreadsPerBlock>>>(d_a,d_Resultado);

  //DEVOLVERLOS AL HOST
  cudaMemcpy(h_Resultado,d_Resultado,sizeof(float),cudaMemcpyDeviceToHost);

  //IMPRIMIRLOS
  cout<<"EL resultado es ="<<h_Resultado[0]<<endl;

  //LIBERAR MEMORIA
  cudaFree(d_a);  cudaFree(d_Resultado);

  return 0;
}
