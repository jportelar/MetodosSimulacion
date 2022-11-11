#include <iostream>
#include <fstream>
using namespace std;

#define Lx 16
#define Nx 8
const int Mx=(Lx+Nx-1)/Nx;

//--------------- KERNELS ----------------
__global__ void SumeDeATres(float *d_a){
  int Left,Center,Right;
  Center=blockIdx.x*blockDim.x+threadIdx.x;
  Right=(Center+1)%Lx;   Left=(Center+Lx-1)%Lx;
  float suma=d_a[Left]+d_a[Center]+d_a[Right];
  //Sincronizar
  __syncthreads();
  //Escribir
  d_a[Center]=suma;
}

int main(void){
  int ix;
  //DECLARAR LAS MATRICES
  float h_a[Lx];
  float *d_a;  cudaMalloc((void**) &d_a,Lx*sizeof(float));

  //INICIALIZAR LOS DATOS
  //Cargarlos en el Host
  for(ix=0;ix<Lx;ix++){
    h_a[ix]=ix;
  }

  //IMPRIMIRLOS
  for(ix=0;ix<Lx;ix++)
    cout<<h_a[ix]<<" ";
  cout<<endl;

  //Enviarlos al Device
  cudaMemcpy(d_a,h_a,Lx*sizeof(float),cudaMemcpyHostToDevice);
  //PROCESAR EN LA TARJETA GRAFICA
  dim3 ThreadsPerBlock(Nx,1,1);
  dim3 BlocksPerGrid(Mx,1,1);
  SumeDeATres<<<BlocksPerGrid,ThreadsPerBlock>>>(d_a);

  //DEVOLVERLOS AL HOST
  cudaMemcpy(h_a,d_a,Lx*sizeof(float),cudaMemcpyDeviceToHost);

  //IMPRIMIRLOS
  for(ix=0;ix<Lx;ix++)
    cout<<h_a[ix]<<" ";
  cout<<endl;

  //LIBERAR MEMORIA
  cudaFree(d_a);

  return 0;
}
