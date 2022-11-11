#include <iostream>
#include <fstream>
using namespace std;

#define Lx 16
#define Nx 8
const int Mx=(Lx+Nx-1)/Nx;

//--------------- KERNELS ----------------
__global__ void SumeDeATres(float *d_a){
  extern __shared__ float Shared[];
  //cargar los datos a la memoria compartida
  int ix; ix=blockIdx.x*blockDim.x+threadIdx.x;
  int is; is=threadIdx.x;
  //Cargar la memoria compartida
  switch(is){
  case 0           : Shared[is+1]=d_a[ix]; Shared[0]=d_a[(ix-1+Lx)%Lx]; break; 
  case Nx-1        : Shared[is+1]=d_a[ix]; Shared[is+2]=d_a[(ix+1)%Lx]; break;
  default          : Shared[is+1]=d_a[ix];
  }
  //Sincronizar
  __syncthreads();
  //Hacer la suma y escribir
  int Left=is,   Center=is+1,  Right=is+2;
  d_a[ix]=Shared[Left]+Shared[Center]+Shared[Right];
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
  SumeDeATres<<<BlocksPerGrid,ThreadsPerBlock,(Nx+2)*sizeof(float)>>>(d_a);

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
