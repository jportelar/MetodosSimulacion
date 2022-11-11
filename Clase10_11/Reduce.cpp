#include <iostream>
#include <cmath>
using namespace std;

#define Lx 16
#define Nx 8
const int Mx=(Lx+Nx-1)/Nx;

//--------------- KERNELS ----------------
__global__ void ReduceThreads(float * d_a,float * d_ReduceResult){
  //Set global index
  int ix=threadIdx.x+blockIdx.x*blockDim.x;
  //Define and load shared memory
  extern __shared__ unsigned int temp[];
  temp[threadIdx.x]=d_a[ix];
  __syncthreads();
  //Reduce loop
  for(unsigned int s=blockDim.x/2;s>0;s>>=1){
    if(threadIdx.x<s)
      temp[threadIdx.x]+=temp[threadIdx.x+s];
    __syncthreads();
  }
  //Write to d_ReduceResult
  switch(threadIdx.x){
  case 0: d_ReduceResult[blockIdx.x]=temp[0]; break;
  default: ;
  }
}

int main(){
  //DECLARAR LAS VARIABLES
  int ix,bx;
  //DECLARAR LAS MATRICES
  //En el Host
  float h_a[Lx],h_ReduceResult[Mx],h_TotalResult[1];
  //En el Device
  float *d_a;             cudaMalloc((void**) &d_a,Lx*sizeof(float));
  float *d_ReduceResult;  cudaMalloc((void**) &d_ReduceResult,Mx*sizeof(float));
  float *d_TotalResult;   cudaMalloc((void**) &d_TotalResult,sizeof(float));

  //INICIALIZAR LOS DATOS
  //Cargarlos en el Host
  for(ix=0;ix<Lx;ix++) h_a[ix]=ix;
  for(bx=0;bx<Mx;bx++) h_ReduceResult[bx]=0;
  h_TotalResult[0]=0;
  //Enviarlos al Device
  cudaMemcpy(d_a,h_a,Lx*sizeof(float),cudaMemcpyHostToDevice);
  cudaMemcpy(d_ReduceResult,h_ReduceResult,Mx*sizeof(float),cudaMemcpyHostToDevice);
  cudaMemcpy(d_TotalResult,h_TotalResult,Mx*sizeof(float),cudaMemcpyHostToDevice);
  //Imprimirlos
  for(ix=0;ix<Lx;ix++) cout<<h_a[ix]<<" "; cout<<endl;
  
  //SUMAR DENTRO DE CADA BLOQUE
  //Procesar en la tarjeta gráfica
  dim3 ThreadsPerBlock1(Nx,1,1);
  dim3 BlocksPerGrid1(Mx,1,1);
  ReduceThreads<<<BlocksPerGrid1,ThreadsPerBlock1,Nx*sizeof(float)>>>(d_a,d_ReduceResult);
  //Imprimir los resultados parciales
  cudaMemcpy(h_ReduceResult,d_ReduceResult,Mx*sizeof(float),cudaMemcpyDeviceToHost);
  for(bx=0;bx<Mx;bx++) cout<<h_ReduceResult[bx]<<" "; cout<<endl;
  
  //SUMAR LAS SUMAS PARCIALES DE LOS BLOQUES EN UN SOLO DATO
  //Procesar en la tarjeta gráfica
  dim3 ThreadsPerBlock2(Mx,1,1);
  dim3 BlocksPerGrid2(1,1,1);
  ReduceThreads<<<BlocksPerGrid2,ThreadsPerBlock2,Mx*sizeof(float)>>>(d_ReduceResult,d_TotalResult);
  //Imprimir los resultados parciales
  cudaMemcpy(h_TotalResult,d_TotalResult,sizeof(float),cudaMemcpyDeviceToHost);
  cout<<h_TotalResult[0]<<endl;

  //LIBERAR MEMORIA
  cudaFree(d_a);  cudaFree(d_ReduceResult);  cudaFree(d_TotalResult);

  return 0;
}
