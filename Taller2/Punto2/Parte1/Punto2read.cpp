#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

const double slope=-5.2;
const double intercept=-174.5;

void Divide(const char * MyFile){
  ifstream infile; int n=0; 
  double ix, iy;
  infile.open(MyFile);
  while(infile>>ix>>iy){
    cout<<ix<<" "<<iy<<endl;
  } 
  infile.close();
}
 
int main(void){
  Divide("LongitudDeOnda20.dat");  
  return 0;
}
