#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

int imax = 100000000;
const double gam = 1.0, m = 1.0, dt = 0.1;
const double X0 = 0.0, Y0 = 0.0, U0 = 0.0, V0 = 0.0;
const double c = 1.0;

double Fr(){
  double cc = sqrt(3.0*c/(2.0*dt));
  return (rand() / (double)RAND_MAX*2 - 1.0) * cc;
}

class calc{
public:
    double x, y, u, v, r;
  void output(int t){
    char filename[20];
    sprintf(filename,"work1_%d.dat", imax);
    FILE *fout;
    fout = fopen(filename, "a");
    if(fout == NULL) exit(1);
    fprintf(fout, "%.10f %.10f %10f %10f %10f\n", t * dt, x, y, u, v);
    fclose(fout);
  };
  calc(){
    x = X0, y = Y0, u = U0, v = V0, r = 0.0;
  };
  void evo(int t){
//    output(t);
    x += dt * u;
    y += dt * v;
    u += dt * (-gam*u + Fr()) / m;
    v += dt * (-gam*v + Fr()) / m;
    maxr();
  };
  void maxr(){
    double rr = sqrt(x*x+y*y);
    if(r < rr) r = rr;
  };
  void distance(){
    char filename[20];
    sprintf(filename,"work1.dat");
    FILE *fout;
    fout = fopen(filename, "a");
    if(fout == NULL) exit(1);
    fprintf(fout, "%d %.10f %.10f\n", imax, sqrt(x*x+y*y), r);
    fclose(fout);
  };
};

int main(){
  calc A;
  for(int i = 0; i <= imax; i++){
    A.evo(i);
    if(i == imax){
      A.distance();
    }
  }
  return 0;
}
