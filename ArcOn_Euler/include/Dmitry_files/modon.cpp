#include <iostream>
#include <iomanip>

#include <typeinfo>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <sstream>
// very important unless you like to type a lot
using namespace std;



double Newton_root(double gamma_in,double kappa_in);

double a,c,lambda;

double modon(double x, double y,double kappa,double gamma)
{

  
  
  double A,B,C,phi,r;
  r = sqrt(pow(x,2.0) + pow(y,2.0));
  

  A = c*a/yn(1,kappa);
  B = c*a*(1.0 + pow(kappa/gamma,2.0));
  C = -pow(kappa/gamma,2.0)*(c*a)/jn(1,gamma);

  if(r > a)
    phi = A*yn(1,kappa*r/a)*(x/r);
  else
    phi = B*(x/a) + C*jn(1,gamma * r/a)*(x/r);

  return phi;
}


int main (int argc, char** argv)
{

  double beta;
  // double c = double(*argv[0]);
  // //double lambda = std::stof (*argv[1]);
  // double lambda = double(*argv[1]);
  // double beta = double(*argv[2]);
  cout << argv[3] <<endl;
  istringstream s(argv[3]);

  if (!(s >> c))
    cerr << "Invalid number " << argv[3] << '\n';
  
  istringstream ss(argv[1]);

  if (!(ss >> lambda))
    cerr << "Invalid number " << argv[1] << '\n';

  istringstream s2(argv[2]);

  if (!(s2 >> beta))
    cerr << "Invalid number " << argv[2] << '\n';
  
  cout << beta <<" "<<lambda/beta<<endl;

  //set the radius a to 1
  double a = .1;
  //we get to pick c at will 
  //c = 50;
  
  //set kappa based on c
  double kappa = a*sqrt(lambda*beta)/c;
  double gamma;

  gamma = Newton_root(M_PI,kappa); //will find the 1st root

  a = .1;
  //c = 1.0;
  // have to pick some x and y values 

  double x = .3;
  double y = .2;


  cout << "modon: "<<kappa<<" "<<gamma<<" "<<modon(.3,.1,kappa,gamma)<<endl;
  

  double phi, n;
  phi = modon(.3,.1,kappa,gamma);
  n =  (lambda/c)*(phi - c*x);
  // cout << "modon: "<<kappa<<" "<<gamma<<" "<<modon(.1,.1,kappa,gamma)<<endl;
  // cout << "modon: "<<kappa<<" "<<gamma<<" "<<modon(.1,1,kappa,gamma)<<endl;

}

// double get_gamma(){
//   return
// }
//you probably have better ways of writting a root finder

double Newton_root(double gamma_in,double kappa){
  
  double atol = 1.0;
  int iter = 0;
  int max_iter = 300;
  double f;
  double J;
  double gamma = gamma_in;

  while (atol > 1e-8 && iter < max_iter)
    {
      // f = (yn(2,kappa)/(kappa*yn(1,kappa)))  + jn(2,gamma)/(gamma*jn(1,gamma));
      // J = (-4.0 + pow(gamma,2)*(1.0 + pow(j0(gamma)/j1(gamma),2.0) )) /pow(3.0,gamma);
      f = (1.0/pow(gamma,2.0))*(yn(2,kappa)*gamma*jn(1,gamma)+kappa*yn(1,kappa)*jn(2,gamma));
      J = -(1.0/pow(gamma,2.0))*(kappa*jn(3.0,gamma)*yn(1,kappa) + gamma*jn(2,gamma)*yn(2,kappa));
      
      // f = (yn(2,kappa)*gamma*jn(1,gamma)+kappa*yn(1,kappa)*jn(2,gamma));
      // J = .5*kappa*(jn(1,gamma) - jn(3,gamma))*yn(1,kappa) + gamma*jn(0,gamma)*yn(2,kappa);
     

      //cout << iter<<": "<<gamma<<" "<<f/J <<endl;
      atol = fabs(f);
      gamma = gamma - 1.0*f/J;
      //cout <<gamma<<" "<<atol <<endl;
      iter++;
      
    }
  //output<< (x_out - x_in)/fabs(x_in)<<endl;
  cout << "atol "<<atol<<" "<<iter<<endl;
  return gamma;
}
