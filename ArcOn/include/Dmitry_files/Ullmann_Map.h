template <int dim>
double arcOn<dim>::UllmannMap(double x, double Lx, double y, double Ly, double x_sol, double eps, double m )
{
  /* Set the constant parameters of the map */
  /* Set these in Run.h :: double Lx = ,  double Ly = ; */
  /* (x0,y0) = bottom lhs of box */
  /* R => number of orbits -> L = 2*pi*R*( # of orbits ) */
  /* k = Chirikov–Taylor constant */
  /* q0 = characteristic value of q */
  /* q = (# toroidal transit) / (# of poloidal transits )*/
  /* int max_orbits */
  /* M_PI = pi const */ 

  int count = 0;
  bool hit_divert = false;
  bool inSOL;
  double q,qmax;

  //int max_orbit = 100;
  
  double L = 0.0;
  //double eps = .5;
  double aa = -.01;
  //double m = 3.0;
  double l = 10.0;
  double R = 85;
  double q0 = 3.0;
  double b = 68.0;
  double a= 60.0;
  
  double nu = 2.0;
  
  int max_orbit = 1000;	
  
 

  // let's consider a region from psi = .7 to 1.3, if rho_s = .2 this
  // corresponds to last closed flux surface +- rho_s * Lx/2
  double width = eps*30.; //a crude estimate
  double offset = x_sol * .2;
  
  x =b - 3*width +  (3*(b-a)+3*width)*(x/Lx); // the total region size is then about 26 to about 35 cm 
  //x = b - width + 3*width*(x/Lx);
  y = y*(2.0*M_PI/Ly);
  //
  //q = q0*pow(x/a,2.0)/(1.0-pow(1.0-x/a,nu+1.0));  

  double x_new;
  double y_new;
  double x_new2;

  double C = ((2.0*m*l*pow(a,2.0))/(R*q0*pow(b,2.0)))*eps;
  x_new = x;
  y_new = y;
  
  qmax = q0*pow(b/a,2.0); 
  while(count < max_orbit and not hit_divert){
    // x_new = x;
    // y_new = y;
    x_new = x_new/(1-aa*sin(y_new));
    //q = q0*pow((x_new/a),2.0);

    q = q0* pow(x_new/a,2.0)/(1.0-double(x_new<a)*pow(1.0-x_new/a,nu+1.0)); 
    if (q >qmax) 
      q = qmax;
    // q = q+ double(x_new>b)*q0*pow(b/a,2.0)/(1.0-pow(1.0-b/a,nu+1.0)); 
    
    C = ((2*m*l*pow(a,2.0))/(R*q*pow(b,2.0)))*eps;
    y_new =  (y_new+ 2*M_PI/q + aa*cos(y_new));
    y_new = fmod(y_new,2*M_PI);
    
    x_new2 = Ullmann_newton_root(x_new,y_new,b,C,m);
    
    //output<< (-1.0*x_new+ x_new2 +(m*b*C)/(m-1.0) * pow(x_new2/b, m-1.0) * sin(m*y_new))<<endl;

    //output<< "old: " <<(-1.0*x_new+ x_new +(m*b*C)/(m-1.0) * pow(x_new/b, m-1.0) * sin(m*y_new))<<endl;

    //chi = (-x_new + x_out +(m*b*C)/(m-1)*(x_out/b)**(m-1) *np.sin(m*y_new))**2
    //x_new2 = (newton_krylov(func,x_new));
    
    //q = q0*pow(x_new2/a,2.0)/(1.0-w(1.0-x_new2/a,nu+1.0));  
    //C = ((2*m*l*pow(a,2.0))/(R*q*pow(b,2.0)))*eps;
    
    y_new = (y_new - C*pow(x_new2/b , m-2) * cos(m*y_new));
    y_new = fmod(y_new,2*M_PI);
    x_new = x_new2;
    //output <<x_new<<endl;
    hit_divert = (x_new > b or x>1.2*b or x_new <0);// or (x_new <  and x < b);
    count++;

    if (!hit_divert) {
      L = L + 2.0*M_PI*q*R;
    }
    else{
      L = L + 2.0*M_PI*qmax*R;
    }

    if(count == max_orbit){
      L = -1;
    }
    
  }
  //output<<L<<endl;
  //if (L == max_orbit) {L = 10.0;}
  return L;
}

