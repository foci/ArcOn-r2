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
  double aa = -.02;
  //double m = 3.0;
  double l = 10.0;
  double R = 85;
  double q0 = 3.0;
  double b = 68.0;
  double a= 60.0;
  
  double nu = 2.0;
  
  int max_orbit = 1000;	
  
  //has to be consistent with Ullmann map that is expressed in cm and rad in x and y respectively, so rho_s = .2 and Ly = 2pi*50
  double rho_s = .2; 
  y = y*(2.0*M_PI/Ly); //remap to rad for ullmann map

  //x  = x/rho_s; //remap to cm for ullmann amp
  //
  //q = q0*pow(x/a,2.0)/(1.0-pow(1.0-x/a,nu+1.0));  

  double x_new;
  double y_new;
  double x_new2;
  
  while(count < max_orbit and not hit_divert){
    // x_new = x;
    // y_new = y;
    x_new = x_new/(1-aa*sin(y_new));
    //q = q0*pow((x_new/a),2.0);
    double C = ((2.0*m*l*pow(a,2.0))/(R*q*pow(b,2.0)))*eps;
    /* q = q0* pow(x_new/a,2.0)/(1.0-double(x_new<a)*pow(1.0-x_new/a,nu+1.0));  */
    /* if (q >qmax)  */
    /*   q = qmax; */
    // q = q+ double(x_new>b)*q0*pow(b/a,2.0)/(1.0-pow(1.0-b/a,nu+1.0)); 
   
    C = ((2*m*l*pow(a,2.0))/(R*q*pow(b,2.0)))*eps;
   

    y_new =  (y_new+ 2*M_PI/q + aa*cos(y_new));
    y_new = fmod(y_new,2*M_PI);
    
    x_new2 = Ullmann_newton_root(x_new,y_new,b,C,m);
    
   
    
    y_new = (y_new - C*pow(x_new2/b , m-2) * cos(m*y_new));
    y_new = fmod(y_new,2*M_PI);
    x_new = x_new2;
    //output <<x_new<<endl;
    hit_divert = (x_new > b);// or x>1.2*b or x_new <0);// or (x_new <  and x < b);
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

