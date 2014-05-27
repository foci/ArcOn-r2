/*These are patch functions for small things in the code*/
template<typename T>
T sum(const std::vector<T>& v)
{
  T s(0);
  for(unsigned int i=0;i<v.size();i++) s+=v[i];
  return s;
}

template<typename T>
void zero(std::vector<T>& v)
{
  for(unsigned int i=0;i<v.size();i++) v[i]=T(0);
}


template<typename T>
inline std::ostream& operator<<( std::ostream& os, const std::vector<T>& v)
{
  for(unsigned int i=0; i<v.size(); i++) {
    os<<" "<<v[i];
  }
  return os;
}

template<typename NUMBER> std::pair<NUMBER,NUMBER> minmax(const Vector<NUMBER>& V){
  NUMBER Vmin=std::numeric_limits<NUMBER>::max();
  NUMBER Vmax=-std::numeric_limits<NUMBER>::max();
  unsigned int sz = V.size();
  for(unsigned int i = 0; i<sz; i++){
    if( V(i) < Vmin ) Vmin = V(i);
    if( V(i) > Vmax ) Vmax = V(i);
  }
  return std::make_pair(Vmin,Vmax);
}
//___________________________________________



/* double Newton_root(double gamma_in,double kappa_in); */
  
/* double Newton_root(double gamma_in,double kappa){ */
    
/*   double atol = 1.0; */
/*   int iter = 0; */
/*   int max_iter = 300; */
/*   double f; */
/*   double J; */
/*   double gamma = gamma_in; */
    
/*   while (atol > 1e-8 && iter < max_iter) */
/*     { */
/* 	// f = (yn(2,kappa)/(kappa*yn(1,kappa)))  + jn(2,gamma)/(gamma*jn(1,gamma)); */
/* 	// J = (-4.0 + pow(gamma,2)*(1.0 + pow(j0(gamma)/j1(gamma),2.0) )) /pow(3.0,gamma); */
/* 	f = (1.0/pow(gamma,2.0))*(yn(2,kappa)*gamma*jn(1,gamma)+kappa*yn(1,kappa)*jn(2,gamma)); */
/* 	J = -(1.0/pow(gamma,2.0))*(kappa*jn(3.0,gamma)*yn(1,kappa) + gamma*jn(2,gamma)*yn(2,kappa)); */
	
/* 	// f = (yn(2,kappa)*gamma*jn(1,gamma)+kappa*yn(1,kappa)*jn(2,gamma)); */
/* 	// J = .5*kappa*(jn(1,gamma) - jn(3,gamma))*yn(1,kappa) + gamma*jn(0,gamma)*yn(2,kappa); */
	
	
/* 	//cout << iter<<": "<<gamma<<" "<<f/J <<std::endl; */
/* 	atol = fabs(f); */
/* 	gamma = gamma - 1.0*f/J; */
/* 	//cout <<gamma<<" "<<atol <<std::endl; */
/* 	iter++; */
	
/*     } */
/*   //output<< (x_out - x_in)/fabs(x_in)<<std::endl; */
/*   pcout << "atol "<<atol<<" "<<iter<<std::endl; */
/*   return gamma; */
/* } */
 

/* double modon(double x, double y,double kappa,double gamma) */
/* { */
  
/*   double A,B,C,phi,r; */
/*   r = sqrt(pow(x,2.0) + pow(y,2.0)); */
  

/*   A = c*a/yn(1,kappa); */
/*   B = c*a*(1.0 + pow(kappa/gamma,2.0)); */
/*   C = -pow(kappa/gamma,2.0)*(c*a)/jn(1,gamma); */

/*   if(r > a) */
/*     phi = A*yn(1,kappa*r/a)*(x/r); */
/*   else */
/*     phi = B*(x/a) + C*jn(1,gamma * r/a)*(x/r); */

/*   return phi; */

/* } */
