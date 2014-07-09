template <int dim>
double arcOn<dim>::newton_root(double gamma, double kappa)
{

  double atol = 1.0;
  int iter = 0;
  int max_iter = 300;
  double f;
  double J;
    
  while (atol > 1e-8 && iter < max_iter)
    {
      // f = (yn(2,kappa)/(kappa*yn(1,kappa)))  + jn(2,gamma)/(gamma*jn(1,gamma));
      // J = (-4.0 + std::pow(gamma,2)*(1.0 + std::pow(j0(gamma)/j1(gamma),2.0) )) /std::pow(3.0,gamma);
      f = (1.0/std::pow(gamma,2.0))*(yn(2,kappa)*gamma*jn(1,gamma)+kappa*yn(1,kappa)*jn(2,gamma));
      J = -(1.0/std::pow(gamma,2.0))*(kappa*jn(3.0,gamma)*yn(1,kappa) + gamma*jn(2,gamma)*yn(2,kappa));
	
      // f = (yn(2,kappa)*gamma*jn(1,gamma)+kappa*yn(1,kappa)*jn(2,gamma));
      // J = .5*kappa*(jn(1,gamma) - jn(3,gamma))*yn(1,kappa) + gamma*jn(0,gamma)*yn(2,kappa);
	
	
      //cout << iter<<": "<<gamma<<" "<<f/J <<std::endl;
      atol = fabs(f);
      gamma = gamma - 1.0*f/J;
      //cout <<gamma<<" "<<atol <<std::endl;
      iter++;
	
    }
  //output<< (x_out - x_in)/fabs(x_in)<<std::endl;
  //pcout << "atol "<<atol<<" "<<iter<<std::endl;
  return gamma;

}
