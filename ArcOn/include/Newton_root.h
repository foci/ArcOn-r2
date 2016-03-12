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
      f = (1.0/std::pow(gamma,2.0))*(yn(2,kappa)*gamma*jn(1,gamma)+kappa*yn(1,kappa)*jn(2,gamma));
      J = -(1.0/std::pow(gamma,2.0))*(kappa*jn(3.0,gamma)*yn(1,kappa) + gamma*jn(2,gamma)*yn(2,kappa));
	
      atol = fabs(f);
      gamma = gamma - 1.0*f/J;
      iter++;
    }

  return gamma;

}
