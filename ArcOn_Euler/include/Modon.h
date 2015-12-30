template <int dim>
double arcOn<dim>::modon(double x, double y,double kappa,double gamma, double c, double a)
{
  
  double A,B,C,phi,r;
  r = std::sqrt(std::pow(x,2.0) + std::pow(y,2.0)) + 1e-12;
    //std::abs(std::numeric_limits<double>::epsilon());
  

  A = c*a/yn(1,kappa);
  B = c*a*(1.0 + std::pow(kappa/gamma,2.0));
  C = -std::pow(kappa/gamma,2.0)*(c*a)/jn(1,gamma);

  if(r > a)
    phi = A*yn(1,kappa*r/a)*(x/r);
  else
    phi = B*(x/a) + C*jn(1,gamma * r/a)*(x/r);

  return phi;

  //pcout << "phi = " << phi << std::endl;

}
