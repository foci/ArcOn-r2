template <int dim>
double arcOn<dim>::Ullmann_newton_root(double x_in, double y_in, double b, double C, double m){
  double x_out = x_in;
  double atol = 1.0;
  int iter = 0;
  int max_iter = 300;
  double f;
  double J;

  while (atol > .0001 && iter < max_iter)
    {
      f = (-1.0*x_in + x_out +(m*b*C)/(m-1.0) * pow(x_out/b, m-1.0) * sin(m*y_in));
      J = 1.0 + (m*C) * pow(x_out/b,m-2.0) * sin(m*y_in);
      atol = fabs(f);
      x_out = x_out - 1.0*f/J;
      iter++;
      
    }
  //output<< (x_out - x_in)/fabs(x_in)<<endl;
  return x_out;
}

