template <int dim>
double arcOn<dim>::Ullman_newton_root2(double x_in,double y_in,double C,double m, double b){
  double y_out = y_in;
  double atol = 1.0;
  int iter = 0;
  int max_iter = 1000;
  double f;
  double J;
  
  while (atol > .0001 && iter < max_iter)
    {
      f = (-y_in + y_out - C*std::pow((x_in/b),(m-2)) *std::cos(m*y_out));
      J = 1.0 + C*std::pow((x_in/b),(m-2))*m*std::sin(m*y_out);
      atol = fabs(f);
      y_out = y_out - 1.0*f/J;
      //std::cout<< y_out<<std::endl;
      iter++;
      
    }
  //std::cout<< (y_out - y_in)<<std::endl;
  return y_out;
}
