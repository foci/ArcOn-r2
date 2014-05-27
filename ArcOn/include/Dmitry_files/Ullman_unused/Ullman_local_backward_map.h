template <int dim>
Point<2> arcOn<dim>::ullman_go_back_local (const Point<2> &in,double eps,double eps_max)
{
  double x0 = in(0);
  double y0 = in(1);
  double aa = -.01;
  double q0 = 2.0;
  double mu = 2.0;
  double a= 60.0;
  double m = 2;
  double b = 68.0;
  double l = 10.0;
  //double eps = .010;
  double R_0 = 85.0;
  
  double C = ((2*m*l*a*a)/(R_0*q0*b*b))*eps;
  double y_old,xx,x_old,y_old2,x_old2,q;
  y_old = Ullman_newton_root2(x0,y0+0.0*(2.0*M_PI/q0),C,m,b);//;
  //std::cout<<y_old<<std::endl;
  x_old = x0 + (m*b*C)/(m-1)*std::pow((x0/b),(m-1)) *std::sin(m*y_old);
  xx = x_old/a;
  
  q = q0*(xx*xx)/(1-(1+2*xx*xx)*std::pow((1-xx),(mu+1)));
  y_old2 = Ullman_newton_root(x_old,y_old,q,aa);
  y_old2 = std::fmod(y_old2+0.0*M_PI,2.0*M_PI);
  x_old2 = x_old*(1.0 - aa*std::sin(y_old2));

  // if (x_old2 > b)
  //   x_old2 = (3.*b+2.*x0)/5.0;
  //   return Point<2>(x0,y_old2);

  double dx = x_old2-x0;
  double dy = y_old2-y0;
  
  // dx = dx*std::exp(-1.0*std::pow(eps-eps_max,2)/.01);
  // dy = dy*std::exp(-1.0*std::pow(eps-eps_max,2)/.01);
  dx = dx*std::exp(-1.0*std::pow(68.0 -x0,2)/25.);
  dy = dy*std::exp(-1.0*std::pow(68.0 -x0,2)/25.);
  
  x_old2 = x0+dx;
  y_old2 = y0+dy;

  return Point<2>(x_old2,y_old2);
}
