/* ---------------------------------------------------------------------
 * $Id: step-1.cc 31349 2013-10-20 19:07:06Z maier $
 *
 * Copyright (C) 1999 - 2013 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------

 */

// @sect3{Include files}

// The most fundamental class in the library is the Triangulation class, which
// is declared here:
#include <deal.II/grid/tria.h>
// We need the following two includes for loops over cells and/or faces:
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
// Here are some functions to generate standard grids:
#include <deal.II/grid/grid_generator.h>
// We would like to use boundaries which are not straight lines, so we import
// some classes which predefine some boundary descriptions:
#include <deal.II/grid/tria_boundary_lib.h>
// Output of grids in various graphics formats:
#include <deal.II/grid/grid_out.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria_boundary_lib.h>
#include <grid/grid_generator.h>
#include <grid/grid_refinement.h>
#include <grid/grid_in.h>
#include <grid/intergrid_map.h>
#include <grid/grid_tools.h>
#include <grid/filtered_iterator.h>

// This is needed for C++ output:
#include <fstream>
// And this for the declarations of the `sqrt' and `fabs' functions:
#include <cmath>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>     


#define _USE_MATH_DEFINES

#include <math.h>
// The final step in importing deal.II is this: All deal.II functions and
// classes are in a namespace <code>dealii</code>, to make sure they don't
// clash with symbols from other libraries you may want to use in conjunction
// with deal.II. One could use these functions and classes by prefixing every
// use of these names by <code>dealii::</code>, but that would quickly become
// cumbersome and annoying. Rather, we simply import the entire deal.II
// namespace for general use:
using namespace dealii;


double newton_root(double x_in,double y_in,double q,double aa);
double newton_root2(double x_in,double y_in,double C,double m, double b);

Point<2> go_back (const Point<2> &in,double eps);
Point<2> go_back_local (const Point<2> &in,double eps,double eps_max);
// @sect3{Creating the first mesh}


Point<2> stretch (const Point<2> &in)
{
  return Point<2>(5.0*in(0),
                  50.0*in(1));
}

Point<2> return_map (const Point<2> &in)
{
  bool edge = (std::abs(in(0)- 68.)/68. < .01);
  // //bool edge = (in(0) == 68.0);
  //if (edge)

  double loc = 2.0;
  double eps_max = .2;
  double center = 68.0;

  double eps = eps_max*std::exp(-1.0*std::pow((in(0)- 68.0)/loc,2));
  // double yoffset = std::exp(-1.0*std::pow((in(0)- 68.)/3.0,2));
  
  Point<2> offset = go_back_local(Point<2>(68.0,0),eps,eps_max);
  Point<2> sol_curve = go_back(Point<2>(68.0,in(1)-offset(1)-0.35*std::sin(2.0*in(1)-M_PI/3.0)),eps);
  
  //Point<2> temp  go_back_local(Point<2>(in(0),in(1)-offset(1)+.25*(eps/eps_max)*std::cos(.030*(in(0)-center))*std::sin(2.0*in(1))),eps,eps_max); 
  //offset(1) = std::fmod(offset(1)+M_PI,M_PI);
  //Point<2> temp =  go_back_local(Point<2>(in(0),in(1)-offset(1)),eps,eps_max); 
  double dx = sol_curve(0) - in(0);
  double dy = sol_curve(1) - in(1);
  
  dx = dx*std::exp(-1.0*std::pow(68.0 -in(0),2)/1.);
  dy = dy*std::exp(-1.0*std::pow(68.0 -in(0),2)/1.); // 3bt
  Point<2> temp = in+Point<2>(dx,dy);
  
  
  if (in(1) ==0 || in(1) == 2.0*M_PI)
    temp(1) = in(1);
  return temp;
  // else
  //   return  Point<2>(in(0),in(1));/

}

// double yoffset(double x_in,eps,eps_max){
//   go_back_local(Point<2>(in(0),
// 				1.0*in(1)+yoffset(in(0))),eps,.02);

// }

void laplace_grid ()
{
 
  double rmin = 38.0;
  double a= 60.0;
  double b = 68.0;
  double rmax = 98.0;
 
  const unsigned int dim = 2;
  
  Triangulation< dim >	tria;

  Point<dim> corner1;
  corner1[0]=rmin;
  corner1[1]=0;
  Point<dim> corner2;
  corner2[0]=rmax;
  corner2[1]=2.0*M_PI;
  
  
  
  GridGenerator::hyper_rectangle(tria,corner1,corner2);

  //boundary indicators?
  tria.begin_active ()->face (0)->set_boundary_indicator (1);
  tria.begin_active ()->face (1)->set_boundary_indicator (2);
  
  tria.begin_active ()->face (2)->set_boundary_indicator (3);
  tria.begin_active ()->face (3)->set_boundary_indicator (4);


  tria.refine_global (2);
  
  

 
  std::map< unsigned int, Point< 2 > > new_points;
  

  //allocate the map
  
  const unsigned int N = 16;

  double a1 = -.01;

 

  Triangulation<2>::active_cell_iterator
    cell= tria.begin_active(),
    endc = tria.end();
 

  const RefinementCase< dim > ref_case;// = RefinementCase< dim>::isotropic_refinement();
  
  //ref_case =  RefinementCase< dim >::cut_axis(0);	
 
  
  int j = 0;
  while (j < 3)
    {
      cell = tria.begin_active();
      endc = tria.end(); 
      for ( ; cell != endc; ++cell){
  	for (unsigned int v=0; v < GeometryInfo<2>::vertices_per_cell;++v) {
  	  double x0 = cell->vertex(v)[0];
  	  double y0 = cell->vertex(v)[1];
  	  if (x0<rmax-(rmax-rmin)/4.0 and x0 > rmin+(rmax-rmin)/4.0) {
	    cell->set_refine_flag (ref_case.cut_axis(1));
	  }
  	}
      }
      tria.execute_coarsening_and_refinement ();
      j++;
    }
                                       

  

  double ymin_r =100;
  double ymax_r = -100;
  double ymin_l =100.0;
  double ymax_l = -100;
  double xmin_t = 100;
  double xmin_b = 100;
  double x0,y0;
  
  
  GridTools::transform(&return_map, tria);
  
  
  j = 0;
  while (j < 2)
    {
      cell = tria.begin_active();
      endc = tria.end(); 
      for ( ; cell != endc; ++cell){
  	for (unsigned int v=0; v < GeometryInfo<2>::vertices_per_cell;++v) {
  	  double x0 = cell->vertex(v)[0];
  	  double y0 = cell->vertex(v)[1];
	  if (x0 < rmax-(rmax-rmin)/4.0 and x0 > rmin+(rmax-rmin)/4.0)
	    cell->set_refine_flag (ref_case.cut_axis(0));
  	}
      }
      tria.execute_coarsening_and_refinement ();
      j++;
    }

  GridTools::transform(&stretch, tria);
  std::cout << "\n\tThere are " <<  tria.n_active_cells()<<std::endl;

  //tria.execute_coarsening_and_refinement ();
  //std::ofstream out ("grid-lap.eps");
  std::ofstream out ("Ullman_mesh.inp");
  GridOut grid_out;
  // grid_out.write_eps (tria, out);
  grid_out.write_ucd (tria, out);


}

// @sect3{Creating the second mesh}



double newton_root(double x_in,double y_in,double q,double aa){
  double y_out = y_in;
  double atol = 1.0;
  int iter = 0;
  int max_iter = 300;
  double f;
  double J;

  while (atol > .0001 && iter < max_iter)
    {
      f = -y_in + y_out + (2.0*M_PI/q) + aa*std::cos(y_out); 
      J = 1.0 - 2.0*aa*std::sin(y_out)*std::cos(y_out);
      atol = fabs(f);
      y_out = y_out - 1.0*f/J;
      //std::cout<< y_out<<std::endl;
      iter++;
      
    }
  //std::cout<< (y_out - y_in)<<std::endl;
  return y_out;
}

double newton_root2(double x_in,double y_in,double C,double m, double b){
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
Point<2> go_back_local (const Point<2> &in,double eps,double eps_max)
{
  double x0 = in(0);
  double y0 = in(1);

  

  double aa = -.02;
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
  y_old = newton_root2(x0,y0+0.0*(2.0*M_PI/q0),C,m,b);//;
  //std::cout<<y_old<<std::endl;
  x_old = x0 + (m*b*C)/(m-1)*std::pow((x0/b),(m-1)) *std::sin(m*y_old);
  xx = x_old/a;
  
  q = q0*(xx*xx)/(1-(1+2*xx*xx)*std::pow((1-xx),(mu+1)));
  y_old2 = newton_root(x_old,y_old,q,aa);
  y_old2 = std::fmod(y_old2+0.0*M_PI,2.0*M_PI);
  x_old2 = x_old*(1.0 - aa*std::sin(y_old2));


  // if (x_old2 >= .999*b && x0 <= b){
  //   double overstep = (x_old2-b);
  //   double start = 1.*(b-x0);
  //   x_old2 = b-1.*std::exp(-1.0*std::pow(1./start,2));
  // }
  // else if (x_old2 > b)
  //   x_old2 = x0;
 

  double dx = x_old2-x0;
  double dy = y_old2-y0;
  
  
  // dx = dx*std::exp(-1.0*std::pow(eps-eps_max,2)/.01);
  // dy = dy*std::exp(-1.0*std::pow(eps-eps_max,2)/.01);
  dx = dx*std::exp(-1.0*std::pow(68.0 -x0,2)/5.);
  dy = dy*std::exp(-1.0*std::pow(68.0 -x0,2)/5.);
  
  x_old2 = x0+dx;
  y_old2 = y0+dy;


  return Point<2>(x_old2,y_old2);
}

Point<2> go_back (const Point<2> &in,double eps)
{
  double x0 = in(0);
  double y0 = in(1);

  

  double aa = -.02;
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
  y_old = newton_root2(x0,y0+0.0*(2.0*M_PI/q0),C,m,b);//;
  //std::cout<<y_old<<std::endl;
  x_old = x0 + (m*b*C)/(m-1)*std::pow((x0/b),(m-1)) *std::sin(m*y_old);
  xx = x_old/a;
  
  q = q0*(xx*xx)/(1-(1+2*xx*xx)*std::pow((1-xx),(mu+1)));
  y_old2 = newton_root(x_old,y_old,q,aa);
  y_old2 = std::fmod(y_old2+0.0*M_PI,2.0*M_PI);
  x_old2 = x_old*(1.0 - aa*std::sin(y_old2));

  // if (x_old2 > b)
  //   return Point<2>(x0,y_old2);

  return Point<2>(x_old2,y_old2);
}

int main ()
{

  laplace_grid ();
}
