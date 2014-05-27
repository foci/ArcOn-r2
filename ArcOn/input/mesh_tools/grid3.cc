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
Point<2> go_back (const Point<2> &in);
// @sect3{Creating the first mesh}


Point<2> stretch (const Point<2> &in)
{
  return Point<2>(in(0),
                  10*in(1));
}
void laplace_grid ()
{
 
  Triangulation< 2 >	tria;
  Point<2> corner1;
  corner1[0]=40.0;
  corner1[1]=0;
  Point<2> corner2;
  corner2[0]=50.0;
  corner2[1]=2.0*M_PI;
  
  const unsigned int dim = 2;
  GridGenerator::hyper_rectangle(tria,corner1,corner2);
  // tria.begin_active ()->face (2)->set_boundary_indicator (1);
  // tria.begin_active ()->face (3)->set_boundary_indicator (1);
  
  tria.refine_global (5);
  std::map< unsigned int, Point< 2 > > new_points;

  //allocate the map
  
  const unsigned int N = 16;

  double a1 = -.01;



  // double nextx(const double x,const double y) const
  // {
  //   return x/(1.0 - a1*std::sin(y));
  // }
 
  for (Triangulation<2>::cell_iterator cell2 = tria.begin();
       cell2 != tria.end(); ++cell2) {
    for (unsigned int f=0; f<GeometryInfo<2>::faces_per_cell; ++f) {
      if (cell2->face(f)->at_boundary() && (cell2->face(f)->center()[0] == 50.0)) 
  	cell2->face(f)->set_boundary_indicator(1);
      if (cell2->face(f)->at_boundary() &&  cell2->face(f)->center()[1] == 2.0*M_PI) 
  	cell2->face(f)->set_boundary_indicator(4);
      if (cell2->face(f)->center()[1] == 0.0) 
  	cell2->face(f)->set_boundary_indicator(5);
      if ((cell2->face(f)->center()[0] == 40.0)) 
  	cell2->face(f)->set_boundary_indicator(2);
      }}

  // //catch the corners
  // cell = tria.begin_active();
  // endc = tria.end(); 
  // for (Triangulation<2>::cell_iterator cell2 = tria.begin();
  //      cell2 != tria.end(); ++cell2) {
  //   for (unsigned int v=0; v<GeometryInfo<2>::vertices_per_cell; ++v) {
  //     if ((cell2->[0] == 40.0)) 
      
  //   }
  // }
  
  
  
  Triangulation<2>::active_cell_iterator
    cell = tria.begin_active(),
    endc = tria.end(); 
  //if (cell->face(f)->at_boundary() && cell->face(f)-> boundary_indicator() ==1)  

  double ymin_r =100;
  double ymax_r = -100;
  double ymin_l =100.0;
  double ymax_l = -100;
  double xmin_t = 100;
  double x0,y0;
  for ( ; cell != endc; ++cell){
    if (cell->at_boundary()){
	for (unsigned int f = 0; f < GeometryInfo<2>::faces_per_cell;++f)	{
	  for (unsigned int v=0; v < GeometryInfo<2>::vertices_per_face;++v) {
	    bool right = cell->face(f)-> boundary_indicator() ==1;
	    bool left = cell->face(f)-> boundary_indicator() ==2;
	    bool top = cell->face(f)-> boundary_indicator() ==4;
	    bool bottom = cell->face(f)-> boundary_indicator() ==5;
	    
	    x0 = cell->face(f)->vertex(v)[0];
	    y0 = cell->face(f)->vertex(v)[1];
	    

	    if (cell->face(f)->at_boundary() && right)  {
	      unsigned int v_n = cell->face(f)->vertex_index(v);
	      new_points[v_n] = go_back(Point<2> (x0,y0));
	      
	      
	      
	      if (new_points[v_n](1)< ymin_r)
		ymin_r = new_points[v_n](1);
	      if (new_points[v_n](1)> ymax_r)
		ymax_r = new_points[v_n](1);
	
	      
	    }
	    if (top)  {
	      Point<2> temp = go_back(Point<2> (x0,y0));				       
	      if (temp(0)< xmin_t)
	    	xmin_t = temp(0);
	    }
	    

}}}}  
  
  std::cout<<"xmin_t: "<<xmin_t<<std::endl;
  cell = tria.begin_active();
  endc = tria.end(); 
  
  for ( ; cell != endc; ++cell){
    if (cell->at_boundary()){
      for (unsigned int f = 0; f < GeometryInfo<2>::faces_per_cell;++f)	{
	for (unsigned int v=0; v < GeometryInfo<2>::vertices_per_face;++v) {
	  bool right = cell->face(f)-> boundary_indicator() ==1;
	  bool left = cell->face(f)-> boundary_indicator() ==2;
	  bool top = cell->face(f)-> boundary_indicator() ==4;
	  bool bottom = cell->face(f)-> boundary_indicator() ==5;
	  
	  unsigned int v_n = cell->face(f)->vertex_index(v);
	  x0 = cell->face(f)->vertex(v)[0];
	  y0 = cell->face(f)->vertex(v)[1];
	  
	  // if (x0 == 40.0 && (y0 == 0 || y0 == 2.0*M_PI)) // 
	  //   new_points[v_n] = Point<2> (x0+10.0,y0+ymin_r);
	  
	  if (cell->face(f)->at_boundary())  {
	    
	    //std::cout<< v_n<<std::endl;
	    
	    
	    if (left)
	      new_points[v_n] = Point<2> (x0,y0+ymin_r);

	    if (top){
	      new_points[v_n] = go_back(Point<2> (x0,y0)); 
	      new_points[v_n] = Point<2> (new_points[v_n](0)+(40.0-xmin_t)*std::exp(-(x0-40.0)/5.0),ymax_r+0.0);
	      //new_points[v_n] = Point<2> (x0,ymax_r);000.0*std::exp(-x0/4.0)
	      ///new_points[v_n] = Point<2> (new_points[v_n](0)+0.0*std::exp(-(x0-40.0)/5.0),ymax_r);
	    }
	    if (bottom){
	      new_points[v_n] = go_back(Point<2> (x0,y0)); 
	      new_points[v_n] = Point<2> (new_points[v_n](0)+(40.0-xmin_t)*std::exp(-(x0-40.0)/5.0),ymin_r);
	      
	    }
	    
	     
	  }}}}
  }  
  
  // cell = tria.begin_active();
  // endc = tria.end(); 
  // for ( ; cell != endc; ++cell){
  //   for (unsigned int v=0; v < GeometryInfo<2>::vertices_per_cell;++v) {
  //     unsigned int v_n = cell->vertex_index(v);
  //     x0 = cell->vertex(v)[0];
  //     y0 = cell->vertex(v)[1];
	  
  //     if (x0 == 40.0 && (y0 == 0 || y0 == 2.0*M_PI))
  // 	new_points[v_n] = Point<2> (x0,y0+ymin_r);
      
  //   }
  // }

  GridTools::transform(&stretch, tria);
  
  GridTools::laplace_transform	(new_points,tria);
  
  std::ofstream out ("grid-lap.eps");
  GridOut grid_out;
  grid_out.write_eps (tria, out);

  // std::ofstream out2 ("transformed.eps");
  // GridOut grid_out2;
  // grid_out.write_msh (tria, out2);

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
Point<2> go_back (const Point<2> &in)
{
  double x0 = in(0);
  double y0 = in(1);

  double aa = -.02;
  double q0 = 3.0;
  double mu = 2.0;
  double a= 40.0;
  double m = 3;
  double b = 50.0;
  double l = 10.0;
  double eps = .1;
  double R_0 = 90.0;
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

  return Point<2>(x_old2,y_old2);
}

int main ()
{

  laplace_grid ();
}
