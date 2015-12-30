template <int dim>
void arcOn<dim>::ullman_mesh()
{
  double rmin = 58.0;
  double a= 60.0;
  double b = 68.0;
  double rmax = 78.0;

  /* Triangulation< 2 >	tria; */
  /* Triangulation< 2 >	tria_open; */
  /* Triangulation< 2 >	tria_all; */
  /* Triangulation< 2 >	tria_win; */

  Point<2> corner1;
  corner1[0]=rmin;
  corner1[1]=0;
  Point<2> corner2;
  corner2[0]=rmax;
  corner2[1]=2.0*M_PI;
  
  GridGenerator::hyper_rectangle(triangulation,corner1,corner2);
  triangulation.refine_global (5);
 
  std::map< unsigned int, Point< 2 > > new_points;
  
  const unsigned int N = 16;

  double a1 = -.01;

  Triangulation<2>::active_cell_iterator
    cell= triangulation.begin_active(),
    endc = triangulation.end();

  double ymin_r =100.0;
  double ymax_r = -100.0;
  double ymin_l =100.0;
  double ymax_l = -100.0;
  double xmin_t = 100.0;
  double xmin_b = 100.0;
  double x0,y0;
  
  GridTools::transform(&arcOn<dim>::ullman_return_map, triangulation);
  //tria.refine_global (3);
  
  std::ofstream out ("grid-lap.eps");
  GridOut grid_out;
  grid_out.write_eps (tria_win, out);

}
