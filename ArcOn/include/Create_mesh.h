/* Make a mesh */
template <int dim>
void arcOn<dim>::create_mesh()
{

  std::vector<unsigned int> subdivisions (dim, 0);

  subdivisions[0] = 1; 
  subdivisions[1] = 1;

  double epspi = 0.0;
  double pi = 3.1415926535897932384626433832795;

  const Point<dim> lb = Point<dim>(0.0,0.0);//(5.0*38.0,0.0);
  const Point<dim> rt = Point<dim>(100.0,20*pi);//(5.0*98.0,2.0*M_PI*50.0);

  GridGenerator::subdivided_hyper_rectangle(triangulation,subdivisions,lb,rt);

  // Let's set the periodic conditions
  // Can do this using the colorize flag in GridGen above
  for (typename Triangulation<dim>::cell_iterator cell = triangulation.begin();
       cell != triangulation.end(); ++cell) {
    for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f) {
      if (cell->face(f)->at_boundary() && (cell->face(f)->center()[1] == 0.0)) {
  	//set bottom boundary
  	cell->face(f)->set_boundary_indicator(1);
      }
      if (cell->face(f)->at_boundary() && (cell->face(f)->center()[1] == 20*pi)) {

  	//set top boundayr
  	cell->face(f)->set_boundary_indicator(2);
      }
      
      if (cell->face(f)->at_boundary() && (cell->face(f)->center()[0] == 0.0)) {
      	//set left boundary
      	cell->face(f)->set_boundary_indicator(3);
      }
      if (cell->face(f)->at_boundary() && (cell->face(f)->center()[0] == 100.0)) {
      	//set right boundary
      	cell->face(f)->set_boundary_indicator(4);
      }
    }
  }

  //Y-periodicity
  std::vector<GridTools::PeriodicFacePair<typename parallel::distributed::Triangulation<dim>::cell_iterator> >
    periodicity_vector;
  
  GridTools::collect_periodic_faces( 	triangulation,	
					1 /* b_id1 */, 
					2 /* b_id2 */,
					1 /* direction */,	
					periodicity_vector );
  
  triangulation.add_periodicity(periodicity_vector);

  std::vector<GridTools::PeriodicFacePair<typename parallel::distributed::Triangulation<dim>::cell_iterator> >
    periodicity_vector2;

  GridTools::collect_periodic_faces( 	triangulation,
					3, //b_id1
					4, //b_id2
					0, //direction
					periodicity_vector2
					);
  
  triangulation.add_periodicity(periodicity_vector2);

  triangulation.refine_global (refinements);

  for (typename Triangulation<dim>::cell_iterator cell = triangulation.begin();
       cell != triangulation.end(); ++cell) {
    for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f) {

      

    }
  }

}
