/* This is the base class template */
template <int dim>
void arcOn<dim>::run()
{
  /* //GridGenerator::subdivided_hyper_cube(triangulation,41,0,140); */

  /* std::vector<unsigned int> subdivisions (dim, 0); */
  /* //subdivisions[0] = 222; */
  /* //subdivisions[1] = 111; */

  /* //subdivisions[0] = 231; */
  /* //subdivisions[1] = 153; */

  /* //base guy */
  /* //subdivisions[0] = 77; */
  /* //subdivisions[1] = 77; */

  /* //subdivisions[0] = 15; */
  /* //subdivisions[1] = 15; */

  /* subdivisions[0] = 1;  */
  /* subdivisions[1] = 1; */

  /* double epspi = 0.0; */
  /* double pi = 1.0 + epspi ; //2.0*3.1415926535897932384626433832795; */

  /* const Point<dim> lb = Point<dim>(5.0*38.0,0.0); */
  /* const Point<dim> rt = Point<dim>(5.0*98.0,2.0*M_PI*50.0); */

  /* //const Point<dim> lb = Point<dim>(epspi,epspi); */
  /* //const Point<dim> rt = Point<dim>(pi,pi); */

  /* //const Point<dim> lb = Point<dim>(-64.0,-192); */
  /* //const Point<dim> rt = Point<dim>(64.0,192); */

  /* //const Point<dim> lb = Point<dim>(-30.0,-20.0); */
  /* //const Point<dim> rt = Point<dim>(30.0,20.0); */

  /* /\* const Point<dim> lb = Point<dim>(-5,-10); *\/ */
  /* /\* const Point<dim> rt = Point<dim>(20,10); *\/ */
    /* //GridGenerator::hyper_rectangle(triangulation,lb,rt); */

  /* //Good one--------- */
  /* /\* GridIn<dim> grid_in; *\/ */
  /* /\* grid_in.attach_triangulation (triangulation); *\/ */
  /* /\* std::ifstream input_file("../../input/Ullman_mesh.inp"); *\/ */
  /* /\* grid_in.read_ucd (input_file); *\/ */
  /* //Good one--------- */
  
  /* //Good one */
  /* GridGenerator::subdivided_hyper_rectangle(triangulation,subdivisions,lb,rt); */

 
  /* // Let's set the periodic conditions */
  /* // Can do this using the colorize flag in GridGen above */
  /* for (typename Triangulation<dim>::cell_iterator cell = triangulation.begin(); */
  /*      cell != triangulation.end(); ++cell) { */
  /*   for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f) { */
  /*     if (cell->face(f)->at_boundary() && (cell->face(f)->center()[1] == 0.0)) { */
  /* 	//set bottom boundary */
  /* 	//pcout << "Left boundary cell index = " << cell->index() << std::endl; */
  /* 	//pcout << "Left boundary face index = " << f << std::endl; */
  /* 	cell->face(f)->set_boundary_indicator(1); */
  /*     } */
  /*     if (cell->face(f)->at_boundary() && (cell->face(f)->center()[1] == 2.0*M_PI*50.0)) { */
  /* 	//set top boundayr */
  /* 	//pcout << "Right boundary cell index = " << cell->index() << std::endl; */
  /* 	//pcout << "Right boundary face index = " << f << std::endl; */
  /* 	cell->face(f)->set_boundary_indicator(2); */
  /*     } */
  /*     /\* if (cell->face(f)->at_boundary() && (cell->face(f)->center()[0] == 290.0)) { *\/ */
  /*     /\* 	//set top boundary *\/ */
  /*     /\* 	//pcout << "Right boundary cell index = " << cell->index() << std::endl; *\/ */
  /*     /\* 	//pcout << "Right boundary face index = " << f << std::endl; *\/ */
  /*     /\* 	cell->face(f)->set_boundary_indicator(15); *\/ */
  /*     /\* } *\/ */
      
      
  /*     if (cell->face(f)->at_boundary() && (cell->face(f)->center()[0] == 5.0*38.0)) { */
  /*     	//set left boundary */
  /*     	//pcout << "Left boundary cell index = " << cell->index() << std::endl; */
  /*     	//pcout << "Left boundary face index = " << f << std::endl; */
  /*     	cell->face(f)->set_boundary_indicator(3); */
  /*     } */
  /*     if (cell->face(f)->at_boundary() && (cell->face(f)->center()[0] == 5.0*98.0)) { */
  /*     	//set right boundary */
  /*     	//pcout << "Right boundary cell index = " << cell->index() << std::endl; */
  /*     	//pcout << "Right boundary face index = " << f << std::endl; */
  /*     	cell->face(f)->set_boundary_indicator(4); */
  /*     } */

  create_mesh();

  setup_system();

  matrixmapper();

  create_dg_periodicity();
 
  InitialValues<dim> IV;

  /* for (unsigned int component=0; component< alphadim; ++component){ */
  /*   subdomain_solution[component].compress();} */

  assemble_system();

}
