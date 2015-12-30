/*Recreate the BC data*/
template <int dim>
void arcOn<dim>::recreate_boundary_data()
{

  // Let's set the periodic conditions
  // Can do this using the colorize flag in GridGen above
  for (typename Triangulation<dim>::cell_iterator cell = triangulation.begin();
       cell != triangulation.end(); ++cell) {
    for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f) {
      if (cell->face(f)->at_boundary() && (cell->face(f)->center()[1] == 0.0)) {
  	//set bottom boundary
  	//pcout << "Left boundary cell index = " << cell->index() << std::endl;
  	//pcout << "Left boundary face index = " << f << std::endl;
  	cell->face(f)->set_boundary_indicator(1);
      }
      if (cell->face(f)->at_boundary() && (cell->face(f)->center()[1] == 2.0*M_PI*50.0)) {
  	//set top boundayr
  	//pcout << "Right boundary cell index = " << cell->index() << std::endl;
  	//pcout << "Right boundary face index = " << f << std::endl;
  	cell->face(f)->set_boundary_indicator(2);
      }
      /* if (cell->face(f)->at_boundary() && (cell->face(f)->center()[0] == 290.0)) { */
      /* 	//set top boundary */
      /* 	//pcout << "Right boundary cell index = " << cell->index() << std::endl; */
      /* 	//pcout << "Right boundary face index = " << f << std::endl; */
      /* 	cell->face(f)->set_boundary_indicator(15); */
      /* } */
      
      
      if (cell->face(f)->at_boundary() && (cell->face(f)->center()[0] == 5.0*38.0)) {
      	//set left boundary
      	//pcout << "Left boundary cell index = " << cell->index() << std::endl;
      	//pcout << "Left boundary face index = " << f << std::endl;
      	cell->face(f)->set_boundary_indicator(3);
      }
      if (cell->face(f)->at_boundary() && (cell->face(f)->center()[0] == 5.0*98.0)) {
      	//set right boundary
      	//pcout << "Right boundary cell index = " << cell->index() << std::endl;
      	//pcout << "Right boundary face index = " << f << std::endl;
      	cell->face(f)->set_boundary_indicator(4);
      }


    }
  }
  
  //Y-periodicity
  /* std::vector<GridTools::PeriodicFacePair<typename parallel::distributed::Triangulation<dim>::cell_iterator> > */
  /*   periodicity_vector; */
  
  /* GridTools::collect_periodic_faces */
  /*   ( 	triangulation,	1 /\* b_id1 *\/, 2 /\* b_id2 *\/,1 /\* direction *\/,	periodicity_vector ); */
  
  /* triangulation.add_periodicity(periodicity_vector); */
  /* //X-periodicity */
  /* //std::vector<GridTools::PeriodicFacePair<typename parallel::distributed::Triangulation<dim>::cell_iterator> > */
  /* //  periodicity_vector2 = */
  /* std::vector<GridTools::PeriodicFacePair<typename parallel::distributed::Triangulation<dim>::cell_iterator> > */
  /*   periodicity_vector2; */

  /* GridTools::collect_periodic_faces 	( 	triangulation, */
  /* 						3, //b_id1 */
  /* 						4, //b_id2 */
  /* 						0, //direction */
  /* 						periodicity_vector2 */
  /* 						); */

  /* triangulation.add_periodicity(periodicity_vector2); */

}
