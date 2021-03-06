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
  	cell->face(f)->set_boundary_indicator(1);
      }
      if (cell->face(f)->at_boundary() && (cell->face(f)->center()[1] == 2.0*M_PI*50.0)) {
  	//set top boundayr
  	cell->face(f)->set_boundary_indicator(2);
      }
      
      if (cell->face(f)->at_boundary() && (cell->face(f)->center()[0] == 5.0*38.0)) {
      	//set left boundary
      	cell->face(f)->set_boundary_indicator(3);
      }
      if (cell->face(f)->at_boundary() && (cell->face(f)->center()[0] == 5.0*98.0)) {
      	//set right boundary
      	cell->face(f)->set_boundary_indicator(4);
      }
    }
  }  
}
