/* This is the mapping class between finite element mass matrices */
template <int dim>
void arcOn<dim>::matrixmapper()
{

  std::vector< FEValues<dim>* > hp_fe_values;
  std::vector< FEFaceValues<dim>* > hp_fe_values_face;

  for (unsigned int component=0; component< alphadim; ++component){

    UpdateFlags updateflags=  update_values | update_gradients | update_quadrature_points | update_JxW_values;

    hp_fe_values.push_back (new FEValues<dim>(*(fe_collection[component]), *(quadrature_collection[component]), updateflags));
    hp_fe_values_face.push_back (new FEFaceValues<dim>(*(fe_collection[component]), 
						       *(face_quadrature_collection[component]),  
						       updateflags | update_normal_vectors));

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler[component]->begin_active(), endc = dof_handler[component]->end();
    for( ;
	 cell!=endc;
	 ++cell )
      if ( cell->is_locally_owned()  )
	{

	  hp_fe_values[component]->reinit (cell);
	  unsigned int CO = cell->index();

	  const FEValues<dim> &fe_values = hp_fe_values[component]->get_present_fe_values ();
	  const Quadrature<dim>& quadrature_formula = fe_values.get_quadrature();
	  const FiniteElement<dim>& fe = cell->get_fe();
	  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
	  const unsigned int   n_q_points    = quadrature_formula.size();
	  const std::vector<double>& quad_weights = quadrature_formula.get_weights();
	  const std::vector<double> &JxW = fe_values.get_JxW_values ();
	  unsigned int n_comp = fe.n_components();
     
	  FullMatrix<double> MassMatrix(dofs_per_cell,dofs_per_cell);
	  MassMatrix = 0;
	  for (unsigned int i=0; i<dofs_per_cell; ++i){
	    for (unsigned int j=0; j<dofs_per_cell; ++j){
	      for (unsigned int q=0; q<n_q_points; ++q){
		//MassMatrix(i,j)+= fe_values[*(alpha[component])].value(i,q) * fe_values[*(alpha[component])].value(j,q) * JxW[q];
		MassMatrix(i,j)+= fe_values[*(alpha[component])].value(i,q) * fe_values[*(alpha[component])].value(j,q) * quad_weights[q];
	      }
	    }
	  }
     
	  for (unsigned int i=0; i<dofs_per_cell; ++i){
	    for (unsigned int j=0; j<dofs_per_cell; ++j){
	      for (unsigned int q=0; q<n_q_points; ++q){
		//MassMatrix(i,j)+= fe_values[*(sigma[component])].value(i,q) * fe_values[*(sigma[component])].value(j,q) * JxW[q];
		MassMatrix(i,j)+= fe_values[*(sigma[component])].value(i,q) * fe_values[*(sigma[component])].value(j,q) * quad_weights[q];
	      }
	    }
	  }

	  //pcout << "Cell-wise mass matrix" << std::endl;
	  //MassMatrix.print_formatted(std::cout);
	  mapinfo[component][0].localMassMatrix = FullMatrix<double>(dofs_per_cell,dofs_per_cell);
	  mapinfo[component][0].localMassMatrix = MassMatrix;
	  mapinfo[component][0].localInverseMassMatrix = FullMatrix<double>(dofs_per_cell,dofs_per_cell);
	  mapinfo[component][0].localInverseMassMatrix.invert(mapinfo[component][0].localMassMatrix);
	  MassMatrix.diagadd(-1);
	  double TestIdentity_norm = MassMatrix.frobenius_norm();
	  //pcout << "norm = " << TestIdentity_norm << std::endl;
	  if( TestIdentity_norm < MassMatrix.m()*1e-14 ){
	    mapinfo[component][0].MapProject=false;
	  } else {
	    mapinfo[component][0].MapProject=true;
	  }
	  
	}
  }
  
}
