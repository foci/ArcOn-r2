/* Let us assemble the stifffness matrix just once, and precondition it just once .... */
template <int dim>
void arcOn<dim>::assemble_cont_stiffness(SolutionVector& subdomain_solution, double delta_t) { 

  (void) delta_t;

  std::vector< FEValues<dim>* > hp_fe_values;
  std::vector< FEFaceValues<dim>* > hp_fe_values_face;
  std::vector< FEFaceValues<dim>* > hp_fe_values_neigh_face;

  double current_time = 1.0;
  for (unsigned int component=0; component< alphadim; ++component){

    UpdateFlags updateflags=  update_values | update_gradients 
      | update_quadrature_points | update_JxW_values;
    
    hp_fe_values.push_back (new FEValues<dim>(*(tfe_collection[component]), 
					      *(tquadrature_collection[component]), updateflags));
    hp_fe_values_face.push_back (new FEFaceValues<dim>(*(tfe_collection[component]), 
						       *(tface_quadrature_collection[component]),  
						       updateflags | update_normal_vectors ));
    hp_fe_values_neigh_face.push_back (new FEFaceValues<dim>(*(tfe_collection[component]), 
							     *(tface_quadrature_collection[component]),  
							     updateflags | update_normal_vectors ));

    if (component == 2){

      cont_poisson_matrix = 0.0;
      typename DoFHandler<dim>::active_cell_iterator 
	cell = tdof_handler[component]->begin_active(), 
	endc = tdof_handler[component]->end();
      for(;
	  cell!=endc;
	  ++cell	)
	if (cell->is_locally_owned()  )
	  {
	    int CO = cell->index();
	    hp_fe_values[component]->reinit (cell);
	    const FEValues<dim> &fe_values = hp_fe_values[component]->get_present_fe_values ();
	    const Quadrature<dim>& quadrature_formula = fe_values.get_quadrature();
	    const FiniteElement<dim>& fe = cell->get_fe();
	    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
	    const unsigned int   n_q_points    = quadrature_formula.size();

	    FullMatrix<double> A_matrix( dofs_per_cell,dofs_per_cell);

	    std::vector<unsigned int> local_dof_indices (dofs_per_cell);

	    A_matrix = 0.0;
	      
	    cell->get_dof_indices (local_dof_indices);

	    const std::vector<double> &JxW = fe_values.get_JxW_values ();
	      
	    for (unsigned int i=0; i<dofs_per_cell; ++i){
	      for (unsigned int j=0; j<dofs_per_cell; ++j){
		for (unsigned int q=0; q<n_q_points; ++q){
		      
		  A_matrix(i,j) += (fe_values[*(alpha[component])].gradient(i,q)) * (fe_values[*(alpha[component])].gradient(j,q)) * JxW[q];

		}
		//		cont_poisson_matrix.add( local_dof_indices[ i ],
		//		    local_dof_indices[ j ],
		//		    A_matrix(i,j));
	      }
	    }
	    
	    telliptic_constraints.distribute_local_to_global(A_matrix,
							     local_dof_indices,
							     cont_poisson_matrix);

	  }

      cont_poisson_matrix.compress (VectorOperation::insert);
      cont_preconditioner.initialize(cont_poisson_matrix.block(0,0),
				PETScWrappers::PreconditionBoomerAMG::AdditionalData(true));
      
    }
  }
  for (unsigned int component=0; component< alphadim; ++component){  
    delete hp_fe_values[component];
    delete hp_fe_values_face[component];
    delete hp_fe_values_neigh_face[component];
  }
}

