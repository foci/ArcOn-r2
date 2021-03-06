
template <int dim>
void arcOn<dim>::revert_density(SolutionVector& subdomain_solution, double delta_t, double current_time, SolutionVector& naive_revert_output)
{

  std::vector< FEValues<dim>* > hp_fe_values;
  std::vector< FEFaceValues<dim>* > hp_fe_values_face;
  
  for (unsigned int component=0; component< alphadim; ++component){
    UpdateFlags updateflags=  update_values | update_gradients | 
      update_quadrature_points | update_JxW_values;
    hp_fe_values.push_back (new FEValues<dim>(*(fe_collection[component]), 
					      *(quadrature_collection[component]), updateflags));
    hp_fe_values_face.push_back (new FEFaceValues<dim>(*(fe_collection[component]), 
						       *(face_quadrature_collection[component]),  
						       updateflags | update_normal_vectors ));
  }
  
  for (unsigned int component=0; component< alphadim-2; ++component){

    std::vector<double> transport_alphas(alphadim,0);
    Functionals<dim> functionals;

    /* We can using ramping to avoid shocking the system */
    double ramp = 0.5; // ramp time in seconds
    double dramp  = 1.0; //std::tanh( 2.0*( current_time/ramp ) ); 
    
    typename DoFHandler<dim>::active_cell_iterator 
      cell = dof_handler[component]->begin_active(), 
      endc = dof_handler[component]->end();
    for(;
	cell!=endc;
	++cell	)
      if ( cell->is_locally_owned() )
	{
	  int CO = cell->index();
	  double pi = 3.1415926535897932384626433;

	  hp_fe_values[component]->reinit (cell);
	  const FEValues<dim> &fe_values = hp_fe_values[component]->get_present_fe_values ();
	  const Quadrature<dim>& quadrature_formula = fe_values.get_quadrature(); 
	  const std::vector<Point<dim> >& quadrature_point = fe_values.get_quadrature_points();
	  const FiniteElement<dim>& fe = cell->get_fe();
	  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
	  const unsigned int   n_q_points    = quadrature_formula.size();

	  FullMatrix<double> n_from_chi (alphadim, dofs_per_cell);
	  FullMatrix<double> exp_MassMatrix (dofs_per_cell, dofs_per_cell);	

	  std::vector<unsigned int> local_dof_indices (dofs_per_cell);

	  for(unsigned int k=0;k<alphadim;k++){
	    syncopated[k] =  std::vector<double>(n_q_points);
	  }

	  //cell->get_dof_indices (local_dof_indices);
	  const std::vector<double> &JxW = fe_values.get_JxW_values ();
	  double JxWsum = 0;
	  for (unsigned int q=0; q<n_q_points; ++q){
	    JxWsum += JxW[q];
	  }

	  fe_values[*(alpha[component])].get_function_values(subdomain_solution[component],
							     syncopated[component]);

	  Tensor<1,alphadim> syncopated_alphas;
	  n_from_chi = 0.0;
	  exp_MassMatrix = 0.0;

	
	  if (component == 0){
	    double alph = 1e-4; //5e-4;
	    double gam = 1.0; //2e-2;
	    double sqrt2pi = 2.506628275;
	    double sigs = 0.002*290.0;
	    for (unsigned int q=0; q<n_q_points; ++q){
	      for (unsigned int i=0; i<dofs_per_cell; ++i){

		n_from_chi(component,i) += std::exp(syncopated[0][q]) * std::exp(fe_values[*(alpha[component])].value(i,q)*JxW[q]) ;


	      }
	    }

	    for (unsigned int q=0; q<n_q_points; ++q){
              for (unsigned int i=0; i<dofs_per_cell; ++i){
		for (unsigned int j=0; j<dofs_per_cell; ++j){
		
		  exp_MassMatrix(i,j) += std::exp(fe_values[*(alpha[component])].value(i,q))*std::exp(fe_values[*(alpha[component])].value(j,q))*JxW[q] ;

		}
              }
            }

	  }

	  FullMatrix<double> localInverseExpMassMatrix(dofs_per_cell,dofs_per_cell);
	  localInverseExpMassMatrix.invert(exp_MassMatrix);
	  
	  cell->get_dof_indices (local_dof_indices);

	  Vector<double> temp_update[alphadim]; 
	  temp_update[component].reinit(dofs_per_cell); 
	  for (unsigned int i=0; i<dofs_per_cell; ++i){
	    temp_update[component](i) = n_from_chi(component,i);
	  }
	  Vector<double> projected(dofs_per_cell); 

	  localInverseExpMassMatrix.vmult(projected,temp_update[component]);
	  projected /= JxWsum;
	  parahyp_constraints[component].distribute_local_to_global (projected,
								     local_dof_indices,
								     naive_revert_output[component]);	  
	  
	}

	  naive_revert_output[component].compress(VectorOperation::add);

  }

  for (unsigned int component=0; component< alphadim; ++component){
    delete hp_fe_values[component];
    delete hp_fe_values_face[component];
  }
}
