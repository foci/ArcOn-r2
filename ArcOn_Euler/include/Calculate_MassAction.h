/* This template function computes all terms with "no derivatives" in the evolution equations.  We might think of these as sources terms (if they do not depend on the primal variables) and reaction terms (if they are functions of the primal variables) .... */

template <int dim>
void arcOn<dim>::Calculate_MassAction_Explicit(SolutionVector& subdomain_solution, double delta_t, double current_time, SolutionVector& mass_action_term)
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
  

    /* We can using ramping to avoid shocking the system */
  double dramp;
  if (Time_ramp > 1.0e-5){
    dramp  = std::tanh( 2.0*( current_time/Time_ramp ) );
  }
  else{
    dramp = 1.0;
  }

  double alph = alpha_parameter*dramp; //1e-4 + dramp*alph_0;
  double volt_b = -bias_parameter ;
  
  for (unsigned int component=0; component< alphadim-1; ++component){

    std::vector<double> transport_alphas(alphadim,0);
    Functionals<dim> functionals;


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

	  FullMatrix<double> mass_action (alphadim, dofs_per_cell);
	
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
	  mass_action = 0.0;

	  double ramp = std::tanh(current_time/100.0);
	
	  if (component == 0){
	    double gam = 1.0; //2e-2;
	    double sqrt2pi = 2.506628275;
	    double sigs = 0.002*290.0;
	    for (unsigned int q=0; q<n_q_points; ++q){
	      for (unsigned int i=0; i<dofs_per_cell; ++i){

		mass_action(component,i) -= ( -7.0*alph*(0.6 + 0.01*quadrature_point[q][0])*std::exp(-std::pow(quadrature_point[q][0]-37.0,2.0)/144.0)/std::exp(syncopated[0][q]) + alph*(0.6 + 0.01*quadrature_point[q][0]))* fe_values[*(alpha[component])].value(i,q)*JxW[q];

	      }
	    }
	  }
	  if (component == 1){
	    fe_values[*(alpha[2])].get_function_values(subdomain_solution[2],syncopated[2]);
	    int q_cent = 0.43764*n_q_points; //Where the bias plate is centered --> center bump function here 
	    double some_constant = 0.5;
	    double phi_b;     //bias voltage function defined by tanh funtion  (corresponds to 2nd plate in Helimak)
	    for (unsigned int q=0; q<n_q_points; ++q){
	      for (unsigned int i=0; i<dofs_per_cell; ++i){

		phi_b = dramp*volt_b*(std::tanh((quadrature_point[q][0] - 20.0)/0.1) + std::tanh((40.0 - quadrature_point[q][0])/0.1))/2;
		mass_action(component,i) += alph * (0.6 + 0.01*quadrature_point[q][0])*fe_values[*(alpha[component])].value(i,q)* (syncopated[2][q] - phi_b ) *JxW[q];

	      }
	    }
	  }
	  
	  cell->get_dof_indices (local_dof_indices);

	  Vector<double> temp_update[alphadim]; 
	  temp_update[component].reinit(dofs_per_cell); 
	  for (unsigned int i=0; i<dofs_per_cell; ++i){
	    temp_update[component](i) = mass_action(component,i);
	  }
	  Vector<double> projected(dofs_per_cell); 
	  mapinfo[component][0].localInverseMassMatrix.vmult(projected,temp_update[component]);
	  projected /= JxWsum;
	  parahyp_constraints[component].distribute_local_to_global (projected,
								     local_dof_indices,
								     mass_action_term[component]);	  
	  
	}

    mass_action_term[component].compress(VectorOperation::insert);

  }

  for (unsigned int component=0; component< alphadim; ++component){
    delete hp_fe_values[component];
    delete hp_fe_values_face[component];
  }
}
