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
  
  for (unsigned int component=0; component< alphadim-1; ++component){

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

	  //double rseed = Utilities::generate_normal_random_number(0.0,0.25); 	
	
	  if (component == 0){
	    double alph = 1e-4; //5e-4;
	    double gam = 1.0; //2e-2;
	    double sqrt2pi = 2.506628275;
	    double sigs = 0.002*290.0;
	    for (unsigned int q=0; q<n_q_points; ++q){
	      for (unsigned int i=0; i<dofs_per_cell; ++i){
		//mass_action(component,i) += alph*(fe_values[*(alpha[component])].value(pinfo[component].alpha_dof_index[i],q) * (syncopated[0][q]-syncopated[2][q]) * JxW[q]);
		/* if ( quadrature_point[q][0] <= 0.0 ){ */
		/*   mass_action(component,i) += 0.001*std::abs(std::sin(0.1*quadrature_point[q][1]))*fe_values[*(alpha[component])].value(i,q)*JxW[q]; */
		/* } */

		//if(syncopated[component][q]>0.0){

		//good one
		mass_action(component,i) += ( 3.0*alph*std::exp(-std::pow(quadrature_point[q][0]-75.0,2.0)/200.0)/std::exp(syncopated[0][q]) - alph ) * fe_values[*(alpha[component])].value(i,q)*JxW[q] ;

		//mass_action(component,i) += ( 3.0*alph*(1.0+std::abs(std::sin(0.08*quadrature_point[q][1]))*std::exp(-0.5*std::pow(quadrature_point[q][0],2.0)/std::pow(sigs,2.0))) - alph ) * fe_values[*(alpha[component])].value(i,q)*JxW[q] ;
		  /* if ( quadrature_point[q][0] > 200.0 ){ */
		    
		  /*   mass_action(component,i) -= 2.0 * alph * fe_values[*(alpha[component])].value(i,q)*JxW[q] ; */

		  /* } */

		  //&& syncopated[0][q] > 2.0*glob_min_ribbon_density
		//	else{
		//mass_action(component,i) += 0.0;}

		//mass_action(component,i) += difs(0) * (fe_values[*(alpha[component])].gradient(i,q)) * (fe_values[*(alpha[component])].gradient(i,q))  * JxW[q];

		  //mass_action(component,i) += std::exp( std::abs(quadrature_point[q][1])/(2.0*pi) ) * ( std::sin(quadrature_point[q][1]*pi) / (quadrature_point[q][1]*pi) * sin(quadrature_point[q][0]*pi) / (quadrature_point[q][0]*pi) )  *   (fe_values[*(alpha[component])].value(pinfo[component].alpha_dof_index[i],q) * JxW[q]);        

		  //mass_action(component,i) += dramp *(0.00101 + 0.001 *( std::cos(0.5*pi*quadrature_point[q][1]) ) ) 
		  //*   (fe_values[*(alpha[component])].value(i,q) * JxW[q]);

		  //mass_action(component,i) += dramp *(0.00095 + 0.001 *( std::cos(0.157*quadrature_point[q][1])+std::sin( pi*quadrature_point[q][1] ) ) ) 
		  // *   (fe_values[*(alpha[component])].value(i,q) * JxW[q]);    
		  
		  /* double xpt = quadrature_point[q][0]; */
		  /* double ypt = quadrature_point[q][1]; */

		  /* double ga = 0.1; */
		  /* double gad = 0.001; */
		  /* double gac = 0.001; */
		  

		  /* /\* mass_action(component,i) += ( -0.5*(std::pow(xpt,2.0)+std::pow(ypt,2.0))*std::exp(-current_time*( (std::pow(xpt,2.0)+std::pow(ypt,2.0) )/2.0 ))  *\/ */
		  /* /\*   + 0.001*xpt*current_time*std::sin(xpt)*std::cos(ypt)*std::exp(-current_time*( (std::pow(xpt,2.0)+std::pow(ypt,2.0) )/2.0 ))   *\/ */
		  /* /\*   - 0.001*ypt*current_time*std::cos(xpt)*std::sin(ypt)*std::exp(-current_time*( (std::pow(xpt,2.0)+std::pow(ypt,2.0) )/2.0 ))   *\/ */
		  /* /\*   + (current_time-std::pow(current_time,2.0)*std::pow(xpt,2.0))*std::exp(-current_time*( (std::pow(xpt,2.0)+std::pow(ypt,2.0) )/2.0 )) *\/ */
		  /* /\* 				*   (fe_values[*(alpha[component])].value(i,q) * JxW[q]) ); *\/ */

		  /* //pcout << "current_time = " << current_time << std::endl; */

		  /* mass_action(component,i) += ( 1.0 + (current_time)*xpt*std::cos(xpt)*std::sin(ypt)*ga */
		  /* 				-(current_time)*ypt*std::sin(xpt)*std::cos(ypt)*ga  */
		  /* 				+ 2.0*gad*current_time - gad*current_time*std::pow(xpt,2.0) - gad * current_time * std::pow(ypt,2.0)  */
		  /* 				- ga*xpt*std::pow(std::cosh(xpt),-1.0)*std::pow(std::cosh(ypt),-1.0)*std::tanh(ypt)*std::exp(std::pow(std::cosh(xpt),-1.0) * std::pow(std::cosh(ypt),-1.0)) )  */
		  /*   //+ 0.02*current_time - */
		  /*   //0.01*(current_time+1.0)*std::pow(xpt,2.0)-0.01*(current_time+1.0)*std::pow(ypt,2.0)) */
		  /*   *std::exp(-(std::pow(xpt,2.0)+std::pow(ypt,2.0))/2.0 ) */
		  /*   *(fe_values[*(alpha[component])].value(i,q) * JxW[q]) ; */
		  
 
		  //} 
	      }
	    }
	  }
	  //for speed
	  /* if (component == 1){ */
	  /*   fe_values[*(alpha[2])].get_function_values(subdomain_solution[2],syncopated[2]); */
	  /*   double alph = 1e-4; //5e-4;  */
	  /*   for (unsigned int q=0; q<n_q_points; ++q){ */
	  /*     for (unsigned int i=0; i<dofs_per_cell; ++i){ */
	  /* 	// mass_action(component,i) += alph*(fe_values[*(alpha[component])].value(pinfo[component].alpha_dof_index[i],q)  */

	  /* 	//mass_action(component,i) += alph*(fe_values[*(alpha[component])].value(i,q)           */
	  /* 	//				  *  (syncopated[2][q]) * JxW[q]); */
		
	  /* 	  double xpt = quadrature_point[q][0]; */
	  /* 	  double ypt = quadrature_point[q][1]; */

	  /* 	  double ga = 0.1; */
	  /* 	  double gad = 0.001; */

	  /* 	  /\* mass_action(component,i) += ( -0.5*(std::pow(xpt,2.0)+std::pow(ypt,2.0))*std::exp(-current_time*( (std::pow(xpt,2.0)+std::pow(ypt,2.0) )/2.0 ))  *\/ */
	  /* 	  /\*   + 0.001*xpt*current_time*std::sin(xpt)*std::cos(ypt)*std::exp(-current_time*( (std::pow(xpt,2.0)+std::pow(ypt,2.0) )/2.0 ))   *\/ */
	  /* 	  /\*   - 0.001*ypt*current_time*std::cos(xpt)*std::sin(ypt)*std::exp(-current_time*( (std::pow(xpt,2.0)+std::pow(ypt,2.0) )/2.0 ))   *\/ */
	  /* 	  /\*   + (current_time-std::pow(current_time,2.0)*std::pow(xpt,2.0))*std::exp(-current_time*( (std::pow(xpt,2.0)+std::pow(ypt,2.0) )/2.0 )) *\/ */
	  /* 	  /\* 				*   (fe_values[*(alpha[component])].value(i,q) * JxW[q]) ); *\/ */

	  /* 	  //pcout << "current_time = " << current_time << std::endl; */

	  /* 	  /\* mass_action(component,i) += ( ( 1.0 + (current_time)*xpt*std::cos(xpt)*std::sin(ypt)*ga *\/ */
	  /* 	  /\* 				-(current_time)*ypt*std::sin(xpt)*std::cos(ypt)*ga  *\/ */
	  /* 	  /\* 				+ 2.0*gad*current_time - gad*current_time*std::pow(xpt,2.0) - gad * current_time * std::pow(ypt,2.0) ) *\/ */
						
	  /* 	  /\*   //+ 0.02*current_time - *\/ */
	  /* 	  /\*   //0.01*(current_time+1.0)*std::pow(xpt,2.0)-0.01*(current_time+1.0)*std::pow(ypt,2.0)) *\/ */
	  /* 	  /\* 				*std::exp(-(std::pow(xpt,2.0)+std::pow(ypt,2.0))/2.0 )  *\/ */
	  /* 	  /\* 				- ga*std::pow(std::cosh(xpt),-1.0)*std::pow(std::cosh(ypt),-1.0)*std::tanh(ypt)*std::exp(std::pow(std::cosh(xpt),-1.0) * std::pow(std::cosh(ypt),-1.0)) ) *\/ */
	  /* 	  /\*   *(fe_values[*(alpha[component])].value(i,q) * JxW[q]) ; *\/ */
		  
 
	  /*     } */
	  /*   } */
	  /* } */

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
