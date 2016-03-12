/* The template function solves the Poisson equation using a global linear solve .... */
template <int dim>
void arcOn<dim>::calc_poisson(SolutionVector& subdomain_solution, double delta_t) { 

  (void) delta_t;

  std::vector< FEValues<dim>* > hp_fe_values;
  std::vector< FEFaceValues<dim>* > hp_fe_values_face;
  double current_time = 1.0;

  for (unsigned int component=0; component< alphadim; ++component){
    UpdateFlags updateflags=  update_values | update_gradients 
      | update_quadrature_points | update_JxW_values;
      
    hp_fe_values.push_back (new FEValues<dim>(*(fe_collection[component]),
					      *(quadrature_collection[component]), updateflags));
    hp_fe_values_face.push_back (new FEFaceValues<dim>(*(fe_collection[component]), 
						       *(face_quadrature_collection[component]),  
						       updateflags | update_normal_vectors ));
      
    if (component == 2){

      poisson_rhs.block(0) = 0.0;
      typename DoFHandler<dim>::active_cell_iterator 
	cell = dof_handler[component]->begin_active(), 
	endc = dof_handler[component]->end();
	
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
	      
	    Vector<double>       cell_rhs (dofs_per_cell);
	    Vector<double>       cell_rhs_boundary (dofs_per_cell);
	    Vector<double>       cell_rhs_composite (dofs_per_cell);
	      
	    std::vector<unsigned int> local_dof_indices (dofs_per_cell);
	    std::vector<unsigned int> neigh_local_dof_indices (dofs_per_cell);

	    cell_rhs_composite = 0.0;
	    cell_rhs = 0.0;
	      
	    for(unsigned int k=0;k<alphadim;k++){
	      prev_soln_alpha[k] =  std::vector<double>(n_q_points);
	      prev_soln_sigma[k] =  std::vector<Tensor<1,dim> >(n_q_points);
	    }
	      
	    cell->get_dof_indices (local_dof_indices);

	    const std::vector<double> &JxW = fe_values.get_JxW_values ();
	      
	    fe_values[*(alpha[1])].get_function_values(subdomain_solution[1],prev_soln_alpha[1]);
	    fe_values[*(sigma[1])].get_function_values(subdomain_solution[1],prev_soln_sigma[1]);

	    std::vector<unsigned int> index_container2 (dofs_per_cell);

	    for (unsigned int i=0; i<dofs_per_cell; ++i){
	      for (unsigned int q=0; q<n_q_points; ++q){
		
		cell_rhs(i) -= (fe_values[*(alpha[component])].value(i,q))* (prev_soln_alpha[1][q]) * JxW[q] ;
		
	      }
	    }

	    //This is for nonperiodic stuff
	    cell_rhs_boundary = 0.0;
	    for(unsigned int face_num=0; face_num < GeometryInfo<dim>::faces_per_cell; ++face_num){
	      typename DoFHandler<dim>::face_iterator face=cell->face(face_num);
	      
	      if ( face->at_boundary()  && face->boundary_indicator() == 0 )
	    	{
	    	  hp_fe_values_face[component]->reinit (cell,face_num);
	    	  const FEFaceValues<dim>& fe_values_face = hp_fe_values_face[component]->get_present_fe_values ();
	    	  const Quadrature<dim-1>& face_quadrature_formula = fe_values_face.get_quadrature();
	    	  const unsigned int n_q_points_face = face_quadrature_formula.size();
		    
	    	  const std::vector<double> &JxW_face = fe_values_face.get_JxW_values ();
	    	  const std::vector<Point<dim> > &normals = fe_values_face.get_normal_vectors ();
		    
	    	  for(unsigned int k=0;k<alphadim;k++){
	    	    prev_soln_alpha_face[k] =  std::vector<double>(n_q_points_face);
	    	    prev_soln_sigma_face[k] =  std::vector<Tensor<1,dim> >(n_q_points_face);
	    	  }

	    	  fe_values_face[*(alpha[component])].get_function_values(subdomain_solution[component],
	    	  							  prev_soln_alpha_face[component]);
	    	  fe_values_face[*(sigma[component])].get_function_values(subdomain_solution[component],
	    	  							  prev_soln_sigma_face[component]);
	    	  fe_values_face[*(alpha[1])].get_function_values(subdomain_solution[1],
	    							  prev_soln_alpha_face[1]);
	    	  fe_values_face[*(sigma[1])].get_function_values(subdomain_solution[1],
	    							  prev_soln_sigma_face[1]);

	    	  std::vector<double> alphas_boundary(n_q_points_face);
		  std::vector<double> Salphas_boundary(n_q_points_face);
	    	  unsigned char boundary_index = face->boundary_indicator();
	    	  const WallBoundaryValues<dim>& wbv = WBV[boundary_index];
	    	  const std::vector< Point<dim> > &quadrature_points = fe_values_face.get_quadrature_points();

	    	  alphas_boundary = prev_soln_alpha_face[component];
	    	  wbv.value_list( quadrature_points, alphas_boundary, component, current_time);

		  const unsigned int normal1 = GeometryInfo<dim>::unit_normal_direction[face_num];
		  
		  double pen1 = cell->extent_in_direction(normal1);
		  
		  double sigma_p1 = degree*( degree + 1.0 ) / pen1;
		  double sigma_p2 = degree*( degree + 1.0 ) / pen1;
		  
		  double sigma_penalty = sigma_prefactor*0.5*(sigma_p1+sigma_p2);

	    	  for (unsigned int i=0; i<dofs_per_cell; ++i){
	    	    for (unsigned int q=0; q<n_q_points_face; ++q){

		      cell_rhs_boundary(i) +=  ( - sigma_penalty * ( alphas_boundary[q]- prev_soln_alpha_face[component][q] )
						 * fe_values_face[*(alpha[component])].value(i,q) ) * JxW_face[q];

	    	    }
	    	  }
		  

		  
	    	}
	      else if ( face->at_boundary()  && face->boundary_indicator() == 15 )
	    	{
	    	  hp_fe_values_face[component]->reinit (cell,face_num);
	    	  const FEFaceValues<dim>& fe_values_face = hp_fe_values_face[component]->get_present_fe_values ();
	    	  const Quadrature<dim-1>& face_quadrature_formula = fe_values_face.get_quadrature();
	    	  const unsigned int n_q_points_face = face_quadrature_formula.size();
		    
	    	  const std::vector<double> &JxW_face = fe_values_face.get_JxW_values ();
	    	  const std::vector<Point<dim> > &normals = fe_values_face.get_normal_vectors ();
		    
	    	  for(unsigned int k=0;k<alphadim;k++){
	    	    prev_soln_alpha_face[k] =  std::vector<double>(n_q_points_face);
	    	    prev_soln_sigma_face[k] =  std::vector<Tensor<1,dim> >(n_q_points_face);
	    	  }

	    	  fe_values_face[*(alpha[component])].get_function_values(subdomain_solution[component],
	    	  							  prev_soln_alpha_face[component]);
	    	  fe_values_face[*(sigma[component])].get_function_values(subdomain_solution[component],
	    	  							  prev_soln_sigma_face[component]);
	    	  fe_values_face[*(alpha[1])].get_function_values(subdomain_solution[1],
	    							  prev_soln_alpha_face[1]);
	    	  fe_values_face[*(sigma[1])].get_function_values(subdomain_solution[1],
	    							  prev_soln_sigma_face[1]);

	    	  std::vector<double> alphas_boundary(n_q_points_face);
		  std::vector<double> Salphas_boundary(n_q_points_face);

	    	  unsigned char boundary_index = face->boundary_indicator();
	    	  const WallBoundaryValues<dim>& wbv = WBV[boundary_index];
	    	  const std::vector< Point<dim> > &quadrature_points = fe_values_face.get_quadrature_points();

	    	  alphas_boundary = prev_soln_alpha_face[component];
	    	  wbv.value_list2( quadrature_points, alphas_boundary, component, current_time);

		  const unsigned int normal1 = GeometryInfo<dim>::unit_normal_direction[face_num];
		  
		  double pen1 = cell->extent_in_direction(normal1);
		  
		  double sigma_p1 = degree*( degree + 1.0 ) / pen1;
		  double sigma_p2 = degree*( degree + 1.0 ) / pen1;
		  
		  double sigma_penalty = sigma_prefactor*0.5*(sigma_p1+sigma_p2);

	    	  for (unsigned int i=0; i<dofs_per_cell; ++i){
	    	    for (unsigned int q=0; q<n_q_points_face; ++q){

		      cell_rhs_boundary(i) +=  ( - sigma_penalty * ( alphas_boundary[q]- prev_soln_alpha_face[component][q] )
						 * fe_values_face[*(alpha[component])].value(i,q) ) * JxW_face[q];

	    	    }
	    	  }
	    	}
	      else{ //The is the face integration step
	      }

	      for (unsigned int i=0; i<dofs_per_cell; ++i){ 
		cell_rhs_composite(i) = cell_rhs_boundary(i)+cell_rhs(i);}

	      elliptic_constraints.distribute_local_to_global(cell_rhs_composite,
							      local_dof_indices,
							      poisson_rhs);
	    }
	  }
      
    
      poisson_rhs.compress (VectorOperation::add);
    
      SolverControl solver_control (n_iters, conv_threshold, true, true);
      
      PETScWrappers::SolverGMRES solver(solver_control, mpi_communicator);
      
      solver.solve (poisson_matrix.block(0,0), 
		    subdomain_solution[component].block(0), 
		    poisson_rhs.block(0),
      		    preconditioner);
      
      pcout << "\033[1;37m   Solved in " << solver_control.last_step()
	    << " iterations " << std::endl;
    }
  }

  for (unsigned int component=0; component< alphadim; ++component){  
    delete hp_fe_values[component];
    delete hp_fe_values_face[component];
  }
}

