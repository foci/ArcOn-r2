/* The template function solves the Poisson equation using a global linear solve .... */
template <int dim>
void arcOn<dim>::calc_poisson_cont(SolutionVector& subdomain_solution, double delta_t) { 

  (void) delta_t;

  std::vector< FEValues<dim>* > hp_fe_values;
  double current_time = 1.0;

  //pcout << "here1?" << std::endl;
    
  for (unsigned int k=0; k< alphadim; ++k){

    cont_global[k] = 0.0;
    cont_output1[k] = 0.0;

    //naive_subdomain_solution[k]= subdomain_solution[k];
    FETools::interpolate (*(dof_handler[k]),  subdomain_solution[k],
			  *(tdof_handler[k]), cont_output1[k]);

    //    cont_global[k].update_ghost_values();
    cont_global[k] = cont_output1[k];

  }

  for (unsigned int component=0; component< alphadim; ++component){
    UpdateFlags updateflags=  update_values | update_gradients 
      | update_quadrature_points | update_JxW_values;
    
    hp_fe_values.push_back (new FEValues<dim>(*(tfe_collection[component]),
					      *(tquadrature_collection[component]), updateflags));
 
    if (component == 2){

      cont_poisson_rhs = 0.0;

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
	      
	    Vector<double>       cell_rhs (dofs_per_cell);
	      
	    std::vector<unsigned int> local_dof_indices (dofs_per_cell);
	    std::vector<unsigned int> neigh_local_dof_indices (dofs_per_cell);

	    cell_rhs = 0.0;
	     
	    //pcout << "here2?" << std::endl;

	    for(unsigned int k=0;k<alphadim;k++){
	      prev_soln_alpha[k] =  std::vector<double>(n_q_points);
	      prev_soln_sigma[k] =  std::vector<Tensor<1,dim> >(n_q_points);
	    }
	      
	    cell->get_dof_indices (local_dof_indices);

	    const std::vector<double> &JxW = fe_values.get_JxW_values ();

	    fe_values[*(alpha[1])].get_function_values(cont_global[1],prev_soln_alpha[1]);
	    fe_values[*(sigma[1])].get_function_values(cont_global[1],prev_soln_sigma[1]);

	    std::vector<unsigned int> index_container2 (dofs_per_cell);

	    for (unsigned int i=0; i<dofs_per_cell; ++i){
	      for (unsigned int q=0; q<n_q_points; ++q){
		
		cell_rhs(i) -= (fe_values[*(alpha[component])].value(i,q))* (prev_soln_alpha[1][q]) * JxW[q] ;
		//cell_rhs(i) -= 0.0010*(fe_values[*(alpha[component])].value(i,q)) * JxW[q] ;
		
	      }
	    }

	    //cell_rhs = -0.0010;

	    telliptic_constraints.distribute_local_to_global(cell_rhs,
							     local_dof_indices,
							     cont_poisson_rhs);
	  }
    
      cont_poisson_rhs.compress (VectorOperation::add);
      //cont_poisson_matrix.compress (VectorOperation::insert); 
      //pcout << "here4?" << std::endl;
      
      SolverControl solver_control (n_iters, conv_threshold, true, true);

      //pcout << "here5?" << std::endl;
      
      PETScWrappers::SolverGMRES solver(solver_control, mpi_communicator);
      //PETScWrappers::SolverCG solver(solver_control, mpi_communicator);

      solver.solve (cont_poisson_matrix.block(0,0), 
		    cont_global[component].block(0), 
		    cont_poisson_rhs.block(0),
		    cont_preconditioner);
      
      pcout << "\033[1;37m   Solved in " << solver_control.last_step()
	    << " iterations " << std::endl;

      //pcout << "here1?" << std::endl;

      cont_output1[component].block(0) = cont_global[component].block(0);

      telliptic_constraints.distribute (cont_output1[component]);

      cont_global[component].block(0) = cont_output1[component].block(0);

      //telliptic_constraints.distribute (cont_global[component]);

      //pcout << "here2?" << std::endl;

 //for (unsigned int k=0; k< alphadim; ++k){
     
    //cont_output1[k] = cont_global[k];

  //  cont_global[2].update_ghost_values();

      //pcout << "here3?" << std::endl;
  
      FETools::interpolate (*(tdof_handler[component]), cont_global[component],
			    *(dof_handler[component]),  naive_subdomain_solution[component]);
      
  //  naive_subdomain_solution[2].update_ghost_values();

      //pcout << "here4?" << std::endl;
  
      subdomain_solution[component].block(0) = naive_subdomain_solution[component].block(0);
      
    }
  }
  
 
  //  }
  
  for (unsigned int component=0; component< alphadim; ++component){  
    delete hp_fe_values[component];
  }
}

