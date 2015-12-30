/* Error Analysis of the solution */
template <int dim>
void arcOn<dim>::cmcmv(SolutionVector& subdomain_solution)
{
  std::vector< FEValues<dim>* > hp_fe_values;
  
  for (unsigned int component=0; component< alphadim; ++component){
    
    if(component == 0){
  
    UpdateFlags updateflags=  update_values | 
      update_gradients | update_quadrature_points | update_JxW_values;
    
    hp_fe_values.push_back (new FEValues<dim>(*(tfe_collection[component]), 
					      *(tquadrature_collection[component]), 
					      updateflags));
    InitialValues<dim> IV;
    double total_densityl = 0.0;
    double cml = 0.0;
    double cmvl = 0.0;
    total_volume = 0.0;
    typename DoFHandler<dim>::active_cell_iterator 
      cell = dof_handler[component]->begin_active(), 
      endc = dof_handler[component]->end();
    for(;
	cell!=endc;
	++cell )
      if (cell->is_locally_owned() )
	{
	  hp_fe_values[component]->reinit (cell);
	  unsigned int CO = cell->index();
	  const FEValues<dim> &fe_values = hp_fe_values[component]->get_present_fe_values ();
	  const Quadrature<dim>& quadrature_formula = fe_values.get_quadrature(); 
          const std::vector<Point<dim> >& quadrature_point = fe_values.get_quadrature_points();
	  const unsigned int   dofs_per_cell = fe_values.dofs_per_cell;
	  const unsigned int   n_q_points    = quadrature_formula.size();

	  const std::vector<double>& quad_weights = quadrature_formula.get_weights();

	  const FiniteElement<dim>& fe = cell->get_fe();

	  std::vector<unsigned int> local_dof_indices (dofs_per_cell);

	  cell->get_dof_indices (local_dof_indices);
	  const std::vector<double> &JxW = fe_values.get_JxW_values ();
	  double JxWsum = 0.0;
	  for (unsigned int q=0; q<n_q_points; ++q){
	    JxWsum += JxW[q];
	  }
	  total_volume += JxWsum;
	  

	  for(unsigned int k=0;k<alphadim;k++){
	    prev_soln_alpha[k] =  std::vector<double>(n_q_points);
	    soln_alpha[k] =  std::vector<double>(n_q_points);
	  }
	  
	  fe_values[*(alpha[component])].get_function_values(subdomain_solution[component],prev_soln_alpha[component]);
	  fe_values[*(alpha[component])].get_function_values(subdomain_solution[component],soln_alpha[component]);
	  
	  const std::vector< Point<dim> >& quadrature_points = fe_values.get_quadrature_points();
	  
	  for(unsigned int q=0;q<n_q_points;q++){
	    for (unsigned int i=0; i<dofs_per_cell; ++i){
	        
	      total_densityl += (soln_alpha[component])[q] *  fe_values[*(alpha[component])].value(i,q) * JxW[q];

	      //pcout << "alpha = " << (soln_alpha[component])[q] << ", base = " << fe_values[*(alpha[component])].value(i,q) << ". x = " <<  quadrature_point[q][0] << ", JxW = " << JxW[q] << std::endl;
	      
	      cml +=  (soln_alpha[component])[q] *  fe_values[*(alpha[component])].value(i,q) * quadrature_point[q][0] * JxW[q];
  
	    }
	  }

	}

    //total_density /= total_volume;
    //cml /= total_volume;

    //cml /= (total_density);
    //cmvl = dt;

    MPI_Allreduce(&cml, &cm, Utilities::MPI::n_mpi_processes(mpi_communicator), MPI_DOUBLE, MPI_SUM, mpi_communicator);
    MPI_Allreduce(&cmvl, &cmv, Utilities::MPI::n_mpi_processes(mpi_communicator), MPI_DOUBLE, MPI_SUM, mpi_communicator);
    MPI_Allreduce(&total_densityl, &total_density, Utilities::MPI::n_mpi_processes(mpi_communicator), MPI_DOUBLE, MPI_SUM, mpi_communicator);


    cm /= total_density;
  
    }
    
  }
  
  
  //for (unsigned int component=0; component< alphadim; ++component){
    delete hp_fe_values[0];
    //}
  //return total_volume;
}
 
