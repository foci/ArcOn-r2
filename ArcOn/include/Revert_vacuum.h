
template <int dim>
void arcOn<dim>::revert_vacuum(SolutionVector& subdomain_solution, double delta_t, double current_time)
{

  std::vector< FEValues<dim>* > hp_fe_values;


    UpdateFlags updateflags=  update_values | update_gradients |
      update_quadrature_points | update_JxW_values;
    hp_fe_values.push_back (new FEValues<dim>(*(fe_collection[2]),
  					      *(quadrature_collection[2]), updateflags));
    
    double ltotal_volume = 0.0;
    double local_Avg = 0.0;
  
    /* We can using ramping to avoid shocking the system */
    double ramp = 0.5; // ramp time in seconds
    double dramp  = 1.0; //std::tanh( 2.0*( current_time/ramp ) );


    typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler[2]->begin_active(),
      endc = dof_handler[2]->end();
    for(;
  	cell!=endc;
  	++cell	)
      if ( cell->is_locally_owned() )
  	{
  	  int CO = cell->index();
  	  double pi = 3.1415926535897932384626433;

  	  hp_fe_values[0]->reinit (cell);
  	  const FEValues<dim> &fe_values = hp_fe_values[0]->get_present_fe_values ();
  	  const Quadrature<dim>& quadrature_formula = fe_values.get_quadrature();
  	  const std::vector<Point<dim> >& quadrature_point = fe_values.get_quadrature_points();
  	  const FiniteElement<dim>& fe = cell->get_fe();
  	  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  	  const unsigned int   n_q_points    = quadrature_formula.size();

  	  std::vector<unsigned int> local_dof_indices (dofs_per_cell);

  	  syncopated[2] =  std::vector<double>(n_q_points);
	  
  	  //cell->get_dof_indices (local_dof_indices);
  	  const std::vector<double> &JxW = fe_values.get_JxW_values ();
  	  double JxWsum = 0;
  	  for (unsigned int q=0; q<n_q_points; ++q){
  	    JxWsum += JxW[q];
  	  }
  	  ltotal_volume += JxWsum;

  	  fe_values[*(alpha[2])].get_function_values(subdomain_solution[2],
  						     syncopated[2]);
	  
  	  for (unsigned int q=0; q<n_q_points; ++q){
  	    for (unsigned int i=0; i<dofs_per_cell; ++i){
	      
  	      local_Avg += syncopated[2][q] * (fe_values[*(alpha[2])].value(i,q)) * JxW[q];
	      
  	    }
  	  }
	
  	  /* cell->get_dof_indices (local_dof_indices); */

  	  /* Vector<double> temp_update[alphadim];  */
  	  /* temp_update[component].reinit(dofs_per_cell);  */
  	  /* for (unsigned int i=0; i<dofs_per_cell; ++i){ */
  	  /*   temp_update[component](i) = n_from_chi(component,i); */
  	  /* } */
  	  /* Vector<double> projected(dofs_per_cell);  */

  	  /* localInverseMassMatrix.vmult(projected,temp_update[component]); */
  	  /* projected /= JxWsum; */
  	  /* parahyp_constraints[component].distribute_local_to_global (projected, */
  	  /* 							     local_dof_indices, */
  	  /* 							     naive_revert_output[component]);	   */
	 
  	}

 

    naive_subdomain_solution[2].block(0) = subdomain_solution[2].block(0);
    //local_Avg = naive_subdomain_solution[2].block(0).mean_value();
    
    MPI_Allreduce(&local_Avg, &global_Avg, 1, MPI_DOUBLE, MPI_SUM, mpi_communicator);
    if (current_time == delta_t){
      MPI_Allreduce(&ltotal_volume, &gtotal_volume, 1, MPI_DOUBLE, MPI_SUM, mpi_communicator);
    }
    
    //global_Avg /=  n_mpi_processes;
    global_Avg /= gtotal_volume;
    
    //pcout << "global_AVG = " << global_Avg << ", global_volume = " << gtotal_volume << std::endl;
    
    // naive_subdomain_solution[2].block(0) = subdomain_solution[2].block(0);
    //if (global_Avg < 0.0){
    naive_subdomain_solution[2].block(0).add(-global_Avg);
      //}
      //else{
      //naive_subdomain_solution[2].block(0).add(global_Avg);
      //}
    naive_subdomain_solution[2].compress(VectorOperation::add);
    subdomain_solution[2].block(0) = naive_subdomain_solution[2].block(0);
    
    delete hp_fe_values[0];

}
