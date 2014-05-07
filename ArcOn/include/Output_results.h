/* This template function outputs the state vector in a cycle each time it is called ... in paraview format (this is easily changed) */
template <int dim>
void arcOn<dim>::output_results (const unsigned int cycle_in) 
{
  static int cycle=0;
  (void)cycle_in;

  //This is the interpolation info                              

  /* for(unsigned int k=0; k<alphadim; k++){ */

  /*   naive_subdomain_solution[k]= subdomain_solution[k]; */
  /*   FETools::interpolate (*(dof_handler[k]),  naive_subdomain_solution[k], */
  /* 			  *(tdof_handler[k]), cont_output1[k]); */

  /*   cont_global[k] = cont_output1[k]; */
    
  /* } */

  for(unsigned int k=0; k<alphadim; k++){
    if(k==0){
      
      DataOut<dim,DoFHandler<dim> >  data_out0;
      data_out0.attach_dof_handler (*(dof_handler[k]));
      data_out0.add_data_vector (subdomain_solution[k],
				 rmhd<dim>::component_names (k),
				 DataOut<dim >::type_dof_data,
				 rmhd<dim>::component_interpretation());
	  
      data_out0.build_patches (degree+1);
	  
      const std::string filename = ("density_output/density-" +
				    Utilities::int_to_string (cycle, 5) +
				    "." +
				    Utilities::int_to_string
				    (triangulation.locally_owned_subdomain(), 4) +
				    ".vtu");
      std::ofstream output (filename.c_str());
      data_out0.write_vtu (output);
	  
      if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
	{
	  std::vector<std::string> filenames;
	  for (unsigned int i=0; i<Utilities::MPI::n_mpi_processes(mpi_communicator); ++i)
	    filenames.push_back (std::string("density-") +
				 Utilities::int_to_string (cycle, 5) +
				 "." +
				 Utilities::int_to_string(i, 4) +
				 ".vtu");
	  const std::string
	    pvtu_master_filename = ("density_output/density-" +
				    Utilities::int_to_string (cycle, 5) +
				    ".pvtu");
	  std::ofstream pvtu_master (pvtu_master_filename.c_str());
	  data_out0.write_pvtu_record (pvtu_master, filenames);
	      
	}
    }
      
    if(k==1){
	  
      //revert_output[0] = std::exp(subdomain_solution[0]);
      DataOut<dim,DoFHandler<dim> >  data_out1;
      data_out1.attach_dof_handler (*(dof_handler[k]));
      data_out1.add_data_vector (subdomain_solution[k],
				 rmhd<dim>::component_names (k),
				 DataOut<dim >::type_dof_data,
				 rmhd<dim>::component_interpretation());
	  
      data_out1.build_patches (degree+1);
	  
      const std::string filename = ("vorticity_output/vorticity-" +
				    Utilities::int_to_string (cycle, 5) +
				    "." +
				    Utilities::int_to_string
				    (triangulation.locally_owned_subdomain(), 4) +
				    ".vtu");
      std::ofstream output (filename.c_str());
      data_out1.write_vtu (output);
	  
      if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
	{
	  std::vector<std::string> filenames;
	  for (unsigned int i=0; i<Utilities::MPI::n_mpi_processes(mpi_communicator); ++i)
	    filenames.push_back (std::string("vorticity-") +
				 Utilities::int_to_string (cycle, 5) +
				 "." +
				 Utilities::int_to_string(i, 4) +
				 ".vtu");
	  const std::string
	    pvtu_master_filename = ("vorticity_output/vorticity-" +
				    Utilities::int_to_string (cycle, 5) +
				    ".pvtu");
	  std::ofstream pvtu_master (pvtu_master_filename.c_str());
	  data_out1.write_pvtu_record (pvtu_master, filenames);

	}
    }
	
    if(k==2){
	  
      DataOut<dim,DoFHandler<dim> >  data_out2;
      data_out2.attach_dof_handler (*(dof_handler[k]));
      data_out2.add_data_vector (subdomain_solution[k],
				 rmhd<dim>::component_names (k),
				 DataOut<dim >::type_dof_data,
				 rmhd<dim>::component_interpretation());
	  
      data_out2.build_patches (degree+1);
	  
      const std::string filename = ("potential_output/potential-" +
				    Utilities::int_to_string (cycle, 5) +
				    "." +
				    Utilities::int_to_string
				    (triangulation.locally_owned_subdomain(), 4) +
				    ".vtu");
      std::ofstream output (filename.c_str());
      data_out2.write_vtu (output);
	  
      if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
	{
	  std::vector<std::string> filenames;
	  for (unsigned int i=0; i<Utilities::MPI::n_mpi_processes(mpi_communicator); ++i)
	    filenames.push_back (std::string("potential-") +
				 Utilities::int_to_string (cycle, 5) +
				 "." +
				 Utilities::int_to_string(i, 4) +
				 ".vtu");
	  const std::string
	    pvtu_master_filename = ("potential_output/potential-" +
				    Utilities::int_to_string (cycle, 5) +
				    ".pvtu");
	  std::ofstream pvtu_master (pvtu_master_filename.c_str());
	  data_out2.write_pvtu_record (pvtu_master, filenames);

	}
    }
  }

  cycle++;

}


