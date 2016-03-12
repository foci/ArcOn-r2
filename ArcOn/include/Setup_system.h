/* Setup the system */
template <int dim>
void arcOn<dim>::setup_system()
{  

  num_blocks = 2;
  
  subdomain_solution.resize(alphadim);
  naive_subdomain_solution.resize(alphadim);
  init_solution.resize(alphadim);
  exact_solution.resize(alphadim);
  solution_diff.resize(alphadim);
  initial_condition.resize(alphadim);
  locally_owned_dofs.resize(alphadim);
  locally_relevant_dofs.resize(alphadim);
  L2_error_interpolant.resize(alphadim);
  naive_L2_error_interpolant.resize(alphadim);
  L2_error_method.resize(alphadim);
  naive_L2_error_method.resize(alphadim);
  L2_interpolate_active.resize(alphadim);
  interpolation_error.resize(alphadim);
  naive_interpolation_error.resize(alphadim);

  subdomain_solution_holder.resize(alphadim);
  interpolated_solution_holder.resize(alphadim);

  pinfo.resize(alphadim);
  mapinfo.resize(alphadim);
  gmax_val_loc.resize(alphadim);
  for (unsigned int component=0; component< alphadim; ++component){
    mapinfo[component].resize(1);
    gmax_val_loc[component] = 0.0;
  }
  
  difs = Vector<double>(alphadim);
  difs(0) = Mass_dif;
  difs(1) = Vort_dif;
  difs(2) = 0.0;
  
  local_solution.resize(alphadim);

  ilocal_solution.resize(alphadim);
  interpolate_base.resize(alphadim);
  interpolate_active.resize(alphadim);
  cont_output1.resize(alphadim);
  cont_global.resize(alphadim);
  revert_output.resize(alphadim);
  naive_revert_output.resize(alphadim);
  proj_solution.resize(alphadim);
  ilocally_owned_dofs.resize(alphadim);
  ilocally_relevant_dofs.resize(alphadim);

  tlocally_owned_dofs.resize(alphadim);
  tlocally_relevant_dofs.resize(alphadim);

  llocal_solution.resize(alphadim);
  linterpolate_active.resize(alphadim);
  llocally_owned_dofs.resize(alphadim);
  llocally_relevant_dofs.resize(alphadim);


  RK_solution.resize(RK_stage+1);
  naive_RK_solution.resize(RK_stage+1);
  naive_RK_solution_temp.resize(RK_stage+1);
  RK_div_flux.resize(RK_stage+1);
  naive_RK_div_flux.resize(RK_stage+1);
  naive_RK_div_flux_temp.resize(RK_stage+1);
  RK_MassAction.resize(RK_stage+1);
  naive_RK_MassAction.resize(RK_stage+1);
  naive_RK_MassAction_temp.resize(RK_stage+1);

  div_flux_integrated.resize(alphadim);
  naive_div_flux_integrated.resize(alphadim);

  MassAction_integrated.resize(alphadim);
  naive_MassAction_integrated.resize(alphadim);

  for(unsigned int s=0; s<RK_stage+1; s++){
    RK_solution[s].resize(alphadim);
    naive_RK_solution[s].resize(alphadim);
    naive_RK_solution_temp[s].resize(alphadim);


    RK_div_flux[s].resize(alphadim);
    naive_RK_div_flux[s].resize(alphadim);
    naive_RK_div_flux_temp[s].resize(alphadim);

    RK_MassAction[s].resize(alphadim);
    naive_RK_MassAction[s].resize(alphadim);
    naive_RK_MassAction_temp[s].resize(alphadim);
  }

  
  for (unsigned int component=0; component< alphadim; ++component){
 
    /* Do this for each component */
 
    std::vector<unsigned int> fesystem_sub_blocks (dim+1,0);
    fesystem_sub_blocks[dim-1] = 1;
    fesystem_sub_blocks[dim]   = 1;

    /////////////////////////////////////
    l_soln.resize(num_blocks);
    i_soln.resize(num_blocks); 
    r_soln.resize(num_blocks); 
    ////////////////////////////////////

    std::vector<unsigned int> ifesystem_sub_blocks (dim+1,0);
    ifesystem_sub_blocks[dim-1] = 1;
    ifesystem_sub_blocks[dim]   = 1;

    std::vector<unsigned int> tfesystem_sub_blocks (dim+1,0);
    tfesystem_sub_blocks[dim-1] = 1;
    tfesystem_sub_blocks[dim]   = 1;
    
    dof_handler[component]->distribute_dofs(*(fe_collection[component]));
    idof_handler[component]->distribute_dofs(*(ife_collection[component]));
    tdof_handler[component]->distribute_dofs(*(tfe_collection[component]));
    //ldof_handler[component]->distribute_dofs(*(lfe_collection[component]));

    DoFRenumbering::component_wise (*(dof_handler[component]), fesystem_sub_blocks );
    DoFRenumbering::component_wise (*(idof_handler[component]), ifesystem_sub_blocks );
    DoFRenumbering::component_wise (*(tdof_handler[component]), tfesystem_sub_blocks );

    // Let's do the indexset gymnastics now
    std::vector<unsigned int> fesystem_dofs_per_block (num_blocks);
    DoFTools::count_dofs_per_block (*(dof_handler[component]), fesystem_dofs_per_block,
				    fesystem_sub_blocks);

    std::vector<unsigned int> ifesystem_dofs_per_block (num_blocks);
    DoFTools::count_dofs_per_block (*(idof_handler[component]), ifesystem_dofs_per_block,
				    ifesystem_sub_blocks);

    std::vector<unsigned int> tfesystem_dofs_per_block (num_blocks);
    DoFTools::count_dofs_per_block (*(tdof_handler[component]), tfesystem_dofs_per_block,
				    tfesystem_sub_blocks);
    
    const unsigned int n_alpha = fesystem_dofs_per_block[0],
      n_sigma = fesystem_dofs_per_block[1];

    const unsigned int n_ialpha = ifesystem_dofs_per_block[0],
      n_isigma = ifesystem_dofs_per_block[1];

    const unsigned int n_talpha = tfesystem_dofs_per_block[0],
      n_tsigma = tfesystem_dofs_per_block[1];

    std::vector<IndexSet> fesystem_partitioning, fesystem_relevant_partitioning;
    std::vector<IndexSet> ifesystem_partitioning, ifesystem_relevant_partitioning;
    std::vector<IndexSet> tfesystem_partitioning, tfesystem_relevant_partitioning;

    locally_owned_dofs[component] = dof_handler[component]->locally_owned_dofs();
    ilocally_owned_dofs[component] = idof_handler[component]->locally_owned_dofs();
    tlocally_owned_dofs[component] = tdof_handler[component]->locally_owned_dofs();
    //llocally_owned_dofs[component] = ldof_handler[component]->locally_owned_dofs();

    fesystem_partitioning.push_back(locally_owned_dofs[component]
				    .get_view(0,n_alpha));
    fesystem_partitioning.push_back(locally_owned_dofs[component]
				    .get_view(n_alpha,n_alpha+n_sigma));

    ifesystem_partitioning.push_back(ilocally_owned_dofs[component]
				     .get_view(0,n_ialpha));
    ifesystem_partitioning.push_back(ilocally_owned_dofs[component]
				     .get_view(n_ialpha,n_ialpha+n_isigma));

    tfesystem_partitioning.push_back(tlocally_owned_dofs[component]
				     .get_view(0,n_talpha));
    tfesystem_partitioning.push_back(tlocally_owned_dofs[component]
				     .get_view(n_talpha,n_talpha+n_tsigma));
    
    DoFTools::extract_locally_relevant_dofs (*(dof_handler[component]),
					     locally_relevant_dofs[component]);
    DoFTools::extract_locally_relevant_dofs (*(idof_handler[component]),
					     ilocally_relevant_dofs[component]);
    DoFTools::extract_locally_relevant_dofs (*(tdof_handler[component]),
					     tlocally_relevant_dofs[component]);
    //DoFTools::extract_locally_relevant_dofs (*(ldof_handler[component]),
    //					     llocally_relevant_dofs[component]);

    fesystem_relevant_partitioning
      .push_back(locally_relevant_dofs[component].get_view(0,n_alpha));
    fesystem_relevant_partitioning
      .push_back(locally_relevant_dofs[component].get_view(n_alpha,n_alpha+n_sigma));
    ifesystem_relevant_partitioning
      .push_back(ilocally_relevant_dofs[component].get_view(0,n_ialpha));
    ifesystem_relevant_partitioning
      .push_back(ilocally_relevant_dofs[component].get_view(n_ialpha,n_ialpha+n_isigma));
    tfesystem_relevant_partitioning
      .push_back(tlocally_relevant_dofs[component].get_view(0,n_talpha));
    tfesystem_relevant_partitioning
      .push_back(tlocally_relevant_dofs[component].get_view(n_talpha,n_talpha+n_tsigma));


    std::vector< unsigned int > num_blocks1;
    std::vector< unsigned int > num_blocks2;
    std::vector< unsigned int > num_blocks3;
    
    num_blocks1.resize(num_blocks,0);
    num_blocks2.resize(num_blocks,0);
    num_blocks3.resize(num_blocks,0);

    subdomain_solution[component].reinit(num_blocks1,mpi_communicator,num_blocks2);
    init_solution[component].reinit(num_blocks1,mpi_communicator,num_blocks2); 
    initial_condition[component].reinit(num_blocks1,mpi_communicator,num_blocks2); 
    exact_solution[component].reinit(num_blocks1,mpi_communicator,num_blocks2);  
    solution_diff[component].reinit(num_blocks1,mpi_communicator,num_blocks2);  
    div_flux_integrated[component].reinit(num_blocks1,mpi_communicator,num_blocks2);
    MassAction_integrated[component].reinit(num_blocks1,mpi_communicator,num_blocks2);

    subdomain_solution_holder[component].reinit(num_blocks1,mpi_communicator,num_blocks2);
    interpolated_solution_holder[component].reinit(num_blocks1,mpi_communicator,num_blocks2);

    naive_subdomain_solution[component].reinit(num_blocks1,mpi_communicator,num_blocks2);
    naive_div_flux_integrated[component].reinit(num_blocks1,mpi_communicator,num_blocks2); 
    naive_MassAction_integrated[component].reinit(num_blocks1,mpi_communicator,num_blocks2); 
    poisson_rhs.reinit(num_blocks1,mpi_communicator,num_blocks2);
    cont_poisson_rhs.reinit(num_blocks1,mpi_communicator,num_blocks2);
    local_solution[component].reinit(num_blocks1); 
    ilocal_solution[component].reinit(num_blocks1); 

    interpolate_base[component].reinit(num_blocks1,mpi_communicator,num_blocks2); 
    cont_output1[component].reinit(num_blocks1,mpi_communicator,num_blocks2); 
    cont_global[component].reinit(num_blocks1,mpi_communicator,num_blocks2); 
    revert_output[component].reinit(num_blocks1,mpi_communicator,num_blocks1);
    naive_revert_output[component].reinit(num_blocks1,mpi_communicator,num_blocks2);
    interpolate_active[component].reinit(num_blocks1,mpi_communicator,num_blocks2); 
    proj_solution[component].reinit(num_blocks1,mpi_communicator,num_blocks2);

    L2_error_interpolant[component].reinit(num_blocks3,mpi_communicator,num_blocks3); 
    naive_L2_error_interpolant[component].reinit(num_blocks3,mpi_communicator,num_blocks3); 
    L2_error_method[component].reinit(num_blocks3,mpi_communicator,num_blocks3); 
    naive_L2_error_method[component].reinit(num_blocks3,mpi_communicator,num_blocks3); 
    L2_interpolate_active[component].reinit(num_blocks3,mpi_communicator,num_blocks3);

    naive_interpolation_error[component].reinit(num_blocks2,mpi_communicator,num_blocks2);
    interpolation_error[component].reinit(num_blocks2,mpi_communicator,num_blocks2);


    for (unsigned int bl=0; bl< num_blocks; ++bl){
      
      subdomain_solution[component].block(bl).reinit(fesystem_partitioning[bl], 
						     fesystem_relevant_partitioning[bl],
						     mpi_communicator);
      init_solution[component].block(bl).reinit(fesystem_partitioning[bl], 
						fesystem_relevant_partitioning[bl],
						mpi_communicator );
      initial_condition[component].block(bl).reinit(fesystem_partitioning[bl], 
						    fesystem_relevant_partitioning[bl],
						    mpi_communicator);
      exact_solution[component].block(bl).reinit(fesystem_partitioning[bl], 
						 fesystem_relevant_partitioning[bl],
						 mpi_communicator);
      solution_diff[component].block(bl).reinit(fesystem_partitioning[bl], 
						fesystem_relevant_partitioning[bl],
						mpi_communicator);
      div_flux_integrated[component].block(bl).reinit(fesystem_partitioning[bl], 
						      fesystem_relevant_partitioning[bl],
						      mpi_communicator);
      MassAction_integrated[component].block(bl).reinit(fesystem_partitioning[bl], 
							fesystem_relevant_partitioning[bl],
							mpi_communicator);
      subdomain_solution_holder[component].block(bl).reinit(fesystem_partitioning[bl], 
							    fesystem_relevant_partitioning[bl],
							    mpi_communicator);
      interpolated_solution_holder[component].block(bl).reinit(fesystem_partitioning[bl], 
							       fesystem_relevant_partitioning[bl],
							       mpi_communicator);
      
      naive_subdomain_solution[component].block(bl).reinit(fesystem_partitioning[bl],
							   mpi_communicator);
      naive_div_flux_integrated[component].block(bl).reinit(fesystem_partitioning[bl],
							    mpi_communicator);
      naive_MassAction_integrated[component].block(bl).reinit(fesystem_partitioning[bl],
							      mpi_communicator);
      poisson_rhs.block(bl).reinit(fesystem_partitioning[bl],mpi_communicator);

      cont_poisson_rhs.block(bl).reinit( tfesystem_partitioning[bl],mpi_communicator );

      interpolate_base[component].block(bl).reinit(fesystem_partitioning[bl],
						   mpi_communicator);
      cont_output1[component].block(bl).reinit(tfesystem_partitioning[bl],
					       mpi_communicator);
      cont_global[component].block(bl).reinit(tfesystem_partitioning[bl], 
					      tfesystem_relevant_partitioning[bl],
					      mpi_communicator);
      revert_output[component].block(bl).reinit(fesystem_partitioning[bl], 
						fesystem_relevant_partitioning[bl],
						mpi_communicator);
      naive_revert_output[component].block(bl).reinit(fesystem_partitioning[bl],
						      mpi_communicator);

      interpolate_active[component].block(bl).reinit(ifesystem_partitioning[bl],
						     mpi_communicator);
      proj_solution[component].block(bl).reinit(fesystem_partitioning[bl],
						mpi_communicator);
      
      L2_error_interpolant[component].block(bl).reinit(tfesystem_partitioning[bl], 
						       tfesystem_relevant_partitioning[bl],
						       mpi_communicator);
      naive_L2_error_interpolant[component].block(bl).reinit(tfesystem_partitioning[bl],
							     mpi_communicator);
      L2_error_method[component].block(bl).reinit(tfesystem_partitioning[bl], 
						  tfesystem_relevant_partitioning[bl],
						  mpi_communicator);
      naive_L2_error_method[component].block(bl).reinit(tfesystem_partitioning[bl],
							mpi_communicator);
      L2_interpolate_active[component].block(bl).reinit(tfesystem_partitioning[bl],
							mpi_communicator);
      
      interpolation_error[component].block(bl).reinit(fesystem_partitioning[bl],
						      fesystem_relevant_partitioning[bl],
						      mpi_communicator);
      naive_interpolation_error[component].block(bl).reinit(fesystem_partitioning[bl],
							    mpi_communicator);


    }

    subdomain_solution[component].collect_sizes();
    init_solution[component].collect_sizes();
    initial_condition[component].collect_sizes();
    exact_solution[component].collect_sizes();
    solution_diff[component].collect_sizes();
    div_flux_integrated[component].collect_sizes();
    MassAction_integrated[component].collect_sizes();

    naive_subdomain_solution[component].collect_sizes();
    naive_div_flux_integrated[component].collect_sizes();
    naive_MassAction_integrated[component].collect_sizes();
    poisson_rhs.collect_sizes();
    cont_poisson_rhs.collect_sizes();

    subdomain_solution_holder[component].collect_sizes();
    interpolated_solution_holder[component].collect_sizes();

    interpolate_base[component].collect_sizes();
    cont_output1[component].collect_sizes();
    cont_global[component].collect_sizes();
    revert_output[component].collect_sizes();
    naive_revert_output[component].collect_sizes();
    interpolate_active[component].collect_sizes();
    proj_solution[component].collect_sizes();

    L2_error_interpolant[component].collect_sizes();
    naive_L2_error_interpolant[component].collect_sizes();
    L2_error_method[component].collect_sizes();
    naive_L2_error_method[component].collect_sizes();
    L2_interpolate_active[component].collect_sizes();

    naive_interpolation_error[component].collect_sizes();
    interpolation_error[component].collect_sizes();


    for(unsigned int s=0; s<RK_stage+1; s++){
      
      RK_solution[s][component].reinit(num_blocks1,mpi_communicator,num_blocks2);
      RK_div_flux[s][component].reinit(num_blocks1,mpi_communicator,num_blocks2);
      RK_MassAction[s][component].reinit(num_blocks1,mpi_communicator,num_blocks2);
      
      naive_RK_solution[s][component].reinit(num_blocks1,mpi_communicator,num_blocks2);
      naive_RK_solution_temp[s][component].reinit(num_blocks1,mpi_communicator,num_blocks2);
      naive_RK_div_flux[s][component].reinit(num_blocks1,mpi_communicator,num_blocks2);
      naive_RK_div_flux_temp[s][component].reinit(num_blocks1,mpi_communicator,num_blocks2);
      naive_RK_MassAction[s][component].reinit(num_blocks1,mpi_communicator,num_blocks2);
      naive_RK_MassAction_temp[s][component].reinit(num_blocks1,mpi_communicator,num_blocks2);

      for (unsigned int bl=0; bl< num_blocks; ++bl){
	
	RK_solution[s][component].block(bl).reinit(fesystem_partitioning[bl], 
						   fesystem_relevant_partitioning[bl],
						   mpi_communicator);
	RK_div_flux[s][component].block(bl).reinit(fesystem_partitioning[bl], 
						   fesystem_relevant_partitioning[bl],
						   mpi_communicator);
	RK_MassAction[s][component].block(bl).reinit(fesystem_partitioning[bl], 
						     fesystem_relevant_partitioning[bl],
						     mpi_communicator);

	naive_RK_solution[s][component].block(bl).reinit(fesystem_partitioning[bl],
							 mpi_communicator);
	naive_RK_solution_temp[s][component].block(bl).reinit(fesystem_partitioning[bl],
							      mpi_communicator);
	naive_RK_div_flux[s][component].block(bl).reinit(fesystem_partitioning[bl],
							 mpi_communicator);
	naive_RK_div_flux_temp[s][component].block(bl).reinit(fesystem_partitioning[bl],
							      mpi_communicator);
	naive_RK_MassAction[s][component].block(bl).reinit(fesystem_partitioning[bl],
							   mpi_communicator);
	naive_RK_MassAction_temp[s][component].block(bl).reinit(fesystem_partitioning[bl],
								mpi_communicator);
	
      }

      RK_solution[s][component].collect_sizes();
      RK_div_flux[s][component].collect_sizes();
      RK_MassAction[s][component].collect_sizes();

      naive_RK_solution[s][component].collect_sizes();
      naive_RK_solution_temp[s][component].collect_sizes();
      naive_RK_div_flux[s][component].collect_sizes();
      naive_RK_div_flux_temp[s][component].collect_sizes();
      naive_RK_MassAction[s][component].collect_sizes();
      naive_RK_MassAction_temp[s][component].collect_sizes();
      
    }

    naive_L2_error_interpolant[component] = L2_error_interpolant[component];
    naive_L2_error_method[component] = L2_error_method[component];
    naive_subdomain_solution[component] = subdomain_solution[component];
    
    poisson_matrix.reinit(num_blocks,num_blocks);
    cont_poisson_matrix.reinit(num_blocks,num_blocks);

    std::vector< unsigned int> n_block_dofs;
    n_block_dofs.resize(num_blocks,0);
    n_block_dofs[0] = n_alpha;
    n_block_dofs[1] = n_sigma;
    
    poisson_matrix.reinit(num_blocks,num_blocks);
    
    for (unsigned int bl=0; bl< num_blocks; ++bl){
      for (unsigned int bl2=0; bl2< num_blocks; ++bl2){
	poisson_matrix.block(bl,bl2).reinit (mpi_communicator,
					     fesystem_partitioning[bl].size(),
					     fesystem_partitioning[bl2].size(),
					     fesystem_partitioning[bl].n_elements(),
					     fesystem_partitioning[bl2].n_elements(),
					     0,
					     false);
      }
    }
    
    poisson_matrix.collect_sizes();
    
    cont_poisson_matrix.reinit(num_blocks,num_blocks); 
    
    for (unsigned int bl=0; bl< num_blocks; ++bl){
      for (unsigned int bl2=0; bl2< num_blocks; ++bl2){
	cont_poisson_matrix.block(bl,bl2).reinit (mpi_communicator,
						  tfesystem_partitioning[bl].size(),
						  tfesystem_partitioning[bl2].size(),
						  tfesystem_partitioning[bl].n_elements(),
						  tfesystem_partitioning[bl2].n_elements(),
						  0,
						  false);
      }
    }
    
    cont_poisson_matrix.collect_sizes();
    
    if(component == 0){
      density_constraints.clear ();
      density_constraints.reinit ( locally_relevant_dofs[0] );
      density_constraints.close ();
    }
    if(component == 1){
      vorticity_constraints.clear ();
      vorticity_constraints.reinit ( locally_relevant_dofs[1] );
      vorticity_constraints.close ();
    }
    if(component == 2){
      potential_constraints.clear ();
      potential_constraints.reinit ( locally_relevant_dofs[2] );
      potential_constraints.close ();
    }
    if(component == 0){
      tdensity_constraints.clear ();
      tdensity_constraints.reinit ( locally_relevant_dofs[0] );
      tdensity_constraints.close ();
    }
    if(component == 1){
      tvorticity_constraints.clear ();
      tvorticity_constraints.reinit ( locally_relevant_dofs[1] );
      tvorticity_constraints.close ();
    }
    if(component == 2){
      tpotential_constraints.clear ();
      tpotential_constraints.reinit ( locally_relevant_dofs[2] );
      tpotential_constraints.close ();
    }

    if(component == 2){
      elliptic_constraints.clear ();
      elliptic_constraints.reinit ( locally_relevant_dofs[2] );
      elliptic_constraints.close ();


      //from here tellipt
      telliptic_constraints.clear ();
      telliptic_constraints.reinit ( tlocally_relevant_dofs[2] );

      DoFTools::make_periodicity_constraints(*tdof_handler[2],
      					     /*b_id*/ 1,
      					     /*b_id*/ 2,
      					     /*direction*/ 1,
      					     telliptic_constraints);
      DoFTools::make_periodicity_constraints(*tdof_handler[2],
      					     /*b_id*/ 3,
      					     /*b_id*/ 4,
      					     /*direction*/ 0,
      					     telliptic_constraints);
      telliptic_constraints.close ();
    }
  }
}
