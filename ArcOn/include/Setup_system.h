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
    // mapinfo[component].resize(triangulation.n_global_active_cells());
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

  parahyp_constraints.resize(alphadim);
  tparahyp_constraints.resize(alphadim);
  
  //alpha_mask.resize(alphadim);
  
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
    //ll_soln.resize(num_blocks); 
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
      
      subdomain_solution[component].block(bl).reinit(mpi_communicator, 
						     fesystem_partitioning[bl], 
						     fesystem_relevant_partitioning[bl] );
      init_solution[component].block(bl).reinit(mpi_communicator, 
						fesystem_partitioning[bl], 
						fesystem_relevant_partitioning[bl] );
      initial_condition[component].block(bl).reinit(mpi_communicator, 
						    fesystem_partitioning[bl], 
						    fesystem_relevant_partitioning[bl] );
      exact_solution[component].block(bl).reinit(mpi_communicator, 
						 fesystem_partitioning[bl], 
						 fesystem_relevant_partitioning[bl] );
      solution_diff[component].block(bl).reinit(mpi_communicator, 
						fesystem_partitioning[bl], 
						fesystem_relevant_partitioning[bl] );
      div_flux_integrated[component].block(bl).reinit(mpi_communicator, 
						      fesystem_partitioning[bl], 
						      fesystem_relevant_partitioning[bl] );
      MassAction_integrated[component].block(bl).reinit(mpi_communicator, 
							fesystem_partitioning[bl], 
							fesystem_relevant_partitioning[bl] );

      subdomain_solution_holder[component].block(bl).reinit(mpi_communicator, 
							    fesystem_partitioning[bl], 
							    fesystem_relevant_partitioning[bl] );
      interpolated_solution_holder[component].block(bl).reinit(mpi_communicator, 
							       fesystem_partitioning[bl], 
							       fesystem_relevant_partitioning[bl] );
      
      naive_subdomain_solution[component].block(bl).reinit(mpi_communicator, 
							   fesystem_partitioning[bl] );
      naive_div_flux_integrated[component].block(bl).reinit(mpi_communicator, 
							    fesystem_partitioning[bl] );
      naive_MassAction_integrated[component].block(bl).reinit(mpi_communicator, 
							      fesystem_partitioning[bl] );
      poisson_rhs.block(bl).reinit(mpi_communicator, fesystem_partitioning[bl] );

      cont_poisson_rhs.block(bl).reinit(mpi_communicator, tfesystem_partitioning[bl] );

      interpolate_base[component].block(bl).reinit(mpi_communicator, 
						   fesystem_partitioning[bl] );
      cont_output1[component].block(bl).reinit(mpi_communicator,
                                                   tfesystem_partitioning[bl] );
      cont_global[component].block(bl).reinit(mpi_communicator, 
						     tfesystem_partitioning[bl], 
						     tfesystem_relevant_partitioning[bl] );
      revert_output[component].block(bl).reinit(mpi_communicator, 
						     fesystem_partitioning[bl], 
						     fesystem_relevant_partitioning[bl] );
      naive_revert_output[component].block(bl).reinit(mpi_communicator, 
						      fesystem_partitioning[bl] );

      /* cont_output1[component].block(bl).reinit(mpi_communicator, */
      /*                                                tfesystem_partitioning[bl], */
      /*                                                tfesystem_relevant_partitioning[bl] ); */
      interpolate_active[component].block(bl).reinit(mpi_communicator, 
						     ifesystem_partitioning[bl] );
      proj_solution[component].block(bl).reinit(mpi_communicator, 
						     fesystem_partitioning[bl] );
      
      L2_error_interpolant[component].block(bl).reinit(mpi_communicator, 
						       tfesystem_partitioning[bl], 
						       tfesystem_relevant_partitioning[bl] );
      naive_L2_error_interpolant[component].block(bl).reinit(mpi_communicator, 
							     tfesystem_partitioning[bl] );
      L2_error_method[component].block(bl).reinit(mpi_communicator, 
						  tfesystem_partitioning[bl], 
						  tfesystem_relevant_partitioning[bl] );
      naive_L2_error_method[component].block(bl).reinit(mpi_communicator, 
							tfesystem_partitioning[bl] );
      L2_interpolate_active[component].block(bl).reinit(mpi_communicator,
                                                        tfesystem_partitioning[bl] );
      
      interpolation_error[component].block(bl).reinit(mpi_communicator,
						      fesystem_partitioning[bl],
						      fesystem_relevant_partitioning[bl] );
      naive_interpolation_error[component].block(bl).reinit(mpi_communicator,
							    fesystem_partitioning[bl] );


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
	
	RK_solution[s][component].block(bl).reinit(mpi_communicator, 
						   fesystem_partitioning[bl], 
						   fesystem_relevant_partitioning[bl] );
	RK_div_flux[s][component].block(bl).reinit(mpi_communicator, 
						   fesystem_partitioning[bl], 
						   fesystem_relevant_partitioning[bl] );
	RK_MassAction[s][component].block(bl).reinit(mpi_communicator, 
						     fesystem_partitioning[bl], 
						     fesystem_relevant_partitioning[bl] );

	naive_RK_solution[s][component].block(bl).reinit(mpi_communicator, 
							 fesystem_partitioning[bl] );
	naive_RK_solution_temp[s][component].block(bl).reinit(mpi_communicator, 
							 fesystem_partitioning[bl] );
	naive_RK_div_flux[s][component].block(bl).reinit(mpi_communicator, 
							 fesystem_partitioning[bl] );
	naive_RK_div_flux_temp[s][component].block(bl).reinit(mpi_communicator, 
							 fesystem_partitioning[bl] );
	naive_RK_MassAction[s][component].block(bl).reinit(mpi_communicator, 
							   fesystem_partitioning[bl] );
	naive_RK_MassAction_temp[s][component].block(bl).reinit(mpi_communicator, 
							   fesystem_partitioning[bl] );
	
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
    
    /* for (unsigned int bl=0; bl< num_blocks; ++bl){ */
          
    /*   local_solution[component].block(bl) = subdomain_solution[component].block(bl); */
    /*   ilocal_solution[component].block(bl) = interpolate_active[component].block(bl); */
    /*   //llocal_solution[component] = linterpolate_active[component]; */
    
    /*   l_soln[bl] =  local_solution[component].block(bl); */
    /*   i_soln[bl] =  ilocal_solution[component].block(bl);  */
    /*   r_soln[bl] =  local_solution[component].block(bl); */

    /* } */

    //*alpha_mask[component] = fe_collection[component]->component_mask(*alpha[component]);

    parahyp_constraints[component].clear ();
    parahyp_constraints[component].reinit ( locally_relevant_dofs[component]);
    /* DoFTools::make_periodicity_constraints(*dof_handler[component], */
    /* 					   /\*b_id*\/ 1, */
    /* 					   /\*b_id*\/ 2, */
    /* 					   /\*direction*\/ 0, */
    /* 					   parahyp_constraints[component]); */
    parahyp_constraints[component].close ();

    tparahyp_constraints[component].clear ();
    tparahyp_constraints[component].reinit ( tlocally_relevant_dofs[component] );
    /* DoFTools::make_periodicity_constraints(*dof_handler[component], */
    /* 					   /\*b_id*\/ 1, */
    /* 					   /\*b_id*\/ 2, */
    /* 					   /\*direction*\/ 1, */
    /* 					   parahyp_constraints[component]); */
    tparahyp_constraints[component].close ();

    poisson_matrix.reinit(num_blocks,num_blocks);
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

 
    
    if(component == 2){
      elliptic_constraints.clear ();
      elliptic_constraints.reinit ( locally_relevant_dofs[2] );
      /* DoFTools::make_periodicity_constraints(*dof_handler[2], */
      /* 					     /\*b_id*\/ 1, */
      /* 					     /\*b_id*\/ 2, */
      /* 					     /\*direction*\/ 0, */
      /* 					     elliptic_constraints); */
      elliptic_constraints.close ();

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

    /* CompressedSimpleSparsityPattern csp (tlocally_relevant_dofs[2]); */
    /* DoFTools::make_sparsity_pattern (*(tdof_handler[2]), csp, */
    /* 				     constraints, false); */
    /* SparsityTools::distribute_sparsity_pattern (csp, */
    /* 						*(tdof_handler[2]).n_locally_owned_dofs_per_processor(), */
    /* 						mpi_communicator, */
    /* 						locally_relevant_dofs); */
    /* system_matrix.reinit (locally_owned_dofs, */
    /* 			  locally_owned_dofs, */
    /* 			  csp, */
    /* 			  mpi_communicator); */
    
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
    
  }
}

  //max_active_cells = triangulation.n_active_cells()* ((unsigned int) std::pow(2.0, (double)dim*(max_refined-init_mesh)));
  
  /* std::pout << "Number of components=" << alphadim << std::endl; */
  /* std::pout << " Total dimension =" << (dim+1)*alphadim << std::endl; */
  /* std::pout << "\n Number of active cells:       " << triangulation.n_active_cells() << "   Number of degrees of freedom for alpha[0]: " << dof_handler[0]->n_dofs() << std::endl; */
  /* std::pout << "\n Number of active cells:       " << triangulation.n_active_cells() << "   Number of degrees of freedom for alpha[1]: " << dof_handler[1]->n_dofs() << std::endl; */
  
  //num_vertex = triangulation.n_vertices();

  //std::cout << "Number of vertices suggested:       " << num_vertex << std::endl;

  /* for (unsigned int component=0; component< alphadim; ++component){ */
  /*   subdomain_solution[component].compress();} */

  //initialize_massmatrix();
  // matrixmapper();

  /*  sum_alpha_dofs = 0; */
  /*   global_sum_alpha_dofs = 0; */

  /*   std::vector<pInfo> pinfo(alphadim); */
      
  /*   std::vector< FEValues<dim>* > hp_fe_values; */
  /*   std::vector< FEFaceValues<dim>* > hp_fe_values_face; */
  /*   std::vector< FEFaceValues<dim>* > hp_fe_values_neigh_face; */

  /*   for (unsigned int component=0; component< alphadim; ++component){ */
      
  /*     UpdateFlags updateflags=  update_values | update_gradients | update_quadrature_points | update_JxW_values; */
	
  /*     hp_fe_values.push_back (new FEValues<dim>(*(fe_collection[component]), *(quadrature_collection[component]), updateflags)); */
  /*     hp_fe_values_face.push_back (new FEFaceValues<dim>(*(fe_collection[component]), *(face_quadrature_collection[component]),  updateflags | update_normal_vectors )); */
  /*     hp_fe_values_neigh_face.push_back (new FEFaceValues<dim>(*(fe_collection[component]), *(face_quadrature_collection[component]),  updateflags | update_normal_vectors )); */

  /*     sum_alpha_dofs = 0; */
  /*     //sum_vertices = 0;                                                                                                                                                                                                                 */
  /*     num_vert = triangulation.n_vertices(); */
 
  /*     //Let's test how many dofs are stored on this core */
  /*     typename DoFHandler<dim>::active_cell_iterator cell = dof_handler[component]->begin_active(), endc = dof_handler[component]->end(); */
  /*     for(; */
  /* 	cell!=endc; */
  /* 	++cell	) */
  /*       if (cell->is_locally_owned()  ) */
  /* 	{ */
	      
  /* 	  if(component == 2) { */
  /* 	    hp_fe_values[component]->reinit (cell); */
  /* 	    const FiniteElement<dim>& fe = cell->get_fe(); */
  /* 	    const std::string fe_name = fe.get_name(); */
  /* 	    pinfo[component] = pInfoFind(fe.get_name(),component); */

  /* 	    const unsigned int alpha_dofs_per_cell = pinfo[component].alpha_dofs_per_cell; */
  /* 	    sum_alpha_dofs = sum_alpha_dofs + alpha_dofs_per_cell; */
	      
  /* 	    for(unsigned int vi = 0; vi < GeometryInfo<dim>::vertices_per_cell; ++vi){ */
	      
  /* 	      int v_dex = cell->vertex_index(vi); */
  /* 	      //std::cout << "vdex = " << v_dex << ", " << vi << ", " <<  GeometryInfo<dim>::vertices_per_cell << std::endl; */
  /* 	      local_vertices.push_back(v_dex); */

  /* 	    } */
  /* 	    //} */
  /* 	  } */
  /* 	} */
  /*   } */

  /*   std::vector< int >::iterator r , w ; */
  /*   std::set< int > tmpset ; */
  
  /*   for( r = local_vertices.begin() , w = local_vertices.begin() ; r != local_vertices.end() ; ++r ) */
  /*     { */
  /*       if( tmpset.insert( *r ).second ) */
  /* 	{ */
  /* 	  *w++ = *r ; */
  /* 	} */
  /*     } */
  /*   lvsr = local_vertices.size(); */

  /*   lvsr_size = *(std::max_element(local_vertices.begin(), local_vertices.end())) + 1; */

  /*   //std::cout << "size = " << lvsr_size << ", num_vert" << std::endl; */
  
  /*   local_vertices.erase( w , local_vertices.end() ); */

  /* /\*   int maxv_local = std::max(local_vertices); *\/ */
  /* /\*   int minv_local = std::min(local_vertices); *\/ */

  /* /\*   std::cout << "Max, min = " << maxv_local << ", " << minv_local << std:endl; *\/ */

  /*   lvs = local_vertices.size(); */
  
  /*   //lvs = std::ceil(num_vertex/n_mpi_processes); */

  /*   for (unsigned int component=0; component< alphadim; ++component){      */
  /*     delete hp_fe_values[component]; */
  /*     delete hp_fe_values_face[component]; */
  /*     delete hp_fe_values_neigh_face[component]; */
  /*   } */

  //scalar_dof_handler.distribute_dofs(scalar_fe);
  //scalar_locally_owned_dofs = scalar_dof_handler.locally_owned_dofs();
  //DoFTools::extract_locally_relevant_dofs (scalar_dof_handler,
  //                                       scalar_locally_relevant_dofs);
  //elliptic_constraints.clear ();
  //elliptic_constraints.reinit (scalar_locally_relevant_dofs);
  /* if (component == 0){ */
  /* alpha_mask = fe_collection[component]->component_mask(*(alpha[component])); */
  /* std::cout << "alpha mask = " << alpha_mask << std::endl; */
  /* VectorTools::interpolate_boundary_values (*(dof_handler[component]), */
  /*                                        0, */
  /*                                        ZeroFunction<dim>(), */
  /*                                        parahyp_constraints[component], */
  /*                                        alpha_mask); */
  /* } */

  /* if (component == 1){ */
  /* alpha_mask1 = fe_collection[component]->component_mask(*(alpha[component])); */
  /* std::cout << "alpha mask = " << alpha_mask << std::endl; */
  /* VectorTools::interpolate_boundary_values (*(dof_handler[component]), */
  /*                                        0, */
  /*                                        ZeroFunction<dim>(), */
  /*                                        parahyp_constraints[component], */
  /*                                        alpha_mask1); */
  //}                                                                                                                                                                                                                                       
  //elliptic_constraints.close ();


  /* MPI_Allreduce(&sum_alpha_dofs, &global_sum_alpha_dofs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); */


  /* //poisson_rhs.reinit(mpi_communicator, scalar_locally_owned_dofs ); */
  /* //poisson_rhs.reinit(mpi_communicator, global_sum_alpha_dofs, sum_alpha_dofs, true ); */

  /* redistributed_solution.reinit(mpi_communicator, global_sum_alpha_dofs, sum_alpha_dofs, true ); */

  /* //std::cout << "WHYWHYWHY? = " << global_sum_alpha_dofs << ", and " << sum_alpha_dofs << std::endl; */




  /* /\* poisson_matrix.reinit( mpi_communicator, *\/ */
  /* /\* 			 global_sum_alpha_dofs, *\/ */
  /* /\* 			 global_sum_alpha_dofs, *\/ */
  /* /\* 			 sum_alpha_dofs, *\/ */
  /* /\* 			 sum_alpha_dofs, *\/ */
  /* /\* 			 0, *\/ */
  /* /\* 			 true  *\/ */
  /* /\* 			 ) ; *\/	 */

  /* //MPI_Allreduce(&lvs, &num_vert, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); */

  /* //  std::cout << " number of proces = " << n_mpi_processes <<", number local = " << std::endl; */

  /* const unsigned int g = 1; */

  /* maxf.reinit(mpi_communicator, */
  /* 	      num_vertex, */
  /* 	      n_mpi_processes, */
  /* 	      lvs, */
  /* 	      g, */
  /* 	      0, */
  /* 	      true */
  /* 	      ) ; */

  /* minf.reinit(mpi_communicator, */
  /* 	      num_vertex, */
  /* 	      n_mpi_processes, */
  /* 	      lvs, */
  /* 	      g, */
  /* 	      0, */
  /* 	      false */
  /* 	      ) ; */


  //minf.reinit(mpi_communicator,num_vertex,lvs,true);


  //std::cout << "Number of vertices actual local:    " << lvs << std::endl;
  //std::cout << "Number of vertices actual:       " << num_vert << std::endl;

  /* for (unsigned int component=0; component< alphadim; ++component){ */
  /*   subdomain_solution[component].compress();} */

  //unsigned int procs = Utilities::MPI::n_mpi_processes(mpi_communicator);
  //unsigned int proc = Utilities::MPI::this_mpi_process(mpi_communicator);
  /* std::vector<unsigned int> COLO_max(procs); */
  /* std::fill(COLO_max.begin(), COLO_max.end(), 0); */
  
  /* std::vector< FEValues<dim>* > hp_fe_values; */
  
  /* UpdateFlags updateflags=  update_values | update_gradients  */
  /*   | update_quadrature_points | update_JxW_values; */
  
  /* hp_fe_values.push_back (new FEValues<dim>(*(fe_collection[0]),  */
  /* 					    *(quadrature_collection[0]), updateflags)); */
  /* //int COLO_min = 100000000; */
  /* //int COLO = 0; */
  /* typename DoFHandler<dim>::active_cell_iterator  */
  /*   cell = dof_handler[0]->begin_active(),  */
  /*   endc = dof_handler[0]->end(); */
  /* for(; */
  /*     cell!=endc; */
  /*     ++cell	) */
  /*   if (cell->is_ghost()  ) */
  /*     { */


  /* 	unsigned int CO = cell->index(); */
  /* 	std::cout << "CO = " << CO <<  ", this mpi_proc = " << Utilities::MPI::this_mpi_process(mpi_communicator)  << std::endl; */
  /* 	//COLO = COLO+1; */


  /* 	COLO_max[proc] = std::max(CO,COLO_max[proc]); */
  /* 	//COLO_min = std::min(CO,COLO_min); */

  /*     } */

  /* delete hp_fe_values[0]; */
   
 
  /* //  COLO = COLO_max; */
  /* std::cout << "COLO = " << COLO_max[proc]   << ", this mpi_proc = " << proc  << std::endl; */
  /* std::cout << "mpi_procs = " << Utilities::MPI::n_mpi_processes(mpi_communicator)  << std::endl; */
  /* //unsigned int procs = Utilities::MPI::n_mpi_processes(mpi_communicator); */
  /* std::cout << "procs = " << procs  << std::endl; */
  //matpatch.resize(procs);
  //for (unsigned int proc=0; proc< procs; ++proc){
  // mapinfo[component].resize(triangulation.n_global_active_cells());
  //matpatch[proc].resize(COLO_max[proc]);
  //matpatch.resize(triangulation.n_global_active_cells());
  //}

  /* std::vector< FEValues<dim>* > hp_fe_values; */
  
  /* UpdateFlags updateflags=  update_values | update_gradients  */
  /*   | update_quadrature_points | update_JxW_values; */
  
  /* hp_fe_values.push_back (new FEValues<dim>(*(fe_collection[2]),  */
  /* 					    *(quadrature_collection[2]), updateflags)); */
  /* //int COLO_min = 100000000; */
  /* //int COLO = 0; */
  /* typename DoFHandler<dim>::cell_iterator  */
  /*   cell = dof_handler[2]->begin_active(),  */
  /*   endc = dof_handler[2]->end(); */
  /* for(; */
  /*     cell!=endc; */
  /*     ++cell	) */
  /*   if (cell->is_locally_owned() ) */
  /*     { */
  /* 	hp_fe_values[component]->reinit (cell); */
  /* 	const FEValues<dim> &fe_values = hp_fe_values[component]->get_present_fe_values (); */
  /* 	const FiniteElement<dim>& fe = cell->get_fe(); */
  /* 	const unsigned int   dofs_per_cell = fe.dofs_per_cell; */
  /* 	glob_dof_ind = dofs_per_cell; */
  /* 	matpatch[CO].localM12 = FullMatrix<double>(dofs_per_cell,dofs_per_cell); */
	  
  /*     } */
  /*   else{ */


  /* 	unsigned int CO = cell->index(); */
  /* 	/\* hp_fe_values[0]->reinit (cell); *\/ */
  /* 	/\* const FiniteElement<dim>& fe = cell->get_fe(); *\/ */
  /* 	/\* const unsigned int   dofs_per_cell = fe.dofs_per_cell; *\/ */
  /* 	//	unsigned int ghostcell_index = cell->neighbor_index(face_num); */
  /* 	//unsigned int proc = Utilities::MPI::this_mpi_process(mpi_communicator);	 */
  /* 	//std::cout << "CO = " << CO <<  ", this mpi_proc = " << Utilities::MPI::this_mpi_process(mpi_communicator)  << std::endl; */
  /* 	//COLO = COLO+1; */
	
  /* 	matpatch[CO].localM12 = FullMatrix<double>(glob_dof_ind,glob_dof_ind); */
	
  /* 	//COLO_max[proc] = std::max(CO,COLO_max[proc]); */
  /* 	//COLO_min = std::min(CO,COLO_min); */
	
  /*   } */

  /* delete hp_fe_values[0]; */

 
  //std::cout << "mpi_procs = " << Utilities::MPI::n_mpi_processes(mpi_communicator)  << std::endl;



