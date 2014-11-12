/* This assembles the basic outer loop of the solver */
template <int dim>
void arcOn<dim>::assemble_system()
{

  for (unsigned int component=0; component< alphadim; ++component){
    naive_subdomain_solution[component] = 0.0;
    subdomain_solution[component] = naive_subdomain_solution[component];
    cont_output1[component] = 0.0;
    cont_global[component] = cont_output1[component];
  }
  
  if(loading_solution == true){

    for (unsigned int component=0; component< alphadim; ++component){

      parallel::distributed::SolutionTransfer<dim,PETScWrappers::MPI::BlockVector> soltrans(*dof_handler[component]);
      
      if(component==0){
	std::string filename = "dat0";
	triangulation.load(filename.c_str());
	soltrans.deserialize(naive_subdomain_solution[component]);
      }
      if(component==1){
	std::string filename = "dat1";
	triangulation.load(filename.c_str());
	soltrans.deserialize(naive_subdomain_solution[component]);
      }
      if(component==2){
	std::string filename = "dat2";
	triangulation.load(filename.c_str());
	soltrans.deserialize(naive_subdomain_solution[component]);
      }
      
      subdomain_solution[component] = naive_subdomain_solution[component];
      
    }
     
    init_solution = subdomain_solution;
    load_top(L2_error_interpolant);
  }
  
  else{

  init_flag = 1;
  load_initial_conditions(subdomain_solution,init_flag);
  //init_flag = 2;
  //load_initial_conditions(subdomain_solution,init_flag);
  //load_top(L2_error_interpolant);

  init_solution = subdomain_solution;

  /* for(unsigned int k=0; k<alphadim; k++){ */
  /*   naive_subdomain_solution[k]= subdomain_solution[k]; */
  /*   FETools::interpolate (*(dof_handler[k]),  naive_subdomain_solution[k], */
  /* 			  *(tdof_handler[k]), cont_output1[k]); */
    
  /*   cont_global[k] = cont_output1[k]; */
  /* } */
  

  }


  bool dyn_adaptation = false;

  if (!dyn_adaptation){
    
    static_periodicity_map(subdomain_solution,0.0);
    if (solver_type == 0)
      {assemble_stiffness(subdomain_solution,0.0);}
    else
      {assemble_cont_stiffness(subdomain_solution,0.0);}
    
    //poisson_matrix.block(0,0).vmult(subdomain_solution[1].block(0),subdomain_solution[2].block(0));
    
    //glob_min_ribbon_density = 0.0;
    
    periodicity_map(subdomain_solution,0.0);
    assemble_sigma(subdomain_solution,0.0,0.0);
    periodicity_map(subdomain_solution,0.0);
    
    //calc_poisson(subdomain_solution,0.0);
    //periodicity_map(subdomain_solution,0.0);
  }

  else {

    static_periodicity_map(subdomain_solution,0.0);
    periodicity_map(subdomain_solution,0.0);

    //Vector<double> subdomain_solution_holder (dof_handler.size());
    subdomain_solution_holder = subdomain_solution;
    
    Vector<float> alpha_measure(triangulation.n_active_cells());
    
    unsigned int qud_pts = (unsigned int)std::ceil(3.0*(double)degree/2.0 - 1.0/2.0);
    
    std::vector< bool > mask_component(dim+1);
    
    mask_component[0] = true; 
    mask_component[1] = false; 
    mask_component[2] = false; 
      
    KellyErrorEstimator<dim>::estimate (*dof_handler[0],
					QGauss<dim-1>(qud_pts),
					typename FunctionMap<dim>::type(),
					subdomain_solution[0],
					alpha_measure
					//mask_component  //Use component mask = 1,0,0 fe_collection[1].
					);   
      
      
    parallel::distributed::GridRefinement::
      refine_and_coarsen_fixed_number (triangulation,
				       alpha_measure,
				       0.3, 0.03);
      
    triangulation.prepare_coarsening_and_refinement();
      
    //Vector of soltrans
      
    parallel::distributed::SolutionTransfer<dim,PETScWrappers::MPI::BlockVector> soltrans0(*dof_handler[0]);
    parallel::distributed::SolutionTransfer<dim,PETScWrappers::MPI::BlockVector> soltrans1(*dof_handler[1]);
    parallel::distributed::SolutionTransfer<dim,PETScWrappers::MPI::BlockVector> soltrans2(*dof_handler[2]);
      
    soltrans0.prepare_for_coarsening_and_refinement(subdomain_solution_holder[0]);
    soltrans1.prepare_for_coarsening_and_refinement(subdomain_solution_holder[1]);
    soltrans2.prepare_for_coarsening_and_refinement(subdomain_solution_holder[2]);
      
    triangulation.execute_coarsening_and_refinement();
      
    //dof_handler[component]->distribute_dofs(*(fe_collection[component]));
    // run setup

    recreate_boundary_data();

    setup_system();
    matrixmapper();
    create_dg_periodicity();
      
    for (unsigned int component=0; component< alphadim; ++component){
	
      naive_subdomain_solution[component] = subdomain_solution_holder[component];
	
      //dof_handler[component]->clear();
	
    }
      
    soltrans0.interpolate(naive_subdomain_solution[0]);
    soltrans1.interpolate(naive_subdomain_solution[1]);
    soltrans2.interpolate(naive_subdomain_solution[2]);
      
    subdomain_solution = naive_subdomain_solution;
 
    static_periodicity_map(subdomain_solution,0.0);
    if (solver_type == 0)
      {
	assemble_stiffness(subdomain_solution,0.0);
      }
    else
      {
	assemble_cont_stiffness(subdomain_solution,0.0);
      }
    periodicity_map(subdomain_solution,0.0);
    assemble_sigma(subdomain_solution,0.0,0.0);
    periodicity_map(subdomain_solution,0.0);
      
  }


  fast_react = false;
  fast_dif = false;
  gvisc_temp1 = 0.0;
  gvisc_temp2 = 0.0;

  if (fstep < 1){ fstep = 0;}
  
  //revert_density(subdomain_solution, 0, 0, naive_revert_output);
  //revert_output=naive_revert_output;
  output_results(0);

  naive_subdomain_solution[2] = 0.0;
  subdomain_solution[2] = naive_subdomain_solution[2];
  
  pcout << "\n\tThere are " <<  triangulation.n_global_active_cells() << 
    " finite element cells. There are " << n_mpi_processes << 
    " processors in use." << std::endl;
  pcout << "\n\tThere are " << alphadim << 
    " components, with degree \'" << degree << "\' basis polynomials, making " 
	<< dof_handler[0]->n_dofs()*alphadim << " total degrees of freedom." 
	<< std::endl;
  pcout << "\033[1;32m\n __________________________________________________________________________________________________________" 
	<< std::endl;


  pcout << "\033[1;32m\n\n ____/\033[1;35marcOn initialized\033[1;32m\\____ "  << std::endl;
  pcout << "\033[1;34m ___________________________ "  << std::endl;

  double current_time_s = 0;
  double current_time = 0;
  int stepper = 0;
  int out_count = 0;
 
  for(unsigned int step =0; step<nstep;++step){
    
    if (step > 0){
      //dt = h_min_dist / ((2.0*degree+1.0)*(CFL_bound) );
      dt = CFL_bound;
      //pcout << "h_min_dist2= " << h_min_dist << ", CFL_bound2 = " << CFL_bound << std::endl;
      //pcout << "dt (" << step << " ) = " << h_min_dist / ((2.0*degree+1.0)*(CFL_bound)) << std::endl;
      //pcout << "Max Jump in domain =" << gmax_jump <<std::endl;
      pcout << "Timestep = " << CFL_bound <<  ", Step(" << stepper << ") = " << current_time_s << std::endl;
    }

    current_time_s = current_time_s + dt;

    if (current_time_s >= hardend) break;

    stepper = stepper + 1;

    if ( (step+1) % (modulus) == 0){
      //stepper = stepper + 1;
      pcout << "\033[1;32m\n  Time (" << stepper << ") = " << current_time_s << "; dt(" << stepper << ") = " << dt << ";" << std::endl;
      pcout << "\033[1;34m ___________________________\n" << std::endl;
    }

    //out_count = out_count + 1;
    //pcout << "timestep = " << out_count << std::endl;

    //clock_t start_clock = clock();
    
    /*Step over reaction/conv-convdiff modes if linear*/
    if (RK_order < 3){

      for(unsigned int fast_nstep = 0; fast_nstep<=fstep; fast_nstep++){
	
	current_time = current_time_s+fast_nstep*fast_dt;
	
	/*Step over fast modes*/

	if (fast_dif == true){
	  
 	  calc_convdiff(subdomain_solution,fast_dt,current_time); 
 
	}
	
	if (fast_react == true){
	  
 	  calc_reaction(subdomain_solution,fast_dt,current_time); 
 
	}	

      }

      /*Compose with slow modes*/

      if (fast_dif == false){

	calc_convdiff(subdomain_solution,dt,current_time_s);

      }

      if (fast_react == false){

	//periodicity_map(subdomain_solution,dt);
	calc_reaction(subdomain_solution,dt,current_time_s);

      }      

    }
    

    /*Step over reaction/convdiff modes if super-linear*/
    if (RK_order > 2){

      /*If any fast modes occur*/
      if (fast_react == true || fast_dif == true){ 

	/*Calculate first fast mode of Strang Splitting*/

	for(unsigned int fast_nstep = 0; fast_nstep<=fstep/2; fast_nstep++){
	  
	  current_time = current_time_s+fast_nstep*fast_dt;
	  
	  /*Step over fast modes*/
	  
	  if (fast_dif == true){
	    
	    calc_convdiff(subdomain_solution,fast_dt,current_time);

	  }
	  
	  if (fast_react == true){
	    
	    calc_reaction(subdomain_solution,fast_dt,current_time);  

	  }	
	  
	}
	
	/*Compose with slow modes*/
	
	if (fast_dif == false){

	  //calc_convdiff(subdomain_solution,dt,current_time_s);

	}
	
	if (fast_react == false){
	  
	  calc_reaction(subdomain_solution,dt,current_time_s);

	}
	
	/*Compose with second Strang Splitting*/
	
	for(unsigned int fast_nstep = 0; fast_nstep<=fstep/2; fast_nstep++){
	  
	  current_time = current_time_s+fast_nstep*fast_dt;
	  
	  /*Step over fast modes*/
	  
	  if (fast_dif == true){
	    
	    calc_convdiff(subdomain_solution,fast_dt,current_time);

	  }
	  
	  if (fast_react == true){
	    
	    calc_reaction(subdomain_solution,fast_dt,current_time);  

	  }	
	}
      }
	      
      if (fast_react == false && fast_dif == false){

      	calc_reaction(subdomain_solution,dt/2.0,current_time_s);
	
      	calc_convdiff(subdomain_solution,dt,current_time_s);

      	calc_reaction(subdomain_solution,dt/2.0,current_time_s);

      }
    }

    periodicity_map(subdomain_solution,dt);

    if (solver_type == 0)
      {
	calc_poisson(subdomain_solution,dt);
      }
    else
      {
	calc_poisson_cont(subdomain_solution,dt);
      }
    
    revert_vacuum(subdomain_solution,dt,current_time_s);

    periodicity_map(subdomain_solution,dt);

    assemble_sigma(subdomain_solution,current_time_s,dt);

    periodicity_map(subdomain_solution,dt);


    //assemble_sigma(subdomain_solution,current_time_s,dt);

    //periodicity_map(subdomain_solution,dt);
    //------------------?


    /* compute_l2_error(subdomain_solution,current_time_s); */
    /* //compute_l2_interpolation(subdomain_solution); */

    /* double l2_err = 0.0; */
    /* double linf_err = 0.0; */
    /* //double l2_interp = 0.0; */
    /* l2_err = L2_error_method[1].block(0).l2_norm(); */
    /* linf_err = L2_error_method[1].block(0).linfty_norm(); */
    /* //l2_interp = interpolation_error[2].block(0).l2_norm(); */
    /* if (l2_err != 0.0){ */
    /*   pcout <<  "L2-error of component[" << 1 << "] = "  << std::setw(5) << std::setprecision(15)  << l2_err << std::endl; */
    /*   pcout <<  "Linf-error of component[" << 1 << "] = "  << std::setw(5) << std::setprecision(15)  << linf_err << std::endl; */
    /*   //pcout <<  "L2 interpolation error of component[" << 2 << "] = "  << std::setw(5) << std::setprecision(15)  << l2_interp << std::endl; */
    /* } */

    //cmcmv(subdomain_solution);

    //cm = dt;
    //cmv = cm /dt;

    //solid
    //pcout << std::setprecision(9) << "cm(" << stepper+1 << ") = " << cm << "; t(" << stepper+1 << ") = " << current_time_s << ";" << std::endl;
    //old
    //pcout << "cmv(" << stepper+1 << ") = " << cmv << "; tt(" << stepper+1 << ") = " << current_time_s << ";" << std::endl;

    cm = 0.0;
    cmv = 0.0;
    
    /* std::ofstream cmfile; */
    /* cmfile.open ("cm.txt", std::ios::out | std::ios::app ); */
    /* cmfile << "cm(" << stepper << ") = " << cm << "; dt(" << stepper << ") = " << dt << ";\r\n"; */
    /* cmfile.close(); */ 
    /* std::ofstream cmvfile; */
    /* cmfile.open ("cmv.txt", std::ios::out | std::ios::app ); */
    /* cmfile << "cmv(" << stepper << ") = " << cmv << "; dt(" << stepper << ") = " << dt << ";\r\n"; */
    /* cmvfile.close(); */


    if ( (step+1) % modulus == 0){
      
      //revert_density(subdomain_solution, dt, current_time_s, naive_revert_output);
      //revert_output=naive_revert_output;
      
      output_results(step+1);
      
      if ( (step+1) % 5*modulus == 0){

      for (unsigned int component=0; component< alphadim; ++component){

	if(saving_solution == true){

	  parallel::distributed::SolutionTransfer<dim,PETScWrappers::MPI::BlockVector> soltrans(*dof_handler[component]);
	  
	  soltrans.prepare_serialization (subdomain_solution[component]);
	  
	  if(component==0){
	    std::string filename = "dat0";
	    triangulation.save(filename.c_str());
	  }
	  if(component==1){
	    std::string filename = "dat1";
	    triangulation.save(filename.c_str());
	  }
	  if(component==2){
	    std::string filename = "dat2";
	    triangulation.save(filename.c_str());
	  }
	  
	}

      }
      }

    }
       
  }
}
