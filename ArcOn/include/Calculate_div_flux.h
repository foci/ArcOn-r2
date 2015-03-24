/* This is the template function for computing the diffusive aand convective evolutionary (time-dependent) fluxes */
template <int dim>
void arcOn<dim>::calculate_div_flux(SolutionVector& substep_solution, double delta_t, double current_time, SolutionVector& div_flux_term) 
{ 

  std::vector< FEValues<dim>* > hp_fe_values;
  std::vector< FEFaceValues<dim>* > hp_fe_values_face;
  std::vector< FEFaceValues<dim>* > hp_fe_values_neigh_face;

  //double sigma_penalty = 5.0
  double max_phi_temp = 0.0;
  double yloc_temp;

  double CFL_temp = 0.0;
  double CFL_approx = 100.0;
  double visc_temp1 = -1000.0;
  double visc_temp2 = 1000.0;
  double h_min_loc = 1.0e6;

  double max_jump_temp = 0.0;

  double art_max = 0.0;
  std::vector< double > max_val_loc;
  max_val_loc.resize(alphadim);

  double pi = 3.1415926535897932384626433832795;

  bool strong = true;
  bool semistrong = false;
  bool quasistrong = false;
  double kappa;

  //  double Time_ramp = 1000.0; // ramp time in seconds                                                                                                                                                                                        
  double dramp;
  if (Time_ramp > 1.0e-5){
    dramp  = std::tanh( 2.0*( current_time/Time_ramp ) );
  }
  else{
    dramp = 1.0;
  }

  double beta = beta_parameter*dramp; //5e-4 + dramp*beta_0;
  double temp = temp_parameter; 
 
  //double local_min_ribbon_density = 0.0;

  int s_flag = 0;
  double beta_pen_coeff;
  double sig_pen_coeff = bpen;
  //LDG_beta;
  if (current_time>10.0*delta_t){
    beta_pen_coeff = 0.00;}
  else{  
    beta_pen_coeff = 0.0;
  }
  
  for (unsigned int component=0; component< alphadim; ++component){
    
    UpdateFlags updateflags=  update_values | update_gradients 
      | update_quadrature_points | update_JxW_values;
    hp_fe_values.push_back (new FEValues<dim>(*(fe_collection[component]), 
					      *(quadrature_collection[component]), updateflags | update_hessians));
    hp_fe_values_face.push_back (new FEFaceValues<dim>(*(fe_collection[component]), 
						       *(face_quadrature_collection[component]),  
						       updateflags | update_normal_vectors  | update_hessians));
    hp_fe_values_neigh_face.push_back (new FEFaceValues<dim>(*(fe_collection[component]), 
							     *(face_quadrature_collection[component]),  
							     updateflags | update_normal_vectors | update_hessians));

    /* local_solution[component] = substep_solution[component]; */
    /* for (unsigned int k=0; k< alphadim; ++k){ */
    /*   if(k!=component){ */
    /* 	local_solution[k] = substep_solution[k]; */
    /*   } */
    /* } */

  }
    
  //std::cout << artificial_visc << " , " << e1 << std::endl; 

 

  for (unsigned int component=0; component< alphadim-1; ++component){

    // The artificial viscosity algorithm needs to be adapted for constraint 
    // matrix consistency (damn periodic BCs!)
    if (artificial_visc == true){ 
      //for (unsigned int bl=0; bl< num_blocks; ++bl){
      //l_soln[bl] = subdomain_solution[component].block(bl);
      //i_soln[bl] = 0.0;
      //r_soln[bl] = 0.0;
      //for (unsigned int k=0; k< alphadim; ++k){
      //interpolate_active[k].reinit(  mpi_communicator,ilocally_owned_dofs[component],ilocally_relevant_dofs[component] );
      //interpolate_active[k].compress();
      //subdomain_solution[k].compress();

      //This is the interpolation info
      /* interpolate_base[component] = subdomain_solution[component]; */
      /* FETools::interpolate (*(dof_handler[component]), interpolate_base[component], */
      /* 			    *(idof_handler[component]), interpolate_active[component]); */
      /* FETools::interpolate (*(idof_handler[component]), interpolate_active[component], */
      /* 			    *(dof_handler[component]), proj_solution[component]); */
      
      
      max_val_loc[component] = 0.0;


  //double visc_median = (gvisc_temp1 - gvisc_temp2) / 2.0;

      if (component == 0){ 
	kappa = kappa_density; //0.1; //+ 4.0*gvisc_median + 1e-8; //0.3; 1.9
	s0 = s0_density; //0.4; //1.6 and 1.5 ubove work 2
	e1 = e1_density; //e1 = 3e-1;
      }
      else{
	kappa = kappa_vorticity; //0.9; //+ 4.0*gvisc_median + 1e-8; //0.3; 1.9
	s0 = s0_vorticity; //1.0; //1.6 and 1.5 ubove work 2
	e1 = e1_vorticity; //e1 = 0.0; //1e-1;
      }
      /* FETools::project_dg (*(dof_handler[component]), interpolate_base[component],  */
      /* 		      *(idof_handler[component]), interpolate_active[component]); */
      /* FETools::project_dg (*(idof_handler[component]), interpolate_active[component],  */
      /* 		      *(dof_handler[component]), proj_solution[component]); */
      //}
    }
 
     //}
      
    //std::vector<double> transport_alphas(alphadim,0);
    //Functionals<dim> functionals;

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
          const std::vector<Point<dim> >& quadrature_point = fe_values.get_quadrature_points();

	  FullMatrix<double> interior_div (alphadim, dofs_per_cell);
	  FullMatrix<double> convection_flux (alphadim, dofs_per_cell);
	  FullMatrix<double> diffusion_flux (alphadim, dofs_per_cell);
	  FullMatrix<double> convection_int (alphadim, dofs_per_cell);

	  std::vector<unsigned int> local_dof_indices (dofs_per_cell);

	  for(unsigned int k=0;k<alphadim;k++){
	    prev_soln_alpha[k] =  std::vector<double>(n_q_points);
	    prev_soln_sigma[k] =  std::vector<Tensor<1,dim> >(n_q_points);
	  }

	  const std::vector<double> &JxW = fe_values.get_JxW_values ();
	  double JxWsum = 0;
	  for (unsigned int q=0; q<n_q_points; ++q){
	    JxWsum += JxW[q];
	  }

	  interior_div = 0.0;
	  convection_int = 0.0;

	  fe_values[*(alpha[component])].get_function_values(subdomain_solution[component],
							     prev_soln_alpha[component]);
	  fe_values[*(sigma[component])].get_function_values(subdomain_solution[component],
							     prev_soln_sigma[component]);

          //fe_values[*(sigma[component])].get_function_values(local_solution[component],prev_soln_sigma[component]);

	  for(unsigned int k=0;k<alphadim;k++){ 
	    if(k!=component){
	      fe_values[*(alpha[k])].get_function_values(subdomain_solution[k],prev_soln_alpha[k]);
	      fe_values[*(sigma[k])].get_function_values(subdomain_solution[k],prev_soln_sigma[k]);
	    }
	  }

	  //Vector<double> eps_smooth(alphadim-1);
	  FullMatrix<double> eps_smooth(alphadim-1, n_q_points);
	  eps_smooth = 0.0;

	  //double kappa = 0.0969; // This is the minor tuning parameter

	  //clean version
	  if (artificial_visc == true ){
	    for (unsigned int q=0; q<n_q_points; ++q){
	      if ( std::abs(prev_soln_alpha[component][q]) > 1e-10 ){
		if(component < 2){
		  max_val_loc[component] = std::max( std::abs(prev_soln_alpha[component][q]), max_val_loc[component] );
		  lnSe = std::max(std::abs(prev_soln_sigma[component][q][0]),std::abs(prev_soln_sigma[component][q][1]))/gmax_val_loc[component];

		  visc_temp1 = std::max(visc_temp1, lnSe );
		  visc_temp2 = std::min(visc_temp2, lnSe );
		  e0 = e1; 

		  if ( lnSe >= s0 + kappa){
		    eps_smooth(component,q) = e0;  
		  }
		  else if ( lnSe >= s0 - kappa &&  lnSe < s0 + kappa){
		    eps_smooth(component,q) = e0*std::sin( pi * ( (lnSe - (s0 - kappa)) / (2.0*kappa)) ); //(e0/2.0) * ( 1 + std::sin( pi*(lnSe-s0) / (2.0*kappa) )  );
		  }
		  else{
		  }
		  art_max = std::max(eps_smooth(component,q),art_max);
		}
	      }
	    }
	  }


/* 	  if (artificial_visc == true ){ */
/* 	    //interpolation stuff */
/* 	    //iprev_soln_alpha[component] =  std::vector<double>(n_q_points); */
/* 	    //fe_values[*(ialpha[component])].get_function_values(proj_solution[component],iprev_soln_alpha[component]); */

/* 	    //fe_values[*(sigma[component])].get_function_values(proj_solution[component],iprev_soln_alpha[component]); */
	    
/* 	    //	  if(component < 2){ */
/* 	    //e0 = 0.05; //1e-2*(cell->diameter())/(max_degree); */
/* 	    //notes = (cell->diameter())/(degree); */
/* 	    //s0 = 1/(std::pow(degree,4.0)); */

/* 	    for (unsigned int q=0; q<n_q_points; ++q){ */

/* 	      o_flow = 0.0;  */
/* 	      l2_base = 0.0; */
	      
/* 	      //o_flow  = std::pow(prev_soln_alpha[component][q]-iprev_soln_alpha[component][q],2); */
/* 	      //o_flow  = std::pow(iprev_soln_alpha[component][q],2); */
/* 	      //l2_base = std::pow(prev_soln_alpha[component][q],2); */

/* 	      if ( std::abs(prev_soln_alpha[component][q]) > 1e-10 ){ */

/* 		//lnSe = std::abs(std::log10(o_flow/l2_base)); */
/* 		//lnSe = std::pow(degree,std::abs(o_flow -l2_base)); */
/* 		if(component < 2){ */
		  
/* 		  max_val_loc[component] = std::max( std::abs(prev_soln_alpha[component][q]), max_val_loc[component] ); */
/* 		  //gmax_val_loc[component] = std::max( gmax_val_loc[component], max_val_loc[component] ); */
/* 		  lnSe = std::max(std::abs(prev_soln_sigma[component][q][0]),std::abs(prev_soln_sigma[component][q][1]))/gmax_val_loc[component]; */
		  
/* 		  //std::cout << " first = " << lnSe << ", " << o_flow << ", " << l2_base <<  std::endl; */
/* 		  visc_temp1 = std::max(visc_temp1, lnSe ); */
/* 		  visc_temp2 = std::min(visc_temp2, lnSe ); */
/* 		  //if ( lnSe != 0.0) std::cout << lnSe <<std::endl; */
/* 		  e0 = e1; //\*std::abs(lnSe); */
		  
		  
/* 		  if ( lnSe >= s0 + kappa){ */
/* 		    eps_smooth(component,q) = e0;   */
/* 		    //std::cout << "t1 = " << current_time << ", lnSe = " << lnSe << ", eps = " << eps_smooth(component,q) << ", l2_base = " << l2_base << ", oflow = " << o_flow << std::endl; */
/* 		  } */
/* 		  else if ( lnSe >= s0 - kappa &&  lnSe < s0 + kappa){ */
		    
/* 		    eps_smooth(component,q) = e0*std::sin( pi * ( (lnSe - (s0 - kappa)) / (2.0*kappa)) ); //(e0/2.0) * ( 1 + std::sin( pi*(lnSe-s0) / (2.0*kappa) )  ); */
		    
/* 		    //std::cout << ", eps = " << eps_smooth(component,q) << ", k = " << kappa << ", component = " << component << std::endl; */
/* 		  } */
/* 		  else{ */
/* 		  } */
/* 		  art_max = std::max(eps_smooth(component,q),art_max); */
/* 		  //if ( eps_smooth(component,q) > e0){pcout << "wtf" << std::endl;} */
/* 		} */
/* 	      } */
/* 	    } */
	    
/* 	    //lnSe = std::log10(o_flow/l2_base); */
	    
/* 	    //	    Vector<double> eps_smooth(alphadim-1); */
/* 	    //      eps_smooth = 0.0; */
	    
/* 	    //	  co = co + 1; */
	    
/* /\* 	    if(component < 2){ *\/ */
/* /\* 	      if ( lnSe > s0 + kappa){ *\/ */
/* /\* 		eps_smooth(component) = e0;   *\/ */
/* /\* 		//	      std::cout << "t1 = " << current_time << ", lnSe = " << lnSe << ", eps = " << eps_smooth << ", k= " << kappa << ", guy = " << component << ", co = " << co << std::endl; *\/ */
/* /\* 	      } *\/ */
/* /\* 	      else if ( lnSe <= s0+kappa && lnSe >= s0-kappa ){ *\/ */
/* /\* 		eps_smooth(component) = e0; //(e0/2.0) * ( 1 + std::sin( pi*(lnSe-s0) / (2.0*kappa) )  ); *\/ */
/* /\* 		//std::cout << "t2 = " << current_time << ", lnSe = " << lnSe << ", eps = " << eps_smooth << ", k = " << kappa << ", guy = " << component << ", co = " << co << std::endl; *\/ */
/* /\* 	      } *\/ */
/* /\* 	      else{ *\/ */
/* /\* 	      } *\/ */
/* /\* 	    } *\/ */

/* 	  } */

	    //eps_smooth = 0.02;
	  //	  std::cout << "e0  = " << e0 << ", diamter = " << cell->diameter() << std::endl;
	  //std::cout << "num  = " << o_flow << std::endl;

	  //std::cout << "den  = " << l2_base << std::endl;
	  //std::cout << "s0-kappa  = " << s0-kappa << ", s0+kappa= " << s0+kappa << std::endl;
	  //std::cout << "e0  = " << e0 << ", diamter = " << cell->diameter() << std::endl;

	  //std::cout << "lnSe  = " << lnSe << std::endl;
	  //std::cout << "eps_smooth  = " << eps_smooth << std::endl;
	  //if(component < 2){  
	  for (unsigned int q=0; q<n_q_points; ++q){
	    /* for (unsigned int k=0; k<alphadim; ++k){ */
	    /*   transport_alphas[k] = prev_soln_alpha[k][q]; */
	    /* } */
	    /* functionals.CE_Transport(1.0,transport_alphas); */
	    /* const TableBase<3,double> D = functionals.Ficks(); */
	    for (unsigned int i=0; i<dofs_per_cell; ++i){
	      interior_div(component,i) += (difs(component)+eps_smooth(component,q)) * (fe_values[*(alpha[component])].gradient(i,q)) * (prev_soln_sigma[component][q]) * JxW[q];

	      if(component == 0){
		//this is the chi stuff
		interior_div(component,i) += (difs(component)+eps_smooth(component,q))*(prev_soln_sigma[component][q])* (prev_soln_sigma[component][q])*(fe_values[*(alpha[component])].value(i,q))* JxW[q];

		/* if ( quadrature_point[q][0] > 225.0 &&  quadrature_point[q][0] <= 260.0  ){ */
		/*   // && quadrature_point[q][1] > 0.0 && quadrature_point[q][1] < 5.0 */
		/*   local_min_ribbon_density = std::min(local_min_ribbon_density, prev_soln_alpha[0][q]); */
		/* } */
	      }

	      //	      std::cout << "what is going on1? = " <<  difs(component)+eps_smooth(component,q) << std::endl;
	      
	    }
	  }
	  
	  convection_int = 0.0;
	  if(component == 0){
	    for (unsigned int q=0; q<n_q_points; ++q){
	      for (unsigned int i=0; i<dofs_per_cell; ++i){
		  
		if ( s_flag == 0 ){
		  convection_int(component,i) +=  (fe_values[*(alpha[component])].value(i,q)) *
		    ( ( ( - prev_soln_sigma[2][q])[0] *  (prev_soln_sigma[component][q])[1]
			+  (prev_soln_sigma[2][q])[1] *  (prev_soln_sigma[component][q])[0] )) * JxW[q];
		
		  convection_int(component,i) -=  2*temp*(prev_soln_sigma[0][q])[1] * fe_values[*(alpha[component])].value(i,q)  * JxW[q];

		  convection_int(component,i) -=  2*(prev_soln_sigma[2][q])[1] * fe_values[*(alpha[component])].value(i,q)  * JxW[q]; 
                //convection_int(component,i) +=  (fe_values[*(alpha[component])].value(i,q)) *
		//( - (prev_soln_sigma[2][q])[1] * 1.0/32.5 * 10.0 * std::exp(-quadrature_point[q][0]/32.5) )  * JxW[q];
		}

		else if ( s_flag == 1 ){//strong
		  
                  convection_int(component,i) +=  0.5 * ( (fe_values[*(alpha[component])].gradient(i,q))[1]
		  					  * ( (prev_soln_sigma[2][q])[0] * prev_soln_alpha[component][q]) * JxW[q]
		  					  - (fe_values[*(alpha[component])].gradient(i,q))[0]
		  					  * ( (prev_soln_sigma[2][q])[1] *  prev_soln_alpha[component][q]) * JxW[q] );
		  
                  convection_int(component,i) +=  0.5 * (  (fe_values[*(alpha[component])].gradient(i,q))[0]
		  					   * ( (prev_soln_sigma[component][q])[1] * prev_soln_alpha[2][q]) * JxW[q]
		  					   - (fe_values[*(alpha[component])].gradient(i,q))[1]
		  					   * ( (prev_soln_sigma[component][q])[0] *  prev_soln_alpha[2][q]) * JxW[q] );

		  convection_int(component,i) +=  2 * temp * (fe_values[*(alpha[component])].gradient(i,q))[1] * prev_soln_alpha[0][q] * JxW[q];

		  convection_int(component,i) -=  2*(prev_soln_sigma[2][q])[1] * fe_values[*(alpha[component])].value(i,q)  * JxW[q];

		  /* convection_int(component,i) +=  0.5 * ( (fe_values[*(alpha[component])].gradient(i,q))[1]  */
		  /* 					  * ( 0.0 * prev_soln_alpha[component][q]) * JxW[q]  */
		  /* 					  - (fe_values[*(alpha[component])].gradient(i,q))[0]  */
		  /* 					  * ( 0.3 *  prev_soln_alpha[component][q]) * JxW[q] ); */
		  
                  /* convection_int(component,i) +=  0.5 * (  (fe_values[*(alpha[component])].gradient(i,q))[0]  */
		  /* 					   * ( (prev_soln_sigma[component][q])[1] * prev_soln_alpha[2][q]) * JxW[q]  */
		  /* 					   - (fe_values[*(alpha[component])].gradient(i,q))[1]  */
		  /* 					   * ( (prev_soln_sigma[component][q])[0] *  prev_soln_alpha[2][q]) * JxW[q] ); */
		  
                }
		
                else {
		  
                  convection_int(component,i) +=  (fe_values[*(alpha[component])].gradient(i,q))[1] 
		    * ( (prev_soln_sigma[2][q])[0] * prev_soln_alpha[component][q]) * JxW[q] 
		    -  (fe_values[*(alpha[component])].gradient(i,q))[0] 
		    * ( (prev_soln_sigma[2][q])[1] *  prev_soln_alpha[component][q]) * JxW[q];
		  
                }

	      }
	    }
	  }
	  
	  if(component == 1){
	    for (unsigned int q=0; q<n_q_points; ++q){
	      //const TableBase<3,double> D = functionals.Ficks();
	      for (unsigned int i=0; i<dofs_per_cell; ++i){
		//double beta = 0.04; //5e-4; //1e-6;
		
                if ( s_flag == 0 ){
		  //this is the good one
		  convection_int(component,i) += (fe_values[*(alpha[component])].value(i,q))
		    * ( ( ( -prev_soln_sigma[2][q])[0]
			  * (prev_soln_sigma[component][q])[1]
			  + (prev_soln_sigma[2][q])[1]
			  *  (prev_soln_sigma[component][q])[0] ))* JxW[q];
		  
		//Pressure type term -- good one (remember to remove beta and background)
		convection_int(component,i) -=  beta*(prev_soln_sigma[0][q])[1] 
		  * fe_values[*(alpha[component])].value(i,q)  * JxW[q];
		/* convection_int(component,i) -=  beta*(prev_soln_sigma[0][q])[1]*std::pow(prev_soln_alpha[0][q],-1.0) * fe_values[*(alpha[component])].value(pinfo[component].alpha_dof_index[i],q) * JxW[q]; */

		}

		else if ( s_flag == 1 ){

                  convection_int(component,i) +=  0.5 * ( (fe_values[*(alpha[component])].gradient(i,q))[1] 
							  * ( (prev_soln_sigma[2][q])[0] * prev_soln_alpha[component][q]) * JxW[q] 
							  -  (fe_values[*(alpha[component])].gradient(i,q))[0] 
							  * ( (prev_soln_sigma[2][q])[1] *  prev_soln_alpha[component][q]) * JxW[q] );
		  
                  convection_int(component,i) +=  0.5 * (  (fe_values[*(alpha[component])].gradient(i,q))[0] 
							   * ( (prev_soln_sigma[component][q])[1] * prev_soln_alpha[2][q]) * JxW[q] 
							   -    (fe_values[*(alpha[component])].gradient(i,q))[1] 
							   * ( (prev_soln_sigma[component][q])[0] *  prev_soln_alpha[2][q]) * JxW[q] );
          
		  //pressure term
		  convection_int(component,i) +=  beta * (fe_values[*(alpha[component])].gradient(i,q))[1] * prev_soln_alpha[0][q] * JxW[q];
		  
		
		  
		}
		else{
                  convection_int(component,i) +=  (fe_values[*(alpha[component])].gradient(i,q))[1] 
		    * ( (prev_soln_sigma[2][q])[0] * prev_soln_alpha[component][q]) * JxW[q] 
		    -  (fe_values[*(alpha[component])].gradient(i,q))[0] 
		    * ( (prev_soln_sigma[2][q])[1] *  prev_soln_alpha[component][q]) * JxW[q];

                }

		double sig0x = std::abs((prev_soln_sigma[2][q])[0]) ; //* fe_values[*(alpha[0])].value(i,q) * JxW[q];
		double sig0y = ( std::abs((prev_soln_sigma[2][q])[1])+std::abs(beta*(prev_soln_sigma[0][q])[1]) ) ; //* fe_values[*(alpha[1])].value(i,q) * JxW[q];

                CFL_temp = std::max(CFL_temp,std::max( sig0x, sig0y ) );

	      }
	    }
	  }
	
	  diffusion_flux = 0;
	  convection_flux = 0;
	  for(unsigned int face_num=0; face_num < GeometryInfo<dim>::faces_per_cell; ++face_num){
	    typename DoFHandler<dim>::face_iterator face=cell->face(face_num);
	    //Point<dim> LDG_beta = Point<dim>(0.0,0.0);
	    //Normal boundaries
	    if ( face->at_boundary() && (face->boundary_indicator() == 0) )
	      {
		hp_fe_values_face[component]->reinit (cell,face_num);
		const FEFaceValues<dim>& fe_values_face = hp_fe_values_face[component]->get_present_fe_values ();
		const Quadrature<dim-1>& face_quadrature_formula = fe_values_face.get_quadrature();
		const unsigned int n_q_points_face = face_quadrature_formula.size();

		const std::vector<double> &JxW_face = fe_values_face.get_JxW_values ();
		const std::vector<Point<dim> > &normals = fe_values_face.get_normal_vectors ();   

		//for(unsigned int k=0;k<alphadim;k++){
		prev_soln_alpha_face[component] =  std::vector<double>(n_q_points_face);
		prev_soln_sigma_face[component] =  std::vector<Tensor<1,dim> >(n_q_points_face);
		prev_soln_alpha_face[0] =  std::vector<double>(n_q_points_face);
		prev_soln_sigma_face[0] =  std::vector<Tensor<1,dim> >(n_q_points_face);
		prev_soln_alpha_face[2] =  std::vector<double>(n_q_points_face);
		prev_soln_sigma_face[2] =  std::vector<Tensor<1,dim> >(n_q_points_face);
		  //}
		fe_values_face[*(alpha[component])].get_function_values(subdomain_solution[component],prev_soln_alpha_face[component]);
		fe_values_face[*(sigma[component])].get_function_values(subdomain_solution[component],prev_soln_sigma_face[component]);

		if (component !=0) {
		  fe_values_face[*(alpha[2])].get_function_values(subdomain_solution[2],prev_soln_alpha_face[2]);
		  fe_values_face[*(sigma[2])].get_function_values(subdomain_solution[2],prev_soln_sigma_face[2]);
		}
		if (component !=2) {
		  fe_values_face[*(alpha[0])].get_function_values(subdomain_solution[0],prev_soln_alpha_face[0]);
		  fe_values_face[*(sigma[0])].get_function_values(subdomain_solution[0],prev_soln_sigma_face[0]);
		}
		
		/* if (component == 0){ */
		/*   fe_values_face[*(alpha[2])].get_function_laplacians(subdomain_solution[2], */
		/* 						 prev_soln_alpha_face[1]); */
		/* } */


		std::vector<double> alphas_boundary(n_q_points_face);
		std::vector<double> alphas_boundary2(n_q_points_face);
		std::vector<double> alphas_boundary0(n_q_points_face);
		//std::vector<double> lap_boundary(n_q_points_face);
		
		std::vector< Tensor< 1, dim > > sigmas_boundary(n_q_points_face);
		std::vector< Tensor< 1, dim > > sigmas_boundary2(n_q_points_face);
		unsigned char boundary_index = face->boundary_indicator();
		const WallBoundaryValues<dim>& wbv = WBV[boundary_index];
		const std::vector< Point<dim> > &quadrature_points = fe_values_face.get_quadrature_points();

		alphas_boundary = prev_soln_alpha_face[component];
		alphas_boundary2 = prev_soln_alpha_face[2];
		alphas_boundary0 = prev_soln_alpha_face[0];
		sigmas_boundary = prev_soln_sigma_face[component];
		sigmas_boundary2 = prev_soln_sigma_face[2];

		//lap_boundary = prev_soln_alpha_face[component];

		/* if (component == 0){ */

		/*   wbv.value_list( quadrature_points, alphas_boundary, component, current_time); */
		/*   wbv.value_list( quadrature_points, alphas_boundary2, 2, current_time); */
		/*   wbv.value_list( quadrature_points, alphas_boundary0, 0, current_time); */
		/*   wbv.gradient_list( quadrature_points, sigmas_boundary, normals, component, current_time); */
		/*   wbv.gradient_list( quadrature_points, sigmas_boundary2, normals, 2, current_time); */

		/* } */
		//else{

		  wbv.value_list( quadrature_points, alphas_boundary, component, current_time);
		  wbv.value_list( quadrature_points, alphas_boundary2, 2, current_time);
		  wbv.value_list( quadrature_points, alphas_boundary0, 0, current_time);
		  wbv.gradient_list( quadrature_points, sigmas_boundary, normals, component, current_time);
		  wbv.gradient_list( quadrature_points, sigmas_boundary2, normals, 2, current_time);

		  //}
		
		FullMatrix<double> eps_smooth_face(alphadim-1, n_q_points_face);
		eps_smooth_face = 0.0;

		if (artificial_visc == true ){
		  for (unsigned int q=0; q<n_q_points_face; ++q){
		    if ( std::abs(prev_soln_alpha_face[component][q]) > 1e-10 ){
		      if(component < 2){
			max_val_loc[component] = std::max( std::abs((prev_soln_alpha_face[component][q]+alphas_boundary[q])/2.0 ), max_val_loc[component] ); 
			lnSe = std::max(std::abs((prev_soln_sigma_face[component][q][0] + sigmas_boundary[q][0])/2.0 ),
					std::abs((prev_soln_sigma_face[component][q][1] + sigmas_boundary[q][1])/2.0 ))/gmax_val_loc[component];
			
			visc_temp1 = std::max(visc_temp1, lnSe );
			visc_temp2 = std::min(visc_temp2, lnSe );
			e0 = e1;
			if ( lnSe >= s0 + kappa){
			  eps_smooth_face(component,q) = e0;  
			}
			else if ( lnSe >= s0 - kappa &&  lnSe < s0 + kappa){
			  eps_smooth_face(component,q) = e0*std::sin( pi * ( (lnSe - (s0 - kappa)) / (2.0*kappa)) ) / 2.0;
			}
			else{
			}
			art_max = std::max(eps_smooth_face(component,q),art_max);
		      }
		    }
		  }

		/* if (artificial_visc == true ){ */

		/*   //iprev_soln_alpha_face[component] =  std::vector<double>(n_q_points_face); */
		/*   //fe_values_face[*(ialpha[component])].get_function_values(proj_solution[component],iprev_soln_alpha_face[component]); */
		  
		/*   for (unsigned int q=0; q<n_q_points_face; ++q){ */
		    
		/*     //o_flow = 0.0; */
		/*     //l2_base = 0.0; */

		/*     //o_flow  = std::pow(iprev_soln_alpha_face[component][q],2); */
		/*     //l2_base = std::pow(prev_soln_alpha_face[component][q],2); */

		/*     if ( std::abs(prev_soln_alpha_face[component][q]) > 1e-10 ){ */

		/*       if(component < 2){ */
			
		/* 	//lnSe = std::pow(degree,std::abs(o_flow -l2_base)); */
			
		/* 	//		      double damn = std::max(std::pow(prev_soln_sigma[component][q][0],2),std::pow(prev_soln_sigma[component][q][1],2)); */
		/* 	//		      lnSe = std::max(std::abs(prev_soln_sigma[component][q][0]),std::abs(prev_soln_sigma[component][q][1])); */
		/* 	//lnSe = std::max(std::abs(prev_soln_sigma[component][q][0]/prev_soln_alpha[component][q]),std::abs(prev_soln_sigma[component][q][1]/prev_soln_alpha[component][q])); */
			
			
		/* 	max_val_loc[component] = std::max( std::abs(prev_soln_alpha_face[component][q]), max_val_loc[component] );  */
		/* 	//gmax_val_loc[component] = std::max( gmax_val_loc[component], max_val_loc[component] ); */
		/* 	lnSe = std::max(std::abs(prev_soln_sigma_face[component][q][0]),std::abs(prev_soln_sigma_face[component][q][1]))/gmax_val_loc[component]; */
			
		/* 	//lnSe = std::pow(degree,damn); */
			
		/* 	//pcout << " second = " << lnSe << ", " << o_flow << ", " << l2_base <<  std::endl; */
		/* 	visc_temp1 = std::max(visc_temp1, lnSe ); */
		/* 	visc_temp2 = std::min(visc_temp2, lnSe ); */
		/* 	e0 = e1; */
			
			
		/* 	if ( lnSe >= s0 + kappa){ */
		/* 	  eps_smooth_face(component,q) = e0;   */
		/* 	  //std::cout << "t1 = " << current_time << ", lnSe = " << lnSe << ", eps = " << eps_smooth(component,q) << ", l2_base = " << l2_base << ", oflow = " << o_flow << std::endl; */
		/* 	} */
		/* 	else if ( lnSe >= s0 - kappa &&  lnSe < s0 + kappa){ */
			  
		/* 	  //eps_smooth(component,q) = e0*std::sin( pi * (lnSe/(2.0*(s0+kappa))) ); //(e0/2.0) * ( 1 + std::sin( pi*(lnSe-s0) / (2.0*kappa) )  ); */
		/* 	  eps_smooth_face(component,q) = e0*std::sin( pi * ( (lnSe - (s0 - kappa)) / (2.0*kappa)) ) / 2.0; */
		/* 	  //std::cout << ", eps = " << eps_smooth(component,q) << ", k = " << kappa << ", component = " << component << std::endl; */
		/* 	} */
		/* 	else{ */
		/* 	} */
		/* 	art_max = std::max(eps_smooth_face(component,q),art_max); */
		/* 	//if ( eps_smooth_face(component,q) > e0){pcout << "wtf" << std::endl;} */
		/*       } */
		/*     } */
		    
		/*     /\* o_flow  = std::pow(prev_soln_alpha_face[component][q]-iprev_soln_alpha_face[component][q],2); *\/ */
		/*     /\* l2_base = std::pow(prev_soln_alpha_face[component][q],2); *\/ */

		/*     /\* if (o_flow > 1e-13 && l2_base > 1e-13){ *\/ */
		      
		/*     /\*   lnSe = std::log10(o_flow/l2_base); *\/ */
		/*     /\*   e0 = e1; //\\*std::abs(lnSe); *\/ */
		/*     /\*   //pcout << " second = " << lnSe << ", " << o_flow << ", " << l2_base <<  std::endl; *\/ */
		  		      
		/*     /\*   //if(component < 2){ *\/ */
		/*     /\* 	if ( lnSe > kappa){ *\/ */
		/*     /\* 	  eps_smooth_face(component,q) = e0;   *\/ */
		/*     /\* 	  //	      std::cout << "t1 = " << current_time << ", lnSe = " << lnSe << ", eps = " << eps_smooth << ", k= " << kappa << ", guy = " << component << ", co = " << co << std::endl; *\/ */
		/*     /\* 	} *\/ */
		/*     /\* 	else if ( lnSe <= kappa && lnSe >= -kappa ){ *\/ */
		/*     /\* 	  eps_smooth_face(component,q) =  e0*std::tanh( (lnSe+kappa) / kappa ); //(e0/2.0) * ( 1 + std::sin( pi*(lnSe-s0) / (2.0*kappa) )  ); *\/ */
		/*     /\* 	  //std::cout << "t2 = " << current_time << ", lnSe = " << lnSe << ", eps = " << eps_smooth << ", k = " << kappa << ", guy = " << component << ", co = " << co << std::endl; *\/ */
		/*     /\* 	} *\/ */
		/*     /\* 	else{ *\/ */
		/*     /\* 	} *\/ */
		/*     /\* 	//} *\/ */
		/*     /\* } *\/ */
		/*   } */


/* 		  } */
		  
/* 		  lnSe = std::log10(o_flow/l2_base); */
/* 		  //eps_smooth = 0.0; */
		  
/* 		  if(component < 2){ */
/* 		    if ( lnSe > s0 + kappa){ */
/* 		      //eps_smooth(component) = e0; */
/* 		      //		    std::cout << "tb1 = " << current_time << ", lnSe = " << lnSe << ", eps = " << eps_smooth << ", k= " << kappa << ", guy = " << component << ", co = " << co << std::endl;                                                         */
/* 		    } */
/* 		    else if ( lnSe <= s0+kappa && lnSe >= s0-kappa ){ */
/* 		      //eps_smooth(component) = e0; //(e0/2.0) * ( 1 + std::sin( pi*(lnSe-s0) / (2.0*kappa) )  ); \ */
/* 		      //std::cout << "tb2 = " << current_time << ", lnSe = " << lnSe << ", eps = " << eps_smooth << ", k = " << kappa << ", guy = " << component << ", co = " << co << std::endl;                          */
/* 		    } */
/* 		    else{ */
/* 		    } */
/* 		  } */
		}
		

		for (unsigned int q=0; q<fe_values_face.n_quadrature_points; ++q){
		  const Point<dim>& normal = normals[q]; 
		  //const TableBase<3,double> D = functionals.Ficks();
		  for (unsigned int i=0; i< dofs_per_cell; ++i){
		    //if( !checksupport || fe.has_support_on_face( pinfo[component].alpha_dof_index[i], face_num ) ){
		    //for(unsigned int l=0; l<dim; ++l){
		    if (component < 2){ 
		      
		      //diffusion_flux(component,i) += (D( TableIndices<3>(component,l,l) ) + eps_test ) *( (prev_soln_sigma_face[component][q])[l] )* fe_values_face[*(alpha[component])].value(pinfo[component].alpha_dof_index[i],q) * normal(l) * JxW_face[q];
		      diffusion_flux(component,i) += ( difs(component) + eps_smooth_face(component,q) )* ( sigmas_boundary[q] ) * (fe_values_face[*(alpha[component])].value(i,q)) * normals[q] * JxW_face[q];
		      //std::cout << "what is going on2? = " <<  difs(component)+eps_smooth_face(component,q) << std::endl;

		      if (component == 0){

			  double norm21 = sigmas_boundary2[q][0] * normal(1);
			  double norm20 = sigmas_boundary2[q][1] * normal(0);
			  double avg_sca0 =  0.5 * ( alphas_boundary[q] + prev_soln_alpha_face[component][q] ) ;
			  
			  double norm01 =  sigmas_boundary[q][0] * normal(1);
			  double norm00 =  sigmas_boundary[q][1] * normal(0);
			  double avg_sca2 =  0.5 * ( alphas_boundary2[q] + prev_soln_alpha_face[2][q]) ;
			
                        if ( s_flag == 1){
			  
                          convection_flux(component,i) -=  0.5 * ( fe_values_face[*(alpha[component])].value(i,q)) * avg_sca0 * (norm21) * JxW_face[q];
			  
                          convection_flux(component,i) -=  0.5 * ( fe_values_face[*(alpha[component])].value(i,q)) * avg_sca2 * (norm00) * JxW_face[q];
			  
                          convection_flux(component,i) +=  0.5 * ( fe_values_face[*(alpha[component])].value(i,q)) * avg_sca0 * (norm20) * JxW_face[q];

                          convection_flux(component,i) +=  0.5 * ( fe_values_face[*(alpha[component])].value(i,q)) * avg_sca2 * (norm01) * JxW_face[q];
			  
			  
                        }
			
        		else if ( s_flag == 2){
			  
                          convection_flux(component,i) -=  ( fe_values_face[*(alpha[component])].value(i,q)) * avg_sca0 * (norm21) * JxW_face[q];
			  
                          convection_flux(component,i) +=  ( fe_values_face[*(alpha[component])].value(i,q)) * avg_sca0 * (norm20) * JxW_face[q];
			  
                        }
                      }

                     if (component == 1){

		       //double beta = 0.04; //5e-4;//6e-3; //1e-6;                                                                                                                                                                                
		       
		       double norm21 = sigmas_boundary2[q][0] * normal(1);
		       double norm20 = sigmas_boundary2[q][1] * normal(0);
		       double avg_sca1 =  0.5 * (alphas_boundary[q] + prev_soln_alpha_face[component][q]);
		       
		       double norm01 =  sigmas_boundary[q][0] * normal(1);
		       double norm00 =  sigmas_boundary[q][1] * normal(0);
		       double avg_sca2 =  0.5 * ( alphas_boundary2[q] + prev_soln_alpha_face[2][q]);
		       
		       if ( s_flag == 1 ){
			 
			 convection_flux(component,i) -=  0.5 * ( fe_values_face[*(alpha[component])].value(i,q)) * avg_sca1 * (norm21) * JxW_face[q];
			 
			 convection_flux(component,i) -=  0.5 * ( fe_values_face[*(alpha[component])].value(i,q)) * avg_sca2 * (norm00) * JxW_face[q];
			 
			 convection_flux(component,i) +=  0.5 * ( fe_values_face[*(alpha[component])].value(i,q)) * avg_sca1 * (norm20) * JxW_face[q];
			 
			 convection_flux(component,i) +=  0.5 * ( fe_values_face[*(alpha[component])].value(i,q)) * avg_sca2 * (norm01) * JxW_face[q];
			 
			 convection_flux(component,i) -= beta*alphas_boundary0[q] * normal(1) * ( fe_values_face[*(alpha[component])].value(i,q)) *  JxW_face[q];
			 
		       }
		       
		       if ( s_flag == 2){
			 
			 convection_flux(component,i) -=  ( fe_values_face[*(alpha[component])].value(i,q)) * avg_sca1 * (norm21) * JxW_face[q];
			 
			 convection_flux(component,i) +=  ( fe_values_face[*(alpha[component])].value(i,q)) * avg_sca1 * (norm20) * JxW_face[q];
			 
			 
		       }
		       
		     }
		     
		     //std::cout << "dif2 = " << D( TableIndices<3>(component,l,l) ) << std::endl;
		     //}
		      //}
		     
		    }
		  }
		}
	      }
	    else if ( face->at_boundary() && (face->boundary_indicator() == 15) )
	      {
		hp_fe_values_face[component]->reinit (cell,face_num);
		const FEFaceValues<dim>& fe_values_face = hp_fe_values_face[component]->get_present_fe_values ();
		const Quadrature<dim-1>& face_quadrature_formula = fe_values_face.get_quadrature();
		const unsigned int n_q_points_face = face_quadrature_formula.size();

		const std::vector<double> &JxW_face = fe_values_face.get_JxW_values ();
		const std::vector<Point<dim> > &normals = fe_values_face.get_normal_vectors ();   

		//for(unsigned int k=0;k<alphadim;k++){
		prev_soln_alpha_face[component] =  std::vector<double>(n_q_points_face);
		prev_soln_sigma_face[component] =  std::vector<Tensor<1,dim> >(n_q_points_face);
		prev_soln_alpha_face[0] =  std::vector<double>(n_q_points_face);
		prev_soln_sigma_face[0] =  std::vector<Tensor<1,dim> >(n_q_points_face);
		prev_soln_alpha_face[2] =  std::vector<double>(n_q_points_face);
		prev_soln_sigma_face[2] =  std::vector<Tensor<1,dim> >(n_q_points_face);
		  //}
		fe_values_face[*(alpha[component])].get_function_values(subdomain_solution[component],prev_soln_alpha_face[component]);
		fe_values_face[*(sigma[component])].get_function_values(subdomain_solution[component],prev_soln_sigma_face[component]);

		if (component !=0) {
		  fe_values_face[*(alpha[2])].get_function_values(subdomain_solution[2],prev_soln_alpha_face[2]);
		  fe_values_face[*(sigma[2])].get_function_values(subdomain_solution[2],prev_soln_sigma_face[2]);
		}
		if (component !=2) {
		  fe_values_face[*(alpha[0])].get_function_values(subdomain_solution[0],prev_soln_alpha_face[0]);
		  fe_values_face[*(sigma[0])].get_function_values(subdomain_solution[0],prev_soln_sigma_face[0]);
		}
		
		/* if (component == 0){ */
		/*   fe_values_face[*(alpha[2])].get_function_laplacians(subdomain_solution[2], */
		/* 						 prev_soln_alpha_face[1]); */
		/* } */


		std::vector<double> alphas_boundary(n_q_points_face);
		std::vector<double> alphas_boundary2(n_q_points_face);
		std::vector<double> alphas_boundary0(n_q_points_face);
		//std::vector<double> lap_boundary(n_q_points_face);
		
		std::vector< Tensor< 1, dim > > sigmas_boundary(n_q_points_face);
		std::vector< Tensor< 1, dim > > sigmas_boundary2(n_q_points_face);
		unsigned char boundary_index = face->boundary_indicator();
		const WallBoundaryValues<dim>& wbv = WBV[boundary_index];
		const std::vector< Point<dim> > &quadrature_points = fe_values_face.get_quadrature_points();

		alphas_boundary = prev_soln_alpha_face[component];
		alphas_boundary2 = prev_soln_alpha_face[2];
		alphas_boundary0 = prev_soln_alpha_face[0];
		sigmas_boundary = prev_soln_sigma_face[component];
		sigmas_boundary2 = prev_soln_sigma_face[2];

		//lap_boundary = prev_soln_alpha_face[component];

		/* if (component == 0){ */

		/*   wbv.value_list( quadrature_points, alphas_boundary, component, current_time); */
		/*   wbv.value_list( quadrature_points, alphas_boundary2, 2, current_time); */
		/*   wbv.value_list( quadrature_points, alphas_boundary0, 0, current_time); */
		/*   wbv.gradient_list( quadrature_points, sigmas_boundary, normals, component, current_time); */
		/*   wbv.gradient_list( quadrature_points, sigmas_boundary2, normals, 2, current_time); */

		/* } */
		//else{

		  wbv.value_list2( quadrature_points, alphas_boundary, component, current_time);
		  wbv.value_list2( quadrature_points, alphas_boundary2, 2, current_time);
		  wbv.value_list2( quadrature_points, alphas_boundary0, 0, current_time);
		  wbv.gradient_list2( quadrature_points, sigmas_boundary, normals, component, current_time);
		  wbv.gradient_list2( quadrature_points, sigmas_boundary2, normals, 2, current_time);

		  //}
		
		FullMatrix<double> eps_smooth_face(alphadim-1, n_q_points_face);
		eps_smooth_face = 0.0;

		if (artificial_visc == true ){
		  for (unsigned int q=0; q<n_q_points_face; ++q){
		    if ( std::abs(prev_soln_alpha_face[component][q]) > 1e-10 ){
		      if(component < 2){
			max_val_loc[component] = std::max( std::abs( (prev_soln_alpha_face[component][q] + alphas_boundary[q])/2.0), max_val_loc[component] ); 
			lnSe = std::max(std::abs((prev_soln_sigma_face[component][q][0]+sigmas_boundary[q][0])/2.0),
					std::abs((prev_soln_sigma_face[component][q][1]+sigmas_boundary[q][1])/2.0))/gmax_val_loc[component];
			visc_temp1 = std::max(visc_temp1, lnSe );
			visc_temp2 = std::min(visc_temp2, lnSe );
			e0 = e1;
			
			if ( lnSe >= s0 + kappa){
			  eps_smooth_face(component,q) = e0;  
			}
			else if ( lnSe >= s0 - kappa &&  lnSe < s0 + kappa){
			  eps_smooth_face(component,q) = e0*std::sin( pi * ( (lnSe - (s0 - kappa)) / (2.0*kappa)) ) / 2.0;
			}
			else{
			}
			art_max = std::max(eps_smooth_face(component,q),art_max);
		      }
		    }
		  }

		/* if (artificial_visc == true ){ */

		/*   //iprev_soln_alpha_face[component] =  std::vector<double>(n_q_points_face); */
		/*   //fe_values_face[*(ialpha[component])].get_function_values(proj_solution[component],iprev_soln_alpha_face[component]); */
		  
		/*   for (unsigned int q=0; q<n_q_points_face; ++q){ */
		    
		/*     //o_flow = 0.0; */
		/*     //l2_base = 0.0; */

		/*     //o_flow  = std::pow(iprev_soln_alpha_face[component][q],2); */
		/*     //l2_base = std::pow(prev_soln_alpha_face[component][q],2); */

		/*     if ( std::abs(prev_soln_alpha_face[component][q]) > 1e-10 ){ */

		/*       if(component < 2){ */
			
		/* 	//lnSe = std::pow(degree,std::abs(o_flow -l2_base)); */
			
		/* 	//		      double damn = std::max(std::pow(prev_soln_sigma[component][q][0],2),std::pow(prev_soln_sigma[component][q][1],2)); */
		/* 	//		      lnSe = std::max(std::abs(prev_soln_sigma[component][q][0]),std::abs(prev_soln_sigma[component][q][1])); */
		/* 	//lnSe = std::max(std::abs(prev_soln_sigma[component][q][0]/prev_soln_alpha[component][q]),std::abs(prev_soln_sigma[component][q][1]/prev_soln_alpha[component][q])); */
			
			
		/* 	max_val_loc[component] = std::max( std::abs(prev_soln_alpha_face[component][q]), max_val_loc[component] );  */
		/* 	//gmax_val_loc[component] = std::max( gmax_val_loc[component], max_val_loc[component] ); */
		/* 	lnSe = std::max(std::abs(prev_soln_sigma_face[component][q][0]),std::abs(prev_soln_sigma_face[component][q][1]))/gmax_val_loc[component]; */
			
		/* 	//lnSe = std::pow(degree,damn); */
			
		/* 	//pcout << " second = " << lnSe << ", " << o_flow << ", " << l2_base <<  std::endl; */
		/* 	visc_temp1 = std::max(visc_temp1, lnSe ); */
		/* 	visc_temp2 = std::min(visc_temp2, lnSe ); */
		/* 	e0 = e1; */
			
			
		/* 	if ( lnSe >= s0 + kappa){ */
		/* 	  eps_smooth_face(component,q) = e0;   */
		/* 	  //std::cout << "t1 = " << current_time << ", lnSe = " << lnSe << ", eps = " << eps_smooth(component,q) << ", l2_base = " << l2_base << ", oflow = " << o_flow << std::endl; */
		/* 	} */
		/* 	else if ( lnSe >= s0 - kappa &&  lnSe < s0 + kappa){ */
			  
		/* 	  //eps_smooth(component,q) = e0*std::sin( pi * (lnSe/(2.0*(s0+kappa))) ); //(e0/2.0) * ( 1 + std::sin( pi*(lnSe-s0) / (2.0*kappa) )  ); */
		/* 	  eps_smooth_face(component,q) = e0*std::sin( pi * ( (lnSe - (s0 - kappa)) / (2.0*kappa)) ) / 2.0; */
		/* 	  //std::cout << ", eps = " << eps_smooth(component,q) << ", k = " << kappa << ", component = " << component << std::endl; */
		/* 	} */
		/* 	else{ */
		/* 	} */
		/* 	art_max = std::max(eps_smooth_face(component,q),art_max); */
		/* 	//if ( eps_smooth_face(component,q) > e0){pcout << "wtf" << std::endl;} */
		/*       } */
		/*     } */
		    
		/*     /\* o_flow  = std::pow(prev_soln_alpha_face[component][q]-iprev_soln_alpha_face[component][q],2); *\/ */
		/*     /\* l2_base = std::pow(prev_soln_alpha_face[component][q],2); *\/ */

		/*     /\* if (o_flow > 1e-13 && l2_base > 1e-13){ *\/ */
		      
		/*     /\*   lnSe = std::log10(o_flow/l2_base); *\/ */
		/*     /\*   e0 = e1; //\\*std::abs(lnSe); *\/ */
		/*     /\*   //pcout << " second = " << lnSe << ", " << o_flow << ", " << l2_base <<  std::endl; *\/ */
		  		      
		/*     /\*   //if(component < 2){ *\/ */
		/*     /\* 	if ( lnSe > kappa){ *\/ */
		/*     /\* 	  eps_smooth_face(component,q) = e0;   *\/ */
		/*     /\* 	  //	      std::cout << "t1 = " << current_time << ", lnSe = " << lnSe << ", eps = " << eps_smooth << ", k= " << kappa << ", guy = " << component << ", co = " << co << std::endl; *\/ */
		/*     /\* 	} *\/ */
		/*     /\* 	else if ( lnSe <= kappa && lnSe >= -kappa ){ *\/ */
		/*     /\* 	  eps_smooth_face(component,q) =  e0*std::tanh( (lnSe+kappa) / kappa ); //(e0/2.0) * ( 1 + std::sin( pi*(lnSe-s0) / (2.0*kappa) )  ); *\/ */
		/*     /\* 	  //std::cout << "t2 = " << current_time << ", lnSe = " << lnSe << ", eps = " << eps_smooth << ", k = " << kappa << ", guy = " << component << ", co = " << co << std::endl; *\/ */
		/*     /\* 	} *\/ */
		/*     /\* 	else{ *\/ */
		/*     /\* 	} *\/ */
		/*     /\* 	//} *\/ */
		/*     /\* } *\/ */
		/*   } */


/* 		  } */
		  
/* 		  lnSe = std::log10(o_flow/l2_base); */
/* 		  //eps_smooth = 0.0; */
		  
/* 		  if(component < 2){ */
/* 		    if ( lnSe > s0 + kappa){ */
/* 		      //eps_smooth(component) = e0; */
/* 		      //		    std::cout << "tb1 = " << current_time << ", lnSe = " << lnSe << ", eps = " << eps_smooth << ", k= " << kappa << ", guy = " << component << ", co = " << co << std::endl;                                                         */
/* 		    } */
/* 		    else if ( lnSe <= s0+kappa && lnSe >= s0-kappa ){ */
/* 		      //eps_smooth(component) = e0; //(e0/2.0) * ( 1 + std::sin( pi*(lnSe-s0) / (2.0*kappa) )  ); \ */
/* 		      //std::cout << "tb2 = " << current_time << ", lnSe = " << lnSe << ", eps = " << eps_smooth << ", k = " << kappa << ", guy = " << component << ", co = " << co << std::endl;                          */
/* 		    } */
/* 		    else{ */
/* 		    } */
/* 		  } */
		}
		

		for (unsigned int q=0; q<fe_values_face.n_quadrature_points; ++q){
		  const Point<dim>& normal = normals[q]; 
		  //const TableBase<3,double> D = functionals.Ficks();
		  for (unsigned int i=0; i< dofs_per_cell; ++i){
		    //if( !checksupport || fe.has_support_on_face( pinfo[component].alpha_dof_index[i], face_num ) ){
		    //for(unsigned int l=0; l<dim; ++l){
		    if (component < 2){ 
		      
		      //diffusion_flux(component,i) += (D( TableIndices<3>(component,l,l) ) + eps_test ) *( (prev_soln_sigma_face[component][q])[l] )* fe_values_face[*(alpha[component])].value(pinfo[component].alpha_dof_index[i],q) * normal(l) * JxW_face[q];
		      diffusion_flux(component,i) += ( difs(component) + eps_smooth_face(component,q) )* ( sigmas_boundary[q] ) * (fe_values_face[*(alpha[component])].value(i,q)) * normals[q] * JxW_face[q];
		      //std::cout << "what is going on2? = " <<  difs(component)+eps_smooth_face(component,q) << std::endl;

		      if (component == 0){

			  double norm21 = sigmas_boundary2[q][0] * normal(1);
			  double norm20 = sigmas_boundary2[q][1] * normal(0);
			  double avg_sca0 =  0.5 * ( alphas_boundary[q] + prev_soln_alpha_face[component][q] ) ;
			  
			  double norm01 =  sigmas_boundary[q][0] * normal(1);
			  double norm00 =  sigmas_boundary[q][1] * normal(0);
			  double avg_sca2 =  0.5 * ( alphas_boundary2[q] + prev_soln_alpha_face[2][q]) ;
			
                        if ( s_flag == 1){
			  
                          convection_flux(component,i) -=  0.5 * ( fe_values_face[*(alpha[component])].value(i,q)) * avg_sca0 * (norm21) * JxW_face[q];
			  
                          convection_flux(component,i) -=  0.5 * ( fe_values_face[*(alpha[component])].value(i,q)) * avg_sca2 * (norm00) * JxW_face[q];
			  
                          convection_flux(component,i) +=  0.5 * ( fe_values_face[*(alpha[component])].value(i,q)) * avg_sca0 * (norm20) * JxW_face[q];

                          convection_flux(component,i) +=  0.5 * ( fe_values_face[*(alpha[component])].value(i,q)) * avg_sca2 * (norm01) * JxW_face[q];
			  
			  
                        }
			
        		else if ( s_flag == 2){
			  
                          convection_flux(component,i) -=  ( fe_values_face[*(alpha[component])].value(i,q)) * avg_sca0 * (norm21) * JxW_face[q];
			  
                          convection_flux(component,i) +=  ( fe_values_face[*(alpha[component])].value(i,q)) * avg_sca0 * (norm20) * JxW_face[q];
			  
                        }
                      }

                     if (component == 1){

		       //double beta = 0.04; //5e-4;//6e-3; //1e-6;                                                                                                                                                                                
		       
		       double norm21 = sigmas_boundary2[q][0] * normal(1);
		       double norm20 = sigmas_boundary2[q][1] * normal(0);
		       double avg_sca1 =  0.5 * (alphas_boundary[q] + prev_soln_alpha_face[component][q]);
		       
		       double norm01 =  sigmas_boundary[q][0] * normal(1);
		       double norm00 =  sigmas_boundary[q][1] * normal(0);
		       double avg_sca2 =  0.5 * ( alphas_boundary2[q] + prev_soln_alpha_face[2][q]);
		       
		       if ( s_flag == 1 ){
			 
			 convection_flux(component,i) -=  0.5 * ( fe_values_face[*(alpha[component])].value(i,q)) * avg_sca1 * (norm21) * JxW_face[q];
			 
			 convection_flux(component,i) -=  0.5 * ( fe_values_face[*(alpha[component])].value(i,q)) * avg_sca2 * (norm00) * JxW_face[q];
			 
			 convection_flux(component,i) +=  0.5 * ( fe_values_face[*(alpha[component])].value(i,q)) * avg_sca1 * (norm20) * JxW_face[q];
			 
			 convection_flux(component,i) +=  0.5 * ( fe_values_face[*(alpha[component])].value(i,q)) * avg_sca2 * (norm01) * JxW_face[q];
			 
			 convection_flux(component,i) -= beta*alphas_boundary0[q] * normal(1) * ( fe_values_face[*(alpha[component])].value(i,q)) *  JxW_face[q];
			 
		       }
		       
		       if ( s_flag == 2){
			 
			 convection_flux(component,i) -=  ( fe_values_face[*(alpha[component])].value(i,q)) * avg_sca1 * (norm21) * JxW_face[q];
			 
			 convection_flux(component,i) +=  ( fe_values_face[*(alpha[component])].value(i,q)) * avg_sca1 * (norm20) * JxW_face[q];
			 
			 
		       }
		       
		     }
		     
		     //std::cout << "dif2 = " << D( TableIndices<3>(component,l,l) ) << std::endl;
		     //}
		      //}
		     
		    }
		  }
		}
	      }
	    // Periodic case
	    else if (  face->at_boundary() && ( face->boundary_indicator() == 3 ||
						face->boundary_indicator() == 4))
	      {

		hp_fe_values_face[component]->reinit(cell, face_num);
		const FEFaceValues<dim>&	   fe_values_face = hp_fe_values_face[component]->get_present_fe_values ();

		const std::vector<double> &JxW_face = fe_values_face.get_JxW_values ();
		const std::vector<Point<dim> > &normals = fe_values_face.get_normal_vectors ();   

		const Quadrature<dim-1>& face_quadrature_formula = fe_values_face.get_quadrature();
		const std::vector<Point<dim> >& quadrature_point = fe_values_face.get_quadrature_points();
		const unsigned int n_q_points_face = face_quadrature_formula.size();

		for(unsigned int k=0;k<alphadim;k++){
		  prev_soln_alpha_face[k] =  std::vector<double>(n_q_points_face);
		  prev_soln_sigma_face[k] = std::vector<Tensor<1,dim> >(n_q_points_face);

		  fe_values_face[*(alpha[component])].get_function_values(subdomain_solution[k],prev_soln_alpha_face[k]);
		  fe_values_face[*(sigma[component])].get_function_values(subdomain_solution[k],prev_soln_sigma_face[k]);
		  
		}

		FullMatrix<double> eps_smooth_face(alphadim-1, n_q_points_face);
		eps_smooth_face = 0.0;

		if (artificial_visc == true ){
		  for (unsigned int q=0; q<n_q_points_face; ++q){
		    if (  std::abs(prev_soln_alpha_face[component][q]) > 1e-10 ){

		      if(component < 2){
			max_val_loc[component] = std::max( (std::abs(prev_soln_alpha_face[component][q] + soln_drawX[component][CO][q] )/2.0 ), max_val_loc[component] ); 
			lnSe = std::max(std::abs( ( prev_soln_sigma_face[component][q][0]+(grad_drawX[component][CO][q])[0])/2.0),
					std::abs( ( prev_soln_sigma_face[component][q][1]+(grad_drawX[component][CO][q])[1])/2.0))/gmax_val_loc[component];		      
			visc_temp1 = std::max(visc_temp1, lnSe );
			visc_temp2 = std::min(visc_temp2, lnSe );
			e0 = e1;
			
			if ( lnSe >= s0 + kappa){
			  eps_smooth_face(component,q) = e0;  
			}
			else if ( lnSe >= s0 - kappa &&  lnSe < s0 + kappa){
			  eps_smooth_face(component,q) = e0*std::sin( pi * ( (lnSe - (s0 - kappa)) / (2.0*kappa)) ) ;	  
			}
			else{
			}
			art_max = std::max(eps_smooth(component,q),art_max);
			if ( eps_smooth_face(component,q) > e0){pcout << "A.D. has problems!!" << std::endl;}
		      }

		    }
		  }
		}

		/* if (artificial_visc == true ){ */
		/*   //iprev_soln_alpha_face[component] =  std::vector<double>(n_q_points_face); */
		/*   //fe_values_face[*(ialpha[component])].get_function_values(proj_solution[component],iprev_soln_alpha_face[component]); */
		
		/*   for (unsigned int q=0; q<n_q_points_face; ++q){ */
		    
		/*     //o_flow = 0.0; */
		/*     //l2_base = 0.0; */
		  
		/*     //o_flow  = std::pow(iprev_soln_alpha[component][q],2); */
		/*     //l2_base = std::pow(prev_soln_alpha[component][q],2); */
		    
		/*     if (  std::abs(prev_soln_alpha_face[component][q]) > 1e-10 ){ */

		/*       if(component < 2){ */
			
		/* 	//lnSe = std::pow(degree,std::abs(o_flow -l2_base)); */
		/* 	//		      double damn = std::max(std::pow(prev_soln_sigma[component][q][0],2),std::pow(prev_soln_sigma[component][q][1],2)); */
		/* 	//lnSe = std::max(std::abs(prev_soln_sigma[component][q][0]/prev_soln_alpha[component][q]),std::abs(prev_soln_sigma[component][q][1]/prev_soln_alpha[component][q])); */
			
		/* 	max_val_loc[component] = std::max( std::abs(prev_soln_alpha_face[component][q]), max_val_loc[component] );  */
		/* 	//gmax_val_loc[component] = std::max( gmax_val_loc[component], max_val_loc[component] ); */
		/* 	lnSe = std::max(std::abs(prev_soln_sigma_face[component][q][0]),std::abs(prev_soln_sigma_face[component][q][1]))/gmax_val_loc[component];		       */
		/* 	//lnSe = std::max(std::abs(prev_soln_sigma[component][q][0]),std::abs(prev_soln_sigma[component][q][1])); */
		/* 	//lnSe = std::pow(degree,damn); */
			
		/* 	visc_temp1 = std::max(visc_temp1, lnSe ); */
		/* 	visc_temp2 = std::min(visc_temp2, lnSe ); */
		/* 	//pcout << " third = " << lnSe << ", " << o_flow << ", " << l2_base <<  std::endl; */
		/* 	e0 = e1; */
			
			
		/* 	if ( lnSe >= s0 + kappa){ */
		/* 	  eps_smooth_face(component,q) = e0;   */
		/* 	  //std::cout << "t1 = " << current_time << ", lnSe = " << lnSe << ", eps = " << eps_smooth(component,q) << ", l2_base = " << l2_base << ", oflow = " << o_flow << std::endl; */
		/* 	} */
		/* 	else if ( lnSe >= s0 - kappa &&  lnSe < s0 + kappa){ */
			  
		/* 	  //eps_smooth(component,q) = e0*std::sin( pi * (lnSe/(2.0*(s0+kappa))) ); //(e0/2.0) * ( 1 + std::sin( pi*(lnSe-s0) / (2.0*kappa) )  ); */
		/* 	  eps_smooth_face(component,q) = e0*std::sin( pi * ( (lnSe - (s0 - kappa)) / (2.0*kappa)) ) ;	   */
		/* 	  //std::cout << ", eps = " << eps_smooth(component,q) << ", k = " << kappa << ", component = " << component << std::endl; */
		/* 	} */
		/* 	else{ */
		/* 	} */
		/* 	art_max = std::max(eps_smooth(component,q),art_max); */
		/* 	if ( eps_smooth_face(component,q) > e0){pcout << "wtf" << std::endl;} */
		/*       } */

		/*     } */
	
		/*    /\*  o_flow  = std::pow(prev_soln_alpha_face[component][q]-iprev_soln_alpha_face[component][q],2); *\/ */
		/*   /\*   l2_base = std::pow(prev_soln_alpha_face[component][q],2); *\/ */
		  
		/*   /\*   if (o_flow > 1e-13 && l2_base > 1e-13){ *\/ */

		/*   /\*     lnSe = std::log10(o_flow/l2_base); *\/ */
		/*   /\*     e0 = e1; //\\*std::abs(lnSe); *\/ */
		/*   /\*     //pcout << " third = " << lnSe << ", " << o_flow << ", " << l2_base <<  std::endl; *\/ */
		/*   /\*     //if(component < 2){ *\/ */
		/*   /\*     if ( lnSe > kappa){ *\/ */
		/*   /\* 	eps_smooth_face(component,q) = e0;   *\/ */
		/*   /\* 	//	      std::cout << "t1 = " << current_time << ", lnSe = " << lnSe << ", eps = " << eps_smooth << ", k= " << kappa << ", guy = " << component << ", co = " << co << std::endl; *\/ */
		/*   /\*     } *\/ */
		/*   /\*     else if ( lnSe <= kappa && lnSe >= -kappa ){ *\/ */
		/*   /\* 	eps_smooth_face(component,q) =   e0*std::tanh( (lnSe+kappa) / kappa ); //(e0/2.0) * ( 1 + std::sin( pi*(lnSe-s0) / (2.0*kappa) )  ); *\/ */
		/*   /\* 	//std::cout << "t2 = " << current_time << ", lnSe = " << lnSe << ", eps = " << eps_smooth << ", k = " << kappa << ", guy = " << component << ", co = " << co << std::endl; *\/ */
		/*   /\*     } *\/ */
		/*   /\*     else{ *\/ */
		/*   /\*     } *\/ */
		/*   /\*     //} *\/ */
		/*   /\*   } *\/ */
		/*   /\* } *\/ */
		/*   } */
		/* } */
		//double sigma_penalty = sig_pen_coeff*( degree + 1.0 )*( degree + 1.0 )*( face->measure())/(cell->measure() );
		double sigma_penalty = sig_pen_coeff;
		//std::cout << "sig_pen = " << sigma_penalty << std::endl;
		//std::cout << "area = " << cell->measure() << std::endl;
		//double sigma_penalty = -(10.0+current_time*0.002)*degree*( degree + 1.0 ) * ( face->measure())/(cell->measure() );
		for (unsigned int q=0; q<n_q_points_face; ++q){
		  const Point<dim>& normal = normals[q];
		  for (unsigned int i=0; i<dofs_per_cell; ++i){
		    
		    // std::cout << "quad point periodic cell [ "<< cell->index() << "] of neighbor [" << Cell_match[cell->index()]<< "] = " <<  quadrature_point[q] << std::endl;

		    //if( !checksupport || fe.has_support_on_face( pinfo[component].alpha_dof_index[i], face_num ) ){
		    //for(unsigned int l=0; l<dim; ++l){
		    if (component < 2){
		      
		      /* if (normal(0)>1e-6){ LDG_beta(0) = 0.5/normal(0); LDG_beta(1) = 0.0; */
		      /* } */
		      /* else if (normal(1)>1e-6) {LDG_beta(0)=0.0; LDG_beta(1) = 0.5/normal(1); */
		      /* } */
		      //LDG_beta(0)=0.0; LDG_beta(1)=0.0;

		      //sigma_penalty = 
		      
		      //std::cout << "penx = " <<  sigma_penalty << ", jumpx = " <<  prev_soln_alpha_face[component][q] - soln_drawX[component][CO][q] << std::endl;

		      diffusion_flux(component,i) += ( difs(component)+eps_smooth_face(component,q) ) 
			* (
			   (0.5 * ( (prev_soln_sigma_face[component][q]) + (grad_drawX[component][CO][q]) )
			    * fe_values_face[*(alpha[component])].value(i,q) * normals[q]) 
			   - sigma_penalty * ( ( prev_soln_alpha_face[component][q] - soln_drawX[component][CO][q] )
					       * fe_values_face[*(alpha[component])].value(i,q)  )
			   - beta_pen_coeff * ( ( prev_soln_sigma_face[component][q]*normals[q] - grad_drawX[component][CO][q]*normals[q] )
						* fe_values_face[*(alpha[component])].value(i,q) )
			   ) * JxW_face[q];
		      
		      //std::cout << "what is going on3? = " <<  difs(component)+eps_smooth_face(component,q) << std::endl;

		      if (component == 0){
			
	    		double norm21 = 0.5 * (  (grad_drawX[2][CO][q])[0] + (prev_soln_sigma_face[2][q])[0] ) * normal(1);
	    		double norm20 = 0.5 * (  (grad_drawX[2][CO][q])[1] + (prev_soln_sigma_face[2][q])[1] ) * normal(0);
	    		//double norm21 =  0.0 * normal(1);
			//double norm20 =  0.3 * normal(0);

	    		double avg_sca0 =  0.5 * (prev_soln_alpha_face[component][q] +  soln_drawX[component][CO][q]);
			
	    		double norm01 = 0.5 * (  (grad_drawX[component][CO][q])[0] + (prev_soln_sigma_face[component][q])[0] ) * normal(1);
	    		double norm00 = 0.5 * (  (grad_drawX[component][CO][q])[1] + (prev_soln_sigma_face[component][q])[1] ) * normal(0);
	    		double avg_sca2 =  0.5 * (prev_soln_alpha_face[2][q] + soln_drawX[2][CO][q]);

			//pcout << "prev_soln[" << q << "] = " << prev_soln_alpha_face[2][q] << std::endl;
			//pcout << "draw[" << q << "] = " << soln_draw[2][CO][q] << std::endl;

			//CFL_temp = std::max(CFL_temp,std::abs(avg_sca0) );
			//CFL_temp = std::max(CFL_temp,std::abs(avg_sca2) );
			
	    		if ( s_flag == 1 ){
			  
	    		  convection_flux(component,i) -=  0.5 * ( fe_values_face[*(alpha[component])].value(i,q)) * avg_sca0 * (norm21) * JxW_face[q];
			  
	    		  convection_flux(component,i) -=  0.5 * ( fe_values_face[*(alpha[component])].value(i,q)) * avg_sca2 * (norm00) * JxW_face[q];
			  
	    		  convection_flux(component,i) +=  0.5 * ( fe_values_face[*(alpha[component])].value(i,q)) * avg_sca0 * (norm20) * JxW_face[q];
			  
	    		  convection_flux(component,i) +=  0.5 * ( fe_values_face[*(alpha[component])].value(i,q)) * avg_sca2 * (norm01) * JxW_face[q];

			  
	    		}
			
	    		if ( s_flag == 2 ){
			  
	    		  convection_flux(component,i) -=  ( fe_values_face[*(alpha[component])].value(i,q)) * avg_sca0 * (norm21) * JxW_face[q];
			  
	    		  convection_flux(component,i) +=  ( fe_values_face[*(alpha[component])].value(i,q)) * avg_sca0 * (norm20) * JxW_face[q];
			  
	    		}
	    	      }

	    	      if (component == 1){
			
	    		//double beta = 0.04; //5e-4; //6e-3; //1e-6;
			
	    		double norm21 = 0.5 * (  (grad_drawX[2][CO][q])[0] + (prev_soln_sigma_face[2][q])[0] ) * normal(1);
	    		double norm20 = 0.5 * (  (grad_drawX[2][CO][q])[1] + (prev_soln_sigma_face[2][q])[1] ) * normal(0);
	    		double avg_sca1 =  0.5*(prev_soln_alpha_face[component][q] +  soln_drawX[component][CO][q]);
			
	    		double norm01 = 0.5 * (  (grad_drawX[component][CO][q])[0] + (prev_soln_sigma_face[component][q])[0] ) * normal(1);
	    		double norm00 = 0.5 * (  (grad_drawX[component][CO][q])[1] + (prev_soln_sigma_face[component][q])[1] ) * normal(0);
	    		double avg_sca2 =  0.5*(prev_soln_alpha_face[2][q] + soln_drawX[2][CO][q]);

			double avg_sca0 = 0.5*(prev_soln_alpha_face[0][q] + soln_drawX[0][CO][q]);

	    		/* double norm21 = 0.5 * ( (prev_soln_sigma_neigh_face[2][q])[0] + (prev_soln_sigma_face[2][q])[0] ) * normal(1); */
	    		/* double norm20 = 0.5 * ( (prev_soln_sigma_neigh_face[2][q])[1] + (prev_soln_sigma_face[2][q])[1] ) * normal(0); */
	    		/* double avg_sca1 =  0.5*(prev_soln_alpha_face[component][q] + prev_soln_alpha_neigh_face[component][q]); */
			
	    		/* double norm01 = 0.5 * ( (prev_soln_sigma_neigh_face[component][q])[0] + (prev_soln_sigma_face[component][q])[0] ) * normal(1); */
	    		/* double norm00 = 0.5 * ( (prev_soln_sigma_neigh_face[component][q])[1] + (prev_soln_sigma_face[component][q])[1] ) * normal(0); */
	    		/* double avg_sca2 =  0.5*(prev_soln_alpha_face[2][q] + prev_soln_alpha_neigh_face[2][q]); */
			
	    		if ( s_flag == 1 ){
			  
	    		  convection_flux(component,i) -=  0.5 * ( fe_values_face[*(alpha[component])].value(i,q)) * avg_sca1 * (norm21) * JxW_face[q];
			  
	    		  convection_flux(component,i) -=  0.5 * ( fe_values_face[*(alpha[component])].value(i,q)) * avg_sca2 * (norm00) * JxW_face[q];
			  
	    		  convection_flux(component,i) +=  0.5 * ( fe_values_face[*(alpha[component])].value(i,q)) * avg_sca1 * (norm20) * JxW_face[q];
			  
	    		  convection_flux(component,i) +=  0.5 * ( fe_values_face[*(alpha[component])].value(i,q)) * avg_sca2 * (norm01) * JxW_face[q];

			  convection_flux(component,i) -=  beta * avg_sca0 * normal(1) * ( fe_values_face[*(alpha[component])].value(i,q)) *  JxW_face[q];

			  
	    		}
			
	    		if( s_flag == 2){
			  
	    		  convection_flux(component,i) -=  ( fe_values_face[*(alpha[component])].value(i,q)) * avg_sca1 * (norm21) * JxW_face[q];
			  
	    		  convection_flux(component,i) +=  ( fe_values_face[*(alpha[component])].value(i,q)) * avg_sca1 * (norm20) * JxW_face[q];
	    		}
			
	    	      }

		    }
		  }
		}
	   
	      }
	    // Periodic case
	    else if (  face->at_boundary() && ( face->boundary_indicator() == 1 ||
						face->boundary_indicator() == 2))
	      {

		hp_fe_values_face[component]->reinit(cell, face_num);
		const FEFaceValues<dim>&	   fe_values_face = hp_fe_values_face[component]->get_present_fe_values ();

		const std::vector<double> &JxW_face = fe_values_face.get_JxW_values ();
		const std::vector<Point<dim> > &normals = fe_values_face.get_normal_vectors ();   

		const Quadrature<dim-1>& face_quadrature_formula = fe_values_face.get_quadrature();
		const std::vector<Point<dim> >& quadrature_point = fe_values_face.get_quadrature_points();
		const unsigned int n_q_points_face = face_quadrature_formula.size();

		for(unsigned int k=0;k<alphadim;k++){
		  prev_soln_alpha_face[k] =  std::vector<double>(n_q_points_face);
		  prev_soln_sigma_face[k] = std::vector<Tensor<1,dim> >(n_q_points_face);
		  
		  fe_values_face[*(alpha[component])].get_function_values(subdomain_solution[k],prev_soln_alpha_face[k]);
		  fe_values_face[*(sigma[component])].get_function_values(subdomain_solution[k],prev_soln_sigma_face[k]);
		  
		}

		FullMatrix<double> eps_smooth_face(alphadim-1, n_q_points_face);
		eps_smooth_face = 0.0;

		if (artificial_visc == true ){
		  for (unsigned int q=0; q<n_q_points_face; ++q){
		    if (  std::abs(prev_soln_alpha_face[component][q]) > 1e-10 ){

		      if(component < 2){
			max_val_loc[component] = std::max( (std::abs(prev_soln_alpha_face[component][q] + soln_drawY[component][CO][q] )/2.0 ), max_val_loc[component] ); 
			lnSe = std::max(std::abs( ( prev_soln_sigma_face[component][q][0]+(grad_drawY[component][CO][q])[0])/2.0),
					std::abs( ( prev_soln_sigma_face[component][q][1]+(grad_drawY[component][CO][q])[1])/2.0))/gmax_val_loc[component];		      
			visc_temp1 = std::max(visc_temp1, lnSe );
			visc_temp2 = std::min(visc_temp2, lnSe );
			e0 = e1;
			
			if ( lnSe >= s0 + kappa){
			  eps_smooth_face(component,q) = e0;  
			}
			else if ( lnSe >= s0 - kappa &&  lnSe < s0 + kappa){
			  eps_smooth_face(component,q) = e0*std::sin( pi * ( (lnSe - (s0 - kappa)) / (2.0*kappa)) ) ;	  
			}
			else{
			}
			art_max = std::max(eps_smooth(component,q),art_max);
			if ( eps_smooth_face(component,q) > e0){pcout << "wtf" << std::endl;}
		      }

		    }
		  }
		}

		/* if (artificial_visc == true ){ */
		/*   for (unsigned int q=0; q<n_q_points_face; ++q){ */
		/*     if (  std::abs(prev_soln_alpha_face[component][q]) > 1e-10 ){ */
		/*       if(component < 2){ */
		/* 	max_val_loc[component] = std::max( std::abs(prev_soln_alpha_face[component][q]), max_val_loc[component] );  */
		/* 	lnSe = std::max(std::abs(prev_soln_sigma_face[component][q][0]),std::abs(prev_soln_sigma_face[component][q][1]))/gmax_val_loc[component];		       */
		/* 	visc_temp1 = std::max(visc_temp1, lnSe ); */
		/* 	visc_temp2 = std::min(visc_temp2, lnSe ); */
		/* 	e0 = e1; */

		/* 	if ( lnSe >= s0 + kappa){ */
		/* 	  eps_smooth_face(component,q) = e0;   */
		/* 	} */
		/* 	else if ( lnSe >= s0 - kappa &&  lnSe < s0 + kappa){ */
		/* 	  eps_smooth_face(component,q) = e0*std::sin( pi * ( (lnSe - (s0 - kappa)) / (2.0*kappa)) ) ;	   */
		/* 	} */
		/* 	else{ */
		/* 	} */
		/* 	art_max = std::max(eps_smooth(component,q),art_max); */
		/* 	if ( eps_smooth_face(component,q) > e0){pcout << "wtf" << std::endl;} */
		/*       } */
		/*     } */
		/*   } */
		/* } */

		/* if (artificial_visc == true ){ */
		/*   //iprev_soln_alpha_face[component] =  std::vector<double>(n_q_points_face); */
		/*   //fe_values_face[*(ialpha[component])].get_function_values(proj_solution[component],iprev_soln_alpha_face[component]); */
		
		/*   for (unsigned int q=0; q<n_q_points_face; ++q){ */
		    
		/*     //o_flow = 0.0; */
		/*     //l2_base = 0.0; */
		  
		/*     //o_flow  = std::pow(iprev_soln_alpha[component][q],2); */
		/*     //l2_base = std::pow(prev_soln_alpha[component][q],2); */
		    
		/*     if (  std::abs(prev_soln_alpha_face[component][q]) > 1e-10 ){ */

		/*       if(component < 2){ */
			
		/* 	//lnSe = std::pow(degree,std::abs(o_flow -l2_base)); */
		/* 	//		      double damn = std::max(std::pow(prev_soln_sigma[component][q][0],2),std::pow(prev_soln_sigma[component][q][1],2)); */
		/* 	//lnSe = std::max(std::abs(prev_soln_sigma[component][q][0]/prev_soln_alpha[component][q]),std::abs(prev_soln_sigma[component][q][1]/prev_soln_alpha[component][q])); */
			
		/* 	max_val_loc[component] = std::max( std::abs(prev_soln_alpha_face[component][q]), max_val_loc[component] );  */
		/* 	//gmax_val_loc[component] = std::max( gmax_val_loc[component], max_val_loc[component] ); */
		/* 	lnSe = std::max(std::abs(prev_soln_sigma_face[component][q][0]),std::abs(prev_soln_sigma_face[component][q][1]))/gmax_val_loc[component];		       */
		/* 	//lnSe = std::max(std::abs(prev_soln_sigma[component][q][0]),std::abs(prev_soln_sigma[component][q][1])); */
		/* 	//lnSe = std::pow(degree,damn); */
			
		/* 	visc_temp1 = std::max(visc_temp1, lnSe ); */
		/* 	visc_temp2 = std::min(visc_temp2, lnSe ); */
		/* 	//pcout << " third = " << lnSe << ", " << o_flow << ", " << l2_base <<  std::endl; */
		/* 	e0 = e1; */
			
			
		/* 	if ( lnSe >= s0 + kappa){ */
		/* 	  eps_smooth_face(component,q) = e0;   */
		/* 	  //std::cout << "t1 = " << current_time << ", lnSe = " << lnSe << ", eps = " << eps_smooth(component,q) << ", l2_base = " << l2_base << ", oflow = " << o_flow << std::endl; */
		/* 	} */
		/* 	else if ( lnSe >= s0 - kappa &&  lnSe < s0 + kappa){ */
			  
		/* 	  //eps_smooth(component,q) = e0*std::sin( pi * (lnSe/(2.0*(s0+kappa))) ); //(e0/2.0) * ( 1 + std::sin( pi*(lnSe-s0) / (2.0*kappa) )  ); */
		/* 	  eps_smooth_face(component,q) = e0*std::sin( pi * ( (lnSe - (s0 - kappa)) / (2.0*kappa)) ) ;	   */
		/* 	  //std::cout << ", eps = " << eps_smooth(component,q) << ", k = " << kappa << ", component = " << component << std::endl; */
		/* 	} */
		/* 	else{ */
		/* 	} */
		/* 	art_max = std::max(eps_smooth(component,q),art_max); */
		/* 	if ( eps_smooth_face(component,q) > e0){pcout << "wtf" << std::endl;} */
		/*       } */

		/*     } */
	
		/*    /\*  o_flow  = std::pow(prev_soln_alpha_face[component][q]-iprev_soln_alpha_face[component][q],2); *\/ */
		/*   /\*   l2_base = std::pow(prev_soln_alpha_face[component][q],2); *\/ */
		  
		/*   /\*   if (o_flow > 1e-13 && l2_base > 1e-13){ *\/ */

		/*   /\*     lnSe = std::log10(o_flow/l2_base); *\/ */
		/*   /\*     e0 = e1; //\\*std::abs(lnSe); *\/ */
		/*   /\*     //pcout << " third = " << lnSe << ", " << o_flow << ", " << l2_base <<  std::endl; *\/ */
		/*   /\*     //if(component < 2){ *\/ */
		/*   /\*     if ( lnSe > kappa){ *\/ */
		/*   /\* 	eps_smooth_face(component,q) = e0;   *\/ */
		/*   /\* 	//	      std::cout << "t1 = " << current_time << ", lnSe = " << lnSe << ", eps = " << eps_smooth << ", k= " << kappa << ", guy = " << component << ", co = " << co << std::endl; *\/ */
		/*   /\*     } *\/ */
		/*   /\*     else if ( lnSe <= kappa && lnSe >= -kappa ){ *\/ */
		/*   /\* 	eps_smooth_face(component,q) =   e0*std::tanh( (lnSe+kappa) / kappa ); //(e0/2.0) * ( 1 + std::sin( pi*(lnSe-s0) / (2.0*kappa) )  ); *\/ */
		/*   /\* 	//std::cout << "t2 = " << current_time << ", lnSe = " << lnSe << ", eps = " << eps_smooth << ", k = " << kappa << ", guy = " << component << ", co = " << co << std::endl; *\/ */
		/*   /\*     } *\/ */
		/*   /\*     else{ *\/ */
		/*   /\*     } *\/ */
		/*   /\*     //} *\/ */
		/*   /\*   } *\/ */
		/*   /\* } *\/ */
		/*   } */
		/* } */
		//double sigma_penalty = sig_pen_coeff*( degree + 1.0 )*( degree + 1.0 )*( face->measure())/(cell->measure() );
		double sigma_penalty = sig_pen_coeff;
		//double sigma_penalty = -(10.0+current_time*0.002)*degree*( degree + 1.0 ) * ( face->measure())/(cell->measure() );
		for (unsigned int q=0; q<n_q_points_face; ++q){
		  const Point<dim>& normal = normals[q];
		  for (unsigned int i=0; i<dofs_per_cell; ++i){
		    
		    // std::cout << "quad point periodic cell [ "<< cell->index() << "] of neighbor [" << Cell_match[cell->index()]<< "] = " <<  quadrature_point[q] << std::endl;

		    //if( !checksupport || fe.has_support_on_face( pinfo[component].alpha_dof_index[i], face_num ) ){
		    //for(unsigned int l=0; l<dim; ++l){
		    if (component < 2){
			
		      /* if (normal(0)>1e-6){ LDG_beta(0) = 0.5/normal(0); LDG_beta(1) = 0.0; */
		      /* } */
		      /* else if (normal(1)>1e-6) {LDG_beta(0)=0.0; LDG_beta(1) = 0.5/normal(1); */
		      /* } */
		      //LDG_beta(0)=0.0; LDG_beta(1)=0.0;
		      
		      diffusion_flux(component,i) += ( difs(component)+eps_smooth_face(component,q) ) 
			* (
			   (0.5 * ( (prev_soln_sigma_face[component][q]) + (grad_drawY[component][CO][q]) )
			    * fe_values_face[*(alpha[component])].value(i,q) * normals[q]) 
			   - sigma_penalty * ( ( prev_soln_alpha_face[component][q] - soln_drawY[component][CO][q] )
					       * fe_values_face[*(alpha[component])].value(i,q) ) 
			   - beta_pen_coeff * ( ( prev_soln_sigma_face[component][q]*normals[q] - grad_drawY[component][CO][q]*normals[q] )
						* fe_values_face[*(alpha[component])].value(i,q) )
			   ) * JxW_face[q];
		      
		      //std::cout << "what is going on3? = " <<  difs(component)+eps_smooth_face(component,q) << std::endl;

		      if (component == 0){
			
	    		double norm21 = 0.5 * (  (grad_drawY[2][CO][q])[0] + (prev_soln_sigma_face[2][q])[0] ) * normal(1);
	    		double norm20 = 0.5 * (  (grad_drawY[2][CO][q])[1] + (prev_soln_sigma_face[2][q])[1] ) * normal(0);
	    		//double norm21 =  0.0 * normal(1);
			//double norm20 =  0.3 * normal(0);

	    		double avg_sca0 =  0.5 * (prev_soln_alpha_face[component][q] +  soln_drawY[component][CO][q]);
			
	    		double norm01 = 0.5 * (  (grad_drawY[component][CO][q])[0] + (prev_soln_sigma_face[component][q])[0] ) * normal(1);
	    		double norm00 = 0.5 * (  (grad_drawY[component][CO][q])[1] + (prev_soln_sigma_face[component][q])[1] ) * normal(0);
	    		double avg_sca2 =  0.5 * (prev_soln_alpha_face[2][q] + soln_drawY[2][CO][q]);

			//pcout << "prev_soln[" << q << "] = " << prev_soln_alpha_face[2][q] << std::endl;
			//pcout << "draw[" << q << "] = " << soln_draw[2][CO][q] << std::endl;

			//CFL_temp = std::max(CFL_temp,std::abs(avg_sca0) );
			//CFL_temp = std::max(CFL_temp,std::abs(avg_sca2) );
			
	    		if ( s_flag == 1 ){
			  
	    		  convection_flux(component,i) -=  0.5 * ( fe_values_face[*(alpha[component])].value(i,q)) * avg_sca0 * (norm21) * JxW_face[q];
			  
	    		  convection_flux(component,i) -=  0.5 * ( fe_values_face[*(alpha[component])].value(i,q)) * avg_sca2 * (norm00) * JxW_face[q];
			  
	    		  convection_flux(component,i) +=  0.5 * ( fe_values_face[*(alpha[component])].value(i,q)) * avg_sca0 * (norm20) * JxW_face[q];
			  
	    		  convection_flux(component,i) +=  0.5 * ( fe_values_face[*(alpha[component])].value(i,q)) * avg_sca2 * (norm01) * JxW_face[q];

			  
	    		}
			
	    		if ( s_flag == 2 ){
			  
	    		  convection_flux(component,i) -=  ( fe_values_face[*(alpha[component])].value(i,q)) * avg_sca0 * (norm21) * JxW_face[q];
			  
	    		  convection_flux(component,i) +=  ( fe_values_face[*(alpha[component])].value(i,q)) * avg_sca0 * (norm20) * JxW_face[q];
			  
	    		}
	    	      }

	    	      if (component == 1){
			
	    		//double beta = 0.04; //5e-4; //6e-3; //1e-6;
			
	    		double norm21 = 0.5 * (  (grad_drawY[2][CO][q])[0] + (prev_soln_sigma_face[2][q])[0] ) * normal(1);
	    		double norm20 = 0.5 * (  (grad_drawY[2][CO][q])[1] + (prev_soln_sigma_face[2][q])[1] ) * normal(0);
	    		double avg_sca1 =  0.5*(prev_soln_alpha_face[component][q] +  soln_drawY[component][CO][q]);
			
	    		double norm01 = 0.5 * (  (grad_drawY[component][CO][q])[0] + (prev_soln_sigma_face[component][q])[0] ) * normal(1);
	    		double norm00 = 0.5 * (  (grad_drawY[component][CO][q])[1] + (prev_soln_sigma_face[component][q])[1] ) * normal(0);
	    		double avg_sca2 =  0.5*(prev_soln_alpha_face[2][q] + soln_drawY[2][CO][q]);

			double avg_sca0 = 0.5*(prev_soln_alpha_face[0][q] + soln_drawY[0][CO][q]);

	    		/* double norm21 = 0.5 * ( (prev_soln_sigma_neigh_face[2][q])[0] + (prev_soln_sigma_face[2][q])[0] ) * normal(1); */
	    		/* double norm20 = 0.5 * ( (prev_soln_sigma_neigh_face[2][q])[1] + (prev_soln_sigma_face[2][q])[1] ) * normal(0); */
	    		/* double avg_sca1 =  0.5*(prev_soln_alpha_face[component][q] + prev_soln_alpha_neigh_face[component][q]); */
			
	    		/* double norm01 = 0.5 * ( (prev_soln_sigma_neigh_face[component][q])[0] + (prev_soln_sigma_face[component][q])[0] ) * normal(1); */
	    		/* double norm00 = 0.5 * ( (prev_soln_sigma_neigh_face[component][q])[1] + (prev_soln_sigma_face[component][q])[1] ) * normal(0); */
	    		/* double avg_sca2 =  0.5*(prev_soln_alpha_face[2][q] + prev_soln_alpha_neigh_face[2][q]); */
			
	    		if ( s_flag == 1 ){
			  
	    		  convection_flux(component,i) -=  0.5 * ( fe_values_face[*(alpha[component])].value(i,q)) * avg_sca1 * (norm21) * JxW_face[q];
			  
	    		  convection_flux(component,i) -=  0.5 * ( fe_values_face[*(alpha[component])].value(i,q)) * avg_sca2 * (norm00) * JxW_face[q];
			  
	    		  convection_flux(component,i) +=  0.5 * ( fe_values_face[*(alpha[component])].value(i,q)) * avg_sca1 * (norm20) * JxW_face[q];
			  
	    		  convection_flux(component,i) +=  0.5 * ( fe_values_face[*(alpha[component])].value(i,q)) * avg_sca2 * (norm01) * JxW_face[q];

			  convection_flux(component,i) -=  beta*avg_sca0 * normal(1) * ( fe_values_face[*(alpha[component])].value(i,q)) *  JxW_face[q];

			  
	    		}
			
	    		if( s_flag == 2){
			  
	    		  convection_flux(component,i) -=  ( fe_values_face[*(alpha[component])].value(i,q)) * avg_sca1 * (norm21) * JxW_face[q];
			  
	    		  convection_flux(component,i) +=  ( fe_values_face[*(alpha[component])].value(i,q)) * avg_sca1 * (norm20) * JxW_face[q];
	    		}
			
	    	      }

		    }
		  }
		}
	   
	      }
	    else { //interior edge
	      typename DoFHandler<dim>::cell_iterator neighbor = cell->neighbor(face_num);
	  
	      const unsigned int neighbor2 = cell->neighbor_of_neighbor(face_num);

	      hp_fe_values_face[component]->reinit(cell, face_num);
	      hp_fe_values_neigh_face[component]->reinit(neighbor,neighbor2);

	      const FEFaceValues<dim>& fe_values_face = hp_fe_values_face[component]->get_present_fe_values ();
	      const FEFaceValues<dim>& fe_values_neigh_face = hp_fe_values_neigh_face[component]->get_present_fe_values();

	      const std::vector<double> &JxW_face = fe_values_face.get_JxW_values ();
	      const std::vector<Point<dim> > &normals = fe_values_face.get_normal_vectors ();   

	      const Quadrature<dim-1>& face_quadrature_formula = fe_values_face.get_quadrature();
	      const unsigned int n_q_points_face = face_quadrature_formula.size();


	      for(unsigned int k=0;k<alphadim;k++){
		prev_soln_alpha_face[k] =  std::vector<double>(n_q_points_face);
		prev_soln_sigma_face[k] = std::vector<Tensor<1,dim> >(n_q_points_face);
		
		prev_soln_alpha_neigh_face[k] =  std::vector<double>(n_q_points_face);;
		prev_soln_sigma_neigh_face[k] =  std::vector<Tensor<1,dim> >(n_q_points_face);
	      }

	      fe_values_face[*(alpha[component])].get_function_values(subdomain_solution[component],prev_soln_alpha_face[component]);
	      fe_values_face[*(sigma[component])].get_function_values(subdomain_solution[component],prev_soln_sigma_face[component]);
	      fe_values_neigh_face[*(alpha[component])].get_function_values(subdomain_solution[component],prev_soln_alpha_neigh_face[component]);
	      fe_values_neigh_face[*(sigma[component])].get_function_values(subdomain_solution[component],prev_soln_sigma_neigh_face[component]);

	      fe_values_face[*(alpha[component])].get_function_values(subdomain_solution[2],prev_soln_alpha_face[2]);
	      fe_values_face[*(sigma[2])].get_function_values(subdomain_solution[2],prev_soln_sigma_face[2]);
	      fe_values_neigh_face[*(alpha[component])].get_function_values(subdomain_solution[2],prev_soln_alpha_neigh_face[2]);
	      fe_values_neigh_face[*(sigma[component])].get_function_values(subdomain_solution[2],prev_soln_sigma_neigh_face[2]);

	      if (component !=0){
		fe_values_face[*(alpha[component])].get_function_values(subdomain_solution[0],prev_soln_alpha_face[0]);
		fe_values_face[*(sigma[component])].get_function_values(subdomain_solution[0],prev_soln_sigma_face[0]);
		fe_values_neigh_face[*(alpha[component])].get_function_values(subdomain_solution[0],prev_soln_alpha_neigh_face[0]);
		fe_values_neigh_face[*(sigma[component])].get_function_values(subdomain_solution[0],prev_soln_sigma_neigh_face[0]);
	      }

	      //}

	      FullMatrix<double> eps_smooth_face(alphadim-1, n_q_points_face);
	      eps_smooth_face = 0.0;

	      if (artificial_visc == true ){
		for (unsigned int q=0; q<n_q_points_face; ++q){
		  if ( std::abs(prev_soln_alpha_face[component][q]) > 1e-10 ){
		    if(component < 2){
		      max_val_loc[component] = std::max( (std::abs(prev_soln_alpha_face[component][q] + prev_soln_alpha_neigh_face[component][q] ) /2.0 ), max_val_loc[component] ); 
		      lnSe = std::max(std::abs((prev_soln_sigma_face[component][q][0]+(prev_soln_sigma_neigh_face[component][q])[0] )/2.0),
				      std::abs((prev_soln_sigma_face[component][q][1]+(prev_soln_sigma_neigh_face[component][q])[1] )/2.0))/gmax_val_loc[component];
		      visc_temp1 = std::max(visc_temp1, lnSe );
		      visc_temp2 = std::min(visc_temp2, lnSe );
		      e0 = e1; //*std::abs(lnSe);
		      if ( lnSe >= s0 + kappa){
			eps_smooth_face(component,q) = e0;  
		      }
		      else if ( lnSe >= s0 - kappa &&  lnSe < s0 + kappa){
			eps_smooth_face(component,q) = e0*std::sin( pi * ( (lnSe - (s0 - kappa)) / (2.0*kappa)) ) ;	
		      }
		      else{
		      }
		      art_max = std::max(eps_smooth_face(component,q),art_max);
		      if ( eps_smooth_face(component,q) > e0){pcout << "A.D has problems" << std::endl;}
		    }
		  }
		}
	      }

	      /* if (artificial_visc == true ){ */
	      /* 	//iprev_soln_alpha_face[component] =  std::vector<double>(n_q_points_face); */
	      /* 	//fe_values_face[*(ialpha[component])].get_function_values(proj_solution[component],iprev_soln_alpha_face[component]); */
		
	      /* 	for (unsigned int q=0; q<n_q_points_face; ++q){ */
		    
	      /* 	  //o_flow = 0.0; */
	      /* 	  //l2_base = 0.0; */

	      /* 	  //o_flow  = std::pow(iprev_soln_alpha[component][q],2); */
	      /* 	  //l2_base = std::pow(prev_soln_alpha[component][q],2); */

	      /* 	  if ( std::abs(prev_soln_alpha_face[component][q]) > 1e-10 ){ */

	      /* 	    if(component < 2){ */

	      /* 	    //lnSe = std::abs(std::log10(o_flow/l2_base)); */
	      /* 	    //lnSe = std::pow(degree,std::abs(o_flow -l2_base)); */

	      /* 	      max_val_loc[component] = std::max( std::abs(prev_soln_alpha_face[component][q]), max_val_loc[component] );  */
	      /* 	      //gmax_val_loc[component] = std::max( gmax_val_loc[component], max_val_loc[component] ); */
	      /* 	      lnSe = std::max(std::abs(prev_soln_sigma_face[component][q][0]),std::abs(prev_soln_sigma_face[component][q][1]))/gmax_val_loc[component]; */
	      /* 	      //double damn = std::max(std::pow(prev_soln_sigma[component][q][0],2),std::pow(prev_soln_sigma[component][q][1],2)); */
	      /* 	      //lnSe = std::max(std::abs(prev_soln_sigma[component][q][0]/prev_soln_alpha[component][q]),std::abs(prev_soln_sigma[component][q][1]/prev_soln_alpha[component][q]));		     */
	      /* 	      //lnSe = std::max(std::abs(prev_soln_sigma[component][q][0]),std::abs(prev_soln_sigma[component][q][1])); */
	      /* 	      //lnSe = std::pow(degree,damn); */
	      /* 	      visc_temp1 = std::max(visc_temp1, lnSe ); */
	      /* 	      visc_temp2 = std::min(visc_temp2, lnSe ); */
	      /* 	      //		    pcout << " fourth = " << lnSe << ", " << o_flow << ", " << l2_base <<  std::endl; */
	      /* 	      //if ( lnSe != 0.0) std::cout << lnSe <<std::endl; */
	      /* 	      e0 = e1; //\*std::abs(lnSe); */
		      
		      
	      /* 	      if ( lnSe >= s0 + kappa){ */
	      /* 		eps_smooth_face(component,q) = e0;   */
	      /* 		//std::cout << "t1 = " << current_time << ", lnSe = " << lnSe << ", eps = " << eps_smooth(component,q) << ", l2_base = " << l2_base << ", oflow = " << o_flow << std::endl; */
	      /* 	      } */
	      /* 	      else if ( lnSe >= s0 - kappa &&  lnSe < s0 + kappa){ */
			
	      /* 		//eps_smooth(component,q) = e0*std::sin( pi * (lnSe/(2.0*(s0+kappa))) ); //(e0/2.0) * ( 1 + std::sin( pi*(lnSe-s0) / (2.0*kappa) )  ); */
	      /* 		eps_smooth_face(component,q) = e0*std::sin( pi * ( (lnSe - (s0 - kappa)) / (2.0*kappa)) ) ;	 */
	      /* 		//std::cout << ", eps = " << eps_smooth(component,q) << ", k = " << kappa << ", component = " << component << std::endl; */
	      /* 	      } */
	      /* 	      else{ */
	      /* 	      } */
	      /* 	      art_max = std::max(eps_smooth_face(component,q),art_max); */
	      /* 	      if ( eps_smooth_face(component,q) > e0){pcout << "wtf" << std::endl;} */
	      /* 	    } */
		    
	      /* 	  } */
	      /* 	  /\* o_flow  = std::pow(prev_soln_alpha_face[component][q]-iprev_soln_alpha_face[component][q],2); *\/ */
	      /* 	  /\* l2_base = std::pow(prev_soln_alpha_face[component][q],2); *\/ */
		  
	      /* 	  /\* if (o_flow > 1e-13 && l2_base > 1e-13){ *\/ */

	      /* 	  /\*   lnSe = std::log10(o_flow/l2_base); *\/ */
	      /* 	  /\*   e0 = e1; //\\*std::abs(lnSe); *\/ */

	      /* 	  /\*   //pcout << " fourth = " << lnSe << ", " << o_flow << ", " << l2_base <<  std::endl; *\/ */
		    
	      /* 	  /\*   //if(component < 2){ *\/ */
	      /* 	  /\*     if ( lnSe > kappa){ *\/ */
	      /* 	  /\* 	eps_smooth_face(component,q) = e0;   *\/ */
	      /* 	  /\* 	//	      std::cout << "t1 = " << current_time << ", lnSe = " << lnSe << ", eps = " << eps_smooth << ", k= " << kappa << ", guy = " << component << ", co = " << co << std::endl; *\/ */
	      /* 	  /\*     } *\/ */
	      /* 	  /\*     else if ( lnSe <= kappa && lnSe >= -kappa ){ *\/ */
	      /* 	  /\* 	eps_smooth_face(component,q) =  e0*std::tanh( (lnSe+kappa) / kappa ); //(e0/2.0) * ( 1 + std::sin( pi*(lnSe-s0) / (2.0*kappa) )  ); *\/ */
	      /* 	  /\* 	//std::cout << "t2 = " << current_time << ", lnSe = " << lnSe << ", eps = " << eps_smooth << ", k = " << kappa << ", guy = " << component << ", co = " << co << std::endl; *\/ */
	      /* 	  /\*     } *\/ */
	      /* 	  /\*     else{ *\/ */
	      /* 	  /\*     } *\/ */
	      /* 	  /\*     //} *\/ */
	      /* 	  /\* } *\/ */
	      /* 	} */
	      /* } */

/* 	      if(artificial_visc == true){ */
/* 		iprev_soln_alpha_face[component] =  std::vector<double>(n_q_points_face); */
/* 		fe_values_face[*(ialpha[component])].get_function_values(r_soln[0],iprev_soln_alpha_face[component]); */
		
/* 		o_flow = 0.0; */
/* 		l2_base = 0.0; */
		
/* 		for (unsigned int q=0; q<n_q_points_face; ++q){ */
/* 		  o_flow  += std::pow(prev_soln_alpha_face[component][q]-iprev_soln_alpha_face[component][q],2)*JxW_face[q]; */
/* 		  l2_base += std::pow(prev_soln_alpha_face[component][q],2)*JxW_face[q]; */
/* 		} */
		
/* 		lnSe = std::log10(o_flow/l2_base); */
/* 		//eps_smooth = 0.0; */
		
/* 		if(component < 2){ */
/* 		  if ( lnSe > s0 + kappa){ */
/* 		    //		    eps_smooth(component) = e0; */
/* 		    //std::cout << "t3 = " << current_time << ", lnSe = " << lnSe << ", eps = " << eps_smooth << ", k= " << kappa << ", guy = " << component << ", co = " << co << std::endl; */
/* 		  } */
/* 		  else if ( lnSe <= s0+kappa && lnSe >= s0-kappa ){ */
/* 		    //eps_smooth(component) = e0; //(e0/2.0) * ( 1 + std::sin( pi*(lnSe-s0) / (2.0*kappa) )  );                                                                                                                                            */
/* 		    //std::cout << "t4 = " << current_time << ", lnSe = " << lnSe << ", eps = " << eps_smooth << ", k = " << kappa << ", guy = " << component << ", notes = " << notes << std::endl; */
		    
/* 		  } */
/* 		  else{ */
/* 		  } */
/* 		} */

	      //double sigma_penalty = sig_pen_coeff*( degree+1.0 )*( degree + 1.0 )*( face->measure())/(cell->measure() );
	      double sigma_penalty = sig_pen_coeff;
	      //double sigma_penalty = -(10.0+current_time*0.002)*degree*( degree + 1.0 ) * ( face->measure())/(cell->measure() );
	      //std::cout << "penalty = " << sigma_penalty << std::endl;
	      for (unsigned int q=0; q<n_q_points_face; ++q){
		const Point<dim>& normal = normals[q];
		//const TableBase<3,double> D = functionals.Ficks();
		for (unsigned int i=0; i<dofs_per_cell; ++i){

		  //if( !checksupport || fe.has_support_on_face( pinfo[component].alpha_dof_index[i], face_num ) ){
		  //for(unsigned int l=0; l<dim; ++l){
		  if (component < 2){

		    for (unsigned int l=0; l<2; ++l){
		      for (unsigned int m=0; m<2; ++m){
		    	double temp =  0.5 * ( (prev_soln_sigma_face[component][q])[l] + (prev_soln_sigma_neigh_face[component][q])[m] );
		    	CFL_temp = std::max(CFL_temp,std::abs(temp) );
		      }
		    }

		    /* if (normal(0)>1e-6){ LDG_beta(0) = 0.5/normal(0); LDG_beta(1) = 0.0; */
		    /* } */
		    /* else if (normal(1)>1e-6) {LDG_beta(0)=0.0; LDG_beta(1) = 0.5/normal(1); */
		    /* } */
		    //LDG_beta(0)=0.0; LDG_beta(1)=0.0;

		    max_jump_temp = std::max(std::abs(prev_soln_alpha_face[component][q] - prev_soln_alpha_neigh_face[component][q]),max_jump_temp);
		    //pcout << "pbet = " <<  beta_pen_coeff  << ", jump = " <<    prev_soln_sigma_face[component][q]*normals[q] - prev_soln_sigma_neigh_face[component][q]*normals[q]  << std::endl;
		    
		    //sigma_penalty = 20.0*  prev_soln_sigma_face[component][q]*normals[q] - prev_soln_sigma_neigh_face[component][q]*normals[q];

		    diffusion_flux(component,i) += ( difs(component)+eps_smooth_face(component,q) ) 
		      * (
			 (0.5 * ( (prev_soln_sigma_face[component][q]) + (prev_soln_sigma_neigh_face[component][q]) )
			  * fe_values_face[*(alpha[component])].value(i,q) * normals[q]) 
			 - sigma_penalty * ( ( prev_soln_alpha_face[component][q] - prev_soln_alpha_neigh_face[component][q] ) 
					     * fe_values_face[*(alpha[component])].value(i,q) ) 
			 - beta_pen_coeff * ( ( prev_soln_sigma_face[component][q]*normals[q] - prev_soln_sigma_neigh_face[component][q]*normals[q] ) 
					      * fe_values_face[*(alpha[component])].value(i,q) )
			 ) * JxW_face[q];
		    
		    //std::cout << "what is going on4? = " <<  difs(component)+eps_smooth_face(component,q) << std::endl;

		    if (component == 0){
		      
		      double norm21 = 0.5 * ( (prev_soln_sigma_neigh_face[2][q])[0] + (prev_soln_sigma_face[2][q])[0] ) * normal(1);
		      double norm20 = 0.5 * ( (prev_soln_sigma_neigh_face[2][q])[1] + (prev_soln_sigma_face[2][q])[1] ) * normal(0);
		      /* double norm21 =  0.0 * normal(1); */
		      /* double norm20 =  0.3 * normal(0); */


		      double avg_sca0 =  0.5 * (prev_soln_alpha_face[0][q] + prev_soln_alpha_neigh_face[0][q]);
		      
		      double norm01 = 0.5 * ( (prev_soln_sigma_neigh_face[component][q])[0] + (prev_soln_sigma_face[component][q])[0] ) * normal(1);
		      double norm00 = 0.5 * ( (prev_soln_sigma_neigh_face[component][q])[1] + (prev_soln_sigma_face[component][q])[1] ) * normal(0);
		      double avg_sca2 =  0.5 * (prev_soln_alpha_face[2][q] + prev_soln_alpha_neigh_face[2][q]);
		      
		      //CFL_temp = std::max(CFL_temp,std::abs(avg_sca0) );
		      //CFL_temp = std::max(CFL_temp,std::abs(avg_sca2) );
		      
		      if ( s_flag == 1 ){
			
			convection_flux(component,i) -=  0.5 * ( fe_values_face[*(alpha[component])].value(i,q)) * avg_sca0 * (norm21) * JxW_face[q];
			
			convection_flux(component,i) -=  0.5 * ( fe_values_face[*(alpha[component])].value(i,q)) * avg_sca2 * (norm00) * JxW_face[q];
			
			convection_flux(component,i) +=  0.5 * ( fe_values_face[*(alpha[component])].value(i,q)) * avg_sca0 * (norm20) * JxW_face[q];
			
			convection_flux(component,i) +=  0.5 * ( fe_values_face[*(alpha[component])].value(i,q)) * avg_sca2 * (norm01) * JxW_face[q];
			
		      }
		      
		      else if (s_flag == 2){
			
			convection_flux(component,i) -=  ( fe_values_face[*(alpha[component])].value(i,q)) * avg_sca0 * (norm21) * JxW_face[q];
			
			convection_flux(component,i) +=  ( fe_values_face[*(alpha[component])].value(i,q)) * avg_sca0 * (norm20) * JxW_face[q];
			
		      }
		      
		      
		    }

		    if (component == 1){
		      
		      //double beta = 0.04; //5e-4; //6e-3; //1e-6;                                                                                                                                                                                     		      
		      double norm21 = 0.5 * ( (prev_soln_sigma_neigh_face[2][q])[0] + (prev_soln_sigma_face[2][q])[0] ) * normal(1);
		      double norm20 = 0.5 * ( (prev_soln_sigma_neigh_face[2][q])[1] + (prev_soln_sigma_face[2][q])[1] ) * normal(0);
		      double avg_sca1 =  0.5 * (prev_soln_alpha_face[component][q] + prev_soln_alpha_neigh_face[component][q]);
		       
		      double norm01 = 0.5 * ( (prev_soln_sigma_neigh_face[component][q])[0] + (prev_soln_sigma_face[component][q])[0] ) * normal(1);
		      double norm00 = 0.5 * ( (prev_soln_sigma_neigh_face[component][q])[1] + (prev_soln_sigma_face[component][q])[1] ) * normal(0);
		      double avg_sca2 =  0.5 * (prev_soln_alpha_face[2][q] + prev_soln_alpha_neigh_face[2][q]);

		      double avg_sca0 =  0.5 * (prev_soln_alpha_face[0][q] + prev_soln_alpha_neigh_face[0][q]);
		      
		      //CFL_temp = std::max(CFL_temp,std::abs(avg_sca1) );
		      //CFL_temp = std::max(CFL_temp,std::abs(avg_sca2) );

		      //pcout << "CFL = " << CFL_temp << std::endl;
		      
		      if ( s_flag == 1 ){
			
			convection_flux(component,i) -=  0.5 * ( fe_values_face[*(alpha[component])].value(i,q)) * avg_sca1 * (norm21) * JxW_face[q];
			
			convection_flux(component,i) -=  0.5 * ( fe_values_face[*(alpha[component])].value(i,q)) * avg_sca2 * (norm00) * JxW_face[q];
			
			convection_flux(component,i) +=  0.5 * ( fe_values_face[*(alpha[component])].value(i,q)) * avg_sca1 * (norm20) * JxW_face[q];
			
			convection_flux(component,i) +=  0.5 * ( fe_values_face[*(alpha[component])].value(i,q)) * avg_sca2 * (norm01) * JxW_face[q];

			convection_flux(component,i) -=  beta * avg_sca0  * normal(1) * ( fe_values_face[*(alpha[component])].value(i,q)) *  JxW_face[q];
			
		      }
		      
		      else if ( s_flag == 2 ){
			
			convection_flux(component,i) -=  ( fe_values_face[*(alpha[component])].value(i,q)) * avg_sca1 * (norm21) * JxW_face[q];
			
			convection_flux(component,i) +=  ( fe_values_face[*(alpha[component])].value(i,q)) * avg_sca1 * (norm20) * JxW_face[q];
		      }

		      
		    }

			//			std::cout << "main compo = " <<  ((0.5 * ( (prev_soln_sigma_face[component][q])[l] + (prev_soln_sigma_neigh_face[component][q])[l] )* fe_values_face[*(alpha[component])].value(pinfo[component].alpha_dof_index[i],q) * normal(l))) << std::endl;
			
			//std::cout << "penalty compo = " << sigma_penalty * ( ((prev_soln_alpha_face[component][q])*normal(l)-(prev_soln_alpha_neigh_face[component][q])*normal(l))*fe_values_face[*(alpha[component])].value(pinfo[component].alpha_dof_index[i],q) ) << std::endl << std::endl;

//std::cout << "dif3 = " << D( TableIndices<3>(component,l,l) ) << std::endl;

			//}
			//}

		  }
		}
	      }
	    }
	  }

	  //Vector<double> temp_update[alphadim];
	  //temp_update[component].reinit(alpha_dofs_per_cell);
	  //std::vector<unsigned int> index_container2 (alpha_dofs_per_cell);
	  //for (unsigned int i=0; i<alpha_dofs_per_cell; ++i){
	  // temp_update[component](i) =  delta_t * ( - interior_div(component,i) + diffusion_flux(component,i) + convection_int(component,i) + convection_flux(component,i ) );
	  //index_container2[i] =  local_dof_indices[ pinfo[component].alpha_dof_index[i] ];
	  //}
	  //Vector<double> projected(alpha_dofs_per_cell);
	  //if(pinfo[component].alphaProject){
	  //  pinfo[component].InverseAlphaMassMatrix.vmult(projected,temp_update[component]);
	  //} else {
	  // projected = temp_update[component];
	  //}
	  //projected /= JxWsum;
	  //parahyp_constraints[component].distribute_local_to_global (projected,
	  //						  index_container2,
	  // 						  div_flux_term[component]);
	  
	  h_min_loc = 1.0e6;
	  h_min_loc = std::min(h_min_loc, cell->extent_in_direction(1));
	  h_min_loc = std::min(h_min_loc, cell->extent_in_direction(0));

	  CFL_temp *= CFL_scaling; //16.0;
	  //std::cout << "h_min = " << h_min_loc << ", degree = " << degree << ", CFL_temp = " << CFL_temp << std::endl;
	  double casting =  h_min_loc / ((2.0*degree+1.0)*CFL_temp); 
	  //std::cout << "casting = " << casting << std::endl;
	  //std::cout << "CFL_approx_before = " << CFL_approx << std::endl;
	  CFL_approx = std::min(CFL_approx, casting);
	  //std::cout << "CFL_approx = " << CFL_approx << std::endl;
	  CFL_temp = 0.0;
	  

	  cell->get_dof_indices (local_dof_indices);

	  Vector<double> temp_update[alphadim]; 
	  temp_update[component].reinit(dofs_per_cell); 
	  for (unsigned int i=0; i<dofs_per_cell; ++i){ 
	    if (component !=2 ){ 
	      temp_update[component](i) =  ( - interior_div(component,i) + diffusion_flux(component,i) + convection_int(component,i) + convection_flux(component,i ) ); 
	    }else{ 
	      temp_update[component](i) =  0.0; 
	    } 
	  } 
	  Vector<double> projected(dofs_per_cell); 
	  mapinfo[component][0].localInverseMassMatrix.vmult(projected,temp_update[component]);
	  projected /= JxWsum;
	  parahyp_constraints[component].distribute_local_to_global (projected,
								     local_dof_indices,
								     div_flux_term[component]);
	  /* if(pinfo[component].alphaProject){  */
	  /*   pinfo[component].InverseAlphaMassMatrix.vmult(projected,temp_update[component]);  */
	  /* } else {  */
	  /*   projected = temp_update[component];  */
	  /* }  */
	  /* projected /= JxWsum;  */
	  /* //if (JxWsum <= 1e-12){std::cout << "jxw = " << JxWsum << std::endl;}  */
	  /* for (unsigned int i=0; i<alpha_dofs_per_cell; ++i){  */
	  /*   int global_index = local_dof_indices[ pinfo[component].alpha_dof_index[i] ];  */
	  /*   div_flux_term[component]( global_index ) = projected( i );  */
	  /* }  */
	}

    div_flux_term[component].compress(VectorOperation::insert);
    
  }

  //CFL_temp *= 16.0; //16.0; //h_min_dist; //10.0;
  //CFL_temp = std::max(CFL_temp,0.01);

  //MPI_Allreduce(&CFL_temp, &CFL_bound, Utilities::MPI::n_mpi_processes(mpi_communicator), MPI_DOUBLE, MPI_MAX, mpi_communicator);
  //MPI_Allreduce(&local_min_ribbon_density, &glob_min_ribbon_density, 1, MPI_DOUBLE, MPI_MIN, mpi_communicator);
  MPI_Allreduce(&CFL_approx, &CFL_bound, 1, MPI_DOUBLE, MPI_MIN, mpi_communicator);
  //MPI_Allreduce(&h_min_loc, &h_min_dist, 1, MPI_DOUBLE, MPI_MIN, mpi_communicator);
  //MPI_Allreduce(&max_jump_temp, &gmax_jump, 1, MPI_DOUBLE, MPI_MAX, mpi_communicator);
  if (artificial_visc == true ){
    /* MPI_Allreduce(&visc_temp1, &gvisc_temp1, 1, MPI_DOUBLE, MPI_MAX, mpi_communicator); */
    /* MPI_Allreduce(&visc_temp2, &gvisc_temp2, 1, MPI_DOUBLE, MPI_MIN, mpi_communicator); */
    /* MPI_Allreduce(&art_max, &gart_max, 1, MPI_DOUBLE, MPI_MAX, mpi_communicator); */
    MPI_Allreduce(&max_val_loc[0], &gmax_val_loc[0], 1, MPI_DOUBLE, MPI_MAX, mpi_communicator);
    MPI_Allreduce(&max_val_loc[1], &gmax_val_loc[1], 1, MPI_DOUBLE, MPI_MAX, mpi_communicator);
    
    //double visc_median = (gvisc_temp1 - gvisc_temp2) / 2.0;
    
    if (current_time == dt)  gmax_val_loc[1] = 1.0;
    
    /* MPI_Allreduce(&visc_median, &gvisc_median, 1, MPI_DOUBLE, MPI_MAX, mpi_communicator); */

    /* /\* The Good one for art visc *\/ */
    /* pcout << "Viscosity max = " << gvisc_temp1 << ", min = " << gvisc_temp2 << ", gvisc_median = " << gvisc_median << ", art_max = " << gart_max << ", gmax_val_loc [0] = " << gmax_val_loc[0] << ", gmax_val_loc [1] = " << gmax_val_loc[1] << std::endl; */
  
}

  //pcout << "CFL = " << CFL_bound << ", " << h_min_dist << std::endl; //", " << cell->extent_in_direction(0) << ", " << cell->extent_in_direction(1) << std::endl;

  //pcout << "h_min_dist= " << h_min_dist << ", CFL_bound = " << CFL_bound << std::endl;

  //div_flux_term = naive_subdomain_solution;

  for (unsigned int component=0; component< alphadim; ++component){  

    delete hp_fe_values[component];
    delete hp_fe_values_face[component];
    delete hp_fe_values_neigh_face[component];

  }
}

