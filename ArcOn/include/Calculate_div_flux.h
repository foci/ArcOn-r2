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

  }

  for (unsigned int component=0; component< alphadim-1; ++component){

    if (artificial_visc == true){ 

      max_val_loc[component] = 0.0;

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
    }

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

	  for(unsigned int k=0;k<alphadim;k++){ 
	    if(k!=component){
	      fe_values[*(alpha[k])].get_function_values(subdomain_solution[k],prev_soln_alpha[k]);
	      fe_values[*(sigma[k])].get_function_values(subdomain_solution[k],prev_soln_sigma[k]);
	    }
	  }

	  FullMatrix<double> eps_smooth(alphadim-1, n_q_points);
	  eps_smooth = 0.0;

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

	  for (unsigned int q=0; q<n_q_points; ++q){
	    for (unsigned int i=0; i<dofs_per_cell; ++i){
	      interior_div(component,i) += (difs(component)+eps_smooth(component,q)) * (fe_values[*(alpha[component])].gradient(i,q)) * (prev_soln_sigma[component][q]) * JxW[q];

	      if(component == 0){
		//this is the chi stuff
		interior_div(component,i) += (difs(component)+eps_smooth(component,q))*(prev_soln_sigma[component][q])* (prev_soln_sigma[component][q])*(fe_values[*(alpha[component])].value(i,q))* JxW[q];

	      }
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
	      for (unsigned int i=0; i<dofs_per_cell; ++i){
                if ( s_flag == 0 ){
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
	    if ( face->at_boundary() && (face->boundary_indicator() == 0) )
	      {
		hp_fe_values_face[component]->reinit (cell,face_num);
		const FEFaceValues<dim>& fe_values_face = hp_fe_values_face[component]->get_present_fe_values ();
		const Quadrature<dim-1>& face_quadrature_formula = fe_values_face.get_quadrature();
		const unsigned int n_q_points_face = face_quadrature_formula.size();

		const std::vector<double> &JxW_face = fe_values_face.get_JxW_values ();
		const std::vector<Point<dim> > &normals = fe_values_face.get_normal_vectors ();   

		prev_soln_alpha_face[component] =  std::vector<double>(n_q_points_face);
		prev_soln_sigma_face[component] =  std::vector<Tensor<1,dim> >(n_q_points_face);
		prev_soln_alpha_face[0] =  std::vector<double>(n_q_points_face);
		prev_soln_sigma_face[0] =  std::vector<Tensor<1,dim> >(n_q_points_face);
		prev_soln_alpha_face[2] =  std::vector<double>(n_q_points_face);
		prev_soln_sigma_face[2] =  std::vector<Tensor<1,dim> >(n_q_points_face);

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
		
		std::vector<double> alphas_boundary(n_q_points_face);
		std::vector<double> alphas_boundary2(n_q_points_face);
		std::vector<double> alphas_boundary0(n_q_points_face);
		
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

		wbv.value_list( quadrature_points, alphas_boundary, component, current_time);
		wbv.value_list( quadrature_points, alphas_boundary2, 2, current_time);
		wbv.value_list( quadrature_points, alphas_boundary0, 0, current_time);
		wbv.gradient_list( quadrature_points, sigmas_boundary, normals, component, current_time);
		wbv.gradient_list( quadrature_points, sigmas_boundary2, normals, 2, current_time);
				
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
		}
		

		for (unsigned int q=0; q<fe_values_face.n_quadrature_points; ++q){
		  const Point<dim>& normal = normals[q]; 
		  for (unsigned int i=0; i< dofs_per_cell; ++i){
		    if (component < 2){ 
		      diffusion_flux(component,i) += ( difs(component) + eps_smooth_face(component,q) )* ( sigmas_boundary[q] ) * (fe_values_face[*(alpha[component])].value(i,q)) * normals[q] * JxW_face[q];
		      
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

		prev_soln_alpha_face[component] =  std::vector<double>(n_q_points_face);
		prev_soln_sigma_face[component] =  std::vector<Tensor<1,dim> >(n_q_points_face);
		prev_soln_alpha_face[0] =  std::vector<double>(n_q_points_face);
		prev_soln_sigma_face[0] =  std::vector<Tensor<1,dim> >(n_q_points_face);
		prev_soln_alpha_face[2] =  std::vector<double>(n_q_points_face);
		prev_soln_sigma_face[2] =  std::vector<Tensor<1,dim> >(n_q_points_face);

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
		
		std::vector<double> alphas_boundary(n_q_points_face);
		std::vector<double> alphas_boundary2(n_q_points_face);
		std::vector<double> alphas_boundary0(n_q_points_face);
				
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

		wbv.value_list2( quadrature_points, alphas_boundary, component, current_time);
		wbv.value_list2( quadrature_points, alphas_boundary2, 2, current_time);
		wbv.value_list2( quadrature_points, alphas_boundary0, 0, current_time);
		wbv.gradient_list2( quadrature_points, sigmas_boundary, normals, component, current_time);
		wbv.gradient_list2( quadrature_points, sigmas_boundary2, normals, 2, current_time);
		
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
		}

		for (unsigned int q=0; q<fe_values_face.n_quadrature_points; ++q){
		  const Point<dim>& normal = normals[q]; 
		  for (unsigned int i=0; i< dofs_per_cell; ++i){
		    if (component < 2){ 
		    
		      diffusion_flux(component,i) += ( difs(component) + eps_smooth_face(component,q) )* ( sigmas_boundary[q] ) * (fe_values_face[*(alpha[component])].value(i,q)) * normals[q] * JxW_face[q];
		      
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

		double sigma_penalty = sig_pen_coeff;
		for (unsigned int q=0; q<n_q_points_face; ++q){
		  const Point<dim>& normal = normals[q];
		  for (unsigned int i=0; i<dofs_per_cell; ++i){		    
		    if (component < 2){
		      
		      diffusion_flux(component,i) += ( difs(component)+eps_smooth_face(component,q) ) 
			* (
			   (0.5 * ( (prev_soln_sigma_face[component][q]) + (grad_drawX[component][CO][q]) )
			    * fe_values_face[*(alpha[component])].value(i,q) * normals[q]) 
			   - sigma_penalty * ( ( prev_soln_alpha_face[component][q] - soln_drawX[component][CO][q] )
					       * fe_values_face[*(alpha[component])].value(i,q)  )
			   - beta_pen_coeff * ( ( prev_soln_sigma_face[component][q]*normals[q] - grad_drawX[component][CO][q]*normals[q] )
						* fe_values_face[*(alpha[component])].value(i,q) )
			   ) * JxW_face[q];
		      
		      if (component == 0){
			
	    		double norm21 = 0.5 * (  (grad_drawX[2][CO][q])[0] + (prev_soln_sigma_face[2][q])[0] ) * normal(1);
	    		double norm20 = 0.5 * (  (grad_drawX[2][CO][q])[1] + (prev_soln_sigma_face[2][q])[1] ) * normal(0);

	    		double avg_sca0 =  0.5 * (prev_soln_alpha_face[component][q] +  soln_drawX[component][CO][q]);
			
	    		double norm01 = 0.5 * (  (grad_drawX[component][CO][q])[0] + (prev_soln_sigma_face[component][q])[0] ) * normal(1);
	    		double norm00 = 0.5 * (  (grad_drawX[component][CO][q])[1] + (prev_soln_sigma_face[component][q])[1] ) * normal(0);
	    		double avg_sca2 =  0.5 * (prev_soln_alpha_face[2][q] + soln_drawX[2][CO][q]);

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
			
	    		double norm21 = 0.5 * (  (grad_drawX[2][CO][q])[0] + (prev_soln_sigma_face[2][q])[0] ) * normal(1);
	    		double norm20 = 0.5 * (  (grad_drawX[2][CO][q])[1] + (prev_soln_sigma_face[2][q])[1] ) * normal(0);
	    		double avg_sca1 =  0.5*(prev_soln_alpha_face[component][q] +  soln_drawX[component][CO][q]);
			
	    		double norm01 = 0.5 * (  (grad_drawX[component][CO][q])[0] + (prev_soln_sigma_face[component][q])[0] ) * normal(1);
	    		double norm00 = 0.5 * (  (grad_drawX[component][CO][q])[1] + (prev_soln_sigma_face[component][q])[1] ) * normal(0);
	    		double avg_sca2 =  0.5*(prev_soln_alpha_face[2][q] + soln_drawX[2][CO][q]);

			double avg_sca0 = 0.5*(prev_soln_alpha_face[0][q] + soln_drawX[0][CO][q]);

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

		double sigma_penalty = sig_pen_coeff;
		for (unsigned int q=0; q<n_q_points_face; ++q){
		  const Point<dim>& normal = normals[q];
		  for (unsigned int i=0; i<dofs_per_cell; ++i){
		    
		    if (component < 2){
			
		      diffusion_flux(component,i) += ( difs(component)+eps_smooth_face(component,q) ) 
			* (
			   (0.5 * ( (prev_soln_sigma_face[component][q]) + (grad_drawY[component][CO][q]) )
			    * fe_values_face[*(alpha[component])].value(i,q) * normals[q]) 
			   - sigma_penalty * ( ( prev_soln_alpha_face[component][q] - soln_drawY[component][CO][q] )
					       * fe_values_face[*(alpha[component])].value(i,q) ) 
			   - beta_pen_coeff * ( ( prev_soln_sigma_face[component][q]*normals[q] - grad_drawY[component][CO][q]*normals[q] )
						* fe_values_face[*(alpha[component])].value(i,q) )
			   ) * JxW_face[q];
		      
		      if (component == 0){
			
	    		double norm21 = 0.5 * (  (grad_drawY[2][CO][q])[0] + (prev_soln_sigma_face[2][q])[0] ) * normal(1);
	    		double norm20 = 0.5 * (  (grad_drawY[2][CO][q])[1] + (prev_soln_sigma_face[2][q])[1] ) * normal(0);

	    		double avg_sca0 =  0.5 * (prev_soln_alpha_face[component][q] +  soln_drawY[component][CO][q]);
			
	    		double norm01 = 0.5 * (  (grad_drawY[component][CO][q])[0] + (prev_soln_sigma_face[component][q])[0] ) * normal(1);
	    		double norm00 = 0.5 * (  (grad_drawY[component][CO][q])[1] + (prev_soln_sigma_face[component][q])[1] ) * normal(0);
	    		double avg_sca2 =  0.5 * (prev_soln_alpha_face[2][q] + soln_drawY[2][CO][q]);

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
			
	    		double norm21 = 0.5 * (  (grad_drawY[2][CO][q])[0] + (prev_soln_sigma_face[2][q])[0] ) * normal(1);
	    		double norm20 = 0.5 * (  (grad_drawY[2][CO][q])[1] + (prev_soln_sigma_face[2][q])[1] ) * normal(0);
	    		double avg_sca1 =  0.5*(prev_soln_alpha_face[component][q] +  soln_drawY[component][CO][q]);
			
	    		double norm01 = 0.5 * (  (grad_drawY[component][CO][q])[0] + (prev_soln_sigma_face[component][q])[0] ) * normal(1);
	    		double norm00 = 0.5 * (  (grad_drawY[component][CO][q])[1] + (prev_soln_sigma_face[component][q])[1] ) * normal(0);
	    		double avg_sca2 =  0.5*(prev_soln_alpha_face[2][q] + soln_drawY[2][CO][q]);

			double avg_sca0 = 0.5*(prev_soln_alpha_face[0][q] + soln_drawY[0][CO][q]);

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

	      double sigma_penalty = sig_pen_coeff;
	      for (unsigned int q=0; q<n_q_points_face; ++q){
		const Point<dim>& normal = normals[q];
		for (unsigned int i=0; i<dofs_per_cell; ++i){

		  if (component < 2){

		    for (unsigned int l=0; l<2; ++l){
		      for (unsigned int m=0; m<2; ++m){
		    	double temp =  0.5 * ( (prev_soln_sigma_face[component][q])[l] + (prev_soln_sigma_neigh_face[component][q])[m] );
		    	CFL_temp = std::max(CFL_temp,std::abs(temp) );
		      }
		    }

		    max_jump_temp = std::max(std::abs(prev_soln_alpha_face[component][q] - prev_soln_alpha_neigh_face[component][q]),max_jump_temp);

		    diffusion_flux(component,i) += ( difs(component)+eps_smooth_face(component,q) ) 
		      * (
			 (0.5 * ( (prev_soln_sigma_face[component][q]) + (prev_soln_sigma_neigh_face[component][q]) )
			  * fe_values_face[*(alpha[component])].value(i,q) * normals[q]) 
			 - sigma_penalty * ( ( prev_soln_alpha_face[component][q] - prev_soln_alpha_neigh_face[component][q] ) 
					     * fe_values_face[*(alpha[component])].value(i,q) ) 
			 - beta_pen_coeff * ( ( prev_soln_sigma_face[component][q]*normals[q] - prev_soln_sigma_neigh_face[component][q]*normals[q] ) 
					      * fe_values_face[*(alpha[component])].value(i,q) )
			 ) * JxW_face[q];

		    if (component == 0){
		      
		      double norm21 = 0.5 * ( (prev_soln_sigma_neigh_face[2][q])[0] + (prev_soln_sigma_face[2][q])[0] ) * normal(1);
		      double norm20 = 0.5 * ( (prev_soln_sigma_neigh_face[2][q])[1] + (prev_soln_sigma_face[2][q])[1] ) * normal(0);

		      double avg_sca0 =  0.5 * (prev_soln_alpha_face[0][q] + prev_soln_alpha_neigh_face[0][q]);
		      
		      double norm01 = 0.5 * ( (prev_soln_sigma_neigh_face[component][q])[0] + (prev_soln_sigma_face[component][q])[0] ) * normal(1);
		      double norm00 = 0.5 * ( (prev_soln_sigma_neigh_face[component][q])[1] + (prev_soln_sigma_face[component][q])[1] ) * normal(0);
		      double avg_sca2 =  0.5 * (prev_soln_alpha_face[2][q] + prev_soln_alpha_neigh_face[2][q]);
		      
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
		      
		      double norm21 = 0.5 * ( (prev_soln_sigma_neigh_face[2][q])[0] + (prev_soln_sigma_face[2][q])[0] ) * normal(1);
		      double norm20 = 0.5 * ( (prev_soln_sigma_neigh_face[2][q])[1] + (prev_soln_sigma_face[2][q])[1] ) * normal(0);
		      double avg_sca1 =  0.5 * (prev_soln_alpha_face[component][q] + prev_soln_alpha_neigh_face[component][q]);
		       
		      double norm01 = 0.5 * ( (prev_soln_sigma_neigh_face[component][q])[0] + (prev_soln_sigma_face[component][q])[0] ) * normal(1);
		      double norm00 = 0.5 * ( (prev_soln_sigma_neigh_face[component][q])[1] + (prev_soln_sigma_face[component][q])[1] ) * normal(0);
		      double avg_sca2 =  0.5 * (prev_soln_alpha_face[2][q] + prev_soln_alpha_neigh_face[2][q]);

		      double avg_sca0 =  0.5 * (prev_soln_alpha_face[0][q] + prev_soln_alpha_neigh_face[0][q]);
		      
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
		  }
		}
	      }
	    }
	  }

	  h_min_loc = 1.0e6;
	  h_min_loc = std::min(h_min_loc, cell->extent_in_direction(1));
	  h_min_loc = std::min(h_min_loc, cell->extent_in_direction(0));

	  CFL_temp *= CFL_scaling; //16.0;
	  double casting =  h_min_loc / ((2.0*degree+1.0)*CFL_temp); 

	  CFL_approx = std::min(CFL_approx, casting);
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
	  if (component ==0){
	    density_constraints.distribute_local_to_global (projected,
							    local_dof_indices,
							    div_flux_term[component]); 
	  }
	  else if (component == 1) {
	    vorticity_constraints.distribute_local_to_global (projected,
							      local_dof_indices,
							      div_flux_term[component]);
	  }
	}

    div_flux_term[component].compress(VectorOperation::add);
    
  }

  MPI_Allreduce(&CFL_approx, &CFL_bound, 1, MPI_DOUBLE, MPI_MIN, mpi_communicator);

  if (artificial_visc == true ){

    MPI_Allreduce(&max_val_loc[0], &gmax_val_loc[0], 1, MPI_DOUBLE, MPI_MAX, mpi_communicator);
    MPI_Allreduce(&max_val_loc[1], &gmax_val_loc[1], 1, MPI_DOUBLE, MPI_MAX, mpi_communicator);
    
    if (current_time == dt)  gmax_val_loc[1] = 1.0;
    
  }

  for (unsigned int component=0; component< alphadim; ++component){  

    delete hp_fe_values[component];
    delete hp_fe_values_face[component];
    delete hp_fe_values_neigh_face[component];

  }
}

