/* Let us assemble the stifffness matrix just once, and precondition it just once .... */
template <int dim>
void arcOn<dim>::assemble_stiffness(SolutionVector& subdomain_solution, double delta_t) { 

  (void) delta_t;

  std::vector< FEValues<dim>* > hp_fe_values;
  std::vector< FEFaceValues<dim>* > hp_fe_values_face;
  std::vector< FEFaceValues<dim>* > hp_fe_values_neigh_face;

  double current_time = 1.0;
  for (unsigned int component=0; component< alphadim; ++component){

    UpdateFlags updateflags=  update_values | update_gradients 
      | update_quadrature_points | update_JxW_values;
    
    hp_fe_values.push_back (new FEValues<dim>(*(fe_collection[component]), 
					      *(quadrature_collection[component]), updateflags));
    hp_fe_values_face.push_back (new FEFaceValues<dim>(*(fe_collection[component]), 
						       *(face_quadrature_collection[component]),  
						       updateflags | update_normal_vectors ));
    hp_fe_values_neigh_face.push_back (new FEFaceValues<dim>(*(fe_collection[component]), 
							     *(face_quadrature_collection[component]),  
							     updateflags | update_normal_vectors ));

    if (component == 2){

      poisson_matrix = 0.0;
      int it = 0;
      typename DoFHandler<dim>::active_cell_iterator 
	cell = dof_handler[component]->begin_active(), 
	endc = dof_handler[component]->end();
      for(;
	  cell!=endc;
	  ++cell	)
	if (cell->is_locally_owned()  )
	  {
	    int CO = cell->index();
	    it = it + 1;
	    hp_fe_values[component]->reinit (cell);
	    const FEValues<dim> &fe_values = hp_fe_values[component]->get_present_fe_values ();
	    const Quadrature<dim>& quadrature_formula = fe_values.get_quadrature();
	    const FiniteElement<dim>& fe = cell->get_fe();
	    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
	    const unsigned int   n_q_points    = quadrature_formula.size();

	    FullMatrix<double> M11( dofs_per_cell,dofs_per_cell);
	    FullMatrix<double> M22( dofs_per_cell,dofs_per_cell);
	    
	    FullMatrix<double> M12( dofs_per_cell,dofs_per_cell);
	    FullMatrix<double> M21( dofs_per_cell,dofs_per_cell);

	    FullMatrix<double> A_matrix( dofs_per_cell,dofs_per_cell);
	    FullMatrix<double> boundary_matrix( dofs_per_cell,dofs_per_cell);

	    std::vector<unsigned int> local_dof_indices (dofs_per_cell);
	    std::vector<unsigned int> neigh_local_dof_indices (dofs_per_cell);

	    A_matrix = 0.0;
	      
	    for(unsigned int k=0;k<alphadim;k++){
	      prev_soln_alpha[k] =  std::vector<double>(n_q_points);
	      prev_soln_sigma[k] =  std::vector<Tensor<1,dim> >(n_q_points);
	    }
	      
	    cell->get_dof_indices (local_dof_indices);

	    const std::vector<double> &JxW = fe_values.get_JxW_values ();
	      
	    fe_values[*(alpha[component])].get_function_values(subdomain_solution[component],prev_soln_alpha[component]);
	    fe_values[*(sigma[component])].get_function_values(subdomain_solution[component],prev_soln_sigma[component]);
	    fe_values[*(alpha[1])].get_function_values(subdomain_solution[1],prev_soln_alpha[1]);
	    fe_values[*(sigma[1])].get_function_values(subdomain_solution[1],prev_soln_sigma[1]);

	    std::vector<unsigned int> index_container2 (dofs_per_cell);

	    for (unsigned int i=0; i<dofs_per_cell; ++i){
	      for (unsigned int j=0; j<dofs_per_cell; ++j){
		for (unsigned int q=0; q<n_q_points; ++q){
		      
		  A_matrix(i,j) += (fe_values[*(alpha[component])].gradient(i,q)) * (fe_values[*(alpha[component])].gradient(j,q)) * JxW[q];

		}

		poisson_matrix.add( local_dof_indices[ i ],
				    local_dof_indices[ j ],
				    A_matrix(i,j));
	      }
	    }

	    for(unsigned int face_num=0; face_num < GeometryInfo<dim>::faces_per_cell; ++face_num){
	      typename DoFHandler<dim>::face_iterator face=cell->face(face_num);
	      
	      M11 = 0.0;
	      M22 = 0.0;
	      M12 = 0.0;
	      M21 = 0.0;
	      
	      boundary_matrix = 0.0;

	      if (face->at_boundary()  && (face->boundary_indicator() == 0 || (face->boundary_indicator() == 15 ) ) )
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

		  fe_values_face[*(alpha[component])].get_function_values(subdomain_solution[component],prev_soln_alpha_face[component]);
		  fe_values_face[*(sigma[component])].get_function_values(subdomain_solution[component],prev_soln_sigma_face[component]);

		  unsigned char boundary_index = face->boundary_indicator();
		  const WallBoundaryValues<dim>& wbv = WBV[boundary_index];
		  const std::vector< Point<dim> > &quadrature_points = fe_values_face.get_quadrature_points();

		  const unsigned int normal1 = GeometryInfo<dim>::unit_normal_direction[face_num];
		  
		  double pen1 = cell->extent_in_direction(normal1);
		  
		  double sigma_p1 = degree*( degree + 1.0 ) / pen1;
		  double sigma_p2 = degree*( degree + 1.0 ) / pen1;
		  
		  double sigma_penalty = sigma_prefactor*0.5*(sigma_p1+sigma_p2);

		  for (unsigned int i=0; i<dofs_per_cell; ++i){
		    for (unsigned int j=0; j<dofs_per_cell; ++j){
		      for (unsigned int q=0; q<n_q_points_face; ++q){
			const Point<dim>& normal = normals[q]; 

			boundary_matrix(i,j) += sigma_penalty *  fe_values_face[*(alpha[component])].value(i,q)
                          * fe_values_face[*(alpha[component])].value(j,q) * JxW_face[q]
			  - normal *  (fe_values_face[*(alpha[component])].gradient(i,q))
                          * fe_values_face[*(alpha[component])].value(j,q) * JxW_face[q]
			  - normal *  (fe_values_face[*(alpha[component])].gradient(j,q))
                          * fe_values_face[*(alpha[component])].value(i,q) * JxW_face[q];
			
		      }
		      
		      poisson_matrix.add( local_dof_indices[ i ],
		      			  local_dof_indices[ j ],
		      			  boundary_matrix(i,j));
		      
		    }
		  }
		}
	      else if (  face->at_boundary() && ( face->boundary_indicator() == 3 ||
						  face->boundary_indicator() == 4))
		{
		  neigh_local_dof_indices = periodic_dof_indX[component][CO];
		  
		  hp_fe_values_face[component]->reinit(cell, face_num);
		  
		  const FEFaceValues<dim>&	   fe_values_face = hp_fe_values_face[component]->get_present_fe_values ();
		  
		  const std::vector<double> &JxW_face = fe_values_face.get_JxW_values ();
		  const std::vector<Point<dim> > &normals = fe_values_face.get_normal_vectors ();   
		  
		  const Quadrature<dim-1>& face_quadrature_formula = fe_values_face.get_quadrature();
		  const unsigned int n_q_points_face = face_quadrature_formula.size();
		  
		  const unsigned int normal1 = GeometryInfo<dim>::unit_normal_direction[face_num];

		  double pen1 = cell->extent_in_direction(normal1);
		  
		  double sigma_p1 = degree*( degree + 1.0 ) / pen1;
		  double sigma_p2 = degree*( degree + 1.0 ) / pen1;
		  
		  double sigma_penalty = sigma_prefactor*0.5*(sigma_p1+sigma_p2);
		  double nu;
		  
		  if (elliptic_type == 1){ 
		    nu =  1.0; //SIPG
		  }
		  else if (elliptic_type == -1){
		    nu = -1.0; //NIPG
		  }
		  else {
		    nu =  0.0; //IIPG
		  }
		  
		  for (unsigned int q=0; q<n_q_points_face; ++q){
		    for (unsigned int i=0; i<dofs_per_cell; ++i){
		      for (unsigned int j=0; j<dofs_per_cell; ++j){
			
			const Point<dim>& normal = normals[q]; 
			
			const double vi = (fe_values_face[*(alpha[component])].value(i,q));
			const double dnvi =  (fe_values_face[*(alpha[component])].gradient(i,q)) * normal;
			const double ve = periodic_basisX[component][CO](i,q);
			const double dnve = periodic_gradX[component][CO][q][i] * normal;
			const double ui = (fe_values_face[*(alpha[component])].value(j,q));
			const double dnui =  (fe_values_face[*(alpha[component])].gradient(j,q)) * normal;
			const double ue = periodic_basisX[component][CO](j,q);
			const double dnue = periodic_gradX[component][CO][q][j] * normal;
			
			M11(i,j) +=  JxW_face[q]*(-0.5*dnvi*ui*nu-0.5*dnui*vi+sigma_penalty*ui*vi);
			M22(i,j) +=  JxW_face[q]*( 0.5*dnve*ue*nu+0.5*dnue*ve+sigma_penalty*ue*ve);//
			M12(i,j) +=  JxW_face[q]*( 0.5*dnvi*ue*nu-0.5*dnue*vi-sigma_penalty*vi*ue);
			M21(i,j) +=  JxW_face[q]*(-0.5*dnve*ui*nu+0.5*dnui*ve-sigma_penalty*ui*ve);//
			
		      }
		    }
		  }
		  
		  for (unsigned int i=0; i<dofs_per_cell; ++i){
		    for (unsigned int j=0; j<dofs_per_cell; ++j){
		      
		      int global_index_i = local_dof_indices[ i ]; 
		      int global_index_j = local_dof_indices[ j ];
		      int global_index_j_neigh = neigh_local_dof_indices[ j ];
		      int global_index_i_neigh = neigh_local_dof_indices[ i ];
		      
		      poisson_matrix.add( global_index_i,
					  global_index_j,
					  M11(i,j) );
		      poisson_matrix.add( global_index_i_neigh,
					  global_index_j_neigh,
					  M22(i,j) );
		      poisson_matrix.add( global_index_i,
					  global_index_j_neigh,
					  M12(i,j));
		      poisson_matrix.add( global_index_i_neigh,
					  global_index_j,
					  M21(i,j));

		    }
		  }
		}
	      else if (  face->at_boundary() && ( face->boundary_indicator() == 1 ||
						  face->boundary_indicator() == 2))
		{
		  neigh_local_dof_indices = periodic_dof_indY[component][CO];
		  
		  hp_fe_values_face[component]->reinit(cell, face_num);
		  
		  const FEFaceValues<dim>&	   fe_values_face = hp_fe_values_face[component]->get_present_fe_values ();
		  
		  const std::vector<double> &JxW_face = fe_values_face.get_JxW_values ();
		  const std::vector<Point<dim> > &normals = fe_values_face.get_normal_vectors ();   
		  
		  const Quadrature<dim-1>& face_quadrature_formula = fe_values_face.get_quadrature();
		  const unsigned int n_q_points_face = face_quadrature_formula.size();
		  
		  const unsigned int normal1 = GeometryInfo<dim>::unit_normal_direction[face_num];
		  
		  double pen1 = cell->extent_in_direction(normal1);
		  
		  double sigma_p1 = degree*( degree + 1.0 ) / pen1;
		  double sigma_p2 = degree*( degree + 1.0 ) / pen1;
		  
		  double sigma_penalty = sigma_prefactor*0.5*(sigma_p1+sigma_p2);
		  double nu;
		  
		  if (elliptic_type == 1){ 
		    nu =  1.0; //SIPG
		  }
		  else if (elliptic_type == -1){
		    nu = -1.0; //NIPG
		  }
		  else {
		    nu =  0.0; //IIPG
		  }
		  
		  for (unsigned int q=0; q<n_q_points_face; ++q){
		    for (unsigned int i=0; i<dofs_per_cell; ++i){
		      for (unsigned int j=0; j<dofs_per_cell; ++j){
			
			const Point<dim>& normal = normals[q]; 
			
			const double vi = (fe_values_face[*(alpha[component])].value(i,q));
			const double dnvi =  (fe_values_face[*(alpha[component])].gradient(i,q)) * normal;
			const double ve = periodic_basisY[component][CO](i,q);
			const double dnve = periodic_gradY[component][CO][q][i] * normal;
			const double ui = (fe_values_face[*(alpha[component])].value(j,q));
			const double dnui =  (fe_values_face[*(alpha[component])].gradient(j,q)) * normal;
			const double ue = periodic_basisY[component][CO](j,q);
			const double dnue = periodic_gradY[component][CO][q][j] * normal;
			
			M11(i,j) +=  JxW_face[q]*(-0.5*dnvi*ui*nu-0.5*dnui*vi+sigma_penalty*ui*vi);
			M22(i,j) +=  JxW_face[q]*( 0.5*dnve*ue*nu+0.5*dnue*ve+sigma_penalty*ue*ve);//
			M12(i,j) +=  JxW_face[q]*( 0.5*dnvi*ue*nu-0.5*dnue*vi-sigma_penalty*vi*ue);
			M21(i,j) +=  JxW_face[q]*(-0.5*dnve*ui*nu+0.5*dnui*ve-sigma_penalty*ui*ve);//
			
		      }
		    }
		  }
		  
		  for (unsigned int i=0; i<dofs_per_cell; ++i){
		    for (unsigned int j=0; j<dofs_per_cell; ++j){
		      
		      int global_index_i = local_dof_indices[ i ]; 
		      int global_index_j = local_dof_indices[ j ];
		      int global_index_j_neigh = neigh_local_dof_indices[ j ];
		      int global_index_i_neigh = neigh_local_dof_indices[ i ];
		      
		      poisson_matrix.add( global_index_i,
					  global_index_j,
					  M11(i,j) );
		      poisson_matrix.add( global_index_i_neigh,
					  global_index_j_neigh,
					  M22(i,j) );
		      poisson_matrix.add( global_index_i,
					  global_index_j_neigh,
					  M12(i,j));
		      poisson_matrix.add( global_index_i_neigh,
					  global_index_j,
					  M21(i,j));

		    }
		  }
		}
	      else{ //The is the interior face integration step
		typename DoFHandler<dim>::cell_iterator neighbor = cell->neighbor(face_num);
		const unsigned int neighbor2 = cell->neighbor_of_neighbor(face_num);

		neighbor->get_dof_indices (neigh_local_dof_indices);
		  
		hp_fe_values_face[component]->reinit(cell, face_num);
		hp_fe_values_neigh_face[component]->reinit(neighbor,neighbor2);
		  
		const FEFaceValues<dim>&	   fe_values_face = hp_fe_values_face[component]->get_present_fe_values ();
		const FEFaceValues<dim>&   fe_values_neigh_face = hp_fe_values_neigh_face[component]->get_present_fe_values();
		  
		const std::vector<double> &JxW_face = fe_values_face.get_JxW_values ();
		const std::vector<Point<dim> > &normals = fe_values_face.get_normal_vectors ();   
	      
		const Quadrature<dim-1>& face_quadrature_formula = fe_values_face.get_quadrature();
		const unsigned int n_q_points_face = face_quadrature_formula.size();

		const unsigned int normal1 = GeometryInfo<dim>::unit_normal_direction[face_num];
		const unsigned int normal2 = GeometryInfo<dim>::unit_normal_direction[neighbor2];		

		double pen1 = cell->extent_in_direction(normal1);
		double pen2 = neighbor->extent_in_direction(normal2);

		double sigma_p1 = degree*( degree + 1.0 ) / pen1;
		double sigma_p2 = degree*( degree + 1.0 ) / pen2;

		double sigma_penalty = sigma_prefactor*0.5*(sigma_p1+sigma_p2);

		double nu;

		if (elliptic_type == 1){ 
		  nu =  1.0; //SIPG
		}
		else if (elliptic_type == -1){
		  nu = -1.0; //NIPG
		}
		else {
		  nu =  0.0; //IIPG
		}

		for (unsigned int q=0; q<n_q_points_face; ++q){
		  for (unsigned int i=0; i<dofs_per_cell; ++i){
		    for (unsigned int j=0; j<dofs_per_cell; ++j){
		      
		      const Point<dim>& normal = normals[q]; 

		      const double vi = (fe_values_face[*(alpha[component])].value(i,q));
		      const double dnvi =  (fe_values_face[*(alpha[component])].gradient(i,q)) * normal;
		      const double ve = (fe_values_neigh_face[*(alpha[component])].value(i,q));
		      const double dnve = (fe_values_neigh_face[*(alpha[component])].gradient(i,q)) * normal;
		      const double ui = (fe_values_face[*(alpha[component])].value(j,q));
		      const double dnui =  (fe_values_face[*(alpha[component])].gradient(j,q)) * normal;
		      const double ue = (fe_values_neigh_face[*(alpha[component])].value(j,q));
		      const double dnue = (fe_values_neigh_face[*(alpha[component])].gradient(j,q)) * normal;

		      M11(i,j) +=  JxW_face[q]*(-0.5*dnvi*ui*nu-0.5*dnui*vi+sigma_penalty*ui*vi);
		      M22(i,j) +=  JxW_face[q]*( 0.5*dnve*ue*nu+0.5*dnue*ve+sigma_penalty*ue*ve);//
		      M12(i,j) +=  JxW_face[q]*( 0.5*dnvi*ue*nu-0.5*dnue*vi-sigma_penalty*vi*ue);
		      M21(i,j) +=  JxW_face[q]*(-0.5*dnve*ui*nu+0.5*dnui*ve-sigma_penalty*ui*ve);//

		    }
		  }
		}
		
		for (unsigned int i=0; i<dofs_per_cell; ++i){
		  for (unsigned int j=0; j<dofs_per_cell; ++j){

		    int global_index_i = local_dof_indices[ i ]; 
		    int global_index_j = local_dof_indices[ j ];
		    int global_index_j_neigh = neigh_local_dof_indices[ j ];
		    int global_index_i_neigh = neigh_local_dof_indices[ i ];

		    poisson_matrix.add( global_index_i,
		    			global_index_j,
		    			M11(i,j) );
                    poisson_matrix.add( global_index_i_neigh,
                                        global_index_j_neigh,
                                        M22(i,j) );
		    poisson_matrix.add( global_index_i,
					global_index_j_neigh,
					M12(i,j));
		    poisson_matrix.add( global_index_i_neigh,
                                        global_index_j,
                                        M21(i,j));

		  }
		}
		
	
		
	      }
	    }
	  }
      
      poisson_matrix.compress (VectorOperation::add);
      preconditioner.initialize(poisson_matrix.block(0,0),
				PETScWrappers::PreconditionBoomerAMG::AdditionalData(true));
    }
  }
  for (unsigned int component=0; component< alphadim; ++component){  
    delete hp_fe_values[component];
    delete hp_fe_values_face[component];
    delete hp_fe_values_neigh_face[component];
  }
}

