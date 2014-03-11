/* This assembles the "auxiliary variables," which are the ordered gradients of the state vector */
template <int dim>
void arcOn<dim>::assemble_sigma(SolutionVector& subdomain_solution, double current_time)
{

  std::vector< FEValues<dim>* > hp_fe_values;
  std::vector< FEFaceValues<dim>* > hp_fe_values_face;
  std::vector< FEFaceValues<dim>* > hp_fe_values_neigh_face;

  for (unsigned int component=0; component< alphadim; ++component){
    
    UpdateFlags updateflags=  update_values | update_gradients 
      | update_quadrature_points | update_JxW_values;

    hp_fe_values.push_back (new FEValues<dim>(*(fe_collection[component]), 
					      *(quadrature_collection[component]), updateflags | update_hessians));
    hp_fe_values_face.push_back (new FEFaceValues<dim>(*(fe_collection[component]), 
						       *(face_quadrature_collection[component]),  
						       updateflags | update_normal_vectors | update_hessians));
    hp_fe_values_neigh_face.push_back (new FEFaceValues<dim>(*(fe_collection[component]), 
							     *(face_quadrature_collection[component]),  
							     updateflags | update_normal_vectors | update_hessians));

    naive_subdomain_solution[component] = 0.0;
    //naive_subdomain_solution[component].block(0) =  subdomain_solution[component].block(0);

    std::vector<double> transport_alphas(alphadim,0.0);
    Functionals<dim> functionals;

    /* Assemble sigma */
    typename DoFHandler<dim>::active_cell_iterator 
      cell = dof_handler[component]->begin_active(), 
      endc = dof_handler[component]->end();
    for( ;
	cell!=endc;
	++cell )
      if ( cell->is_locally_owned()  )
	{

	  hp_fe_values[component]->reinit (cell);
	  unsigned int CO = cell->index();

	  const FEValues<dim> &fe_values = hp_fe_values[component]->get_present_fe_values ();
	  const Quadrature<dim>& quadrature_formula = fe_values.get_quadrature(); 
	  const FiniteElement<dim>& fe = cell->get_fe();
	  const unsigned int   n_q_points    = quadrature_formula.size();
	  const unsigned int   dofs_per_cell = fe.dofs_per_cell;

	  Vector<double> sigma_interior[alphadim];
	  sigma_interior[component].reinit(dofs_per_cell);
	  Vector<double> sigma_flux[alphadim];
	  sigma_flux[component].reinit(dofs_per_cell);

	  std::vector<unsigned int> local_dof_indices (dofs_per_cell);
	  cell->get_dof_indices (local_dof_indices);

	  const std::vector<double> &JxW = fe_values.get_JxW_values ();
	  double JxWsum = 0.0;
	  for (unsigned int q=0; q<n_q_points; ++q){
	    JxWsum += JxW[q];
	  }

	  soln_alpha[component] =  std::vector<double>(n_q_points);
	  soln_sigma[component] =  std::vector<Tensor<1,dim> >(n_q_points);

	  fe_values[*(alpha[component])].get_function_values(subdomain_solution[component],
							     soln_alpha[component]);

	  sigma_interior[component]=0.0;
	  for (unsigned int i=0; i<dofs_per_cell; ++i){
	    for (unsigned int q=0; q<n_q_points; ++q){

	      double add_it = 0.0;
	      //if(component == 0){add_it = 0.1;}

	      sigma_interior[component](i) +=  fe_values[*(sigma[component])].divergence(i,q) 
		* (soln_alpha[component][q] + add_it) * JxW[q];
	    }
	  }

	  sigma_flux[component] = 0.0;
	  for(unsigned int face_num=0; face_num < GeometryInfo<dim>::faces_per_cell; ++face_num){
	    double current_face_jump=0.0;
	    typename DoFHandler<dim>::face_iterator face=cell->face(face_num);
	    /* We might be on the domain boundary */
	    if ( face->at_boundary() && (face->boundary_indicator() == 0) ){

	      hp_fe_values_face[component]->reinit (cell,face_num);
	      const FEFaceValues<dim>& fe_values_face =
	    	hp_fe_values_face[component]->get_present_fe_values ();
	      const Quadrature<dim-1>& face_quadrature_formula = fe_values_face.get_quadrature();
	      const unsigned int n_q_points_face = face_quadrature_formula.size();

	      const std::vector<Point<dim> >& quadrature_point =
	    	fe_values_face.get_quadrature_points();

	      const std::vector<double> &JxW_face = fe_values_face.get_JxW_values ();
	      const std::vector<Point<dim> > &normals = fe_values_face.get_normal_vectors ();

	      soln_alpha_face[component] =  std::vector<double>(n_q_points_face);
	      fe_values_face[*(alpha[component])].get_function_values(subdomain_solution[component],
	    							      soln_alpha_face[component]);


	      if (component == 0){
		soln_alpha_face[1] =  std::vector<double>(n_q_points_face);
		fe_values_face[*(alpha[2])].get_function_laplacians(subdomain_solution[2],
								    soln_alpha_face[1]);
	      }


	      std::vector<double> alphas_boundary(n_q_points_face);
	      //std::vector<double> lap_boundary(n_q_points_face);
	      unsigned char boundary_index = face->boundary_indicator();
	      const WallBoundaryValues<dim>& wbv = WBV[boundary_index];
	      alphas_boundary = soln_alpha_face[component];
	      

	      if (component == 0){

		//lap_boundary = soln_alpha_face[1];
		wbv.value_list( quadrature_point, alphas_boundary, component, current_time);
		
	      }
	      else{
		
		wbv.value_list( quadrature_point, alphas_boundary, component, current_time);
		
	      }
	      
	      for (unsigned int i=0; i<dofs_per_cell; ++i){
	    	for (unsigned int q=0; q<fe_values_face.n_quadrature_points; ++q){
		  
		  if (component == 0){
		    
		      sigma_flux[component](i) += 0.5*(soln_alpha_face[component][q] + alphas_boundary[q])
			* ( fe_values_face[*(sigma[component])].value(i,q) * normals[q] ) * JxW_face[q];
		      
		  }
		    
		  
		  
		  else{
		    
		    sigma_flux[component](i) += 0.5*(soln_alpha_face[component][q] + alphas_boundary[q])
		      * ( fe_values_face[*(sigma[component])].value(i,q) * normals[q] ) * JxW_face[q];
		    
		  }
	    	}
	      }
	    }
	    else if ( face->at_boundary() && (face->boundary_indicator() == 15) ){

	      hp_fe_values_face[component]->reinit (cell,face_num);
	      const FEFaceValues<dim>& fe_values_face =
	    	hp_fe_values_face[component]->get_present_fe_values ();
	      const Quadrature<dim-1>& face_quadrature_formula = fe_values_face.get_quadrature();
	      const unsigned int n_q_points_face = face_quadrature_formula.size();

	      const std::vector<Point<dim> >& quadrature_point =
	    	fe_values_face.get_quadrature_points();

	      const std::vector<double> &JxW_face = fe_values_face.get_JxW_values ();
	      const std::vector<Point<dim> > &normals = fe_values_face.get_normal_vectors ();

	      soln_alpha_face[component] =  std::vector<double>(n_q_points_face);
	      fe_values_face[*(alpha[component])].get_function_values(subdomain_solution[component],
	    							      soln_alpha_face[component]);


	      if (component == 0){
		soln_alpha_face[1] =  std::vector<double>(n_q_points_face);
		fe_values_face[*(alpha[2])].get_function_laplacians(subdomain_solution[2],
								    soln_alpha_face[1]);
	      }


	      std::vector<double> alphas_boundary(n_q_points_face);
	      //std::vector<double> lap_boundary(n_q_points_face);
	      unsigned char boundary_index = face->boundary_indicator();
	      const WallBoundaryValues<dim>& wbv = WBV[boundary_index];
	      alphas_boundary = soln_alpha_face[component];

	      if (component == 0){

		//alphas_boundary = glob_min_ribbon_density;
		//alphas_boundary = total_density/ triangulation.n_global_active_cells();
		//lap_boundary = soln_alpha_face[1];
		wbv.value_list2( quadrature_point, alphas_boundary, component, glob_min_ribbon_density);
		//std::cout <<  glob_min_ribbon_density << std::endl;
		
	      }
	      else{
		
		wbv.value_list2( quadrature_point, alphas_boundary, component, current_time);
		
	      }
	      
	      for (unsigned int i=0; i<dofs_per_cell; ++i){

	
	    	for (unsigned int q=0; q<fe_values_face.n_quadrature_points; ++q){
		  
		  if (component == 0){
		    
		    sigma_flux[component](i) += 0.5*(soln_alpha_face[component][q] + alphas_boundary[q])
		      * ( fe_values_face[*(sigma[component])].value(i,q) * normals[q] ) * JxW_face[q]; 
		    //total_density/ triangulation.n_global_active_cells(); //alphas_boundary[q];
		    
		  }
		  else{
		    
		    sigma_flux[component](i) += 0.5*(soln_alpha_face[component][q] + alphas_boundary[q])
		      * ( fe_values_face[*(sigma[component])].value(i,q) * normals[q] ) * JxW_face[q];
		    
		  }
	    	}
	      }
	    }
	    else if (  face->at_boundary() && ( face->boundary_indicator() == 3 ||
					       face->boundary_indicator() == 4)){

	      //	      typename DoFHandler<dim>::cell_iterator neighbor = cell->neighbor(face_num);
	      //const unsigned int neighbor2 = cell->neighbor_of_neighbor(face_num);
	      hp_fe_values_face[component]->reinit(cell, face_num);
	      //hp_fe_values_neigh_face[component]->reinit(neighbor,neighbor2);

	      const FEFaceValues<dim>& fe_values_face = 
		hp_fe_values_face[component]->get_present_fe_values ();
	      //const FEFaceValues<dim>& fe_values_neigh_face = 
	      //hp_fe_values_neigh_face[component]->get_present_fe_values();

	      const std::vector<double> &JxW_face = fe_values_face.get_JxW_values ();
	      const std::vector<Point<dim> > &normals = fe_values_face.get_normal_vectors ();   

	      const Quadrature<dim-1>& face_quadrature_formula = fe_values_face.get_quadrature();
	      const unsigned int n_q_points_face = face_quadrature_formula.size();

	      soln_alpha_face[component] =  std::vector<double>(n_q_points_face);
	      fe_values_face[*(alpha[component])].get_function_values(subdomain_solution[component],
								      soln_alpha_face[component]);

	      

	      //soln_alpha_neigh_face[component] =  std::vector<double>(n_q_points_face);
	      //fe_values_neigh_face[*(alpha[component])].get_function_values(subdomain_solution[component],
	      //							    soln_alpha_neigh_face[component]);

	      for (unsigned int i=0; i<dofs_per_cell; ++i){
		for (unsigned int q=0; q<n_q_points_face; ++q){

		  double add_it = 0.0;
		  //if(component == 0){add_it = 0.1;}

	      	  sigma_flux[component](i) += ( 0.5 * ( soln_alpha_face[component][q] 
							+ soln_drawX[component][CO][q] + 2.0*add_it) 
						* (  fe_values_face[*(sigma[component])].value(i,q) 
						     * normals[q] ) ) * JxW_face[q];
	      	}
	      }


	    }
	    else if (  face->at_boundary() && ( face->boundary_indicator() == 1 ||
						face->boundary_indicator() == 2)){
	      
	      //	      typename DoFHandler<dim>::cell_iterator neighbor = cell->neighbor(face_num);
	      //const unsigned int neighbor2 = cell->neighbor_of_neighbor(face_num);
	      hp_fe_values_face[component]->reinit(cell, face_num);
	      //hp_fe_values_neigh_face[component]->reinit(neighbor,neighbor2);

	      const FEFaceValues<dim>& fe_values_face = 
		hp_fe_values_face[component]->get_present_fe_values ();
	      //const FEFaceValues<dim>& fe_values_neigh_face = 
	      //hp_fe_values_neigh_face[component]->get_present_fe_values();

	      const std::vector<double> &JxW_face = fe_values_face.get_JxW_values ();
	      const std::vector<Point<dim> > &normals = fe_values_face.get_normal_vectors ();   

	      const Quadrature<dim-1>& face_quadrature_formula = fe_values_face.get_quadrature();
	      const unsigned int n_q_points_face = face_quadrature_formula.size();

	      soln_alpha_face[component] =  std::vector<double>(n_q_points_face);
	      fe_values_face[*(alpha[component])].get_function_values(subdomain_solution[component],
								      soln_alpha_face[component]);
	      //soln_alpha_neigh_face[component] =  std::vector<double>(n_q_points_face);
	      //fe_values_neigh_face[*(alpha[component])].get_function_values(subdomain_solution[component],
	      //							    soln_alpha_neigh_face[component]);

	      for (unsigned int i=0; i<dofs_per_cell; ++i){
		for (unsigned int q=0; q<n_q_points_face; ++q){

		  double add_it = 0.0;
		  //if(component == 0){add_it = 0.1;}
		  
	      	  sigma_flux[component](i) += ( 0.5 * ( soln_alpha_face[component][q] 
							+ soln_drawY[component][CO][q] + 2.0*add_it) 
						* (  fe_values_face[*(sigma[component])].value(i,q) 
						     * normals[q] ) ) * JxW_face[q];
	      	}
	      }
	      
	      
	    }
	    /* Else we are on an interior element boundary with neighbors */
	    else {

	      typename DoFHandler<dim>::cell_iterator neighbor = cell->neighbor(face_num);
	      const unsigned int neighbor2 = cell->neighbor_of_neighbor(face_num);
	      hp_fe_values_face[component]->reinit(cell, face_num);
	      hp_fe_values_neigh_face[component]->reinit(neighbor,neighbor2);

	      const FEFaceValues<dim>& fe_values_face = 
		hp_fe_values_face[component]->get_present_fe_values ();
	      const FEFaceValues<dim>& fe_values_neigh_face = 
		hp_fe_values_neigh_face[component]->get_present_fe_values();

	      const std::vector<double> &JxW_face = fe_values_face.get_JxW_values ();
	      const std::vector<Point<dim> > &normals = fe_values_face.get_normal_vectors ();   

	      const Quadrature<dim-1>& face_quadrature_formula = fe_values_face.get_quadrature();
	      const unsigned int n_q_points_face = face_quadrature_formula.size();

	      soln_alpha_face[component] =  std::vector<double>(n_q_points_face);
	      fe_values_face[*(alpha[component])].get_function_values(subdomain_solution[component],
								      soln_alpha_face[component]);
	      soln_alpha_neigh_face[component] =  std::vector<double>(n_q_points_face);
	      fe_values_neigh_face[*(alpha[component])].get_function_values(subdomain_solution[component],
									    soln_alpha_neigh_face[component]);

	      for (unsigned int i=0; i<dofs_per_cell; ++i){
		for (unsigned int q=0; q<n_q_points_face; ++q){
		  double add_it = 0.0;
		  //if(component == 0){add_it = 0.1;}

	      	  sigma_flux[component](i) += ( 0.5 * ( soln_alpha_face[component][q] 
							+ soln_alpha_neigh_face[component][q] + 2.0*add_it) 
						* (  fe_values_face[*(sigma[component])].value(i,q) 
						     * normals[q] ) ) * JxW_face[q];
	      	}
	      }
	    }
	  }

	  cell->get_dof_indices(local_dof_indices);

	  Vector<double> temp_update[alphadim];
	  temp_update[component].reinit(dofs_per_cell);
	  
	  for (unsigned int i=0; i<dofs_per_cell; ++i){
	    temp_update[component](i) = sigma_flux[component](i) - sigma_interior[component](i);
	  }

	  Vector<double> projected(dofs_per_cell);
	  mapinfo[component][0].localInverseMassMatrix.vmult(projected,temp_update[component]);
	  projected /= JxWsum;
	  parahyp_constraints[component].distribute_local_to_global (projected,
								     local_dof_indices,
								     naive_subdomain_solution[component]);
       	}
    
    naive_subdomain_solution[component].compress(VectorOperation::add);
    subdomain_solution[component].block(1) = naive_subdomain_solution[component].block(1);
    
  }
  
  for (unsigned int component=0; component< alphadim; ++component){
    delete hp_fe_values[component];
    delete hp_fe_values_face[component];
    delete hp_fe_values_neigh_face[component];
  }

}
