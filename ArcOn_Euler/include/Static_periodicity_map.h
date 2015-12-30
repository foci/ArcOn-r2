template <int dim>
void arcOn<dim>::static_periodicity_map(SolutionVector& subdomain_solution, double delta_t) { 

  (void) delta_t;

  std::vector< FEValues<dim>*  > fe_vals;
  std::vector< FEFaceValues<dim>*  > fe_vals_face;

  for (unsigned int component=0; component< alphadim; ++component){

    UpdateFlags updateflags=  update_values | update_gradients | update_quadrature_points | update_JxW_values;
    fe_vals.push_back (new FEValues<dim>(*(fe_collection[component]), 
					      *(quadrature_collection[component]), updateflags));
    fe_vals_face.push_back (new FEFaceValues<dim>(*(fe_collection[component]), 
						  *(face_quadrature_collection[component]),  
						  updateflags | update_normal_vectors ));
    
    typename DoFHandler<dim>::active_cell_iterator 
      celli = dof_handler[component]->begin_active(), 
      endci = dof_handler[component]->end();
    
    for (; celli!=endci; ++celli){
      if (!(celli->is_artificial())){

	fe_vals[component]->reinit (celli);
	const FiniteElement<dim>& fe = celli->get_fe();
	const unsigned int   dofs_per_cell = fe.dofs_per_cell;
	std::vector<unsigned int> local_dof_indices (dofs_per_cell);
	celli->get_dof_indices (local_dof_indices);

	for(unsigned int face_numi=0;
	    face_numi < GeometryInfo<dim>::faces_per_cell; ++face_numi){
	  
	  typename DoFHandler<dim>::face_iterator facei = celli->face(face_numi);
	  
	  if ( facei->at_boundary() 
	       && facei->boundary_indicator() != 0 && facei->boundary_indicator() > 2)
	    {
	      
	      fe_vals_face[component]->reinit(celli, face_numi);
	      const FEFaceValues<dim>& fe_values_face =
		fe_vals_face[component]->get_present_fe_values ();
	      
	      const Quadrature<dim-1>& face_quadrature_formula = fe_values_face.get_quadrature();
	      const std::vector<Point<dim> >& quadrature_point = fe_values_face.get_quadrature_points();
	      const unsigned int n_q_points_face = face_quadrature_formula.size();

	      periodic_basisX[component][Cell_matchX[celli->index()]] = FullMatrix<double>(dofs_per_cell,fe_values_face.n_quadrature_points);

	      periodic_gradX[component][Cell_matchX[celli->index()]].resize(n_q_points_face);

	      periodic_dof_indX[component][Cell_matchX[celli->index()]] =  local_dof_indices;
	      
	      /* soln_draw[component][Cell_match[celli->index()]] = std::vector<double>(n_q_points_face); */
	      /* prev_soln_alpha_face[component] =  std::vector<double>(n_q_points_face); */
	      /* fe_values_face[*(alpha[component])].get_function_values(subdomain_solution[component], */
	      /* 							      prev_soln_alpha_face[component]);   */

	      /* grad_draw[component][Cell_match[celli->index()]] = std::vector<Tensor<1,dim> >(n_q_points_face); */
	      /* prev_soln_sigma_face[component] = std::vector<Tensor<1,dim> >(n_q_points_face); */
	      /* fe_values_face[*(sigma[component])].get_function_values(subdomain_solution[component], */
	      /* 							      prev_soln_sigma_face[component]);   */


	      for (unsigned int q=0; q<n_q_points_face; ++q){

		periodic_gradX[component][Cell_matchX[celli->index()]][q]  =  std::vector<Tensor<1,dim> >(dofs_per_cell); 

		for (unsigned int i=0; i<dofs_per_cell; ++i){
		  
		  periodic_basisX[component][Cell_matchX[celli->index()]](i,q) = (fe_values_face[*(alpha[component])].value(i,q));
		  periodic_gradX[component][Cell_matchX[celli->index()]][q][i]  = (fe_values_face[*(alpha[component])].gradient(i,q));
		
		}
		
		//periodic_dof_ind[component][Cell_match[celli->index()]][i] = local_dof_indices[i];
		
	      }

	    }
	  else if ( facei->at_boundary() 
	       && facei->boundary_indicator() != 0 && facei->boundary_indicator() < 3)
	    {
	      
	      fe_vals_face[component]->reinit(celli, face_numi);
	      const FEFaceValues<dim>& fe_values_face =
		fe_vals_face[component]->get_present_fe_values ();
	      
	      const Quadrature<dim-1>& face_quadrature_formula = fe_values_face.get_quadrature();
	      const std::vector<Point<dim> >& quadrature_point = fe_values_face.get_quadrature_points();
	      const unsigned int n_q_points_face = face_quadrature_formula.size();

	      periodic_basisY[component][Cell_matchY[celli->index()]] = FullMatrix<double>(dofs_per_cell,fe_values_face.n_quadrature_points);

	      periodic_gradY[component][Cell_matchY[celli->index()]].resize(n_q_points_face);

	      periodic_dof_indY[component][Cell_matchY[celli->index()]] =  local_dof_indices;
	      
	      for (unsigned int q=0; q<n_q_points_face; ++q){

		periodic_gradY[component][Cell_matchY[celli->index()]][q]  =  std::vector<Tensor<1,dim> >(dofs_per_cell); 

		for (unsigned int i=0; i<dofs_per_cell; ++i){
		  
		  periodic_basisY[component][Cell_matchY[celli->index()]](i,q) = (fe_values_face[*(alpha[component])].value(i,q));
		  periodic_gradY[component][Cell_matchY[celli->index()]][q][i]  = (fe_values_face[*(alpha[component])].gradient(i,q));
		
		}
	      }
	    }


	}
      }
    }
    delete fe_vals[component];
    delete fe_vals_face[component];
  }
}
