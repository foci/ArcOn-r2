template <int dim>
void arcOn<dim>::periodicity_map(SolutionVector& subdomain_solution, double delta_t) { 

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
	for(unsigned int face_numi=0;
	    face_numi < GeometryInfo<dim>::faces_per_cell; ++face_numi){
	
	  typename DoFHandler<dim>::face_iterator facei = celli->face(face_numi);
	
	  if ( facei->at_boundary() 
	       && facei->boundary_indicator() != 0 && facei->boundary_indicator() > 2)
	    {
	      //   && celli != cellj ) {			    
	      fe_vals_face[component]->reinit(celli, face_numi);
	      const FEFaceValues<dim>& fe_values_face =
		fe_vals_face[component]->get_present_fe_values ();
	    
	      const Quadrature<dim-1>& face_quadrature_formula = fe_values_face.get_quadrature();
	      const std::vector<Point<dim> >& quadrature_point = fe_values_face.get_quadrature_points();
	      const unsigned int n_q_points_face = face_quadrature_formula.size();

	      soln_drawX[component][Cell_matchX[celli->index()]] = std::vector<double>(n_q_points_face);
	      prev_soln_alpha_face[component] =  std::vector<double>(n_q_points_face);
	      fe_values_face[*(alpha[component])].get_function_values(subdomain_solution[component],
								      prev_soln_alpha_face[component]);  

	      grad_drawX[component][Cell_matchX[celli->index()]] = std::vector<Tensor<1,dim> >(n_q_points_face);
	      prev_soln_sigma_face[component] = std::vector<Tensor<1,dim> >(n_q_points_face);
	      fe_values_face[*(sigma[component])].get_function_values(subdomain_solution[component],
								      prev_soln_sigma_face[component]);  

	      soln_drawX[component][Cell_matchX[celli->index()]] = prev_soln_alpha_face[component];
	      grad_drawX[component][Cell_matchX[celli->index()]] = prev_soln_sigma_face[component];
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
	      
	      soln_drawY[component][Cell_matchY[celli->index()]] = std::vector<double>(n_q_points_face);
	      prev_soln_alpha_face[component] =  std::vector<double>(n_q_points_face);
	      fe_values_face[*(alpha[component])].get_function_values(subdomain_solution[component],
								      prev_soln_alpha_face[component]);  

	      grad_drawY[component][Cell_matchY[celli->index()]] = std::vector<Tensor<1,dim> >(n_q_points_face);
	      prev_soln_sigma_face[component] = std::vector<Tensor<1,dim> >(n_q_points_face);
	      fe_values_face[*(sigma[component])].get_function_values(subdomain_solution[component],
								      prev_soln_sigma_face[component]);  

	      soln_drawY[component][Cell_matchY[celli->index()]] = prev_soln_alpha_face[component];
	      grad_drawY[component][Cell_matchY[celli->index()]] = prev_soln_sigma_face[component];

	    }
	}
      }
    }
    delete fe_vals[component];
    delete fe_vals_face[component];
  }
}
