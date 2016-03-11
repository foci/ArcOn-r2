/* Error Analysis of the solution */
template <int dim>
void arcOn<dim>::compute_l2_interpolation(SolutionVector& subdomain_solution)
{
  std::vector< FEValues<dim>* > hp_fe_values;
  
  for (unsigned int component=0; component< alphadim; ++component){
    
    UpdateFlags updateflags=  update_values | 
      update_gradients | update_quadrature_points | update_JxW_values;
    
    hp_fe_values.push_back (new FEValues<dim>(*(fe_collection[component]), 
					      *(quadrature_collection[component]), 
					      updateflags));
  /* for (unsigned int component=0; component< alphadim; ++component){ */
  /*   naive_subdomain_solution[component] = subdomain_solution[component];  */
  /*   naive_L2_error_method[component] = 0.0; */

  /*   interpolate_base[component] = subdomain_solution[component]; */
  /*   L2_interpolate_active[component] = 0.0; */

  /*   FETools::interpolate (*(dof_handler[component]), interpolate_base[component] ,  */
  /*   			  *(tdof_handler[component]), L2_interpolate_active[component]); */
  /*   /\* This is not fully implemented yet for modal DG error analysis in parallel *\/ */
  /*   //FETools::project_dg (*(dof_handler[component]), naive_subdomain_solution[component] ,  */
  /*   //			 *(tdof_handler[component]), L2_error_method[component]); */
  /* } */
    
    naive_interpolation_error[component] = 0.0;
    
    InitialValues<dim> IV;
    double total_volume = 0.0;
    typename DoFHandler<dim>::active_cell_iterator 
      cell = dof_handler[component]->begin_active(), 
      endc = dof_handler[component]->end();
    for(;
	cell!=endc;
	++cell )
      if (cell->is_locally_owned() )
	{
	  hp_fe_values[component]->reinit (cell);
	  unsigned int CO = cell->index();
	  const FEValues<dim> &fe_values = hp_fe_values[component]->get_present_fe_values ();
	  const Quadrature<dim>& quadrature_formula = fe_values.get_quadrature(); 
          const std::vector<Point<dim> >& quadrature_point = fe_values.get_quadrature_points();
	  const unsigned int   dofs_per_cell = fe_values.dofs_per_cell;
	  const unsigned int   n_q_points    = quadrature_formula.size();

	  const std::vector<double>& quad_weights = quadrature_formula.get_weights();

	  const FiniteElement<dim>& fe = cell->get_fe();

	  std::vector<unsigned int> local_dof_indices (dofs_per_cell);

	  cell->get_dof_indices (local_dof_indices);
	  const std::vector<double> &JxW = fe_values.get_JxW_values ();
	  double JxWsum = 0.0;
	  for (unsigned int q=0; q<n_q_points; ++q){
	    JxWsum += JxW[q];
	  }
	  total_volume += JxWsum;


	  for(unsigned int k=0;k<alphadim;k++){
	    prev_soln_alpha[k] =  std::vector<double>(n_q_points);
	    soln_alpha[k] =  std::vector<double>(n_q_points);
	  }

	  fe_values[*(alpha[component])].get_function_values(init_solution[component],prev_soln_alpha[component]);
	  //fe_values[*(alpha[component])].get_function_values(L2_error_interpolant[component],soln_alpha[component]);

	  FullMatrix<double> MassMatrix(dofs_per_cell,dofs_per_cell);
	  MassMatrix = 0.0;
	  for (unsigned int i=0; i<dofs_per_cell; ++i){
	    for (unsigned int j=0; j<dofs_per_cell; ++j){
	      for (unsigned int q=0; q<n_q_points; ++q){
	  	MassMatrix(i,j)+= fe_values[*(alpha[component])].value(i,q)
	  	  * fe_values[*(alpha[component])].value(j,q) *  quad_weights[q];
	      }
	    }
	  }
	  for (unsigned int i=0; i<dofs_per_cell; ++i){
	    for (unsigned int j=0; j<dofs_per_cell; ++j){
	      for (unsigned int q=0; q<n_q_points; ++q){
	  	MassMatrix(i,j)+= fe_values[*(sigma[component])].value(i,q)
	  	  * fe_values[*(sigma[component])].value(j,q) *  quad_weights[q];
	      }
	    }
	  }

	  FullMatrix<double> InverseMassMatrix(dofs_per_cell,dofs_per_cell);
	  InverseMassMatrix.invert(MassMatrix);
	  
	  const std::vector< Point<dim> >& quadrature_points = fe_values.get_quadrature_points();
	  Vector<double> alpha_sum[alphadim];
	  alpha_sum[component].reinit(dofs_per_cell);
	  Vector<double> qpvalue(alphadim);
	  
 	  for(unsigned int q=0;q<n_q_points;q++){
	    IV.vector_value(quadrature_points[q],qpvalue);
	    for (unsigned int i=0; i<dofs_per_cell; ++i){
	      //if(component ==2){
	      //	std::cout << "base_val = " << qpvalue(component) << "integrated val = " << (prev_soln_alpha[component])[q]
	      //  *fe_values[*(alpha[component])].value(i,q) << std::endl;}
	      
	      //alpha_sum[component](i) += std::pow( (qpvalue( component ) - (prev_soln_alpha[component])[q] )* fe_values[*(alpha[component])].value(i,q), 1.0) * JxW[q];
	      alpha_sum[component](i) += ( qpvalue( component ) - (prev_soln_alpha[component])[q] ) * JxW[q];
	      
	      //alpha_sum[component](i) += std::pow( ( (soln_alpha[component])[q] - (prev_soln_alpha[component])[q] )* fe_values[*(alpha[component])].value(i,q), 1.0) * JxW[q];


	      /* for (unsigned int j=0; j<dofs_per_cell; ++i){ */
	      /* alpha_sum[component](i) -= (prev_soln_alpha[component])[q] *fe_values[*(alpha[component])].value(i,q) */
	      /* 	*fe_values[*(alpha[component])].value(j,q)* JxW[q]; */
	      /* } */
	      

	      //alpha_sum[component](i) += ( qpvalue( component )*fe_values[*(alpha[component])].value(i,q)
	      //				   - preprojected(i)[q]  ) * JxW[q];

	      // alpha_sum[component](i) += std::pow((prev_soln_alpha[0])[q]-(prev_soln_alpha[component])[q],1.0)
	      //	*fe_values[*(alpha[component])].value(i,q)*JxW[q];
	    }
	  }

	  //Vector<double> projected(dofs_per_cell);
	  //	  InverseMassMatrix.vmult(projected,alpha_sum[component]);
	  //projected /= JxWsum;
	  //alpha_sum[component] /= JxWsum;

	  if (component == 0){
            density_constraints.distribute_local_to_global (projected,
                                                            local_dof_indices,
                                                            naive_interpolation_error[component]);
          }
          else if (component == 1) {
            vorticity_constraints.distribute_local_to_global (projected,
                                                              local_dof_indices,
                                                              naive_interpolation_error[component]);
          }
          else if (component == 2) {
            potential_constraints.distribute_local_to_global (projected,
                                                              local_dof_indices,
                                                              naive_interpolation_error[component]);
          }

	  /* parahyp_constraints[component].distribute_local_to_global (alpha_sum[component], */
	  /* 							     local_dof_indices, */
	  /* 							     naive_interpolation_error[component]); */
	  
	}
    
    naive_interpolation_error[component].compress(VectorOperation::add);
    interpolation_error[component].block(0) = naive_interpolation_error[component].block(0);

  }
  


  for (unsigned int component=0; component< alphadim; ++component){
    delete hp_fe_values[component];
  }
  //return total_volume;
}
 
