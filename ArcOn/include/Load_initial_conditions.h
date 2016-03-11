/* Load the initial Conditions using an L2-projection into the basis */
template <int dim>
void arcOn<dim>::load_initial_conditions(SolutionVector& subdomain_solution, unsigned int init_flag)
{
  std::vector< FEValues<dim>* > hp_fe_values;
  std::vector< double > max_val_loc;
  max_val_loc.resize(alphadim);

  //-------------------------------------------------------------------------------------------------------------------
  //Modon stuff

  double beta;
  double a,c,lambda;

  // double c = double(*argv[0]);                                                                                                                                                                                                            
  // //double lambda = std::stof (*argv[1]);                                                                                                                                                                                                 
  // double lambda = double(*argv[1]);                                                                                                                                                                                                       
  // double beta = double(*argv[2]);                                                                                                                                                                                                         
  /* pcout << std::argv[3] <<std::endl; */
  /* std::istringstream s(std::argv[3]); */

  /* if (!(s >> c)) */
  /*   cerr << "Invalid number " << std::argv[3] << '\n'; */

  /* std::istringstream ss(std::argv[1]); */

  /* if (!(ss >> lambda)) */
  /*   cerr << "Invalid number " << std::argv[1] << '\n'; */

  /* std::istringstream s2(std::argv[2]); */

  /* if (!(s2 >> beta)) */
  /*   cerr << "Invalid number " << std::argv[2] << '\n'; */

  /* pcout << beta <<" "<<lambda/beta<<std::endl; */

  //set the radius a to 1                                                                                                                                                                                                                    
  //Ullmann map
  double x_sol = .3;
  double rho_s = .2;
  double y_length = 314.0;
  double mhmm = 2.0;
  double eps = 0.2;



  //a = .1;
  double xlen = 1.0;
  double scaling = 1.0;
  //lambda = 
  
  beta = 1.0;
  a = .1*xlen;//*((mesh->dx)[0][0]);
  c = scaling*(xlen/1.0);
  lambda = .10/xlen;//;*(
  //we get to pick c at will                                                                                                                                                                                                                 
  //c = 50;                                                                                                                                                                                                                                   //set kappa based on c                                                                                                                                                                                                                     
  double kappa = a*std::sqrt(lambda*beta)/c;
  //pcout << "kappa = " << kappa << std::endl;
  double gamma;
  double pi = 3.1415926535897932384626433832795;


  gamma = newton_root(pi,kappa);   
  //pcout << "gamma = " << gamma;

  //c = 1.0;                                                                                                                                                                                                                                 
  // have to pick some x and y values                                                                                                                                                                                                       

  //double x = .3;
  //double y = .2;

  //pcout << "modon: "<<kappa<<" "<<gamma<<" "<<modon(.3,.1,kappa,gamma)<<std::endl;

  //double phi, n;
  //phi = modon(.3,.1,kappa,gamma);
  //n =  (lambda/c)*(phi - c*x);


  //--------------------------------------------------------------------------------------------------------------------
  for (unsigned int component=0; component< alphadim; ++component){
    
    UpdateFlags updateflags=  update_values | 
      update_gradients | update_quadrature_points | update_JxW_values | update_hessians;
    
    hp_fe_values.push_back (new FEValues<dim>(*(fe_collection[component]), 
					      *(quadrature_collection[component]), 
					      updateflags));
    
    InitialValues<dim> IV;
    double total_volume = 0.0;
    max_val_loc[component] = 1.0;

    if (init_flag == 1){
      naive_subdomain_solution[component] = 0.0;}
    

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

	  if (component != 2 && init_flag == 2){
	    prev_soln_alpha[0] =  std::vector<double>(n_q_points);
	    prev_soln_alpha[1] =  std::vector<double>(n_q_points);}

	  cell->get_dof_indices (local_dof_indices);
	  const std::vector<double> &JxW = fe_values.get_JxW_values ();
	  double JxWsum = 0.0;
	  for (unsigned int q=0; q<n_q_points; ++q){
	    JxWsum += JxW[q];
	  }
	  total_volume += JxWsum;

	  if (component != 2 && init_flag == 2){
	    fe_values[*(alpha[2])].get_function_values(subdomain_solution[2],
						       prev_soln_alpha[0]);
	    fe_values[*(alpha[2])].get_function_laplacians(subdomain_solution[2],
							   prev_soln_alpha[1]);
	  }

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
	      if(component == 0){ 
		//conserved_mass+= scaling * qpvalue( component )
		// *fe_values[*(alpha[component])].value(i,q)*JxW[q]; 
		//if (quadrature_point[q](0) <= 0.0){
		  //alpha_sum[component](i) -= 0.02*quadrature_point[q][0]*fe_values[*(alpha[component])].value(real_index,q)*JxW[q];
		//}
	      }

	      alpha_sum[component](i) += qpvalue( component )
		*fe_values[*(alpha[component])].value(i,q)*JxW[q];

	      /* //Ullmann */
	      /* if (component == 1 && init_flag ==  1){ */
	      /* 	//alpha_sum[component](i) += qpvalue( component ) */
	      /* 	//\*fe_values[*(alpha[component])].value(i,q)*JxW[q]; */
		
	      /* 	double x_in = quadrature_point[q](0)/5.0; */
	      /* 	double y_in =  quadrature_point[q](1); */
	      
	      /* 	//pcout << "x_pre = " << x_in << ", y_pre = " << y_in << std::endl; */

	      /* 	double Ulmap = UllmannMap(x_in, 1.0, y_in, y_length, x_sol, eps, mhmm); */
	      /* 	if ( UllmannMap(x_in, 1.0, y_in, y_length, x_sol, eps, mhmm) < 0.0 ) {Ulmap = 0.0;} */

	      /* 	alpha_sum[component](i) +=  (2.0*rho_s / Ulmap) */
	      /* 	  *fe_values[*(alpha[component])].value(i,q)*JxW[q]; */

	      /* 	//pcout << "modon = " << modon(quadrature_point[q](0)-5.0, quadrature_point[q](1)-5.0, kappa, gamma,c,a)  << std::endl; */
		
	      /* 	/\* else if  (component == 1){ *\/ */
	      /* 	/\*   alpha_sum[component](i) += modon *\/ */
	      /* 	/\*     *fe_values[*(alpha[component])].value(i,q)*JxW[q]; *\/ */
	      /* 	/\* } *\/ */
	      /* } */



	      /* else if (component == 0 && init_flag == 2){ */

	      /* 	double phi = prev_soln_alpha[0][q]; */

	      /* 	alpha_sum[component](i) +=  (lambda/c)*(phi - c*quadrature_point[q](0))  */
	      /* 	  *fe_values[*(alpha[component])].value(i,q)*JxW[q]; */
		
	      /* } */
	      /* else if (  (component == 1 && init_flag == 2)){ */

	      /* 	double lap_phi = prev_soln_alpha[1][q]; */
	      /* 	alpha_sum[component](i) +=  lap_phi */
	      /* 	  *fe_values[*(alpha[component])].value(i,q)*JxW[q]; */

	      /* } */
		//alpha_sum[component](i) += modon
		//*fe_values[*(alpha[component])].value(i,q)*JxW[q];
	
	      
	      
	      //else { alpha_sum[component](i) += 0.0;}
	      //}
	    }
	    if(component < 2) max_val_loc[component] = std::max( std::abs(qpvalue(component)), max_val_loc[component] ); 
	  }

	  Vector<double> projected(dofs_per_cell);
	  InverseMassMatrix.vmult(projected,alpha_sum[component]);
	  projected /= JxWsum;

	  if (component == 0){
            density_constraints.distribute_local_to_global (projected,
                                                            local_dof_indices,
                                                            naive_subdomain_solution[component]);
          }
          else if (component == 1) {
            vorticity_constraints.distribute_local_to_global (projected,
                                                              local_dof_indices,
                                                              naive_subdomain_solution[component]);
          }
	  else if (component == 2) {
            potential_constraints.distribute_local_to_global (projected,
                                                              local_dof_indices,
                                                              naive_subdomain_solution[component]);
          }	  
	  /* parahyp_constraints[component].distribute_local_to_global (projected, */
	  /* 							     local_dof_indices, */
	  /* 							     naive_subdomain_solution[component]); */
	  
	}

    naive_subdomain_solution[component].compress(VectorOperation::add);
    subdomain_solution[component].block(0) = naive_subdomain_solution[component].block(0);

    /* 	  Vector<double> projected(dofs_per_cell); */
    /* 	  projected = 0.0; */
    /* 	  if (init_flag < 2 && component == 2 ){ */
	    
    /* 	    InverseMassMatrix.vmult(projected,alpha_sum[component]); */
    /* 	    projected /= JxWsum; */
    /* 	    parahyp_constraints[component].distribute_local_to_global (projected, */
    /* 								       local_dof_indices, */
    /* 								       naive_subdomain_solution[component]); */
	    
    /* 	  } */
    /* 	  else if (component ==0 && init_flag == 2){ */
    /* 	      Vector<double> projected(dofs_per_cell); */
    /* 	      InverseMassMatrix.vmult(projected,alpha_sum[component]); */
    /* 	      projected /= JxWsum; */
    /* 	      parahyp_constraints[component].distribute_local_to_global (projected, */
    /* 									 local_dof_indices, */
    /* 									 naive_subdomain_solution[component]); */
    /* 	    } */
    /* 	  else if (component ==1 && init_flag == 2){ */
	    
    /* 	    Vector<double> projected(dofs_per_cell); */
    /* 	    InverseMassMatrix.vmult(projected,alpha_sum[component]); */
    /* 	    projected /= JxWsum; */
    /* 	    parahyp_constraints[component].distribute_local_to_global (projected, */
    /* 	    							       local_dof_indices, */
    /* 	    							       naive_subdomain_solution[component]); */
    /* 	  } */

	  
    /* 	} */
    
    /* if (init_flag < 2 && component == 2 ){ */
    /*   naive_subdomain_solution[component].compress(VectorOperation::insert); */
    /*   subdomain_solution[component].block(0) = naive_subdomain_solution[component].block(0); */
    /* } */

       
    /* if (init_flag == 2 && component !=2 ){ */ 
   /*   naive_subdomain_solution[component].compress(VectorOperation::insert); */
    /*   subdomain_solution[component].block(0) = naive_subdomain_solution[component].block(0); */
    /* } */

  }
  
  MPI_Allreduce(&max_val_loc[0], &gmax_val_loc[0], 1, MPI_DOUBLE, MPI_MAX, mpi_communicator);
  MPI_Allreduce(&max_val_loc[1], &gmax_val_loc[1], 1, MPI_DOUBLE, MPI_MAX, mpi_communicator);

  //pcout << gmax_val_loc[0] << ", gmax_val_loc [1] = " << gmax_val_loc[1] << std::endl;

  for (unsigned int component=0; component< alphadim; ++component){
    delete hp_fe_values[component];
  }
  //return total_volume;
}
 
