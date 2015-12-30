/* This is the BDS slopelimiter ... this needs to be fixed (easy) */
template <int dim>
void arcOn<dim>::slopelimiter( SolutionVector& subdomain_solution, int step )
{

  std::vector< FEValues<dim>* > hp_fe_values;
  //std::vector< FEValues<dim>* > lhp_fe_values;
  

//std::vector< FEFaceValues<dim>* > hp_fe_values_face;
  //std::vector< FEFaceValues<dim>* > hp_fe_values_neigh_face;

  for (unsigned int component=0; component< alphadim; ++component){

    UpdateFlags updateflags=  update_values | update_gradients | update_quadrature_points | update_JxW_values;
    
    hp_fe_values.push_back (new FEValues<dim>(*(fe_collection[component]), *(quadrature_collection[component]), updateflags));
    //lhp_fe_values.push_back (new FEValues<dim>(*(lfe_collection[component]), *(lquadrature_collection[component]), updateflags));
 
   //hp_fe_values_face.push_back (new FEFaceValues<dim>(*(fe_collection[component]), *(face_quadrature_collection[component]),  updateflags | update_normal_vectors ));
    //hp_fe_values_neigh_face.push_back (new FEFaceValues<dim>(*(fe_collection[component]), *(face_quadrature_collection[component]),  updateflags | update_normal_vectors ));

    local_solution[component] = subdomain_solution[component];

    // Let's move into linear spaces

    //    l_soln = local_solution[component];
    //i_soln = 0.0;
    //r_soln = 0.0;
    //FETools::interpolate (*(dof_handler[component]), l_soln, *(idof_handler[component]), i_soln);
    //FETools::interpolate (*(ldof_handler[component]), i_soln, *(dof_handler[component]), r_soln);
 
    //Vector<double> newmax(num_vert);
    //Vector<double> newmin(num_vert);

    //newmin =  9999.0;
    //newmax = -9999.0;



    //std::map<int, double> oldquad;
    //std::vector<double>           syncopated[alphadim];
    //std::vector< Tensor<1,dim> >  syncopated2[alphadim];
    //Vector local_verts;
    //first get the size of the local object

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler[component]->begin_active(), endc = dof_handler[component]->end();
    for(;
	cell!=endc;
	++cell	)
      if (cell->is_locally_owned()  )
	{
	  hp_fe_values[component]->reinit (cell);
	  for(unsigned int vi = 0; vi < GeometryInfo<dim>::vertices_per_cell; ++vi){
	    int v_dex = cell->vertex_index(vi);
	      
	    newmax[v_dex] = -9999.0;
	    newmin[v_dex] =  9999.0; 

	    newmax0[v_dex] = -9999.0;
	    newmin0[v_dex] =  9999.0; 

	    newmax1[v_dex] = -9999.0;
	    newmin1[v_dex] =  9999.0; 
	  }
	}
      else if (cell->is_ghost()  )
	{
	  hp_fe_values[component]->reinit (cell);
	  for(unsigned int vi = 0; vi < GeometryInfo<dim>::vertices_per_cell; ++vi){
	    int v_dex = cell->vertex_index(vi);
	    //std::cout << "vdex = " << v_dex << ", " << vi << ", " <<  GeometryInfo<dim>::vertices_per_cell << std::endl;
	    newmax[v_dex] = -9999.0; 
	    newmin[v_dex] =  9999.0;

	    newmax0[v_dex] = -9999.0; 
	    newmin0[v_dex] =  9999.0;

	    newmax1[v_dex] = -9999.0; 
	    newmin1[v_dex] =  9999.0;
	  }
	}

    //First find the stencil extrema
    //typename DoFHandler<dim>::active_cell_iterator cell = dof_handler[component]->begin_active(), endc = dof_handler[component]->end();
    cell = dof_handler[component]->begin_active(), endc = dof_handler[component]->end();
    for(;
	cell!=endc;
	++cell	)
      if ( cell->is_locally_owned() )
	{

	  hp_fe_values[component]->reinit (cell);
	  const FEValues<dim> &fe_values = hp_fe_values[component]->get_present_fe_values ();
	  const Quadrature<dim>& quadrature_formula = fe_values.get_quadrature(); 
	  const FiniteElement<dim>& fe = cell->get_fe();
	  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
	  const unsigned int   n_q_points    = quadrature_formula.size();
	  
	  const std::string fe_name = fe.get_name();
	  pinfo[component] = pInfoFind(fe.get_name(),component);
	  
	  std::vector<unsigned int> local_dof_indices (dofs_per_cell); 
	  const unsigned int alpha_dofs_per_cell = pinfo[component].alpha_dofs_per_cell;
	  const unsigned int sigma_dofs_per_cell = pinfo[component].sigma_dofs_per_cell;
	  
	  double cell_avg = 0;

	  double cell_avg0 = 0;
	  double cell_avg1 = 0;

	  syncopated[component] =  std::vector<double>(n_q_points);
	  syncopated2[component] =  std::vector<Tensor<1,dim> >(n_q_points);
	  //syncopated2[component] =  std::vector<double>(n_q_points);

	  cell->get_dof_indices (local_dof_indices);
	  //	  int global_index = local_dof_indices[ pinfo[component].alpha_dof_index[0] ];


	  const std::vector<double> &JxW = fe_values.get_JxW_values ();
	  double JxWsum = 0;
	  for (unsigned int q=0; q<n_q_points; ++q){
	    JxWsum += JxW[q];
	  }

	  fe_values[*(alpha[component])].get_function_values(local_solution[component],syncopated[component]);
	  fe_values[*(sigma[component])].get_function_values(local_solution[component],syncopated2[component]);


	  for (unsigned int q=0; q<n_q_points; ++q){
	    for (unsigned int i=0; i<alpha_dofs_per_cell; ++i){
	      cell_avg += (fe_values[*(alpha[component])].value(pinfo[component].alpha_dof_index[i],q) * syncopated[component][q])*JxW[q];
	    }
	  }

	  cell_avg /= JxWsum;

	  for (unsigned int q=0; q<n_q_points; ++q){
	    for (unsigned int i=0; i<alpha_dofs_per_cell; ++i){
	      cell_avg0 += (fe_values[*(alpha[component])].value(pinfo[component].alpha_dof_index[i],q) * (syncopated2[component][q])[0]  )*JxW[q];
	      cell_avg1 += (fe_values[*(alpha[component])].value(pinfo[component].alpha_dof_index[i],q) * (syncopated2[component][q])[1]  )*JxW[q];
	    }
	  }

	  cell_avg0 /= JxWsum;
	  cell_avg1 /= JxWsum;

	  
	  //cell_avg = local_solution[component](global_index);
 
	  //std::cout << "cell_avg 1st = " << cell_avg1 <<std::endl;
	  //std::cout << "cell_avg1 1st = " << cell_avg1 <<std::endl;

	  for(unsigned int vi = 0; vi < GeometryInfo<dim>::vertices_per_cell; ++vi){

	    int v_dex = cell->vertex_index(vi);
	    // int procer = this_mpi_process;

	    if ( cell_avg >=  newmax[v_dex] ){
	      newmax[v_dex] = cell_avg;}
	    
	    if ( cell_avg <  newmin[v_dex] ) {
	      newmin[v_dex] = cell_avg;}

	    if ( cell_avg0 >=  newmax0[v_dex] ){
	      newmax0[v_dex] = cell_avg0;}
	    
	    if ( cell_avg0 <  newmin0[v_dex] ) {
	      newmin0[v_dex] = cell_avg0;}

	    if ( cell_avg1 >=  newmax1[v_dex] ){
	      newmax1[v_dex] = cell_avg1;}
	    
	    if ( cell_avg1 <  newmin1[v_dex] ) {
	      newmin1[v_dex] = cell_avg1;}

	  }
	}
      else if ( cell->is_ghost() )
	{
	  hp_fe_values[component]->reinit (cell);
	  const FEValues<dim> &fe_values = hp_fe_values[component]->get_present_fe_values ();
	  //const FiniteElement<dim>& fe = cell->get_fe();
	  //const unsigned int   dofs_per_cell = fe.dofs_per_cell;

          const Quadrature<dim>& quadrature_formula = fe_values.get_quadrature();
          const FiniteElement<dim>& fe = cell->get_fe();
          const unsigned int   dofs_per_cell = fe.dofs_per_cell;
          const unsigned int   n_q_points    = quadrature_formula.size();

	  
	  const std::string fe_name = fe.get_name();
	  pinfo[component] = pInfoFind(fe.get_name(),component);
	  
	  std::vector<unsigned int> local_dof_indices (dofs_per_cell); 
	  const unsigned int alpha_dofs_per_cell = pinfo[component].alpha_dofs_per_cell;
	  const unsigned int sigma_dofs_per_cell = pinfo[component].sigma_dofs_per_cell;
	  
          //cell->get_dof_indices (local_dof_indices);
          //int global_index = local_dof_indices[ pinfo[component].alpha_dof_index[0] ];

	  syncopated[component] =  std::vector<double>(n_q_points);
	  syncopated2[component] =  std::vector<Tensor<1,dim> >(n_q_points);
          //const std::vector<double> &JxW = fe_values.get_JxW_values (); 

	  fe_values[*(alpha[component])].get_function_values(local_solution[component],syncopated[component]);
	  fe_values[*(sigma[component])].get_function_values(local_solution[component],syncopated2[component]);

	  const std::vector<double> &JxW = fe_values.get_JxW_values ();
          double JxWsum = 0;
          for (unsigned int q=0; q<n_q_points; ++q){
	    JxWsum += JxW[q];
          }

	  double cell_avg = 0;

	  double cell_avg0 = 0;
	  double cell_avg1 = 0;

	  for (unsigned int q=0; q<n_q_points; ++q){
            for (unsigned int i=0; i<alpha_dofs_per_cell; ++i){
              cell_avg += (fe_values[*(alpha[component])].value(pinfo[component].alpha_dof_index[i],q) * syncopated[component][q])*JxW[q];
            }
          }

	  cell_avg /= JxWsum;

	  for (unsigned int q=0; q<n_q_points; ++q){
	    for (unsigned int i=0; i<alpha_dofs_per_cell; ++i){
	      cell_avg0 += (fe_values[*(alpha[component])].value(pinfo[component].alpha_dof_index[i],q) * (syncopated2[component][q])[0]  )*JxW[q];
	      cell_avg1 += (fe_values[*(alpha[component])].value(pinfo[component].alpha_dof_index[i],q) * (syncopated2[component][q])[1]  )*JxW[q];
	    }
	  }

	  cell_avg0 /= JxWsum;
	  cell_avg1 /= JxWsum;


	  //          cell_avg = local_solution[component](global_index);

	  //std::cout << "cell_avg 2nd = " << cell_avg1 <<std::endl;

	  for(unsigned int vi = 0; vi < GeometryInfo<dim>::vertices_per_cell; ++vi){
	    
	    int v_dex = cell->vertex_index(vi);
	    
	    if ( cell_avg >=  newmax[v_dex] ){
	      newmax[v_dex] = cell_avg;}
	    
	    if ( cell_avg <  newmin[v_dex] ) {
	      newmin[v_dex] = cell_avg;}

	    if ( cell_avg0 >=  newmax0[v_dex] ){
	      newmax0[v_dex] = cell_avg0;}
	    
	    if ( cell_avg0 <  newmin0[v_dex] ) {
	      newmin0[v_dex] = cell_avg0;}

	    if ( cell_avg1 >=  newmax1[v_dex] ){
	      newmax1[v_dex] = cell_avg1;}
	    
	    if ( cell_avg1 <  newmin1[v_dex] ) {
	      newmin1[v_dex] = cell_avg1;}


	  }
	}    
    
    // Now we can update the solution
    cell = dof_handler[component]->begin_active();
    endc = dof_handler[component]->end();
    for(;
	cell!=endc;
	++cell	)
      if (cell->is_locally_owned()  )
	{

	  hp_fe_values[component]->reinit (cell);
	  const FEValues<dim> &fe_values = hp_fe_values[component]->get_present_fe_values ();
	  const Quadrature<dim>& quadrature_formula = fe_values.get_quadrature(); 
	  const FiniteElement<dim>& fe = cell->get_fe();
	  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
	  const unsigned int   n_q_points    = quadrature_formula.size();
	  
	  const std::string fe_name = fe.get_name();
	  pinfo[component] = pInfoFind(fe.get_name(),component);
	  
	  std::vector<unsigned int> local_dof_indices (dofs_per_cell); 
	  const unsigned int alpha_dofs_per_cell = pinfo[component].alpha_dofs_per_cell;
	  const unsigned int sigma_dofs_per_cell = pinfo[component].sigma_dofs_per_cell;
	  
	  cell->get_dof_indices (local_dof_indices);
	  //
          syncopated[component] =  std::vector<double>(n_q_points);
	  syncopated2[component] =  std::vector<Tensor<1,dim> >(n_q_points);


	  fe_values[*(alpha[component])].get_function_values(local_solution[component],syncopated[component]);    
	  fe_values[*(sigma[component])].get_function_values(local_solution[component],syncopated2[component]);

	  //          int global_index = local_dof_indices[ pinfo[component].alpha_dof_index[0] ];

	  const std::vector<double> &JxW = fe_values.get_JxW_values ();
          double JxWsum = 0;
          for (unsigned int q=0; q<n_q_points; ++q){
	    JxWsum += JxW[q];
          }

          double cell_avg = 0;
	  //double cell_avg = 0;

	  double cell_avg0 = 0;
	  double cell_avg1 = 0;

          for (unsigned int q=0; q<n_q_points; ++q){
            for (unsigned int i=0; i<alpha_dofs_per_cell; ++i){
              cell_avg += (fe_values[*(alpha[component])].value(pinfo[component].alpha_dof_index[i],q) * syncopated[component][q])*JxW[q];
            }
          }

	  cell_avg /= JxWsum;

	  for (unsigned int q=0; q<n_q_points; ++q){
	    for (unsigned int i=0; i<alpha_dofs_per_cell; ++i){
	      cell_avg0 += (fe_values[*(alpha[component])].value(pinfo[component].alpha_dof_index[i],q) * (syncopated2[component][q])[0]  )*JxW[q];
	      cell_avg1 += (fe_values[*(alpha[component])].value(pinfo[component].alpha_dof_index[i],q) * (syncopated2[component][q])[1]  )*JxW[q];
	    }
	  }

	  cell_avg0 /= JxWsum;
	  cell_avg1 /= JxWsum;


	  //          cell_avg = local_solution[component](global_index);
	  //std::cout << "cell_avg 3rd = " << cell_avg1 <<std::endl;

	  std::vector<double> reset_value(4);
	  std::vector<double> old_val(4);
	  Vector<double> oldquad(n_q_points);
	  std::vector<double> compare(n_q_points);

	  std::vector<double> reset_value0(4);
	  std::vector<double> old_val0(4);
	  Vector<double> oldquad0(n_q_points);
	  std::vector<double> compare0(n_q_points);

	  std::vector<double> reset_value1(4);
	  std::vector<double> old_val1(4);
	  Vector<double> oldquad1(n_q_points);
	  std::vector<double> compare1(n_q_points);
	  
	  // First loop, to reset value
	  //	  for(unsigned int vi = 0; vi < GeometryInfo<dim>::vertices_per_cell; ++vi){

	    
	  //int v_dex = cell->vertex_index(vi);
	  //Point<dim> v_dex_loc = cell->vertex(vi);

	  //std::cout << std::setprecision (15) << "vdex_loc = " << v_dex_loc <<std::endl;

	  //Vector<double> val(3);

	  //VectorTools::point_value( *(dof_handler[component]),
	  //		      local_solution[component],
	  //		      v_dex_loc,
	  //		      val); 

	    
	  //	  double v_dex_val = val(0);
	  //double maf = newmax[v_dex];
	  //double mif = newmin[v_dex];
	  
	  //oldval[v_dex] = val(0);
	  
	  //  reset_value[vi] = std::max( std::min( v_dex_val, maf ), mif );
	  //}
	  //syncopated[component] =  std::vector<double>(n_q_points);      
	  // const std::vector<double> &JxW = fe_values.get_JxW_values ();                                                                                                                                                               
          //fe_values[*(alpha[component])].get_function_values(local_solution[component],syncopated[component]);                                                                                                                              
	  
	  //double JxWsum = 0;				
	  //for (unsigned int q=0; q<n_q_points; ++q){		
	  //JxWsum += JxW[q];	   
	  //}                                                                                                                                                                                                                                  

          // += (fe_values[*(alpha[component])].value(pinfo[component].alpha_dof_in


	  for(unsigned int vi = 0; vi < GeometryInfo<dim>::vertices_per_cell; ++vi){ 
	    int v_dex = cell->vertex_index(vi);
	    Point<dim> v_dex_loc = cell->vertex(vi); 
	    
	    for (unsigned int q=0; q<n_q_points; ++q){
	      oldquad(q) = 0.0;
	      oldquad0(q)= 0.0;
	      oldquad1(q)= 0.0;
	      for (unsigned int i=0; i<alpha_dofs_per_cell; ++i){		
		//oldquad(q) = syncopated[component][q];
		oldquad(q) += (fe_values[*(alpha[component])].value(pinfo[component].alpha_dof_index[i],q) * syncopated[component][q]); 

		oldquad0(q) += (fe_values[*(alpha[component])].value(pinfo[component].alpha_dof_index[i],q) * (syncopated2[component][q])[0] ); 
		oldquad1(q) += (fe_values[*(alpha[component])].value(pinfo[component].alpha_dof_index[i],q) * (syncopated2[component][q])[1] ); 
		
		//std::cout << "Components of soln = " << syncopated[component][q] << std::endl;
		//std::cout << "Components of shape = " << (fe_values[*(alpha[component])].value(pinfo[component].alpha_dof_index[i],q) ) << std::endl;
		
	      }

	      Point<dim> qup = hp_fe_values[component]->quadrature_point(q);
	      
              if (v_dex_loc == qup){

                reset_value[vi] = oldquad[q];
		old_val[vi] = reset_value[vi];

                reset_value0[vi] = oldquad0[q];
		old_val0[vi] = reset_value0[vi];

                reset_value1[vi] = oldquad1[q];
		old_val1[vi] = reset_value1[vi];


		double maf = newmax[v_dex];
		double mif = newmin[v_dex];

		double maf0 = newmax0[v_dex];
		double mif0 = newmin0[v_dex];

		double maf1 = newmax1[v_dex];
		double mif1 = newmin1[v_dex];

		//old_val[vi] = reset_value[vi];                                                                                                                                                                                                  
		//old_val[vi] = reset_value[vi];                                                                                                                                                                                                  
		reset_value[vi]  = std::max( std::min( reset_value[vi], maf ), mif );

		reset_value0[vi] = std::max( std::min( reset_value0[vi], maf0 ), mif0 );

		reset_value1[vi] = std::max( std::min( reset_value1[vi], maf1 ), mif1 );

		//std::cout << "old_quadrature values? = " << oldquad1[q] << std::endl;
	      }
	    }
          }  
	  
	  // Second loop, to redistribute at vertices
	  double sumloc = 0.0;
	  double sumloc0 = 0.0;
	  double sumloc1 = 0.0;

	  //std::cout << "reset_value 1st = " << reset_value1 << std::endl;

	  for(unsigned int vl = 0; vl < GeometryInfo<dim>::vertices_per_cell; ++vl){

	    sumloc  = ( reset_value[0] +reset_value[1] +reset_value[2] +reset_value[3]  ) / 4.0;

	    sumloc0 = ( reset_value0[0]+reset_value0[1]+reset_value0[2]+reset_value0[3] ) / 4.0;
	    sumloc1 = ( reset_value1[0]+reset_value1[1]+reset_value1[2]+reset_value1[3] ) / 4.0;

	    int v_dex;
	    int vsym;

	    //std::cout << "sumloc = " << sumloc << ", cell_avg" << std::endl;
	    /*     if (vl == 1){ */
	    /* 	      vsym = 2; */
	    /* 	      v_dex = cell->vertex_index(vsym); */
	    /* 	    } */
	    /* 	    else if (vl == 2){ */
	    /*               vsym = 3; */
	    /*               v_dex = cell->vertex_index(vsym); */
	    /* 	    } */
	    /* 	    else if (vl == 3){ */
	    /*               vsym = 1; */
	    /*               v_dex = cell->vertex_index(vsym); */
	    /* 	    } */
	    /* 	    else { */
	    /*               v_dex = cell->vertex_index(vl); */

	    /* 	    } */

	    v_dex = cell->vertex_index(vl);

	    double sumdif = ( sumloc - cell_avg ) * 4.0;
	    double sign_dif = 1.0;

	    double sumdif0 = ( sumloc0 - cell_avg0 ) * 4.0;
	    double sign_dif0 = 1.0;

	    double sumdif1 = ( sumloc1 - cell_avg1 ) * 4.0;
	    double sign_dif1 = 1.0;
	    
	    //std::cout << "sumdif = " << sumdif << std::endl;
	    
	    //double dif1 = 1.0;
	    //double dif2 = 1.0;
	    //double dif3 = 1.0;
	    //double dif4 = 1.0;

	    if (sumdif >= 0.0) sign_dif = 1.0;
	    else sign_dif = -1.0;

	    if (sumdif0 >= 0.0) sign_dif0 = 1.0;
	    else sign_dif0 = -1.0;

	    if (sumdif1 >= 0.0) sign_dif1 = 1.0;
	    else sign_dif1 = -1.0;

	  
	    std::vector<double> dif(4,0);
	    std::vector<double> dif0(4,0);
	    std::vector<double> dif1(4,0);
	  
	    dif[0] = (reset_value[0]-cell_avg)*sign_dif;
	    dif[1] = (reset_value[1]-cell_avg)*sign_dif;
	    dif[2] = (reset_value[2]-cell_avg)*sign_dif;
	    dif[3] = (reset_value[3]-cell_avg)*sign_dif;

	    dif0[0] = (reset_value0[0]-cell_avg0)*sign_dif0;
	    dif0[1] = (reset_value0[1]-cell_avg0)*sign_dif0;
	    dif0[2] = (reset_value0[2]-cell_avg0)*sign_dif0;
	    dif0[3] = (reset_value0[3]-cell_avg0)*sign_dif0;

	    dif1[0] = (reset_value1[0]-cell_avg1)*sign_dif1;
	    dif1[1] = (reset_value1[1]-cell_avg1)*sign_dif1;
	    dif1[2] = (reset_value1[2]-cell_avg1)*sign_dif1;
	    dif1[3] = (reset_value1[3]-cell_avg1)*sign_dif1;

	    //std::cout << "dif = " << dif << std::endl;
	  
	    double inc1 = 0.0;
	    if (dif[0]>0.0) inc1 = 1.0;
	    double inc2 = 0.0;
	    if (dif[1]>0.0) inc2 = 1.0;
	    double inc3 = 0.0; 
	    if (dif[2]>0.0) inc3 = 1.0;
	    double inc4 = 0.0; 
	    if (dif[3]>0.0) inc4 = 1.0;
	    double kdp = inc1+inc2+inc3+inc4;

	    double inc10 = 0.0;
	    if (dif0[0]>0.0) inc10 = 1.0;
	    double inc20 = 0.0;
	    if (dif0[1]>0.0) inc20 = 1.0;
	    double inc30 = 0.0; 
	    if (dif0[2]>0.0) inc30 = 1.0;
	    double inc40 = 0.0; 
	    if (dif0[3]>0.0) inc40 = 1.0;
	    double kdp0 = inc10+inc20+inc30+inc40;

	    double inc11 = 0.0;
	    if (dif1[0]>0.0) inc11 = 1.0;
	    double inc21 = 0.0;
	    if (dif1[1]>0.0) inc21 = 1.0;
	    double inc31 = 0.0; 
	    if (dif1[2]>0.0) inc31 = 1.0;
	    double inc41 = 0.0; 
	    if (dif1[3]>0.0) inc41 = 1.0;
	    double kdp1 = inc11+inc21+inc31+inc41;
	  
	    // Need to redistribute now
	    for(unsigned int vj = 0; vj < GeometryInfo<dim>::vertices_per_cell; ++vj){
	      
	      //int v_dex = cell->vertex_index(vj);
	      
	      double div  = std::max(1.0,kdp);
	      double div0 = std::max(1.0,kdp0);
	      double div1 = std::max(1.0,kdp1);
	      
	      //std::cout << "div = " << div << std::endl;
	      double redfac;
	      double redmax;

	      double redfac0;
	      double redmax0;
	      double redfac1;
	      double redmax1;
	      
	      if(dif[vj] > 0.0){
		redfac = sumdif*sign_dif/div;
		kdp = kdp - 1.0;}
	      else{
		redfac = 0.0;}
	      if(sign_dif > 0.0){
		redmax = reset_value[vj] - newmin[v_dex];}
	      else{
		redmax = newmax[v_dex] - reset_value[vj];}
	      redfac = std::min(redfac,redmax);
	      sumdif = sumdif - redfac*sign_dif;
	      reset_value[vj] = reset_value[vj] - redfac*sign_dif;

	      if(dif0[vj] > 0.0){
		redfac0 = sumdif0*sign_dif0/div0;
		kdp0 = kdp0 - 1.0;}
	      else{
		redfac0 = 0.0;}
	      if(sign_dif0 > 0.0){
		redmax0 = reset_value0[vj] - newmin0[v_dex];}
	      else{
		redmax0 = newmax0[v_dex] - reset_value0[vj];}
	      redfac0 = std::min(redfac0,redmax0);
	      sumdif0 = sumdif0 - redfac0*sign_dif0;
	      reset_value0[vj] = reset_value0[vj] - redfac0*sign_dif0;

	      if(dif1[vj] > 0.0){
		redfac1 = sumdif1*sign_dif1/div1;
		kdp1 = kdp1 - 1.0;}
	      else{
		redfac1 = 0.0;}
	      if(sign_dif1 > 0.0){
		redmax1 = reset_value1[vj] - newmin1[v_dex];}
	      else{
		redmax1 = newmax1[v_dex] - reset_value1[vj];}
	      redfac1 = std::min(redfac1,redmax1);
	      sumdif1 = sumdif1 - redfac1*sign_dif1;
	      reset_value1[vj] = reset_value1[vj] - redfac1*sign_dif1;
	      
	    }
	    //std::cout << "reset_value_final = " << reset_value << std::endl;
	  }
	  
	  // We need to restict to the linear part, and interpolate
	  // for (unsigned int i =0; i<alpha_dofs_per_cell; ++i){
	  // int gi0 = local_dof_indices[ pinfo[component].alpha_dofs_index[0] ];
	  
	  //std::vector<double> navel(2);

	  //Point<dim> navel = cell->center();
	  //Point<dim> ref_coord = cell->vertex(0);
	  
	  /* 	  std::cout << "point = " << ref_coord << std::endl; */
	  /* 	  std::cout << "pointx = " << ref_coord(0) << std::endl; */
	  /* 	  std::cout << "pointy = " << ref_coord(1) << std::endl; */

	  /* 	  std::cout << "npoint = " << navel << std::endl; */
	  /* 	  std::cout << "npointx = " << navel(0) << std::endl; */
	  /* 	  std::cout << "npointy = " << navel(1) << std::endl; */	  
	  
	  for (unsigned int i=0; i<alpha_dofs_per_cell; ++i){
	    
	    //	    int gi0 = local_dof_indices[ pinfo[component].alpha_dof_index[0] ];
	    //int gi1 = local_dof_indices[ pinfo[component].alpha_dof_index[2] ];
	    //int gi2 = local_dof_indices[ pinfo[component].alpha_dof_index[6] ];
	    //int gi3 = local_dof_indices[ pinfo[component].alpha_dof_index[8] ];

	    int gi = local_dof_indices[ pinfo[component].alpha_dof_index[i] ];

	    //	    double deltx =  ( reset_value[1] + reset_value[3] - reset_value[2] - reset_value[0] ) /  (  cell->extent_in_direction(0)*2.0  );
	    //double delty =  ( reset_value[2] + reset_value[3] - reset_value[0] - reset_value[1] ) /  (  cell->extent_in_direction(1)*2.0  );

	    /* 	  double sec_termx =  deltx*ref_coord(0)/navel(0); */
	    /* 	  double sec_termy =  delty*ref_coord(1)/navel(1); */

	    /* 	  double t_termx = ( reset_value[2] + reset_value[0] )/ (2.0*navel(0)); */
	    /* 	  double t_termy = ( reset_value[0] + reset_value[1] )/ (2.0*navel(1)); */
 	    
	    //	  double weightx =  subdomain_solution[component]( gi1 ) / deltx;
	    //  double weighty =  subdomain_solution[component]( gi2 ) / delty;

	    //          if (step == 0 && reset_value == old_val){

	    //    weightx =  subdomain_solution[component]( gi1 ) / deltx; 
	    //weighty =  subdomain_solution[component]( gi2 ) / delty;
	    //continue;;

	    //}
	    //else if (step >= 1){
	  
	    //redistributed_solution ( gi0 ) = reset_value[0];
	    //redistributed_solution ( gi1 ) = reset_value[1];  // - sec_termx + t_termx;
	    //redistributed_solution ( gi2 ) = reset_value[2];
	    //redistributed_solution ( gi3 ) = reset_value[3]; //- sec_termy + t_termy;
	    
	    //redistributed_solution (gi1 ) = ( navel(1)-ref_coord(1) + deltx*
	    
	    //std::cout << "delta x = " << cell->extent_in_direction(0) << ", deltay = " << cell->extent_in_direction(0) << std::endl;
	    
	    //	    if ( std::abs(subdomain_solution[component]( gi1 )) < std::abs(redistributed_solution (gi1 )) ){ //&& old_val == reset_value ){
	    //std::cout << std::setprecision (15) << "measure = " << cell->measure() << std::endl;

	    //std::cout << std::setprecision (15) << "prev_slope0 = " <<   subdomain_solution[component]( gi0 ) << ", new_slope0 = " <<   redistributed_solution (gi0 ) << std::endl;
	    //	    std::cout << std::setprecision (15) << "prev_slope0 = " <<   subdomain_solution[component]( gi0 ) << ", new_slope0 = " <<   reset_value[0] << std::endl;
	    //std::cout << std::setprecision (15) << "prev_slope1 = " <<   subdomain_solution[component]( gi1 ) << ", new_slope1 = " <<   reset_value[1] << std::endl;
	    //std::cout << std::setprecision (15) << "prev_slope2 = " <<   subdomain_solution[component]( gi2 ) << ", new_slope2 = " <<   reset_value[2] << std::endl;
	    //std::cout << std::setprecision (15) << "prev_slope3 = " <<   subdomain_solution[component]( gi3 ) << ", new_slope3 = " <<   reset_value[3] << std::endl;
	    //std::cout << std::setprecision (15) << "reset_val = " << reset_value[0] << ", " << reset_value[1] << ", " << reset_value[2] << ", " << reset_value[3] << std::endl;
	    //std::cout << std::setprecision (15) << "oldval = " << old_val[0] << ", " << old_val[1] << ", " << old_val[2] << ", " << old_val[3] << std::endl;

	    //std::cout << std::setprecision (15) << "oldquad = " << oldquad << ", x_extent = " << cell->extent_in_direction(0) << std::endl;
	    //std::cout << std::setprecision (15) << "compare= " << compare << ", x_extent = " << cell->extent_in_direction(0) << std::endl;
	      
	    //}
	    
	    // Now we need to reset the solution in the basis
	    
	    //	    for (unsigned int i=0; i<alpha_dofs_per_cell; ++i){
	    //int global_index = local_dof_indices[ pinfo[component].alpha_dof_index[i] ];
	    //if ( i == 1 ){

 	    double bound = 5.0; 
	    double dx = cell->extent_in_direction(0); 
	    double dy = cell->extent_in_direction(1); 

	    //std::cout << "slope1 = " << std::abs(old_val[0]-old_val[1])/dx << std::endl;
	    //std::cout << "slope2 = " << std::abs(old_val[3]-old_val[2])/dx << std::endl;
	    //std::cout << "slope3 = " << std::abs(old_val[0]-old_val[3])/dy << std::endl;
	    //std::cout << "slope4 = " << std::abs(old_val[1]-old_val[2])/dy << std::endl;



	    if ( (std::abs(old_val[0]-old_val[1])/dx >= bound) || (std::abs(old_val[3]-old_val[2])/dx >= bound) || (std::abs(old_val[0]-old_val[3])/dy >= bound) || (std::abs(old_val[1]-old_val[2])/dy >= bound) ) { 
	      
	      /* 	  double mild = -99.0; */
	      /* 	  for(unsigned int vl = 0; vl < GeometryInfo<dim>::vertices_per_cell; ++vl){ */
	      /* 	  for(unsigned int vj = 0; vj < GeometryInfo<dim>::vertices_per_cell; ++vj){ */
	      /* 	    mild = std::max(std::abs(old_val[vl]-old_val[vj]),mild  ); */
	      /* 	  } */
	      /* 	  } */
	      
	      /* 	  std::cout << "max = " << mild/cell->extent_in_direction(0) << std::endl; */
	      
	      /* 	  std::cout << "extent = " << cell->extent_in_direction(0) << std::endl; */
	      
	      
	      /* 	    	      std::cout << "tell me, component = " << component << std::endl; */
	      /* 		    } */
	      

	      //std::cout << "slope1 = " << std::abs(old_val[0]-old_val[1])/dx << std::endl;                                                                                                                                                    
	      //std::cout << "slope2 = " << std::abs(old_val[3]-old_val[2])/dx << std::endl;                                                                                                                                                    
	      //std::cout << "slope3 = " << std::abs(old_val[0]-old_val[3])/dy << std::endl;                                                                                                                                                    
	      //std::cout << "slope4 = " << std::abs(old_val[1]-old_val[2])/dy << std::endl; 
	      
	      //For higher than linears need to arrange dofs correctly
	      if(i == 0){ 
		//subdomain_solution[component]( gi ) =  reset_value[0]; 
	      } 
	      else if(i == 1){ 
		//subdomain_solution[component]( gi ) =  reset_value[1]; 
	      } 
	      else if(i == 2){ 
		//subdomain_solution[component]( gi ) =  reset_value[2]; 
	      } 
	      else if(i == 3){ 
	      	//subdomain_solution[component]( gi ) =  reset_value[3]; 
	      } 
	      else{ 
	      	//subdomain_solution[component]( gi ) =  0.0; 
	      } 

	    }
	  }
	  
	  for (unsigned int i=0; i<sigma_dofs_per_cell; ++i){
	    int gi0 = local_dof_indices[  pinfo[component].sigma_dof_index[0][i] ];
	    int gi1 = local_dof_indices[  pinfo[component].sigma_dof_index[1][i] ];

	    //std::cout << std::setprecision (15) << "prev_slope1 = " <<   subdomain_solution[component](  local_dof_indices[  pinfo[component].sigma_dof_index[0][0] ] ) << ", i= " <<  i << std::endl; 
	    //std::cout << std::setprecision (15) << "prev_slope2 = " <<   subdomain_solution[component](  local_dof_indices[  pinfo[component].sigma_dof_index[0][1] ] ) << ", i= " <<  i << std::endl; 
	    //std::cout << std::setprecision (15) << "prev_slope3 = " <<   subdomain_solution[component](  local_dof_indices[  pinfo[component].sigma_dof_index[0][max_degree+1] ]) << ", i= " <<   i << std::endl; 
	    //std::cout << std::setprecision (15) << "prev_slope4 = " <<   subdomain_solution[component](  local_dof_indices[  pinfo[component].sigma_dof_index[0][max_degree+2] ]) << ", i= " <<   i << std::endl; 
	    //std::cout << std::setprecision (15) << "reset_val = " << reset_value0[0] << ", " << reset_value0[1] << ", " << reset_value0[2] << ", " << reset_value0[3] << std::endl; 
	    //std::cout << std::setprecision (15) << "oldval = " << old_val0[0] << ", " << old_val0[1] << ", " << old_val0[2] << ", " << old_val0[3] << std::endl; 
	    //std::cout << "gi0 = " << gi0 << " ,gi1 = " << gi1 << std:endl;

	    //int gi00 = local_dof_indices[  pinfo[component].sigma_dof_index[1][0] ]; 
 	    //int gi10 = local_dof_indices[  pinfo[component].sigma_dof_index[1][1] ]; 
 	    //int gi20 = local_dof_indices[  pinfo[component].sigma_dof_index[1][2] ]; 
	    //int gi30 = local_dof_indices[  pinfo[component].sigma_dof_index[1][3] ]; 

	    // 	    std::cout << std::setprecision (15) << "prev_slope0 = " <<   subdomain_solution[component]( gi00 ) << ", new_slope0 = " <<   reset_value1[0] << std::endl; 
 	    //std::cout << std::setprecision (15) << "prev_slope1 = " <<   subdomain_solution[component]( gi10 ) << ", new_slope1 = " <<   reset_value1[1] << std::endl; 
 	    //std::cout << std::setprecision (15) << "prev_slope2 = " <<   subdomain_solution[component]( gi20 ) << ", new_slope2 = " <<   reset_value1[2] << std::endl; 
 	    //std::cout << std::setprecision (15) << "prev_slope3 = " <<   subdomain_solution[component]( gi30 ) << ", new_slope3 = " <<   reset_value1[3] << std::endl; 
 	    //std::cout << std::setprecision (15) << "reset_val = " << reset_value1[0] << ", " << reset_value1[1] << ", " << reset_value1[2] << ", " << reset_value1[3] << std::endl; 
 	    //std::cout << std::setprecision (15) << "oldval = " << old_val1[0] << ", " << old_val1[1] << ", " << old_val1[2] << ", " << old_val1[3] << std::endl; 
	    //std::cout << "i = " << i << " ,gi00 = " << gi00 << " ,gi10 = " << gi10 << std::endl;

	    if( i == 0 ){
	      subdomain_solution[component]( gi0 ) =  reset_value0[0];
	      subdomain_solution[component]( gi1 ) =  reset_value1[0];
	    }
	    else if( i == 1 ){
	      subdomain_solution[component]( gi0 ) =  reset_value0[1];
	      subdomain_solution[component]( gi1 ) =  reset_value1[1];
	    }
	    else if( i == 2 ){
	      subdomain_solution[component]( gi0 ) =  reset_value0[2];
	      subdomain_solution[component]( gi1 ) =  reset_value1[2];
	    }
	    else if( i == 3 ){
	      subdomain_solution[component]( gi0 ) =  reset_value0[3];
	      subdomain_solution[component]( gi1 ) =  reset_value1[3];
	    }
	    else{
	      subdomain_solution[component]( gi0 ) =  0.0;
	      subdomain_solution[component]( gi1 ) =  0.0;
	    }

	    
/* 	    if(gi1 == 0){ */
/* 	      subdomain_solution[component]( gi1 ) =  reset_value1[0]; */
/* 	    } */
	    
/* 	    else if(gi1 == 1){ */
/* 	      subdomain_solution[component]( gi1 ) =  reset_value1[1]; */
/* 	    } */
	    
/* 	    else if(gi1 == max_degree+1){ */
/* 	      subdomain_solution[component]( gi1 ) =  reset_value1[2]; */
/* 	    } */
/* 	    else if(gi1 == max_degree+2){ */
/* 	      subdomain_solution[component]( gi1 ) =  reset_value1[3]; */
/* 	    } */
/* 	    else{ */
/* 	      subdomain_solution[component]( gi1 ) =  0.0; */
/* 	    } */
	    
	  }

	  //}
	  //	      else if ( i == max_degree+1){
	  //	subdomain_solution[component]( global_index ) = redistributed_solution( global_index );
	  //}
	  //else if ( i == max_degree+2){
	  //subdomain_solution[component]( global_index ) = redistributed_solution( global_index );
	  //}

	  //else if (i != 0){
	  //subdomain_solution[component]( global_index ) = 0.0;
	  //}
	  //	    }
	    
	  //}
	}
    //redistributed_solution.compress();
  }

  for (unsigned int component=0; component< alphadim; ++component){  
    delete hp_fe_values[component];
    //    delete lhp_fe_values[component];
    //delete hp_fe_values_face[component];
    //delete hp_fe_values_neigh_face[component];
  }
}














