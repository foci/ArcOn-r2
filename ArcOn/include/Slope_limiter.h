/* This is the BDS slopelimiter ... this needs to be fixed (easy) */
template <int dim>
void arcOn<dim>::slopelimiter( SolutionVector& subdomain_solution, int step )
{

  std::vector< FEValues<dim>* > hp_fe_values;
  for (unsigned int component=0; component< alphadim; ++component){

    UpdateFlags updateflags=  update_values | update_gradients | update_quadrature_points | update_JxW_values;
    
    hp_fe_values.push_back (new FEValues<dim>(*(fe_collection[component]), *(quadrature_collection[component]), updateflags));
    local_solution[component] = subdomain_solution[component];

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

	    newmax[v_dex] = -9999.0; 
	    newmin[v_dex] =  9999.0;

	    newmax0[v_dex] = -9999.0; 
	    newmin0[v_dex] =  9999.0;

	    newmax1[v_dex] = -9999.0; 
	    newmin1[v_dex] =  9999.0;
	  }
	}

    //First find the stencil extrema

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

	  cell->get_dof_indices (local_dof_indices);

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
      else if ( cell->is_ghost() )
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
	  
	  syncopated[component] =  std::vector<double>(n_q_points);
	  syncopated2[component] =  std::vector<Tensor<1,dim> >(n_q_points);

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
          syncopated[component] =  std::vector<double>(n_q_points);
	  syncopated2[component] =  std::vector<Tensor<1,dim> >(n_q_points);


	  fe_values[*(alpha[component])].get_function_values(local_solution[component],syncopated[component]);    
	  fe_values[*(sigma[component])].get_function_values(local_solution[component],syncopated2[component]);

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
	  
	  for(unsigned int vi = 0; vi < GeometryInfo<dim>::vertices_per_cell; ++vi){ 
	    int v_dex = cell->vertex_index(vi);
	    Point<dim> v_dex_loc = cell->vertex(vi); 
	    
	    for (unsigned int q=0; q<n_q_points; ++q){
	      oldquad(q) = 0.0;
	      oldquad0(q)= 0.0;
	      oldquad1(q)= 0.0;
	      for (unsigned int i=0; i<alpha_dofs_per_cell; ++i){		

		oldquad(q) += (fe_values[*(alpha[component])].value(pinfo[component].alpha_dof_index[i],q) * syncopated[component][q]); 

		oldquad0(q) += (fe_values[*(alpha[component])].value(pinfo[component].alpha_dof_index[i],q) * (syncopated2[component][q])[0] ); 
		oldquad1(q) += (fe_values[*(alpha[component])].value(pinfo[component].alpha_dof_index[i],q) * (syncopated2[component][q])[1] ); 
		
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

		reset_value[vi]  = std::max( std::min( reset_value[vi], maf ), mif );

		reset_value0[vi] = std::max( std::min( reset_value0[vi], maf0 ), mif0 );

		reset_value1[vi] = std::max( std::min( reset_value1[vi], maf1 ), mif1 );

	      }
	    }
          }  
	  
	  // Second loop, to redistribute at vertices
	  double sumloc = 0.0;
	  double sumloc0 = 0.0;
	  double sumloc1 = 0.0;

	  for(unsigned int vl = 0; vl < GeometryInfo<dim>::vertices_per_cell; ++vl){

	    sumloc  = ( reset_value[0] +reset_value[1] +reset_value[2] +reset_value[3]  ) / 4.0;

	    sumloc0 = ( reset_value0[0]+reset_value0[1]+reset_value0[2]+reset_value0[3] ) / 4.0;
	    sumloc1 = ( reset_value1[0]+reset_value1[1]+reset_value1[2]+reset_value1[3] ) / 4.0;

	    int v_dex;
	    int vsym;

	    v_dex = cell->vertex_index(vl);

	    double sumdif = ( sumloc - cell_avg ) * 4.0;
	    double sign_dif = 1.0;

	    double sumdif0 = ( sumloc0 - cell_avg0 ) * 4.0;
	    double sign_dif0 = 1.0;

	    double sumdif1 = ( sumloc1 - cell_avg1 ) * 4.0;
	    double sign_dif1 = 1.0;
	    
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
	      
	      double div  = std::max(1.0,kdp);
	      double div0 = std::max(1.0,kdp0);
	      double div1 = std::max(1.0,kdp1);

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
	  }
	  
	  for (unsigned int i=0; i<alpha_dofs_per_cell; ++i){
	    

	    int gi = local_dof_indices[ pinfo[component].alpha_dof_index[i] ];

 	    double bound = 5.0; 
	    double dx = cell->extent_in_direction(0); 
	    double dy = cell->extent_in_direction(1); 

	    if ( (std::abs(old_val[0]-old_val[1])/dx >= bound) || (std::abs(old_val[3]-old_val[2])/dx >= bound) || (std::abs(old_val[0]-old_val[3])/dy >= bound) || (std::abs(old_val[1]-old_val[2])/dy >= bound) ) { 
	      
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

	  }

	}
  }

  for (unsigned int component=0; component< alphadim; ++component){  
    delete hp_fe_values[component];
  }
}














