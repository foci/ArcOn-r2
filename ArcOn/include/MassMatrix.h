/* This compute the mass matrix for basis functions for so-called "primitive" basis functions, i.e. non-coupled basis functions.  This is deprecated at present, see the MatrixMapper template. */
template<int dim>
void arcOn<dim>::initialize_massmatrix()
{

  //std::vector<pInfo>  pinfo(alphadim);

  for(unsigned int index=0; index<alphadim; index++){

    pmap.push_back( pMap() );

    UpdateFlags updateflags =  update_values | update_gradients | update_quadrature_points | update_JxW_values;
    std::cout << "***Initializing list of mass matrices***" << std::endl;
    std::cout << "Mass matrix for component[" << index << "]" << std::endl;
    Triangulation<dim> fake_triangulation;
    GridGenerator::hyper_cube(fake_triangulation);
    typename Triangulation<dim>::active_cell_iterator fake_cell = fake_triangulation.begin_active(); 
    std::cout << "                     Collection Index = " << index << std::endl;
    const FiniteElement<dim>& fe = *(fe_collection[index]);
    const Quadrature<dim>& quadrature_formula = *(quadrature_collection[index]); 
    FEValues<dim> fe_values(fe,quadrature_formula,updateflags);
    fe_values.reinit( fake_cell );
    const std::vector<double>& quad_weights = quadrature_formula.get_weights();
    const unsigned int n_q_points = quadrature_formula.size();
    const std::string fe_name = fe.get_name();
    bool itexists = pInfoExists(fe_name,index);
    unsigned int n_comp = fe.n_components();
    pcout << "Total State Vector components = " << n_comp <<std::endl;

    if(!itexists){
      //pInfo  pinfo;
      //newly created -- need to calculate
      pcout << std::endl;
      pcout << "Found new FESystem: " << fe_name << std::endl << std::endl;
      const unsigned int n_base_elements = fe.n_base_elements();
      for(unsigned int i=0;i<n_base_elements;i++){
	pcout << "Base element " << i << " which is " << fe.base_element(i).get_name() << " has multiplicity " << fe.element_multiplicity(i) << std::endl;
      }
      std::cout << "System to finite element (state vector) component index " << std::endl;
      // component sizes is the "total length" of the state vector
      unsigned int component_sizes[dim+1];
      for(unsigned int l=0; l<dim+1; l++) component_sizes[l]=0;
      std::vector<int> roll(n_comp);
      std::fill(roll.begin(), roll.end(), 0);
      for(unsigned int i=0; i<fe.dofs_per_cell; ++i){
	//component_to_system_index returns a pair, the first entry is component, the second is the component index
	std::pair<int,int> component_i = fe.system_to_component_index(i);
	
	//New we can store the last index number in the component sizes for later use (here we assume contiguity)
	if( component_sizes[ component_i.first ] < (unsigned int)component_i.second ) component_sizes[component_i.first] = component_i.second;

	//Alternatively, if we want to make no contiguity assumptions, we can just roll over the count

	for(unsigned int l=0; l<n_comp; l++){
	  if (component_i.first == l) 
	    roll[l] = roll[l] + 1;
	}
	
	pcout << i << "-->(" << component_i.first << "," << component_i.second << ") ";
	pcout << std::endl;
      }
      for(unsigned int l=0; l<n_comp; l++){
	//pcout << "roll = " << roll[l] << std::endl;
	//pcout << std::endl;
      }
      //Now component_sizes has the last for each component
      pinfo[index].alpha_dofs_per_cell = roll[0];
      pinfo[index].sigma_dofs_per_cell = roll[1];
      /* pinfo[index].alpha_dofs_per_cell = component_sizes[0]+1; */
      /* pinfo[index].sigma_dofs_per_cell = component_sizes[1]+1; */

      pinfo[index].alpha_dof_index = std::vector<int>( pinfo[index].alpha_dofs_per_cell );
      for(unsigned int l=0;l<dim;l++){ pinfo[index].sigma_dof_index[l] = std::vector<int>( pinfo[index].sigma_dofs_per_cell ); }

      pcout << "component to system index " << std::endl;
      for(unsigned int component=0; component<dim+1; component++){
	for(unsigned int index_i=0; index_i < roll[component] ; index_i++){
	  int global_i = fe.component_to_system_index(component,index_i);
	  if (component==0){
	    pinfo[index].alpha_dof_index[ index_i ] = global_i;
	  } else {
	    pinfo[index].sigma_dof_index[ component-1 ][ index_i ]  = global_i;
	  }
	  pcout << "(" << component << ", " << index_i << ")-->" << global_i << " ";  
	}
      }
      pcout << std::endl;
      pcout << "dofs_per_cell = " << fe.dofs_per_cell << " , alpha_dofs_per_cell= " << pinfo[index].alpha_dofs_per_cell << " sigma_dofs_per_cell= " << pinfo[index].sigma_dofs_per_cell << std::endl;
      std::cout << "DOFs for alpha: " << pinfo[index].alpha_dof_index << std::endl;
      for(unsigned int l=0;l<dim;l++){ std::cout << "DOFs for sigma_" << l << " " << pinfo[index].sigma_dof_index[l] << std::endl; }

      /*
	std::cout << "DOF has support on face " << std::endl;
	for(unsigned int face_num=0; face_num < GeometryInfo<dim>::faces_per_cell; ++face_num){ std::cout << "\t" << face_num; }
	std::cout << std::endl; 
	for(unsigned int i=0;i<pinfo.alpha_dofs_per_cell;i++){
	std::cout << "[" << pinfo.alpha_dof_index[i] << "]\t"; 
	for(unsigned int face_num=0; face_num < GeometryInfo<dim>::faces_per_cell; ++face_num){
	if( fe.has_support_on_face( pinfo.alpha_dof_index[i] , face_num ) ){
	std::cout << "X\t"; 
	} else {
	std::cout << "-\t";
	}
	}
	std::cout << std::endl;
	}
	std::cout << std::endl;
      */

      FullMatrix<double> MassMatrix(pinfo[index].alpha_dofs_per_cell,pinfo[index].alpha_dofs_per_cell);
      MassMatrix = 0;
      std::cout << "Alpha Mass Matrix numerical integral over " << n_q_points << " quad points of phi_i * phi_j " << pinfo[index].alpha_dofs_per_cell << "-by-" << pinfo[index].alpha_dofs_per_cell << std::endl; 

      for (unsigned int i=0; i<pinfo[index].alpha_dofs_per_cell; ++i){
	//std::cout << "first i = "  << i << std::endl;
	for (unsigned int j=0; j<pinfo[index].alpha_dofs_per_cell; ++j){
	  //std::cout << "first j = " << j << std::endl;
	  for (unsigned int q=0; q<n_q_points; ++q){
	    //std::cout << "fe_values = " << fe_values[*(alpha[0])].value(pinfo.alpha_dof_index[i],q) * fe_values[*(alpha[0])].value(pinfo.alpha_dof_index[j],q) * quad_weights[q] << ", i = " << i << ", q =" << q << ", j = " << j << std::endl;
	    MassMatrix(i,j)+= fe_values[*(alpha[index])].value(pinfo[index].alpha_dof_index[i],q) * fe_values[*(alpha[index])].value(pinfo[index].alpha_dof_index[j],q) * quad_weights[q];
	  }
	}
      }
      pcout << "Alpha Mass Matrix" << std::endl;
      MassMatrix.print_formatted(std::cout);
      pinfo[index].AlphaMassMatrix = MassMatrix;

      pinfo[index].InverseAlphaMassMatrix = FullMatrix<double>(pinfo[index].alpha_dofs_per_cell,pinfo[index].alpha_dofs_per_cell);
      pinfo[index].InverseAlphaMassMatrix.invert(pinfo[index].AlphaMassMatrix);
      MassMatrix.diagadd(-1);
      double TestIdentity_norm = MassMatrix.frobenius_norm();
      pcout << "alpha norm = " << TestIdentity_norm << std::endl;
      if( TestIdentity_norm < MassMatrix.m()*1e-14 ){
	pinfo[index].alphaProject=false;
      } else {
	pinfo[index].alphaProject=true;
      }
      pcout << "Alpha Inverse Mass Matrix";
      if(!pinfo[index].alphaProject) std::cout << " is close to the Identity. ";
      std::cout << std::endl;
      std::cout << "Alpha Mass Matrix" << std::endl;
      pinfo[index].AlphaMassMatrix.print_formatted(std::cout);
      std::cout << "Inverse Alpha Mass Matrix" << std::endl;
      pinfo[index].InverseAlphaMassMatrix.print_formatted(std::cout);
      std::cout << "Alpha Mass Matrix minus identity" << std::endl;
      MassMatrix.print_formatted(std::cout);



      MassMatrix = FullMatrix<double>(pinfo[index].sigma_dofs_per_cell,pinfo[index].sigma_dofs_per_cell);
      MassMatrix = 0;
      std::cout << "Sigma Mass Matrix numerical integral over " << n_q_points << " quad points of phi_i * phi_j " << pinfo[index].sigma_dofs_per_cell << "-by-" << pinfo[index].sigma_dofs_per_cell << std::endl; 

      for (unsigned int i=0; i<pinfo[index].sigma_dofs_per_cell; ++i){
	for (unsigned int j=0; j<pinfo[index].sigma_dofs_per_cell; ++j){
	  for (unsigned int q=0; q<n_q_points; ++q){
	    MassMatrix(i,j)+= fe_values[*(sigma[index])].value(pinfo[index].sigma_dof_index[0][i],q) * fe_values[*(sigma[index])].value(pinfo[index].sigma_dof_index[0][j],q) * quad_weights[q];
	  }
	}
      }
      //std::cout << "Sigma Mass Matrix" << std::endl;
      //MassMatrix.print_formatted(std::cout);
      pinfo[index].SigmaMassMatrix = MassMatrix;

      pinfo[index].InverseSigmaMassMatrix = FullMatrix<double>(pinfo[index].sigma_dofs_per_cell,pinfo[index].sigma_dofs_per_cell);
      pinfo[index].InverseSigmaMassMatrix.invert(MassMatrix);

      MassMatrix.diagadd(-1);
      TestIdentity_norm = MassMatrix.frobenius_norm();
      std::cout << "sigma norm = " << TestIdentity_norm << std::endl;
      if( TestIdentity_norm < MassMatrix.m()*1e-14 ){
	pinfo[index].sigmaProject=false;
      } else {
	pinfo[index].sigmaProject=true;
      }
      std::cout << "Inverse Sigma Mass Matrix";
      if(!pinfo[index].sigmaProject) std::cout << " is close to the Identity. ";
      std::cout << std::endl;
      //pinfo[index].InverseSigmaMassMatrix.print_formatted(std::cout);
      std::cout << std::endl;

      pmap[index].insert( std::make_pair(fe_name,pinfo[index]) );
    } 
  }
}
