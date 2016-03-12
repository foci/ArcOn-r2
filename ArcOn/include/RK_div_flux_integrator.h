/* Implements the time integrators for diffuse and convective components */
template <int dim>
void arcOn<dim>::RK_div_flux_integrator(SolutionVector& subdomain_solution, double delta_t, double current_time, SolutionVector& div_flux_integrated)
{
  for (unsigned int component=0; component< alphadim; ++component){
    for(unsigned int s=0; s<RK_stage+1; s++){

      naive_RK_solution[s][component] = 0.0;
      naive_RK_solution_temp[s][component] = 0.0;
      naive_RK_div_flux[s][component] = 0.0;
      naive_RK_div_flux_temp[s][component] = 0.0;

    }
  }
  
  RK_solution = naive_RK_solution;
  RK_solution[0] = subdomain_solution;
  
  if (RKtype == 1){

    for(unsigned int s=0; s<RK_stage; s++){
      
      calculate_div_flux( RK_solution[s], delta_t, current_time + delta_t*RK_mu2(s), naive_RK_div_flux[s] );
      periodicity_map( RK_solution[s],delta_t); //wrap in periodicity boolean
      naive_RK_solution[s] =  RK_solution[s];

      for(unsigned int k=0; k<=s; k++){
	for (unsigned int component=0; component< alphadim; ++component){
	  if(RK_alpha[s][k] != 0){
	    
	    naive_RK_solution_temp[k][component] = naive_RK_solution[k][component];
	    naive_RK_solution_temp[k][component] *= RK_alpha[s][k];
	    naive_RK_solution[s+1][component] += naive_RK_solution_temp[k][component];

	  }
	  if(RK_beta[s][k] != 0 ){
	    
	    naive_RK_div_flux_temp[k][component] = naive_RK_div_flux[k][component];
	    naive_RK_div_flux_temp[k][component] *= delta_t*RK_beta[s][k];
	    naive_RK_solution[s+1][component] += naive_RK_div_flux_temp[k][component];

	  }
	}
      }

      for (unsigned int component=0; component< alphadim; ++component){
      	naive_RK_solution[s+1][component].compress(VectorOperation::add);
      }
      
      RK_solution[s+1] = naive_RK_solution[s+1];

    }
  }

  /* This can be easily set for RKC methods (low priority at the moment) */
  if (RKtype == 2){
  }
  
  div_flux_integrated = RK_solution[RK_stage];
  
}

