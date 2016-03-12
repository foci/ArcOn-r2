/* Implements the time integrators for the source and reaction terms */
template <int dim>
void arcOn<dim>::RK_MassAction_integrator(SolutionVector& subdomain_solution, double delta_t, double current_time, SolutionVector& MassAction_integrated )
{
  for (unsigned int component=0; component< alphadim; ++component){
    for(unsigned int s=0; s<RK_stage+1; s++){

      naive_RK_solution[s][component] = 0.0;
      naive_RK_solution_temp[s][component] = 0.0;
      naive_RK_MassAction[s][component] = 0.0;
      naive_RK_MassAction_temp[s][component] = 0.0;

    }
  }

  RK_solution = naive_RK_solution;
  RK_solution[0] = subdomain_solution;

  if (RKtype == 1){

    for(unsigned int s=0; s<RK_stage; s++){

      Calculate_MassAction_Explicit(RK_solution[s], delta_t, current_time+delta_t*RK_mu2(s), naive_RK_MassAction[s]);
      naive_RK_solution[s] =  RK_solution[s];

      for(unsigned int k=0; k<=s; k++){
	for (unsigned int component=0; component< alphadim; ++component){
	  if(RK_alpha[s][k] != 0 ){
	    //naive_RK_solution[s+1][component].add(RK_alpha[s][k],naive_RK_solution[s][component]);
	    naive_RK_solution_temp[k][component] = naive_RK_solution[k][component];
	    naive_RK_solution_temp[k][component] *= RK_alpha[s][k];
	    naive_RK_solution[s+1][component] += naive_RK_solution_temp[k][component];

	  }
	  if(RK_beta[s][k] != 0 ){
	    //naive_RK_solution[s+1][component].add(delta_t*RK_beta[s][k],naive_RK_MassAction[k][component]);
	    naive_RK_MassAction_temp[k][component] = naive_RK_MassAction[k][component];
	    naive_RK_MassAction_temp[k][component] *= delta_t*RK_beta[s][k];
	    naive_RK_solution[s+1][component] += naive_RK_MassAction_temp[k][component];
	  }
	}
      }
      for (unsigned int component=0; component< alphadim; ++component){
	naive_RK_solution[s+1][component].compress(VectorOperation::add);
      }
      RK_solution[s+1] = naive_RK_solution[s+1];
     }
  }

  /* Easily updated for RKC */
  if (RKtype == 2){
    
  }
  
  MassAction_integrated = RK_solution[RK_stage];

}

