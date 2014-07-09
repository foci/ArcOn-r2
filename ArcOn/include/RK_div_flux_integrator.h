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
      
      //if (s=0){naive_RK_solution[s] = 0.0;}

      for(unsigned int k=0; k<=s; k++){
	for (unsigned int component=0; component< alphadim; ++component){
	  if(RK_alpha[s][k] != 0){
	    
	    //naive_RK_solution[s+1][component].add(RK_alpha[s][k],naive_RK_solution[s][component]);
	    //naive_RK_solution[k][component]   *= RK_alpha[s][k];
	    naive_RK_solution_temp[k][component] = naive_RK_solution[k][component];
	    naive_RK_solution_temp[k][component] *= RK_alpha[s][k];
	    naive_RK_solution[s+1][component] += naive_RK_solution_temp[k][component];
	    
	    //pcout << "Timestepping, s = " << s << ", k = " << k << std::endl; 
	    //pcout << "RK-alpha = " << RK_alpha[s][k] << std::endl;
	  }
	  if(RK_beta[s][k] != 0 ){
	    
	    //naive_RK_solution[s+1][component].add(delta_t*RK_beta[s][k],naive_RK_div_flux[k][component]);
	    naive_RK_div_flux_temp[k][component] = naive_RK_div_flux[k][component];
	    naive_RK_div_flux_temp[k][component] *= delta_t*RK_beta[s][k];
	    naive_RK_solution[s+1][component] += naive_RK_div_flux_temp[k][component];
	    //pcout << "Timestepping, s = " << s << ", k = " << k << std::endl; 
	    //pcout << "RK-beta = " << RK_beta[s][k] << std::endl;

	  }
	}
      }

     for (unsigned int component=0; component< alphadim; ++component){
      	naive_RK_solution[s+1][component].compress(VectorOperation::add);
      }
      
      RK_solution[s+1] = naive_RK_solution[s+1];
      //periodicity_map( RK_solution[s+1],delta_t);

    }
  }

  /* This can be easily set for RKC methods (low priority at the moment) */
  if (RKtype == 2){

/*     for(unsigned int s=0; s<RK_stage; s++){ */

/*      for (unsigned int component=0; component< alphadim; ++component){ */
/*        RK_solution[s+1][component].reinit(  mpi_communicator,locally_owned_dofs[component],locally_relevant_dofs[component] ); */

/*        RK_div_flux[s][component].reinit(  mpi_communicator,locally_owned_dofs[component],locally_relevant_dofs[component] ); */

/*        RK_solution[s+1][component].compress(); */
/*        RK_div_flux[s][component].compress(); */
/*      } */

/*       calculate_div_flux( RK_solution[s] , delta_t*RKC_c(s), current_time, RK_div_flux[s] ); */

/*     for (unsigned int component=0; component< alphadim; ++component){ */
/*       RK_solution[s][component].update_ghost_values(); */
/*       RK_solution[s][component].compress(); */
/*       RK_div_flux[s][component].update_ghost_values(); */
/*       RK_div_flux[s][component].compress(); */
/*     } */

/*       for (unsigned int component=0; component< alphadim; ++component){ */

/* 	if(s==0){ */
/* 	  RK_solution[s+1][component].add(1.0,RK_solution[s][component]); */
/* 	  RK_solution[s+1][component].add(RKC_tildemu(1),RK_div_flux[s][component]); */
/* 	} */
/* 	else{ */
/* 	  RK_solution[s+1][component].add(1.0-RKC_mu(s+1)-RKC_nu(s+1),RK_solution[0][component]); */
/* 	  RK_solution[s+1][component].add(RKC_mu(s+1),RK_solution[s][component]); */
/* 	  RK_solution[s+1][component].add(RKC_nu(s+1),RK_solution[s-1][component]); */
/* 	  RK_solution[s+1][component].add(RKC_tildemu(s+1),RK_div_flux[s][component]); */
/* 	  RK_solution[s+1][component].add(RKC_gamma(s+1),RK_div_flux[0][component]); */
/* 	} */
/* 	//RK_solution[s+1][component].print(std::cout); */
/*       } */
/*       //std::cout << s << "\t" << std::flush; */
/*     } */
  }

  div_flux_integrated = RK_solution[RK_stage];

}

