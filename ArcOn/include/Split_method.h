 /* Implement the splitting scheme, for "operator splitting" */
template <int dim>
void arcOn<dim>::calc_convdiff(SolutionVector& subdomain_solution, double delta_t, double current_time) 
{
  for (unsigned int component=0; component< alphadim; ++component){
    naive_div_flux_integrated[component] = 0.0;
  }

  div_flux_integrated =  naive_div_flux_integrated;
  RK_div_flux_integrator( subdomain_solution, delta_t, current_time, div_flux_integrated );
  subdomain_solution  = div_flux_integrated;  
    
}

template <int dim>
void arcOn<dim>::calc_reaction(SolutionVector& subdomain_solution, double delta_t, double current_time)
{
  
  if(method == 1){

    for (unsigned int component=0; component< alphadim; ++component){
      naive_MassAction_integrated[component] = 0.0;
    }

    MassAction_integrated = naive_MassAction_integrated;
    RK_MassAction_integrator( subdomain_solution, delta_t, current_time, MassAction_integrated );
    subdomain_solution = MassAction_integrated;
    
  } 
  
  else { 
    
    std::cout << "This has been turned off due to lack of candidate subsystems" << std::endl;
    
  }
  
}
