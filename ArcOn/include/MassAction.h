/* Compute mass action function, in case we have chemical reactions proper */
template <int dim>
std::vector<double> Functionals<dim>::MassAction(Tensor<1, alphadim>& prev_alpha, double dt) const {

  // Set the law of mass action for each species

  std::vector<double> A(alphadim,0);

  std::vector<double> B(alphadim,0);

  //B[0] =  Kfb[1]*prev_alpha[1] - Kfb[0]*prev_alpha[0];
 
  //B[1] =  Kfb[0]*prev_alpha[0] - Kfb[1]*prev_alpha[1];

  //std::cout << "B[0] = " << B[0] << std::endl;
  //std::cout << "B[1] = " << B[1] << std::endl;
    
  /* A[0] = 1.0/dt * ( std::exp(- Kfb[0] * dt) * (prev_alpha[0] - prev_alpha[1]/(Keq) ) + prev_alpha[1]/(Keq) - prev_alpha[0]);  */
 
  /* A[1] = 1.0/dt * ( std::exp(- Kfb[2] * dt) * (prev_alpha[1] - (prev_alpha[0] * Keq ) ) +  (prev_alpha[0] * Keq ) - prev_alpha[1]); */ 

  //std::cout << "first= " << A[0] << std::endl;
  //std::cout << "second= " << A[1] << std::endl;

  std::vector<double> ReactionQuotient_f(num_reactions,1.0);
  std::vector<double> ReactionQuotient_b(num_reactions,1.0);

  for(unsigned int j=0;j<num_reactions;j++) {
    for(unsigned int i=0;i<alphadim;i++) {
  
      ReactionQuotient_f[j] *= std::pow(prev_alpha[i],nu_f(j,i));
      ReactionQuotient_b[j] *= std::pow(prev_alpha[i],nu_b(j,i));

    }
  }

  for(unsigned int i=0;i<alphadim;i++) {

    for(unsigned int j=0;j<num_reactions;j++) {

      A[i] = A[i] + MM[i] * ( nu_b(j,i) - nu_f(j,i) ) * ( Kfb[2*j]*ReactionQuotient_f[j] -  Kfb[(2*j)+1]*ReactionQuotient_b[j] );

      //std::cout << "A[i] = " << A[i] << std::endl;
	
    }

  }

  return A;
}
