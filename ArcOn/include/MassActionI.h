/* Compute mass action function implicitly */

template <int dim>
std::vector<double> Functionals<dim>::MassActionI(Vector<double>& prev_alpha, double dt) const {

  // Set the law of mass action for each species

  std::vector<double> A(alphadim,0);

  /* A[0] =  Kfb[1]*prev_alpha(1) - Kfb[0]*prev_alpha(0); */
 
  /* A[1] =  Kfb[0]*prev_alpha(0) - Kfb[1]*prev_alpha(1); */


  std::vector<double> ReactionQuotient_f(num_reactions,1.0);
  std::vector<double> ReactionQuotient_b(num_reactions,1.0);
  
  for(unsigned int j=0;j<num_reactions;j++) {
    for(unsigned int i=0;i<alphadim;i++) {
      
      ReactionQuotient_f[j] *= std::pow(prev_alpha(i),nu_f(j,i));
      ReactionQuotient_b[j] *= std::pow(prev_alpha(i),nu_b(j,i));

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


template <int dim>
FullMatrix<double> Functionals<dim>::MassActionI_Jac(Vector<double>& prev_alpha, double dt) const {

  // Set the Jacobians of law of mass action for each species

  FullMatrix<double> J(alphadim,alphadim);

  J(0,0) =  -Kfb[0]; 
  J(0,1) =   Kfb[1];

  J(1,0) =   Kfb[0];
  J(1,1) =  -Kfb[1];

  return J;
}

