/*Computes a barotropic pressure function.*/
/* At present not in use .... */

template <int dim>
double Functionals<dim>::Pressure(Tensor<1, alphadim>& prev_alpha){

  double P = 0;

  for(unsigned int i=0;i<alphadim;i++){
    P += std::pow( prev_alpha[i], gamma[i] );

  }
  return P;
}
