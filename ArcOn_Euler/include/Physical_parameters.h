/* These need to be documented in detail, and diffusion needs to be fixed */
/* Generally, these are constants and functionals having to due with associated */
/* MHD parameters, like viscosity, diffusion, pressure, & kinetic expressions, etc .... */

template<int dim>
class Functionals
{
public:
  Functionals();

  const  FullMatrix<double >&		Binary_constants() const;
  const  TableBase< 3, double >&	Ficks() const;
  void							Transport(double P, std::vector<double>& A);
  void							CE_Transport(double P, std::vector<double>& A);
  std::vector<double>				MassAction(Tensor<1,alphadim>& prev_alpha, double dt) const;
  std::vector<double>				MassAction_exact(Tensor<1,alphadim>& prev_alpha, Tensor<1,alphadim>& init_alpha, double dt) const;
  std::vector<double>				MassActionI(Vector<double>& prev_alpha, double dt) const;
  FullMatrix<double>				MassActionI_Jac(Vector<double>& prev_alpha, double dt) const;
  std::vector<double>				MassActionProd(Tensor<1,alphadim>& alpha) const;
  double							Pressure(Tensor<1,alphadim>& prev_alpha);
  double                                        kappa_sum;
  double                                        Keq_prod;

  std::vector<double> MM;
  std::vector<double> rad;
  std::vector<double> pureredtemp;
  // std::vector<double> Rc;
  std::vector<double> Knud;
  std::vector<double> Kfb;
  std::vector<double> Ea;
  std::vector<double> gamma;
  std::vector<double> Keq;
  Vector<double> kappa_circ;
  Vector<double> sum_nu_f;
  Vector<double> sum_nu_b;

  FullMatrix<double> nu_f; 
  FullMatrix<double> nu_b;
  

  double Homogeneous_Diffusion;
  double vartheta;
  double Ri;
  double Hbf;
  double proton;
  unsigned int Num_reactions;

  bool fixedBinary;
  bool fixedFicks;

private:
  FullMatrix< double > _BinaryConstants; // alphadim x alphadim
  FullMatrix< double > sigmaij;
  FullMatrix< double > bijmax;
  FullMatrix< double > reducedtemp;
  TableBase< 3, double > D; // alphadim x dim x dim 

};

template<int dim>
Functionals<dim>::Functionals() {

  fixedBinary = true;
  fixedFicks  = true;

  _BinaryConstants = FullMatrix< double >( alphadim, alphadim);
  sigmaij = FullMatrix< double >( alphadim, alphadim);
  bijmax  = FullMatrix< double >( alphadim, alphadim);
  reducedtemp = FullMatrix< double >( alphadim, alphadim);
  D  = TableBase< 3, double >( TableIndices<3>(alphadim, dim, dim) );

  nu_f = FullMatrix<double>(num_reactions,alphadim);
  nu_b = FullMatrix<double>(num_reactions,alphadim);

  kappa_circ = Vector<double>(alphadim);
  sum_nu_f = Vector<double>(alphadim);
  sum_nu_b = Vector<double>(alphadim);


  /* Set the molar masses */

  double MolarMass [alphadim] = { 1.0,1.0 }; //

  for(unsigned int i=0;i<alphadim;i++) {

    MM.push_back(MolarMass[i]);

  } 

  /* radii for collisional cross sections */
	
  double radius [alphadim] = {3.35,1.2};

  for(unsigned int i=0;i<alphadim;i++) {

    rad.push_back(radius[i]);

  } 

  /* Reduced temperatures */

  double purereducedtemp [alphadim] = {0.343,0.364};

  for(unsigned int i=0;i<alphadim;i++) {

    pureredtemp.push_back(purereducedtemp[i]);

  } 

  /* Basic parameters */

  Homogeneous_Diffusion        = 0.01;  // sum of constants (this is too high, but set for now)
  vartheta                     = 297;   // constant temp
  Ri                           = 8.314472;     // J/K*mol; ideal gas law
  Hbf                          = 0.0;    // Mass transfer correction, should not need this for quiescent case
  proton                       = 0.8;      // proton concentration

  /* Set the coefficients of the reaction rates */

  double ReactionRates [2*alphadim] = { 0.6,0.4,0,0};         // we use the experimental values here, forward,backward ...


  /* Set backward and forward reaction rates */

  for(unsigned int i=0;i<2*alphadim;i++) {

    Kfb.push_back(ReactionRates[i]); //we use experimental values here at constant Temp

  }

  // double scalfact = 1;
  // note that the first half are forward rxn rates, second half backward rxn rates
  // kf1, kf2, kf3, kf4, kf5,
  // double scal = 314.0;

  /* double Reactor_size [num_reactions] = { 3}; */

  /* for(unsigned int i=0;i<num_reactions;i++) { */

  /*   Rc.push_back(Reactor_size[i]); */

  /* } */

  double Stoichdiff [alphadim] = { 1.0,1.0};

  for(unsigned int i=0;i<alphadim;i++) {

    Knud.push_back(Stoichdiff[i]);

  }

  /* Set the stoichiometric coefficients */		
  /* WARNING!!! Decide which is better here */		
							
  //double Stoich [alphadim] = {1.0,1.0};			
							
  for(unsigned int i=0;i<num_reactions;i++) {		
    for(unsigned int j=0;j<alphadim;j++) {		
							
      nu_f(i,j)= 1.0;					
      nu_b(i,j)= 1.0;					
							
    }							  }							
							
  /* or hardcode the stoichiometric coeffs by hand */	
							
  nu_f(0,0) = 1.0;					
  nu_b(0,0) = 0.0;
  nu_f(0,1) = 0.0;
  nu_b(0,1) = 1.0;
    

  /* Set the reaction activation energies */

  double ActivationEnergy [2*alphadim] = { 1,2,3,4}; // these are not used yet

  for(unsigned int i=0;i<2*alphadim;i++) {

    Ea.push_back(ActivationEnergy[i]);

  }

  /* Set equilibrium constants using reaction rates Kfb */

  for(unsigned int i=0;i<num_reactions;i++) {

    Keq.push_back(ReactionRates[2*i]/ReactionRates[(2*i)+1]); //we use experimental values here at constant Temp

    //std::cout << "Keq[i] = " << Keq[i] << ", i= "<< i << std::endl;

  }

  /* Some Auxiliary constants */

  Keq_prod = 0.0;

    for(unsigned int j=0;j<num_reactions;j++) {
          
      Keq_prod = Keq_prod + std::pow(Keq[j],-1.0);

    }


    sum_nu_f = 0.0;
    sum_nu_b = 0.0;

  for(unsigned int i=0;i<alphadim;i++) {
    for(unsigned int j=0;j<num_reactions;j++) {
      
      sum_nu_f(i) = sum_nu_f(i) + nu_f(j,i);
      sum_nu_b(i) = sum_nu_b(i) + nu_b(j,i);
      
    }
  }


  /* Set the species independent constants kappa_circ for entropy calc*/

  for(unsigned int j=0;j<alphadim;j++) {

    kappa_circ(j) = 1.0;

    //std::cout << "kappa_circ(j) = " << kappa_circ(j) << std::endl;

  }
 
  for(unsigned int j=0;j<alphadim;j++) {
 
    for(unsigned int i=0;i<num_reactions;i++) {

      kappa_circ(j) *=  std::pow( (double) Keq[i],(double)(-1.0/((alphadim*(nu_b(i,j)-nu_f(i,j))))));

      //std::cout << "Keq[i] = " << Keq[i] << ", alphadim = " << alphadim << ", nu_b(i,j) = " << nu_b(i,j) << ", nu_f(i,j) = " << nu_f(i,j) << ", kappa_circ[j] = " << kappa_circ(j) << ", i = " << i << ", j = " << j << ", temp = " << temp << ", temp2 = " << temp2 << std::endl;

    }
  }

  kappa_sum = 0.0;

  for(unsigned int i=0;i<alphadim;i++) {
    kappa_sum = kappa_sum + kappa_circ(i);
  }

  /* Set the adiabtic indices */

  /* double Gammas [alphadim] = { 5/3,5/3}; //Using helium bath reference */

  /* std::vector<double> gamma; */

  /* for(unsigned int i=0;i<2*alphadim;i++) { */

  /*   gamma.push_back(Gammas[i]); */

  /* } */

}

template<int dim>
void Functionals<dim>::Transport(double P, std::vector<double>& A) {
  
  _BinaryConstants = 0;
  D.fill(0);

  for(unsigned int i=0;i<alphadim;i++){
    for(unsigned int j=0;j<alphadim;j++){
      if (fixedBinary) {  // For known binary constants
	_BinaryConstants( i,j ) = Homogeneous_Diffusion;
      } else {   // For functional coefficients
	_BinaryConstants( i,j ) = i*j;
      }
    }
  }


  for(unsigned int i=0;i<alphadim;i++){
    for(unsigned int j=0;j<dim;j++){
      if (fixedFicks) { // For constant in 3D
	if (i<2){
	  D( TableIndices<3>(i,j,j) ) = Homogeneous_Diffusion;}
	else{
	  D( TableIndices<3>(i,j,j) ) = 0.0;}
	  
      } else { // For functional coefficients 
	for(unsigned int k=0;k<dim;k++){
	  D( TableIndices<3>(i,j,k) ) = i*j*k*P*A[i];
	}
      }
    }
  }

}

template<int dim>
void Functionals<dim>::CE_Transport(double P, std::vector<double>& A) {
  _BinaryConstants = 0;
  sigmaij          = 0;
  reducedtemp      = 0;
  bijmax           = 0;
  D.fill(0);

  for(unsigned int i=0;i<alphadim;i++){
    for(unsigned int j=0;j<alphadim;j++){
      bijmax( i,j ) = rad[i]+rad[j];
    }
  }

  for(unsigned int i=0;i<alphadim;i++){
    for(unsigned int j=0;j<alphadim;j++){
      sigmaij( i,j ) = M_PI*bijmax(i,j)*bijmax(i,j);
    }
  }

  double sumAs=0;
  for(unsigned int i=0;i<alphadim;i++){
    sumAs = sumAs + A[i];
  }

  for(unsigned int i=0;i<alphadim;i++){
    for(unsigned int j=0;j<alphadim;j++){
		  
      reducedtemp( i,j ) = (pureredtemp[i] * A[i] + pureredtemp[j] * A[j])/ (A[i]+A[j]);
    }
  }

  for(unsigned int i=0;i<alphadim;i++){
    for(unsigned int j=0;j<alphadim;j++){
      if (fixedBinary) {  // For known binary constants
	_BinaryConstants( i,j ) = 1e-7;
      } else {   // For functional coefficients
	_BinaryConstants( i,j ) = 2.628e-3 / ( sigmaij(i,j)*sigmaij(i,j)*reducedtemp(i,j) );

	//if (A[1] > .08){
	//std::cout << "_BC(i,j)=" << _BinaryConstants( i,j ) << "   sigmaij=" << sigmaij(i,j) << "    reducedtemp(i,j) " << reducedtemp(i,j) << " A[i]=" << A[i] << " A[j]=" << A[j] << " sumAs=" << sumAs <<  " " << i << "," << j <<  std::endl;}
      }
    }
  }

  if (fixedFicks) { // For constant in 3D
    for(unsigned int i=0;i<alphadim;i++){
      for(unsigned int j=0;j<dim;j++){
	if (i<2){
	  D( TableIndices<3>(i,j,j) ) = Homogeneous_Diffusion;}
	else{
	  D( TableIndices<3>(i,j,j) ) = 0.0;}
      }
    }
  } else {

    FullMatrix<double> BMD(alphadim,alphadim);

    double vartheta3 = vartheta*vartheta*vartheta;

    for(unsigned int i=0;i<alphadim;i++){
      for(unsigned int j=0;j<alphadim;j++){
	BMD(i,j) = _BinaryConstants(i,j)/P * std::sqrt(vartheta3 * (MM[i]+MM[j])/(2*MM[i]*MM[j]) );
      }
    }

    double nSum = 0;
    std::vector<double> n(alphadim);
    for(unsigned int i=0; i<alphadim; i++){
      n[i] = A[i]/MM[i];
      nSum += n[i];
    }
    std::vector<double> x(alphadim);
    for(unsigned int i=0; i<alphadim; i++){
      x[i] = n[i] / nSum;
    }

    FullMatrix<double> K(alphadim,alphadim);

    for(unsigned int i=0; i<alphadim; i++){
      std::vector<double> xD(alphadim);
      double xDsum = 0;
      for(unsigned int j=0;j<alphadim;j++){
	xD[j] = x[j] / BMD(i,j);
	xDsum += xD[j];
      }
      for(unsigned int j=0;j<alphadim;j++){
	if( i==j){
	  K(i,j)=0;
	} else {
	  K(i,j) = x[i] / BMD(i,j) + MM[j]/MM[i]*(xDsum-x[i]/BMD(i,i));
	}
      }
    }

    /*
      std::cout << "A=" << A << std::endl;
      std::cout << "n=" << n << " sum=" << nSum << std::endl;
      std::cout << "x=" << x << std::endl;
      std::cout << "BMD matrix" << std::endl;
      BMD.print_formatted(std::cout);
      std::cout << "K matrix" << std::endl;
      K.print_formatted(std::cout);
    */

    FullMatrix<double> Kinv(alphadim,alphadim);
    Kinv.invert(K);

    FullMatrix<double> MMD(alphadim,alphadim);
    Vector<double> MMDsum(alphadim);
    for(unsigned int i=0; i<alphadim; i++){
      for(unsigned int j=0;j<alphadim;j++){
	MMD(i,j) = 1/MM[j]/nSum*(Kinv(i,j)-Kinv(i,i));
	MMDsum(i) += MMD(i,j); 
      }
    }

    //No viscosity terms
    for(unsigned int i=0; i<alphadim; i++){
      for(unsigned int k=0;k<dim;k++){
	D( TableIndices<3>(i,k,k) ) = MMDsum(i);
      }
    }

    /*
      std::cout << "Kinv matrix" << std::endl;
      Kinv.print_formatted(std::cout);
      std::cout << "MMD matrix" << std::endl;
      MMD.print_formatted(std::cout);
    */
    double Dmin=std::numeric_limits<double>::max();
    double Dmax=-std::numeric_limits<double>::max();
    double Dsum=0;
    for(unsigned int i=0;i<alphadim;i++){
      Dsum += MMDsum(i);
      if (Dmin>MMDsum(i)) Dmin=MMDsum(i);
      if (Dmax<MMDsum(i)) Dmax=MMDsum(i);
    }
    //std::cout << "D=" << MMDsum << " min=" << Dmin << " max=" << Dmax << " avg=" << Dsum/alphadim << std::endl;

  }
}		

template<int dim>
const FullMatrix< double >& Functionals<dim>::Binary_constants() const {
  return _BinaryConstants;
}


template<int dim>
const  TableBase< 3, double >&  Functionals<dim>::Ficks() const {
  
  return D;
}
