/* Set up the Runge-Kutta template */
template<int dim>
void arcOn<dim>::setup_RK()
{

  /* Strong Stability Preserving (SSP) RK method */

  RK_alpha = FullMatrix<double>(RK_stage,RK_stage);
  RK_alpha = 0;
  RK_beta  = FullMatrix<double>(RK_stage,RK_stage);
  RK_beta  = 0;
  
  RK_mu1 = FullMatrix<double>(RK_stage,RK_stage);
  RK_mu1 = 0;
  RK_mu2 = Vector<double>(RK_stage);  
  RK_mu2 = 0;

  /* Runge Kutta Chebyshev method */

  RKC_T = Vector<double>(RK_stage+1);
  RKC_U = Vector<double>(RK_stage+1); 
  RKC_Tprime = Vector<double>(RK_stage+1); 
  RKC_Tdprime = Vector<double>(RK_stage+1); 
  RKC_b = Vector<double>(RK_stage+1);
  RKC_a = Vector<double>(RK_stage+1);
  RKC_c = Vector<double>(RK_stage+1);

  RKC_T = 0;
  RKC_U = 0; 
  RKC_Tprime = 0; 
  RKC_Tdprime = 0; 
  RKC_b = 0;
  RKC_a = 0;
  RKC_c = 0;

  RKC_mu = Vector<double>(RK_stage+1);
  RKC_tildemu = Vector<double>(RK_stage+1);
  RKC_nu = Vector<double>(RK_stage+1);
  RKC_gamma = Vector<double>(RK_stage+1);
 
  RKC_mu = 0;
  RKC_tildemu = 0;
  RKC_nu = 0;
  RKC_gamma = 0;
 

  double casum;

  switch (RK_order){
  case 1:
  default:
    {
      switch (RK_stage) {
      case 1:
      default:
	{
	  RK_alpha(0,0) = 1.0;
	  RK_beta(0,0) = 1.0;
	  break;
	}
      }
      break;
    }

  case 2:
    {
      switch (RK_stage) {
      case 2:
      default:
	{ 
	  double nrk = RK_stage;
	  for(unsigned int i=0; i < RK_stage; i++){
	    for(unsigned int j=0; j < RK_stage; j++){
	      
	      if (j == i && i<RK_stage){
		
		RK_alpha(i,j) = 1.0;
		RK_beta(i,j)  = 1.0/(nrk - 1.0);
	      }
	      if (j == 0 && i == RK_stage-1){ 
		
		RK_alpha(i,j) = 1.0/nrk;
	      }
	      if (j == RK_stage - 1 &&  i == RK_stage -1){

		RK_alpha(i,j) = (nrk - 1.0)/nrk;
		RK_beta(i,j)  = 1.0/(nrk);
		
	      }
	    }
	  }

	  break;
	}
      }
      break;
    }
      

  case 3:
    {
      switch (RK_stage) {
      case 3:
      default:
	{
	  RK_alpha(0,0) = 1.0;
	  RK_alpha(1,0) = 3.0/4.0;
	  RK_alpha(1,1) = 1.0/4.0;
	  RK_alpha(2,0) = 1.0/3.0;
	  RK_alpha(2,2) = 2.0/3.0;

	  RK_beta(0,0) = 1.0;
	  RK_beta(1,1) = 1.0/4.0;
	  RK_beta(2,2) = 2.0/3.0;	
	  break;
	}

      case 4:
	{
	  RK_alpha(0,0) = 1.0;
	  RK_alpha(1,1) = 1.0;
	  RK_alpha(2,0) = 2.0/3.0;
	  RK_alpha(2,2) = 1.0/3.0;
	  RK_alpha(3,3) = 1.0;

	  RK_beta(0,0) = 1.0/2.0;
	  RK_beta(1,1) = 1.0/2.0;
	  RK_beta(2,2) = 1.0/6.0;
	  RK_beta(3,3) = 1.0/2.0;
	  break;
	}

      case 5:
	{
	  RK_alpha(0,0) = 1.0;
	  RK_alpha(1,1) = 1.0;
	  RK_alpha(2,0) = 0.355909775063327;
	  RK_alpha(2,2) = 0.644090224936674;
	  RK_alpha(3,0) = 0.367933791638137;
	  RK_alpha(3,3) = 0.632066208361863;
	  RK_alpha(4,2) = 0.237593836598569;
	  RK_alpha(4,4) = 0.762406163401431;

	  RK_beta(0,0) = 0.377268915331368;
	  RK_beta(1,1) = 0.377268915331368;
	  RK_beta(2,2) = 0.242995220537396;
	  RK_beta(3,3) = 0.238458932846290;
	  RK_beta(4,4) = 0.287632146308408;
	  break;
	}
      }
      break;
    }

  case 4:
    {
      switch (RK_stage) {
      case 5:
      default: 
	{
	  RK_alpha(0,0) = 1.0;
	  RK_alpha(1,0) = 0.44437049406734;
	  RK_alpha(1,1) = 0.55562950593266;
	  RK_alpha(2,0) = 0.62010185138540;
	  RK_alpha(2,2) = 0.37989814861460;
	  RK_alpha(3,0) = 0.17807995410773;
	  RK_alpha(3,3) = 0.82192004589227;
	  RK_alpha(4,0) = 0.00683325884039;
	  RK_alpha(4,2) = 0.51723167208978;
	  RK_alpha(4,3) = 0.12759831133288;
	  RK_alpha(4,4) = 0.34833675773694;

	  RK_beta(0,0) = 0.39175222700392;
	  RK_beta(1,1) = 0.36841059262959;
	  RK_beta(2,2) = 0.25189177424738;
	  RK_beta(3,3) = 0.54497475021237;
	  RK_beta(4,3) = 0.08460416338212;
	  RK_beta(4,4) = 0.22600748319395;
	  break;
	}

      case 6:
	{
	  RK_alpha(0,0) = 1.0;
	  RK_alpha(1,0) = 0.30948026455053;
	  RK_alpha(1,1) = 0.69051973544947;
	  RK_alpha(2,0) = 0.54205244285557;
	  RK_alpha(2,2) = 0.45794755714443;
	  RK_alpha(3,0) = 0.35984960863377;
	  RK_alpha(3,3) = 0.64015039136623;
	  RK_alpha(4,4) = 1.0;
	  RK_alpha(5,0) = 0.05776282890116;
	  RK_alpha(5,2) = 0.44216432622405;
	  RK_alpha(5,4) = 0.10115567086469;
	  RK_alpha(5,5) = 0.39891717401009;

	  RK_beta(0,0) = 0.39270746575722;
	  RK_beta(1,1) = 0.30154043149172;
	  RK_beta(2,2) = 0.19997937335132;
	  RK_beta(3,3) = 0.27954483459696;
	  RK_beta(4,4) = 0.43668618869443;
	  RK_beta(5,2) = 0.09150931531680;
	  RK_beta(5,4) = 0.04417328437472;
	  RK_beta(5,5) = 0.14911300530736;
	  break;
	}

      }
      break;
    }
  }

  //Compute the time-dependent parameters

  if(RK_stage == 1 ){RK_mu2(0) = 1.0;}
  
  if(RK_stage > 1  ){

    for(unsigned int k = 0; k <= RK_stage-1; k++){
      for(unsigned int i = 1; i <= RK_stage; i++){
      
  	casum = 0.0;
	
	if(i>1){
	  for(unsigned int l=k+1; l <= i-1; l++){
	    casum = casum + RK_mu1(l-1,k)*RK_alpha(i-1,l);
	  }
	}
	  
	RK_mu1(i-1,k) = RK_beta(i-1,k) + casum;
	
      }
    }
    
    for(unsigned int k=1; k<RK_stage ;k++){

      RK_mu2(k) = 0.0;
      
      for(unsigned int l=0;l<=k; l++){

  	RK_mu2(k) = RK_mu2(k) + RK_mu1(k-1,l);

      }

    }

    double RKC_omega0;
    double RKC_omega1;

    //Compute the RKC version

    RKC_T(0) = 1.0;
    RKC_U(0) = 1.0;
   
    RKC_omega0 = 1.0 + eps_const*std::pow(RK_stage,-2.0);

    RKC_c(0) = 0.0;

    RKC_T(1) = RKC_omega0;
    RKC_U(1) = 2.0*RKC_omega0;

    RKC_Tprime(1) = std::pow(RK_stage,2.0);
    RKC_Tdprime(1) = 1/3 * (std::pow(RK_stage,2.0)*(std::pow(RK_stage,2.0)-1.0));

    for(unsigned int j=0;j<=RK_stage; j++){

      //Chebyshev polynomial stuff

      if(j>1){

	RKC_T(j) = 2.0 * RKC_omega0 * RKC_T(j-1) - RKC_T(j-2);
	RKC_U(j) = 2.0 * RKC_omega0 * RKC_U(j-1) - RKC_U(j-2);
      
      }


      if(j>1){

	RKC_Tprime(j) = j*RKC_U(j-1);
	RKC_Tdprime(j) = (j*(j*RKC_T(j) - RKC_omega0*RKC_U(j-1))/(std::pow(RKC_omega0,2.0)-1.0));
	RKC_b(j) = RKC_Tdprime(j)*std::pow(RKC_Tprime(j),-2.0);

      }
    
      RKC_a(j)  = 1.0 - RKC_b(j)*RKC_T(j);

    }

    if(RK_stage > 1){
      RKC_b(0) = RKC_b(2);
    }

    RKC_b(1) = 1.0/RKC_omega0;

    RKC_omega1 =  RKC_Tprime(RK_stage)/RKC_Tdprime(RK_stage);

    RKC_tildemu(1)= std::pow(RKC_omega0,-1.0)*RKC_omega1;

    for(unsigned int j=2;j<=RK_stage; j++){

      RKC_mu(j) = 2.0*RKC_b(j)*RKC_omega0/RKC_b(j-1);

      RKC_tildemu(j) = 2.0*RKC_b(j)*RKC_omega1/RKC_b(j-1);

      RKC_nu(j) = - RKC_b(j)/RKC_b(j-2);

      RKC_gamma(j) = - RKC_a(j-1)*RKC_tildemu(j);

      RKC_c(j) = (std::pow(j,2.0)-1.0)/(std::pow(RK_stage,2.0)-1.0);

    }

    RKC_c(1) = (RKC_c(2))/ (4.0*RKC_omega0);

    RKC_c(RK_stage) = 0.0;
    RKC_c(RK_stage) = 1.0;

  }
}
