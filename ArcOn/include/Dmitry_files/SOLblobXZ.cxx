/*******************************************************************
 * Bare Bones Bout
 *
 * D. Meyerson 2013
 *******************************************************************/

#include <bout.hxx>
#include <boutmain.hxx>

#include <derivs.hxx>
#include <initialprofiles.hxx>
#include <invert_laplace.hxx>
#include <invert_parderiv.hxx>
#include <invert_laplace_gmres.hxx>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <sstream>
// Evolving variables 
Field3D u, n; //vorticity, density

//Background field
Field2D n0;
//derived variables
Field3D phi,brkt;
int phi_flags;

//other fields
Field3D test1, test2,modon_IC;


//Constrained 
Field3D C_phi;

//other params
BoutReal alpha, nu, mu,gam, beta;


//solver options
bool use_jacobian, use_precon;

//experimental
bool use_constraint;

int MZ;

FieldGroup comms; // Group of variables for communications

double Newton_root(double gamma_in,double kappa_in);

const Field3D mybracket(const Field3D &phi, const Field3D &A);
int jacobian(BoutReal t); // Jacobian-vector multiply
int precon(BoutReal t, BoutReal cj, BoutReal delta); // Preconditioner

//int precon_phi(BoutReal t, BoutReal cj, BoutReal delta);
//int jacobian_constrain(BoutReal t); // Jacobian-vector multiply

double c, lambda,a,r;//kappa;//,gamma;
  // beta =1.0;
  // lambda = 1.0;
  // c = 1.0;
  // a = 1.0;
double modon(double x, double y,double kappa,double gamma);


#define bracket3D(f, g) ( b0xGrad_dot_Grad(f, g) )
#define LapXZ(f)(mesh->g11*D2DX2(f) + mesh->g33*D2DZ2(f))

int physics_init(bool restarting)
{
  // // 2D initial profiles
  // Field2D N0, P0;
  // Vector2D V0;
  // BoutReal v0_multiply;

  // // Read initial conditions

  // mesh->get(N0, "density");
  // mesh->get(P0, "pressure");
  // V0.covariant = false; // Read contravariant components
  // V.covariant = false; // Evolve contravariant components
  // mesh->get(V0, "v");
  // g.covariant = false;
  // mesh->get(g, "g");
  Options *globaloptions = Options::getRoot();
  Options *options = globaloptions->getSection("physics");
  Options *solveropts = globaloptions->getSection("solver");

  OPTION(options, phi_flags, 0);
  OPTION(options, alpha,3e-5);
  OPTION(options, nu, 2e-3);
  //OPTION(options, mu, 0.040);
  OPTION(options, mu, 2e-3);
  OPTION(options, gam, 1e1);
  OPTION(options, beta, 6e-4);

  OPTION(globaloptions,MZ,33);

  OPTION(solveropts,use_precon,false);
  OPTION(solveropts,use_jacobian,true);
  OPTION(solveropts,use_constraint,false);


  // read options
  //phi.setLocation(CELL_CENTRE);
  
  //overide
  //beta = 1e-5;
  //mu = 1e-2;
  //nu = 1e-2;
  n0 = 0.0;

  bout_solve(u, "u");
  comms.add(u);
  //phi = invert_laplace(u, phi_flags);
  static Field2D A = 0.0;
  static Field2D C = 1e-12;
  static Field2D D = 1.0;
  
  phi = invert_laplace(u, phi_flags,&A,&C,&D);
  //Laplacian *lap = Laplacian::create();
  
  bout_solve(n, "n");
  comms.add(n);
  //u.setBoundary("u");
  //n.setBoundary("n");

  //brkt = b0xGrad_dot_Grad(phi, u);

  //dump.add(phi,"phi",1);
  dump.add(brkt,"brkt",1);
  dump.add(test1,"test1",1);
  dump.add(test2,"test2",1);

  if (use_constraint){
    //solver->setPrecon(precon_phi);
    //solver->setJacobian(jacobian_constrain);
    phi.setBoundary("phi");
    //bout_solve(phi,"phi");
  }else
    dump.add(phi,"phi",1);

  comms.add(phi); //super duper important 

  if (use_jacobian)
    solver->setJacobian(jacobian);

  // if (use_precon)
  //   solver->setPrecon(precon);
    
  output.write("use jacobian %i \n",use_jacobian);
  output.write("use precon %i \n",use_precon);
  output.write("DONE WITH PHYSICS_INIT\n");
 

  modon_IC.allocate();
  BoutReal ***m = modon_IC.getData();
  BoutReal ***n_ijk = n.getData();
  double temp;


  BoutReal xlen = (mesh->GlobalNx)*((mesh->dx)[0][0]);
  xlen = 1.0;
  

  a = .1*xlen;//*((mesh->dx)[0][0]);
  c = 1.00*(xlen/1e0);
  lambda = .10/xlen;//;*((mesh->dx)[0][0]);
  
  //a = 1;//*((mesh->dx)[0][0]);
  // c = 1.00*(xlen/1e0);
  // lambda = .0010/xlen;//;*((mesh->dx)[0][0]);
  
  double kappa = a*sqrt(lambda*beta)/c;
  double gamma = Newton_root(2*M_PI,kappa);


  for(int jz=0;jz<mesh->ngz;jz++) 
    for(int jx=0;jx<mesh->ngx;jx++){
      output.write("use jacobian %g, %g \n",kappa,gamma);
      temp = modon((mesh->GlobalX(jx)-.5)


,(mesh->dz*jz/mesh->zlength) - .5,
		   kappa,gamma);
      
      for(int jy=0;jy<mesh->ngy;jy++){
  	m[jx][jy][jz]=temp;
	n_ijk[jx][jy][jz] = (lambda/c)*(temp - c*(mesh->GlobalX(jx)-.5));
      }
    }

  phi = modon_IC;


  u = LapXZ(phi);
  //n = (lambda/c)*(phi - 
  
  return 0;
}



int physics_run(BoutReal t)
{
  // Run communications
  mesh->communicate(comms);
  //phi = invert_laplace(u, phi_flags);
  
  static Field2D A = 0.0;
  static Field2D C = 1e-24;
  static Field2D D = 1.0;
  
  phi = invert_laplace(u, phi_flags,&A,&C,&D);
  //phi = LaplaceGMRES
 
  phi.applyBoundary("dirichlet");
  // Density
  //f = lowPass(f,1);
  //f = lowPass(g,1);
  mesh->communicate(comms);
  //mesh->communicate(phi);
  ddt(u)=0;
  ddt(n)=0;
 

 
  ddt(u) += bracket3D(phi,u);
  //ddt(u) += alpha * phi;
  //ddt(u) += nu * LapXZ(u);
  
  //ddt(u) -= beta * DDY(n); 
  ddt(u) -= beta* DDZ(n);
 
  
  //ddt(n) -= mybracket(phi,n);
  ddt(n)  += bracket3D(phi,n);
 

  
 
  return 0;
}


const Field3D mybracket(const Field3D &phi, const Field3D &A)
{
  Field3D dpdx, dpdy, dpdz;
  Field3D vx, vy, vz;
  Field3D result;

  //output.write("mesh->Bxy = %e\n", (mesh->J*sqrt(mesh->g_22))[2][2]);

#ifdef CHECK
  int msg_pos = msg_stack.push("b0xGrad_dot_Grad( Field3D , Field3D )");
#endif

  // Calculate phi derivatives
  #pragma omp parallel sections
  {
    #pragma omp section
    dpdx = DDX(phi); 
    
    #pragma omp section
    dpdy = DDY(phi);
    
    #pragma omp section
    dpdz = 0;
  }
  
  // Calculate advection velocity
  #pragma omp parallel sections
  {

    #pragma omp section
    vx = mesh->g_23*dpdz - mesh->g_33*dpdy; 
    //vx = mesh->g_22*dpdz - mesh->g_23*dpdy;
    
    #pragma omp section
    vy = mesh->g_33*dpdx - mesh->g_13*dpdz;
      //vy = mesh->g_23*dpdx - mesh->g_12*dpdz;
    
    #pragma omp section
    //vz = mesh->g_12*dpdy - mesh->g_22*dpdx;
    vz =0;
  }


  // Upwind A using these velocities
  
  Field3D ry, rz;
  #pragma omp parallel sections
  {
    #pragma omp section
    result = VDDX(vx, A);
    
    #pragma omp section
    ry = VDDY(vy, A);
    
    #pragma omp section
    rz = VDDZ(vz, A);
  }
  //output.write("mesh->g_22: %g  \n" ,vx[4][4]);
  //result = (ry + rz); //usually  convention
  result = (result + ry) / (mesh->J*sqrt(mesh->g_22));

#ifdef TRACK
  result.name = "b0xGrad_dot_Grad("+phi.name+","+A.name+")";
#endif
#ifdef CHECK
  msg_stack.pop(msg_pos);
#endif
  return result;
}



/* computes Jv, where ddt() is holding v and the system state holds Jv */ 
int jacobian(BoutReal t) {
  mesh->communicate(ddt(u),ddt(n));
  
  static Field2D A = 0.0;
  static Field2D C = 1e-12;
  static Field2D D = 1.0;
  
  ddt(phi) = invert_laplace(ddt(u), phi_flags,&A,&C,&D);
  //ddt(phi) = invert_laplace(ddt(u), phi_flags); 

  mesh->communicate(ddt(phi));

  u=0;
  n=0;

  //u -= mybracket(ddt(phi),ddt(u));
  u += bracket3D(ddt(phi),ddt(u));
  //ddt(u) += alpha * phi;
  u += nu * Delp2(ddt(u));
  //ddt(u) -= beta * DDY(n)/n; 
  //ddt(u) -= beta* Grad_par(n)/n; 
  //u -= Grad_par(ddt(n)); 
  //ddt(u).applyBoundary("dirichlet");
  u -= beta* DDZ(n); 
  //mesh->communicate(comms); no don't do this here
  //.applyBoundary();
  //brkt = VDDY(DDY(phi), n) +  VDDZ(DDZ(phi), n) ;
 
  n += bracket3D(ddt(phi),ddt(n));
  //n -= mybracket(ddt(phi),ddt(n));
  n += mu * Delp2(ddt(n));
  //n -= alpha* n;
  
  n.applyBoundary();
  u.applyBoundary();
  return 0;
}

/* computes P^-1 r, where ddt() holds (-r) and the system state hold P^-1 r*/

int precon(BoutReal t, BoutReal gamma, BoutReal delta) {
  // mesh->communicate(rhscomms);
  mesh->communicate(ddt(n),ddt(u));

  n= 0;
  u =0;

  n += ddt(n);
  // mesh->communicate(n);
  // Ni -= (mesh->Bxy*mesh->Bxy*ddt(Ni) - ddt(rho))/(mesh->Bxy*mesh->Bxy);

  u += ddt(u);
  u -= gamma * Grad_par(ddt(n)); 
 
return 0;
 
}

double modon(double x, double y,double kappa,double gamma)
{
  double A,B,C,phi_out;//,phi,r,a,c;

  BoutReal xlen = (mesh->GlobalNx)*((mesh->dx)[0][0]);
  xlen = 1.0;
  r = sqrt(pow(x,2.0) + pow(y,2.0))+ 1e-12;
  r = xlen*r;
  // a = .50;
  // c = 10;

  A =c*a/yn(1,kappa);
  B = c*a*(1.0 + pow(kappa/gamma,2.0));
  C = -pow(kappa/gamma,2.0)*(c*a)/jn(1,gamma);
  

  if(r > a)
    phi_out = A*yn(1,kappa*r/a)*(x/r);
  else
    phi_out = B*(x/a) + C*jn(1,gamma * r/a)*(x/r);
  

  return xlen* phi_out;
}

double Newton_root(double gamma_in,double kappa){
  
  double atol = 1.0;
  int iter = 0;
  int max_iter = 300;
  double f;
  double J;
  double gamma_out = gamma_in;

  while (atol > 1e-8 && iter < max_iter)
    {
      // f = (yn(2,kappa)/(kappa*yn(1,kappa)))  + jn(2,gamma)/(gamma*jn(1,gamma));
      // J = (-4.0 + pow(gamma,2)*(1.0 + pow(j0(gamma)/j1(gamma),2.0) )) /pow(3.0,gamma);
      f = (1.0/pow(gamma_out,2.0))*(yn(2,kappa)*gamma_out*jn(1,gamma_out)+kappa*yn(1,kappa)*jn(2,gamma_out));
      J = -(1.0/pow(gamma_out,2.0))*(kappa*jn(3.0,gamma_out)*yn(1,kappa) + gamma_out*jn(2,gamma_out)*yn(2,kappa));
      
      // f = (yn(2,kappa)*gamma*jn(1,gamma)+kappa*yn(1,kappa)*jn(2,gamma));
      // J = .5*kappa*(jn(1,gamma) - jn(3,gamma))*yn(1,kappa) + gamma*jn(0,gamma)*yn(2,kappa);
     

      //cout << iter<<": "<<gamma<<" "<<f/J <<endl;
      atol = fabs(f);
      gamma_out = gamma_out - 1.0*f/J;
      //cout <<gamma<<" "<<atol <<endl;
      iter++;
      
    }
  //output<< (x_out - x_in)/fabs(x_in)<<endl;
  //cout << "atol "<<atol<<" "<<iter<<endl;
  return gamma_out;
}
