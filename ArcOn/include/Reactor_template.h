/* Set up the reactor template */
template <int dim>
arcOn<dim>::arcOn(const unsigned int deg, const unsigned int refine,
		  const unsigned int saved_h, const unsigned int saved_p,
		  ParameterHandler& param)
  :
  mpi_communicator (MPI_COMM_WORLD),
  n_mpi_processes (Utilities::MPI::n_mpi_processes(mpi_communicator)),
  this_mpi_process (Utilities::MPI::this_mpi_process(mpi_communicator)),
  degree ( deg ),
  refinements( refine ),
  triangulation(mpi_communicator, Triangulation<dim>::limit_level_difference_at_vertices),
  saved_h(saved_h),
  saved_p(saved_p),
  prm(param),
  pcout (std::cout,
	 (Utilities::MPI::this_mpi_process(mpi_communicator)
	  == 0))
{

  unsigned int qud_pts = (unsigned int)std::ceil(3.0*(double)deg/2.0 - 1.0/2.0); // For artificial diff might need to go to 3p/2+1/2 

  for (unsigned int component=0; component< alphadim; ++component){
    dof_handler.push_back (new DoFHandler<dim>(triangulation));
    alpha.push_back (new FEValuesExtractors::Scalar(0));
    sigma.push_back (new FEValuesExtractors::Vector(1));

    quadrature_collection.push_back ( new QGauss<dim>(qud_pts));
    face_quadrature_collection.push_back ( new QGauss<dim-1>(qud_pts));
    line_quadrature_collection.push_back ( new QGauss<1>(qud_pts));

    fe_collection.push_back ( new FESystem<dim>(FE_DGQArbitraryNodes<dim>( *(line_quadrature_collection[component]) ), 1,
    						FE_DGQArbitraryNodes<dim>( *(line_quadrature_collection[component]) ), dim) );

    idof_handler.push_back (new DoFHandler<dim>(triangulation));
    ialpha.push_back (new FEValuesExtractors::Scalar(0));
    isigma.push_back (new FEValuesExtractors::Vector(1));

    iquadrature_collection.push_back ( new QGauss<dim>(qud_pts-1));
    line_iquadrature_collection.push_back ( new QGauss<1>(qud_pts-1));
    iface_quadrature_collection.push_back ( new QGauss<dim-1>(qud_pts-1));

    ife_collection.push_back ( new FESystem<dim>(FE_DGQArbitraryNodes<dim>( *line_iquadrature_collection[component] ), 1,
    						 FE_DGQArbitraryNodes<dim>( *line_iquadrature_collection[component] ), dim));

    ldof_handler.push_back (new DoFHandler<dim>(triangulation));
    lalpha.push_back (new FEValuesExtractors::Scalar(0));
    lsigma.push_back (new FEValuesExtractors::Vector(1));

    lquadrature_collection.push_back ( new QGauss<dim>(qud_pts));
    line_lquadrature_collection.push_back ( new QGauss<1>(qud_pts));
    lface_quadrature_collection.push_back ( new QGauss<dim-1>(qud_pts));

    lfe_collection.push_back ( new FESystem<dim>(FE_DGQArbitraryNodes<dim>( *line_lquadrature_collection[component] ), 1,
    						 FE_DGQArbitraryNodes<dim>( *line_lquadrature_collection[component] ), dim));

    unsigned int top_deg = deg;
    unsigned int tqud_pts = (unsigned int)std::ceil(3.0*(double)top_deg/2.0);

    tdof_handler.push_back (new DoFHandler<dim>(triangulation));
    talpha.push_back (new FEValuesExtractors::Scalar(0));
    tsigma.push_back (new FEValuesExtractors::Vector(1));

    tquadrature_collection.push_back ( new QGauss<dim>(qud_pts));
    line_tquadrature_collection.push_back ( new QGauss<1>(qud_pts));
    tface_quadrature_collection.push_back ( new QGauss<dim-1>(qud_pts));
 
    tfe_collection.push_back ( new FESystem<dim>(FE_Q<dim>(top_deg), 1, 
						 FE_Q<dim>(top_deg), dim)); 

   }

  for(unsigned int j=0; j<num_boundaries;j++){
    RBV[j] = RobinBoundaryValues<dim>(j);
    WBV[j] = WallBoundaryValues<dim>(j);
  }
  
  prm.enter_subsection("Physics");
  Mass_dif  = prm.get_double("Mass diffusion");
  Vort_dif  = prm.get_double("Vorticity diffusion");
  alpha_parameter  = prm.get_double("alpha");
  beta_parameter  = prm.get_double("beta");
  bias_parameter  = -prm.get_double("bias");
  prm.leave_subsection ();

  prm.enter_subsection ("Time");
  dt  = prm.get_double("Timestep size dt");
  nstep = prm.get_integer("Number of timesteps");
  hardend = prm.get_double("Hard endtime");
  fstep = prm.get_integer("fstep");
  RK_order = prm.get_integer("RK order");
  RK_stage = prm.get_integer("RK stage");
  method = prm.get_integer("Method");
  RKtype = prm.get_integer("RK type");
  eps_const = prm.get_double("eps");
  modulus = prm.get_integer("Output_modulus");
  output_type = prm.get_integer("Output type");
  CFL_scaling = prm.get_double("CFL scaling");
  Time_ramp = prm.get_double("Ramp");
  prm.leave_subsection ();

  prm.enter_subsection("Regularity");
  artificial_visc  = prm.get_bool("Artificial diffusion");
  e1_density  = prm.get_double("epsilon weight density");
  e1_vorticity  = prm.get_double("epsilon weight vorticity");
  s0_density =  prm.get_double("s0 for density");
  s0_vorticity =  prm.get_double("s0 for vorticity");
  kappa_density =   prm.get_double("kappa for density");
  kappa_vorticity =   prm.get_double("kappa for vorticity");
  bpen  = prm.get_double("Brezzi_penalty");
  prm.leave_subsection ();

  fast_dt = dt/fstep;
  if(fstep == 0)fast_dt = 0.0;

  prm.enter_subsection ("Convergence");
  conv_threshold = prm.get_double("Threshold");
  n_iters = prm.get_integer("Iterations");
  fast_conv_threshold = prm.get_double("Fast Threshold");
  fast_n_iters = prm.get_integer("Fast Iterations");
  prm.leave_subsection ();

  prm.enter_subsection ("EllipticSolver");
  solver_type = prm.get_integer("Continuity type");
  elliptic_type = prm.get_integer("Penalty type");
  sigma_prefactor  = prm.get_double("Sigma penalty");
  prm.leave_subsection ();

  prm.enter_subsection ("SavedSolution");
  loading_solution = prm.get_bool("Load");
  saving_solution  = prm.get_bool("Save");
  prm.leave_subsection ();

  checksupport = false; // Set for speed

  cell_fraction = 1 / (std::pow(2.0,(double)dim)-1);

  pcout << "\033[1;32m\n ____________________________________________Solver information____________________________________________" << std::endl;
  pcout << "\033[1;37m\n\tdt = " << dt << ", Number of timesteps = " << nstep << ", Convergence threshold = " << conv_threshold << ", Iterations allowed = " << n_iters << std::endl;
  pcout << "\n\tfast_dt = " << fast_dt << ", fstep = " << fstep << ", FAST Convergence threshold = " << fast_conv_threshold << ", Iterations allowed = " << fast_n_iters << std::endl;
  pcout << "\n\tRunge-Kutta order = " << RK_order << ", Runge-Kutta stage = " << RK_stage << std::endl;

  setup_RK();


}

template <int dim>
arcOn<dim>::~arcOn(){
  
  for (unsigned int component=0; component< alphadim; ++component){
      dof_handler[component]->clear();
      idof_handler[component]->clear();
      ldof_handler[component]->clear();
  }
  
  pcout << "\033[0;37m\n  Closing and cleaning up ...." << std::endl;

}
