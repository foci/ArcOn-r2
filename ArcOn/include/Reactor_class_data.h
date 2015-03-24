/*Set up the reactor class*/

typedef std::vector < PETScWrappers::MPI::BlockVector > SolutionVector;
unsigned int init_flag;

template <int dim>
class arcOn
{
public:
  arcOn(const unsigned int deg, const unsigned int refine, const unsigned int saved_h, const unsigned int saved_p, ParameterHandler& param);
  ~arcOn();

  void run ();
  bool init;
  bool fast_dif;
  bool fast_react;
  bool checksupport;
  bool p_adapt;
  bool h_adapt;
  bool loading_solution;
  bool saving_solution;

private:
  void initialize_massmatrix();
  void create_mesh ();
  void recreate_boundary_data();
  void setup_system();
  void matrixmapper();
  void create_dg_periodicity();
  double newton_root(double gamma, double kappa);
  double modon(double x, double y, double kappa, double gamma, double c, double a);
  double Ullmann_newton_root(double x_in, double y_in, double b, double C, double m);
  double UllmannMap(double x, double Lx, double y, double Ly, double x_sol, double eps, double m );
  void load_initial_conditions(SolutionVector& subdomain_solution, unsigned int init_flag);
  void load_top(SolutionVector& L2_error_interpolant);
  void periodicity_map(SolutionVector& subdomain_solution, double delta_t);
  void static_periodicity_map(SolutionVector& subdomain_solution, double delta_t);
  void setup_RK();

  void assemble_system();
  void assemble_sigma(SolutionVector& subdomain_solution, double current_time, double delta_t);
  void assemble_stiffness(SolutionVector& subdomain_solution, double delta_t);
  void assemble_cont_stiffness(SolutionVector& subdomain_solution, double delta_t);

  void Calculate_MassAction_Explicit( SolutionVector& subdomain_solution, 
				      double delta_t, double current_time,
				      SolutionVector& mass_action_term);  

  void revert_density( SolutionVector& subdomain_solution, 
				      double delta_t, double current_time,
				      SolutionVector& naive_revert_output);  

  void revert_vacuum( SolutionVector& subdomain_solution, 
				      double delta_t, double current_time);

  /* void Calculate_MassAction_Implicit( SolutionVector& solution, double delta_t, SolutionVector& mass_action_term) const;  */
  /* void calculate_mass_action_exact(SolutionVector& init_solution, SolutionVector& solution, double delta_t, SolutionVector& mass_action_term_exact) const; */

  void calc_convdiff(SolutionVector& subdomain_solution,  double delta_t, double current_time); 
  void calc_reaction(SolutionVector& subdomain_solution, double delta_t, double current_time); 

  void slopelimiter(SolutionVector& subdomain_solution, int step);

  void calculate_div_flux(SolutionVector& substep_solution, double delta_t, 
			  double current_time, SolutionVector& div_flux_term); 
  void RK_div_flux_integrator(SolutionVector& subdomain_solution, double delta_t, 
			      double current_time, SolutionVector& div_flux_integrated); 

  void RK_MassAction_integrator(SolutionVector& subdomain_solution, double delta_t, 
				double current_time, SolutionVector& MassAction_integrated); 

  void calc_poisson(SolutionVector& subdomain_solution, double delta_t);

  void calc_poisson_cont(SolutionVector& subdomain_solution, double delta_t);

  void compute_l2_error(SolutionVector& subdomain_solution, double current_time);

  void cmcmv(SolutionVector& subdomain_solution);

  void compute_l2_interpolation(SolutionVector& subdomain_solution);

  void output_results(const unsigned int cycle);

  MPI_Comm                           mpi_communicator;
  const unsigned int n_mpi_processes;
  const unsigned int this_mpi_process;

  unsigned int degree;
  unsigned int refinements;
  unsigned int num_vertex;
  unsigned int num_blocks;

  parallel::distributed::Triangulation<dim>       triangulation;
 
  int  sum_alpha_dofs;
  int  global_sum_alpha_dofs;
  int  sum_vertices;
  int  num_vert;
  int  lvs;
  int  lvsr;
  int  lvsr_size;
  int  saved_h;
  int  saved_p;
  int  glob_dof_ind;

  std::vector<int> local_vertices;

  std::vector< DoFHandler<dim> * >        dof_handler;
  std::vector< FESystem<dim> * >          fe_collection;
  std::vector< QGauss<dim> * >            quadrature_collection;
  std::vector< QGauss<1> * >              line_quadrature_collection;
  std::vector< QGauss<dim-1> * >          face_quadrature_collection;

  /* std::vector< QGaussLobatto<dim> * >            quadrature_collection; */
  /* std::vector< QGaussLobatto<dim-1> * >          face_quadrature_collection;  */

  std::vector< DoFHandler<dim> * >        idof_handler;
  std::vector< FESystem<dim> * >          ife_collection;
  std::vector< QGauss<dim> * >            iquadrature_collection;
  std::vector< QGauss<1> * >              line_iquadrature_collection;
  std::vector< QGauss<dim-1> * >          iface_quadrature_collection;

  /* std::vector< QGaussLobatto<dim> * >            iquadrature_collection;          */
  /* std::vector< QGaussLobatto<dim-1> * >          iface_quadrature_collection;   */

  std::vector< DoFHandler<dim> * >        ldof_handler;
  std::vector< FESystem<dim> * >          lfe_collection;
  std::vector< QGauss<dim> * >            lquadrature_collection;
  std::vector< QGauss<1> * >              line_lquadrature_collection;
  std::vector< QGauss<dim-1> * >          lface_quadrature_collection;

  /* std::vector< QGaussLobatto<dim> * >            lquadrature_collection; */
  /* std::vector< QGaussLobatto<dim-1> * >          lface_quadrature_collection;   */

  //DoFHandler<dim>                        scalar_dof_handler;
  //FE_DGQ<dim>                            scalar_fe;

  std::vector< DoFHandler<dim> * >        tdof_handler;
  std::vector< FESystem<dim> * >          tfe_collection;
  std::vector< QGauss<dim> * >            tquadrature_collection;
  std::vector< QGauss<1> * >              line_tquadrature_collection;
  std::vector< QGauss<dim-1> * >          tface_quadrature_collection;

  /* std::vector< QGaussLobatto<dim> * >            tquadrature_collection; */
  /* std::vector< QGaussLobatto<dim-1> * >          tface_quadrature_collection;   */

  //std::vector< DoFHandler<dim> * >        trapz_dof_handler;
  //std::vector< FESystem<dim> * >          trapz_fe_collection;
  //std::vector< QTrapez<dim> * >           trapz_quad_collection;
  //std::vector< QTrapez<dim-1> * >         trapz_face_quad_collection;

  Vector<double>        cell_L2_error;
  Vector<double>        cell_L2_error_init;

  Vector<double>		  was_refined;	

  Vector<double>                  counter;
  Vector<double>                  counterh;

  std::vector< Vector<double> > l_soln;
  std::vector< Vector<double> > i_soln;
  std::vector< Vector<double> > r_soln;
  std::vector< Vector<double> > ll_soln;

  double weightx;
  double weighty;

  double				  average_entropy_density;
  double				  old_average_entropy_density;
  double				  global_entropic_jump;
  double				  old_global_entropic_jump;

  double                                conserved_mass;
  double                                center_of_mass;
  double				sup_entropystate;

  Vector<double>      max_face_jump;
  double				total_volume;

  std::vector <IndexSet>           locally_owned_dofs;
  std::vector <IndexSet>           locally_relevant_dofs;
  std::vector <IndexSet>           ilocally_owned_dofs;
  std::vector <IndexSet>           ilocally_relevant_dofs;
  std::vector <IndexSet>           llocally_owned_dofs;
  std::vector <IndexSet>           llocally_relevant_dofs;
  std::vector <IndexSet>           tlocally_owned_dofs;
  std::vector <IndexSet>           tlocally_relevant_dofs;

  IndexSet           scalar_locally_owned_dofs;
  IndexSet           scalar_locally_relevant_dofs;


  struct pInfo {
    FullMatrix<double> AlphaMassMatrix;
    FullMatrix<double> SigmaMassMatrix;
    FullMatrix<double> InverseAlphaMassMatrix;
    FullMatrix<double> InverseSigmaMassMatrix;
    bool alphaProject;
    bool sigmaProject;
    unsigned int dofs_per_cell;
    unsigned int alpha_dofs_per_cell;
    unsigned int sigma_dofs_per_cell;
    std::vector<int> alpha_dof_index; 
    std::vector<int> sigma_dof_index[dim];
  };

  struct MapInfo {
    FullMatrix<double> localMassMatrix;
    FullMatrix<double> localInverseMassMatrix;
    bool MapProject;
    unsigned int dofs_per_cell;
  };

  struct MatPatch {
    FullMatrix<double> localM12;
    unsigned int dofs_per_cell;
    std::pair<unsigned int, unsigned int> global_ind;
  };

  typedef std::map<std::string, pInfo > pMap;
  std::vector<pMap> pmap;

  typedef std::map<std::string, MapInfo > MatMap;
  std::vector<MatMap> matmap;

  typedef std::map<std::string, MatPatch > PatchMap;
  std::vector<PatchMap> patchmap;

  bool pInfoExists(const std::string& name, const int index) const;
  const pInfo& pInfoFind(const std::string& name, const int index) const;
  const FullMatrix<double>& findProjection(const std::string& old_name,const std::string& new_name) const;
  std::vector<pInfo> pinfo;

  bool MapInfoExists(const std::string& name, const int index) const;
  const MapInfo& MapInfoFind(const std::string& name, const int index) const;
  const FullMatrix<double>& findMapProjection(const std::string& old_name,const std::string& new_name) const;
  std::vector < std::vector<MapInfo> > mapinfo ;

  bool PatchInfoExists(const std::string& name, const int index) const;
  const MatPatch& PatchInfoFind(const std::string& name, const int index) const;
  const FullMatrix<double>& findPatchProjection(const std::string& old_name,const std::string& new_name) const;
  std::vector<MatPatch> matpatch ;

  std::vector<double>           prev_soln_alpha[alphadim];
  std::vector< Tensor<1,dim> >  prev_soln_sigma[alphadim];
  std::vector<double>           prev_soln_alpha_face[alphadim];
  std::vector< Tensor<1,dim> >  prev_soln_sigma_face[alphadim];
  std::vector<double>           prev_soln_alpha_neigh_face[alphadim];
  std::vector< Tensor<1,dim> >  prev_soln_sigma_neigh_face[alphadim];
  std::vector<double>           iprev_soln_alpha[alphadim];
  std::vector<double>           iprev_soln_alpha_face[alphadim];

  double glob_min_ribbon_density;

  std::vector<double>           soln_alpha[alphadim];
  std::vector< Tensor<1,dim> >  soln_sigma[alphadim];
  std::vector<double>           soln_alpha_face[alphadim];
  std::vector<double>           isoln_alpha_face[alphadim];
  std::vector<double>           soln_alpha_neigh_face[alphadim];


  std::vector<double>           syncopated[alphadim];
  std::vector< Tensor<1,dim> >  syncopated2[alphadim];

  std::map< int, int > Cell_matchX;
  std::vector < std::map< int, std::vector<double> > > soln_drawX;
  std::vector < std::map< int, std::vector< Tensor<1,dim> > > > grad_drawX;

  std::vector <  std::map< int, std::vector<unsigned int> > > periodic_dof_indX;
  std::vector <  std::map< int, FullMatrix<double> > > periodic_basisX;
  std::vector <  std::map< int, std::vector< std::vector< Tensor<1,dim> > > > > periodic_gradX;

  std::map< int, int > Cell_matchY;
  std::vector < std::map< int, std::vector<double> > > soln_drawY;
  std::vector < std::map< int, std::vector< Tensor<1,dim> > > > grad_drawY;

  std::vector <  std::map< int, std::vector<unsigned int> > > periodic_dof_indY;
  std::vector <  std::map< int, FullMatrix<double> > > periodic_basisY;
  std::vector <  std::map< int, std::vector< std::vector< Tensor<1,dim> > > > > periodic_gradY;
  
  typedef std::map< std::pair<std::string,std::string>, FullMatrix<double> > enrichMap;
  enrichMap emap;

  std::vector< FEValuesExtractors::Scalar * > alpha;
  std::vector< FEValuesExtractors::Vector * > sigma;

  //std:: vector< ComponentMask * > alpha_mask;
  ComponentMask  alpha_mask;

  std::vector< FEValuesExtractors::Scalar * > ialpha;
  std::vector< FEValuesExtractors::Vector * > isigma;

  std::vector< FEValuesExtractors::Scalar * > lalpha;
  std::vector< FEValuesExtractors::Vector * > lsigma;

  std::vector< FEValuesExtractors::Scalar * > talpha;
  std::vector< FEValuesExtractors::Vector * > tsigma;


  SolutionVector       subdomain_solution;
  SolutionVector       naive_subdomain_solution;
  SolutionVector       init_solution;
  SolutionVector       exact_solution;
  SolutionVector       solution_diff;
  SolutionVector       initial_condition;
  SolutionVector       interpolate_base;
  SolutionVector       interpolate_active;
  SolutionVector       proj_solution;
  SolutionVector       linterpolate_active;
  SolutionVector       top_interpolant;
  SolutionVector       L2_error_interpolant;
  SolutionVector       naive_L2_error_interpolant;
  SolutionVector       L2_error_method;
  SolutionVector       naive_L2_error_method;
  SolutionVector       L2_interpolate_active;
  SolutionVector       interpolation_error;
  SolutionVector       naive_interpolation_error;

  SolutionVector       subdomain_solution_holder;
  SolutionVector       interpolated_solution_holder;

  SolutionVector       revert_output;
  SolutionVector       naive_revert_output;
  SolutionVector       cont_output1;
  SolutionVector       cont_global;

  SolutionVector div_flux_integrated;
  SolutionVector naive_div_flux_integrated;
  SolutionVector MassAction_integrated;
  SolutionVector naive_MassAction_integrated;

  //PETScWrappers::MPI::SparseMatrix temp_poisson_matrix;
  //PETScWrappers::MPI::Vector temp_poisson_rhs;
  //PETScWrappers::MPI::Vector redistributed_solution;

  std::map<int, double> newmax;
  std::map<int, double> newmin;
  std::map<int, double> oldval;
  
  std::map<int, double> newmax0;
  std::map<int, double> newmin0;
  std::map<int, double> oldval0;
  
  std::map<int, double> newmax1;
  std::map<int, double> newmin1;
  std::map<int, double> oldval1;
    
  PETScWrappers::MPI::BlockSparseMatrix poisson_matrix;
  PETScWrappers::MPI::BlockVector poisson_rhs; 
  PETScWrappers::PreconditionBoomerAMG preconditioner;

  PETScWrappers::MPI::BlockSparseMatrix cont_poisson_matrix;
  PETScWrappers::MPI::BlockVector cont_poisson_rhs; 
  PETScWrappers::PreconditionBoomerAMG cont_preconditioner;

  //PETScWrappers::MPI::Vector poisson_rhs;

  //PETScWrappers::MPI::Vector maxf;
  //PETScWrappers::MPI::Vector minf;
  PETScWrappers::MPI::SparseMatrix maxf;
  PETScWrappers::MPI::SparseMatrix minf;

  std::vector < PETScWrappers::BlockVector > local_solution;
  std::vector < PETScWrappers::BlockVector > ilocal_solution;
  std::vector < PETScWrappers::BlockVector > llocal_solution;

  ConstraintMatrix elliptic_constraints;
  std::vector < ConstraintMatrix > parahyp_constraints;

  ConstraintMatrix telliptic_constraints;
  std::vector < ConstraintMatrix > tparahyp_constraints;

  std::vector  < SolutionVector >   RK_solution;
  std::vector  < SolutionVector >   naive_RK_solution;
  std::vector  < SolutionVector >   naive_RK_solution_temp;
  std::vector  < SolutionVector >   RK_div_flux;
  std::vector  < SolutionVector >   naive_RK_div_flux;
  std::vector  < SolutionVector >   naive_RK_div_flux_temp;
  std::vector  < SolutionVector >   RK_MassAction;
  std::vector  < SolutionVector >   naive_RK_MassAction;
  std::vector  < SolutionVector >   naive_RK_MassAction_temp;
  
 private:

  RobinBoundaryValues<dim> RBV[num_boundaries];
  WallBoundaryValues<dim>  WBV[num_boundaries];

  ParameterHandler&	               prm;
  ConditionalOStream                 pcout;

  static const unsigned int  fulldim = (dim+1)*alphadim;

  double dt; // time discretization
  unsigned int nstep;
  unsigned int fstep;
  double fast_dt;
  double hardend;

  double conv_threshold;
  unsigned int n_iters;
  double fast_conv_threshold;
  unsigned int fast_n_iters;

  unsigned int RK_order;
  unsigned int RK_stage;
  unsigned int method;
  unsigned int RKtype;
  double eps_const;
  unsigned int modulus;
  unsigned int output_type;
  double CFL_scaling;
  double Time_ramp;


  FullMatrix<double> RK_alpha;
  FullMatrix<double> RK_beta;
  FullMatrix<double> RK_mu1;
  Vector<double> RK_mu2;

  Vector<double> RKC_T;
  Vector<double> RKC_U;
  Vector<double> RKC_Tprime;
  Vector<double> RKC_Tdprime;
  Vector<double> RKC_b;
  Vector<double> RKC_a;
  Vector<double> RKC_c;

  Vector<double> RKC_mu;
  Vector<double> RKC_tildemu;
  Vector<double> RKC_nu;
  Vector<double> RKC_gamma;

  bool artificial_visc;
  double e1;
  double e1_density, e1_vorticity, s0_density, s0_vorticity, kappa_density, kappa_vorticity;
  double bpen;
  double e0; //1e-2*(cell->diameter())/(max_degree);     
  double notes;
  double s0;
  double o_flow;
  double l2_base;
  double lnSe;

  unsigned int elliptic_type;
  unsigned int solver_type;
  double sigma_prefactor;

  double cell_fraction;
  double Mass_dif;
  double Vort_dif;
  double alpha_parameter;
  double beta_parameter;
  double bias_parameter;
  double temp_parameter;
  Vector<double> difs;

  double refine_above;
  double coarsen_below;
  double enrich_below;
  double unenrich_above;
  double CFL_bound; 
  double gmax_jump; 
  double gvisc_temp1;
  double gvisc_temp2;
  double gvisc_avg;
  double gvisc_median;
  double gart_max;
  std::vector <double> gmax_val_loc;
  double h_min_dist;
  double gtotal_volume;
  double global_Avg;

  double cm;
  double cmv;
  double total_density;

};

void petsc_RK_vec (SolutionVector& source, SolutionVector& destination) 
{

  destination.resize( source.size() );
  
  for(unsigned int s=0; s< source.size(); s++){
    
    destination[s]= source[s];
    
  }
}

