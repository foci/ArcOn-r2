#include <deal.II/base/table.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parsed_function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/thread_local_storage.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/index_set.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/compressed_simple_sparsity_pattern.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/full_matrix.h>

#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_dgp_monomial.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_bdm.h>
#include <deal.II/fe/fe_dg_vector.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/intergrid_map.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/filtered_iterator.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/fe_field_function.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/fe_field_function.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>

#include <ctime>
#include <cmath>

//Dmitry's modon functions
#include <math.h>
#include <stdio.h>

#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#include <deal.II/lac/petsc_block_sparse_matrix.h>
#include <deal.II/lac/petsc_parallel_block_sparse_matrix.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/parallel_vector.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/compressed_simple_sparsity_pattern.h>
#include <deal.II/lac/block_sparsity_pattern.h>


#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/index_set.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/solution_transfer.h>


inline double elapsed(clock_t start_clock,clock_t end_clock)
{
  return ((double)(end_clock-start_clock))/CLOCKS_PER_SEC;
}

using namespace dealii;

//#include <Multivector.h>

class ParameterReader : public Subscriptor
{
public:
  ParameterReader(ParameterHandler &);
  void read_parameters(const std::string);

private:
  void declare_parameters();
  ParameterHandler &prm;
};

ParameterReader::ParameterReader(ParameterHandler &paramhandler)
  :
  prm(paramhandler)
{}

void ParameterReader::declare_parameters()
{
  prm.enter_subsection ("Mesh");
  {
    prm.declare_entry("Number of spatial refinements", "0",
		      Patterns::Integer(0),
		      "Refinement level from initial macromesh");
    prm.declare_entry("Polynomial Degree", "0",
		      Patterns::Integer(0),
		      "Degree of polynomials to use");
  }
  prm.leave_subsection ();

  prm.enter_subsection ("Physics");
  {
    prm.declare_entry("Mass diffusion", "0",
		      Patterns::Double(0),
		      "Diffusion in density transport");
    prm.declare_entry("Vorticity diffusion", "0",
		      Patterns::Double(0),
		      "Diffusion in vorticity transport");
    prm.declare_entry("alpha", "0",
		      Patterns::Double(0),
		      "alpha parameter");
    prm.declare_entry("beta", "0",
		      Patterns::Double(0),
		      "beta parameter");
    prm.declare_entry("bias", "0",
		      Patterns::Double(0),
		      "bias scaling");
  }
  prm.leave_subsection ();


  prm.enter_subsection ("Time");
  {
    prm.declare_entry("Timestep size dt", "0.1",
		      Patterns::Double(0),
		      "Time step");

    prm.declare_entry("Number of timesteps", "10",
		      Patterns::Integer(0),
		      "Number of (coarse) time steps");

    prm.declare_entry("Hard endtime", "30",
		      Patterns::Double(0),
		      "hard endtime in seconds");

    prm.declare_entry("fstep", "0",
		      Patterns::Integer(0),
		      "Number of (fine) time steps");

    prm.declare_entry("RK order", "2",
		      Patterns::Integer(0),
		      "Runge-Kutta Order");

    prm.declare_entry("RK stage", "2",
		      Patterns::Integer(0),
		      "Number of Runge-Kutta Stages");
    prm.declare_entry("Method", "1",
		      Patterns::Integer(0),
		      "Explicit/Implicit reactions");
    prm.declare_entry("RK type", "1",
		      Patterns::Integer(0),
		      "RKSSP or RKC");
    prm.declare_entry("eps", "0.153846154",
		      Patterns::Double(0),
		      "RKC epsilon");
    prm.declare_entry("Output_modulus", "1",
		      Patterns::Integer(0),
		      "Output solution at every 'modulus'");
    prm.declare_entry("Output type", "0",
		      Patterns::Integer(0),
		      "Output continuous or discontinuous data");
    prm.declare_entry("CFL scaling", "10.0",
                      Patterns::Double(0),
                      "nonlinear CFL constant");
    prm.declare_entry("Ramp", "0.0",
                      Patterns::Double(0),
                      "Time ramp");
  }
  prm.leave_subsection ();

  prm.enter_subsection ("Regularity");
  {
    prm.declare_entry("Artificial diffusion", "false",
		      Patterns::Bool(),
		      "Should artificial diffusion be turned on?");
    prm.declare_entry("epsilon weight density", "0",
		      Patterns::Double(0),
		      "Artificial diffusion = epsilon_weight*local_entropy");
    prm.declare_entry("epsilon weight vorticity", "0",
		      Patterns::Double(0),
		      "Artificial diffusion = epsilon_weight*local_entropy");
    prm.declare_entry("s0 for density", "0",
		      Patterns::Double(0),
		      "Artificial diffusion parameter");
    prm.declare_entry("s0 for vorticity", "0",
		      Patterns::Double(0),
		      "Artificial diffusion parameter");
   prm.declare_entry("kappa for density", "0",
		      Patterns::Double(0),
		      "Artificial diffusion parameter");
   prm.declare_entry("kappa for vorticity", "0",
		      Patterns::Double(0),
		      "Artificial diffusion parameter");
    prm.declare_entry("Brezzi_penalty", "0",
		      Patterns::Double(0),
		      "Coeff for penalty method");
  }
  prm.leave_subsection ();

  prm.enter_subsection ("Convergence");
  {
    prm.declare_entry("Threshold", "1e-12",
		      Patterns::Double(0),
		      "Consider converged if difference has norm below this threshold");
    prm.declare_entry("Iterations", "10",
		      Patterns::Integer(0),
		      "Maximum number of iterations");
    prm.declare_entry("Fast Threshold", "1e-12",
		      Patterns::Double(0),
		      "Consider FAST converged if difference has norm below this threshold");
    prm.declare_entry("Fast Iterations", "10",
		      Patterns::Integer(0),
		      "Maximum number of FAST iterations");
  }
  prm.leave_subsection ();


  prm.enter_subsection ("EllipticSolver");
  {
    prm.declare_entry("Continuity type", "0",
                      Patterns::Integer(0),
                      "Continuous Galerkin (1) or discontinuous Galerkin (0)");
    prm.declare_entry("Penalty type", "1",
                      Patterns::Integer(0),
                      "SIPG = 1, NIPG = -1, IIPG = 0");
    prm.declare_entry("Sigma penalty", "2",
                      Patterns::Double(0),
                      "The penalty for the DG method");

  }
  prm.leave_subsection ();

  prm.enter_subsection ("SavedSolution");
  {
    prm.declare_entry("Load", "true",
		      Patterns::Bool(),
		      "Load Saved 'Exact' Solution");
    prm.declare_entry("Save", "false",
		      Patterns::Bool(),
		      "Save Final solution as 'Exact' Solution");
    prm.declare_entry("h", "5",
		      Patterns::Integer(0),
		      "Saved refinement");
    prm.declare_entry("p", "3",
		      Patterns::Integer(0),
		      "Saved degree");
  }
  prm.leave_subsection ();

  /* 
     prm.enter_subsection ("Output parameters");
     {
     prm.declare_entry("Output file", "solution",
     Patterns::Anything(),
     "Name of the output file (without extension)");
     DataOutInterface<1>::declare_parameters (prm);
     }
     prm.leave_subsection ();
  */
}

void ParameterReader::read_parameters (const std::string parameter_file)
{
  declare_parameters();
  prm.read_input (parameter_file);
}

const unsigned int num_boundaries = 2;
const unsigned int alphadim = 3;
const unsigned int dimension = 2;
const unsigned int num_reactions = 1;

bool fast = false;

#include <Physical_parameters.h>
#include <MassAction.h>
#include <Equation_of_state.h>
#include <Component_parser.h>
#include <Boundary_conditions.h>
#include <Initial_Conditions.h>
#include <Reactor_class_data.h>
#include <Reactor_template.h>
#include <Matrix_projection.h>
#include <Runge_Kutta_integrator.h>
#include <MassMatrix.h>
#include <MatrixMapper.h>
#include <Recreate_boundary_data.h>
#include <Create_mesh.h>
#include <Setup_system.h>
#include <Create_DG_periodicity.h>
#include <Assemble_system.h>
#include <Load_initial_conditions.h>
#include <Load_top.h>
#include <Periodicity_map.h>
#include <Static_periodicity_map.h>
#include <Calculate_div_flux.h>
#include <RK_div_flux_integrator.h>
#include <Calculate_MassAction.h>
#include <Revert_density.h>
#include <Revert_vacuum.h>
#include <Slope_limiter.h>
#include <Assemble_stiffness.h>
#include <Assemble_cont_stiffness.h>
#include <Calculate_Poisson.h>
#include <Calculate_Poisson_Continuous.h>
#include <Compute_l2_error.h>
#include <Compute_cm.h>
#include <Compute_l2_interpolation.h>
#include <RK_MassAction_integrator.h>
#include <Split_method.h>
#include <Assemble_sigma.h>
#include <Output_results.h>
#include <Run.h>
#include <Generic_functions.h>
#include <Newton_root.h>
#include <Modon.h>
#include <Ullmann_Newton_root.h>
#include <Ullmann_Map.h>

int main (int argc, char *argv[]) 
{
  try
    {

      //PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);
      deallog.depth_console (0);
      
      ParameterHandler  prm;
      ParameterReader   param(prm);
      param.read_parameters("../../input/Parameters.txt");

      prm.enter_subsection ("Mesh");
      const unsigned int refine  = prm.get_integer("Number of spatial refinements");
      const unsigned int deg = prm.get_integer("Polynomial Degree");
      prm.leave_subsection ();

      prm.enter_subsection ("SavedSolution");
      const unsigned int saved_p  = prm.get_integer("p");
      const unsigned int saved_h  = prm.get_integer("h");
      prm.leave_subsection ();

      //Assert( (min_refine<=n_refinements) && (n_refinements<=max_refine) , ExcInternalError() );
      //Assert( (min_degree<=initial_degree) && (initial_degree<=max_degree) , ExcInternalError() );

      arcOn<dimension> run_arcOn(deg,refine,saved_h,saved_p,prm);
      run_arcOn.run ();

      //PetscFinalize();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      std::cerr << "Exception on processing: " << std::endl
		<< exc.what() << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      return 1;
    }
  catch (...) 
    {
      std::cerr << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      std::cerr << "Unknown exception!" << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      return 1;
    };

  exit(0);

  return 0;
}


