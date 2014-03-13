#include <base/table.h>
#include <base/parameter_handler.h>
#include <base/parsed_function.h>
#include <base/quadrature_lib.h>
#include <base/function.h>
#include <base/point.h>
#include <base/tensor.h>
#include <base/thread_local_storage.h>
#include <base/utilities.h>
#include <base/conditional_ostream.h>
#include <base/index_set.h>

#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/constraint_matrix.h>
#include <lac/solver_control.h>
#include <lac/solver_gmres.h>
#include <lac/precondition.h>
#include <lac/solver_cg.h>
#include <lac/compressed_simple_sparsity_pattern.h>
#include <lac/sparsity_tools.h>
#include <lac/full_matrix.h>

#include <lac/sparse_matrix.h>
#include <lac/solver_cg.h>

#include <fe/fe_dgq.h>
#include <fe/fe_dgp.h>
#include <fe/fe_dgp_monomial.h>
#include <fe/fe_q.h>
#include <fe/fe_bdm.h>
#include <fe/fe_dg_vector.h>
#include <fe/fe_system.h>
#include <fe/fe_values.h>
#include <fe/fe_tools.h>
#include <fe/mapping_q1.h>

#include <dofs/dof_handler.h>
#include <dofs/dof_renumbering.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>

#include <grid/tria.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria_boundary_lib.h>
#include <grid/grid_generator.h>
#include <grid/grid_refinement.h>
#include <grid/grid_in.h>
#include <grid/intergrid_map.h>
#include <grid/grid_tools.h>
#include <grid/filtered_iterator.h>

#include <numerics/data_out.h>
#include <numerics/error_estimator.h>
#include <numerics/solution_transfer.h>
#include <numerics/fe_field_function.h>
#include <numerics/matrix_tools.h>
#include <numerics/fe_field_function.h>
#include <numerics/vector_tools.h>

#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>

#include <ctime>
#include <cmath>

//Dmitry's modon functions
#include <math.h>
#include <stdio.h>

#include <lac/petsc_parallel_vector.h>
#include <lac/petsc_block_vector.h>
#include <lac/petsc_vector.h>
#include <lac/petsc_parallel_sparse_matrix.h>
#include <lac/petsc_block_sparse_matrix.h>
#include <lac/petsc_parallel_block_sparse_matrix.h>
#include <lac/petsc_solver.h>
#include <lac/petsc_precondition.h>
#include <lac/parallel_vector.h>
#include <lac/sparsity_tools.h>
#include <lac/compressed_simple_sparsity_pattern.h>
#include <lac/block_sparsity_pattern.h>


#include <base/utilities.h>
#include <base/conditional_ostream.h>
#include <base/index_set.h>
#include <distributed/tria.h>
#include <distributed/grid_refinement.h>
#include <distributed/solution_transfer.h>


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
  }
  prm.leave_subsection ();

  prm.enter_subsection ("Regularity");
  {
    prm.declare_entry("Artificial diffusion", "false",
		      Patterns::Bool(),
		      "Should artificial diffusion be turned on?");
    prm.declare_entry("epsilon weight", "0",
		      Patterns::Double(0),
		      "Artificial diffusion = epsilon_weight*local_entropy");
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
#include <Assemble_system.h>
#include <Load_initial_conditions.h>
#include <Load_top.h>
#include <Periodicity_map.h>
#include <Static_periodicity_map.h>
#include <Calculate_div_flux.h>
#include <RK_div_flux_integrator.h>
#include <Calculate_MassAction.h>
#include <Revert_density.h>
#include <Slope_limiter.h>
#include <Assemble_stiffness.h>
#include <Calculate_Poisson.h>
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


int main (int argc, char *argv[]) 
{
  try
    {

      //PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);
      Utilities::System::MPI_InitFinalize mpi_initialization(argc, argv);
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


