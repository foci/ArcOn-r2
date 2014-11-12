// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


// check FETools::interpolate on parallel vector

#include "tests.h"
#include <deal.II/base/utilities.h>
#include <deal.II/base/index_set.h>
#include <deal.II/lac/parallel_vector.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/base/function.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/petsc_parallel_block_vector.h>
#include <deal.II/lac/constraint_matrix.h> 
#include <deal.II/fe/fe_system.h>



void test ()
{

  MPI_Comm mpi_communicator (MPI_COMM_WORLD);

  const unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  const unsigned int n_processes = Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD);

  const unsigned int alphadim = 2;
  const unsigned int dim = 2;
  const unsigned int num_blocks = 2;
  const unsigned int deg = 2;

  parallel::distributed::Triangulation<dim> triangulation(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(3);

  typedef std::vector < PETScWrappers::MPI::BlockVector > SolutionVector;

  SolutionVector       subdomain_solution;
  std::vector< DoFHandler<dim> * >        dof_handler;
  std::vector< FESystem<dim> * >          fe_collection;
  std::vector< QGauss<dim> * >            quadrature_collection;
  std::vector< QGauss<1> * >              line_quadrature_collection;
  std::vector< QGauss<dim-1> * >          face_quadrature_collection;

  std::vector <IndexSet>           locally_owned_dofs;
  std::vector <IndexSet>           locally_relevant_dofs;

  std::vector< FEValuesExtractors::Scalar * > alpha;
  std::vector< FEValuesExtractors::Vector * > sigma;

  unsigned int qud_pts = (unsigned int)std::ceil(3.0*(double)deg/2.0 - 1.0/2.0); 

  for (unsigned int component=0; component< alphadim; ++component){
    dof_handler.push_back (new DoFHandler<dim>(triangulation));
    alpha.push_back (new FEValuesExtractors::Scalar(0));
    sigma.push_back (new FEValuesExtractors::Vector(1));

    quadrature_collection.push_back ( new QGauss<dim>(qud_pts));
    face_quadrature_collection.push_back ( new QGauss<dim-1>(qud_pts));
    line_quadrature_collection.push_back ( new QGauss<1>(qud_pts));

    fe_collection.push_back ( new FESystem<dim>(FE_DGQArbitraryNodes<dim>( *(line_quadrature_collection[component]) ), 1,
    						FE_DGQArbitraryNodes<dim>( *(line_quadrature_collection[component]) ), dim) );
   
   }


  subdomain_solution.resize(alphadim);

  
  for (unsigned int component=0; component< alphadim; ++component){
 
    /* Do this for each component */
 
    std::vector<unsigned int> fesystem_sub_blocks (dim+1,0);
    fesystem_sub_blocks[dim-1] = 1;
    fesystem_sub_blocks[dim]   = 1;

    dof_handler[component]->distribute_dofs(*(fe_collection[component]));
    DoFRenumbering::component_wise (*(dof_handler[component]), fesystem_sub_blocks );

    std::vector<unsigned int> fesystem_dofs_per_block (num_blocks);
    DoFTools::count_dofs_per_block (*(dof_handler[component]), fesystem_dofs_per_block,
				    fesystem_sub_blocks);
    
    const unsigned int n_alpha = fesystem_dofs_per_block[0],
      n_sigma = fesystem_dofs_per_block[1];

    std::vector<IndexSet> fesystem_partitioning, fesystem_relevant_partitioning;

    locally_owned_dofs[component] = dof_handler[component]->locally_owned_dofs();

    fesystem_partitioning.push_back(locally_owned_dofs[component]
				    .get_view(0,n_alpha));
    fesystem_partitioning.push_back(locally_owned_dofs[component]
				    .get_view(n_alpha,n_alpha+n_sigma));

    DoFTools::extract_locally_relevant_dofs (*(dof_handler[component]),
					     locally_relevant_dofs[component]);
    
    fesystem_relevant_partitioning
      .push_back(locally_relevant_dofs[component].get_view(0,n_alpha));
    fesystem_relevant_partitioning
      .push_back(locally_relevant_dofs[component].get_view(n_alpha,n_alpha+n_sigma));


    std::vector< unsigned int > num_blocks1;
    std::vector< unsigned int > num_blocks2;
    std::vector< unsigned int > num_blocks3;
    
    num_blocks1.resize(num_blocks,0);
    num_blocks2.resize(num_blocks,0);
    num_blocks3.resize(num_blocks,0);

    subdomain_solution[component].reinit(num_blocks1,mpi_communicator,num_blocks2);

    for (unsigned int bl=0; bl< num_blocks; ++bl){
      
      subdomain_solution[component].block(bl).reinit(mpi_communicator, 
						     fesystem_partitioning[bl], 
						     fesystem_relevant_partitioning[bl] );
    }

    subdomain_solution[component].collect_sizes();
  }
}



int main (int argc, char **argv)
{
  Utilities::System::MPI_InitFinalize mpi_initialization(argc, argv);

  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      std::ofstream logfile("output");
      deallog.attach(logfile);
      deallog << std::setprecision(4);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      test();
    }
  else
    test();

}
