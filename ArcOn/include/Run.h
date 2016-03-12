/* This is the base class template */
template <int dim>
void arcOn<dim>::run()
{
  create_mesh();

  setup_system();

  matrixmapper();

  create_dg_periodicity();
 
  InitialValues<dim> IV;

  assemble_system();

}
