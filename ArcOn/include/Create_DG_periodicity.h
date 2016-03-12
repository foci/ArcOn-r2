/* Create_dg_periodicity */
template <int dim>
void arcOn<dim>::create_dg_periodicity()
{ 

  soln_drawX.clear();
  grad_drawX.clear();
  periodic_basisX.clear();
  periodic_gradX.clear();
  periodic_dof_indX.clear();
  Cell_matchX.clear();

  soln_drawX.resize(alphadim);
  grad_drawX.resize(alphadim);
  periodic_basisX.resize(alphadim);
  periodic_gradX.resize(alphadim);
  periodic_dof_indX.resize(alphadim);

  typename DoFHandler<dim>::active_cell_iterator 
    celli = dof_handler[0]->begin_active(), 
    endci = dof_handler[0]->end();

  for (; celli!=endci; ++celli){
    typename DoFHandler<dim>::active_cell_iterator 
      cellj = dof_handler[0]->begin_active(), 
      endcj = dof_handler[0]->end();
    for (; cellj!=endcj; ++cellj){
      if (!(celli->is_artificial())){
	if (!(cellj->is_artificial())){
	  for(unsigned int face_numi=0;
	      face_numi < GeometryInfo<dim>::faces_per_cell; ++face_numi){
	    for(unsigned int face_numj=0;
		face_numj < GeometryInfo<dim>::faces_per_cell; ++face_numj){
	  
	      typename DoFHandler<dim>::face_iterator facei = celli->face(face_numi);
	      typename DoFHandler<dim>::face_iterator facej = cellj->face(face_numj);
	  
	      if ( facei->at_boundary() && facej->at_boundary()
		   && facei->center()[0] != facej->center()[0]
		   && facei->center()[1] == facej->center()[1] 
		   && facei->boundary_indicator() != facej->boundary_indicator())
		{
		  Cell_matchX[celli->index()] = cellj->index();
		}
	    }
	  }
	}
      }
    }
  }

  soln_drawY.clear();
  grad_drawY.clear();
  periodic_basisY.clear();
  periodic_gradY.clear();
  periodic_dof_indY.clear();
  Cell_matchY.clear();

  soln_drawY.resize(alphadim);
  grad_drawY.resize(alphadim);
  periodic_basisY.resize(alphadim);
  periodic_gradY.resize(alphadim);
  periodic_dof_indY.resize(alphadim);

  typename DoFHandler<dim>::active_cell_iterator
    cellk = dof_handler[0]->begin_active(), 
    endck = dof_handler[0]->end();
  
  for (; cellk!=endck; ++cellk){
    typename DoFHandler<dim>::active_cell_iterator
      celll = dof_handler[0]->begin_active(), 
      endcl = dof_handler[0]->end();
    for (; celll!=endcl; ++celll){
      if (!(cellk->is_artificial())){
	if (!(celll->is_artificial())){
	  for(unsigned int face_numk=0;
	      face_numk < GeometryInfo<dim>::faces_per_cell; ++face_numk){
	    for(unsigned int face_numl=0;
		face_numl < GeometryInfo<dim>::faces_per_cell; ++face_numl){
	      
	      typename DoFHandler<dim>::face_iterator facek = cellk->face(face_numk);
	      typename DoFHandler<dim>::face_iterator facel = celll->face(face_numl);
	      
	      if ( facek->at_boundary() && facel->at_boundary()
		   && facek->center()[0] == facel->center()[0]
		   && facek->center()[1] != facel->center()[1] 
		   && facek->boundary_indicator() != facel->boundary_indicator())
		{
		  Cell_matchY[cellk->index()] = celll->index();
		}
	    }
	  }
	}
      }
    }
  }

}
