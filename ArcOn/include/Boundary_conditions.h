/* templates for setting non-periodic boundary conditions explicitly*/

template <int dim>
class RobinBoundaryValues
{
 public:
  RobinBoundaryValues(int boundarycomponent=0) : bc(boundarycomponent) {}

  void value_list (const std::vector<Point<dim> > &points,
		   std::vector<double> &values,
		   const unsigned int component,
		   const double time) const;
  void gradient_list (const std::vector< Point< dim > > &points,
		      std::vector< Tensor< 1, dim > > &gradients,
		      const std::vector< Point<dim> >& normals, 
		      const unsigned int component,
		      const double time) const;	 

 private:
  int bc;

};

template <int dim>
void RobinBoundaryValues<dim>::value_list(const std::vector<Point<dim> > &points,
					  std::vector<double> &values,
					  const unsigned int component,
					  const double /*time*/) const
{
  Assert(values.size()==points.size(),
	 ExcDimensionMismatch(values.size(),points.size()));

  for (unsigned int i=0; i<values.size(); ++i)
    {
      if ( ((int)component)==(bc-1)){
	values[i]=.4;
      } else {
	values[i]= ((int)component==4) ? .6 : .1;
      }
    }
}


template <int dim>
void RobinBoundaryValues<dim>::gradient_list (const std::vector< Point< dim > > &  points,
					      std::vector< Tensor< 1, dim > > &  gradients,
					      const std::vector< Point<dim> >& normals, 
					      const unsigned int component,
					      const double /*time*/) const {
  Assert(gradients.size()==points.size(),
	 ExcDimensionMismatch(gradients.size(),points.size()));

  for (unsigned int i=0; i<gradients.size(); ++i)
    {
      Tensor<1,dim> result;
      if ( ((int)component)==(bc-1)){
	result = normals[i];
      } else {
	result = 0;
      }
      gradients[i] = result;
    }
}

template <int dim>
class WallBoundaryValues
{
 public:
  WallBoundaryValues(int boundarycomponent=0) : bc(boundarycomponent) {}

  void value_list (const std::vector<Point<dim> > &points,
		   std::vector<double> &values,
		   const unsigned int component,
		   const double time) const;
  void value_list2 (const std::vector<Point<dim> > &points,
		   std::vector<double> &values,
		   const unsigned int component,
		   const double time) const;
  void gradient_list (const std::vector< Point< dim > > &points,
		      std::vector< Tensor< 1, dim > > &gradients,
		      const std::vector< Point<dim> >& normals, 
		      const unsigned int component,
		      const double time) const;	 
  void gradient_list2 (const std::vector< Point< dim > > &points,
		      std::vector< Tensor< 1, dim > > &gradients,
		      const std::vector< Point<dim> >& normals, 
		      const unsigned int component,
		       const double time) const;

 private:
  int bc;

};

template <int dim>
void WallBoundaryValues<dim>::value_list(const std::vector<Point<dim> > &points,
					 std::vector<double> &values,
					 const unsigned int component,
					 const double time ) const
{
  Assert(values.size()==points.size(),
	 ExcDimensionMismatch(values.size(),points.size()));

  for (unsigned int i=0; i<values.size(); ++i)
    {
      /*
	if ( ((int)component)==(bc-1)){
	values[i]=0;
	} else {
	values[i]= ((int)component==4) ? 0.5 : 0;
	}
      */

      (void) time;

      if(component == 0 ){ values[i] = values[i];}; //0.0;}; //values[i]*0.9999;};(0.1) * ( values[i] - points[i](0) ) ;};
      if(component == 1 ){ values[i] = values[i];}; //values[i];};
      if(component == 2 ){ values[i] = values[i];};  // Force dirichlet vanishing?
      
      /* leave values[i] as passed in */
    }
}
  

template <int dim>
void WallBoundaryValues<dim>::value_list2(const std::vector<Point<dim> > &points,
					 std::vector<double> &values,
					 const unsigned int component,
					 const double rightside ) const
{
  Assert(values.size()==points.size(),
	 ExcDimensionMismatch(values.size(),points.size()));

  for (unsigned int i=0; i<values.size(); ++i)
    {
      /*
	if ( ((int)component)==(bc-1)){
	values[i]=0;
	} else {
	values[i]= ((int)component==4) ? 0.5 : 0;
	}
      */

      (void) time;

      if(component == 0 ){ values[i] = 1.1*rightside;}; //-3.0*1e-4*time;}; //values[i];};//0.0;}; //values[i]*0.9999;};(0.1) * ( values[i] - points[i](0) ) ;};
      if(component == 1 ){ values[i] = 0.0;}; //values[i];};
      if(component == 2 ){ values[i] = values[i];};  // Force dirichlet vanishing?
      
      /* leave values[i] as passed in */
    }
}
  


template <int dim>
void WallBoundaryValues<dim>::gradient_list (const std::vector< Point< dim > > &  points,
					     std::vector< Tensor< 1, dim > > &  gradients,
					     const std::vector< Point<dim> >& /* normals */, 
					     const unsigned int component,
					     const double /*time*/) const {
  Assert(gradients.size()==points.size(),
	 ExcDimensionMismatch(gradients.size(),points.size()));

  //Tensor<1,dim> result(gradients.size());
  std::vector< Tensor< 1, dim > > result(gradients.size());

  for (unsigned int i=0; i<gradients.size(); ++i)
    {
      for (unsigned int j=0; j<dim; ++j){

	if (component == 0){ result[i][j] = 0.0;}; //gradients[i][j];}; //gradients[i][j];}; //do nothing
	if (component == 1){ result[i][j] = 0.0;}; //
	if (component == 2){ result[i][j] = 0.0;}; 

	/*       if ( ((int)component)==(bc-1)){ */
	/* 	result = 0.0; //gradients[i]; //- normals[i]; */
	/*       } else { */
	/* 	result = 0.0; //gradients[i]; //reflective */
	/*       } */
	gradients[i][j] = result[i][j];

      }
    }
}



template <int dim>
void WallBoundaryValues<dim>::gradient_list2 (const std::vector< Point< dim > > &  points,
					     std::vector< Tensor< 1, dim > > &  gradients,
					     const std::vector< Point<dim> >& /* normals */, 
					     const unsigned int component,
					     const double /*time*/) const {
  Assert(gradients.size()==points.size(),
	 ExcDimensionMismatch(gradients.size(),points.size()));

  //Tensor<1,dim> result(gradients.size());
  std::vector< Tensor< 1, dim > > result(gradients.size());

  for (unsigned int i=0; i<gradients.size(); ++i)
    {
      for (unsigned int j=0; j<dim; ++j){

	if (component == 0){ result[i][j] = 0.0;}; //gradients[i][j];}; //do nothing
	if (component == 1){ result[i][j] = 0.0;}; // 
	if (component == 2){ result[i][j] = gradients[i][j];}; 


	/*       if ( ((int)component)==(bc-1)){ */
	/* 	result = 0.0; //gradients[i]; //- normals[i]; */
	/*       } else { */
	/* 	result = 0.0; //gradients[i]; //reflective */
	/*       } */
	gradients[i][j] = result[i][j];
      }
    }
}
