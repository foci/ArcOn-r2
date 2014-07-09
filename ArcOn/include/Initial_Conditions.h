/* Set the initial conditions to load */
template <int dim>
class InitialValues : public Function<dim> 
{
public:
  InitialValues () : Function<dim>( alphadim ), initial_n(alphadim) {

    // A parameter handler
    ParameterHandler prm;

    // Declare a section for the function we need
    prm.enter_subsection("Alphas");
    Functions::ParsedFunction<dim>::declare_parameters(prm, alphadim);
    prm.leave_subsection();


    // Parse an input file.
    if (dim==1) prm.read_input("../../input/ICond1.txt");
    if (dim==2) prm.read_input("../../input/ICond2.txt");
    if (dim==3) prm.read_input("../../input/ICond3.txt");

    // Initialize the ParsedFunction object with the given file
    prm.enter_subsection("Alphas");
    initial_n.parse_parameters(prm);
    prm.leave_subsection();

	N1 = 0.5;	
	N2 = 0.5;
	N3 = 0.5;

	U1 = 0.5;
	U2 = - 0.5;
	U3 = 0.5;
  }

  virtual double value (const Point<dim>   &p,
			const unsigned int  component = 0) const;

  virtual void vector_value (const Point<dim> &p, 
			     Vector<double>   &value) const;

  Functions::ParsedFunction<dim> initial_n;

	double N1;
	double N2;
	double N3;
	double U1;
	double U2;
	double U3;
};


template <int dim>
double
InitialValues<dim>::value (const Point<dim>  &p,
			   const unsigned int component) const 
{
    Vector<double> vals(alphadim);
    vector_value(p,vals);
    return vals(component);
}

template <int dim>
void
InitialValues<dim>::vector_value (const Point<dim> &p,
				  Vector<double>   &values) const 
{
    Vector<double> vals(alphadim);
    initial_n.vector_value(p,vals);

    double 	I1 = 0.0;
    double	I2 = 0.0;
    double	I3 = 0.0;
    
	double factor = 1 / (I3*I1 - ( I1*U1*N2*N3 + I2*U2*N1*N3 + I1*I2*U3*N1*N2 ) ) ;
	(void) factor;

	values(0) = vals(0); //I1*N1*N2*N3    * factor; 
	values(1) = vals(1); //I2*N1*N2*N3    * factor;  
	values(2) = vals(2);

	//values(2) = (vals(0)* vals(1)* vals(2))/1000.0; //I1*I2*N1*N2*N3 * factor;
	/* if( (std::abs(values(0)) < 1e-16) || (std::abs(values(1)) < 1e-16) || (std::abs(values(2)) < 1e-16) ){ */
	/* 	std::cout << "WARNING alphas small at point " << p;  */
        /*         std::cout << " I1=" << I1;  */
        /*         std::cout << " I2=" << I2;  */
        /*         std::cout << " I3=" << I3;  */
        /*         std::cout << " alpha=" << values; */
        /*         std::cout << std::endl; */
	/* } */
}
