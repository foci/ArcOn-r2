/* This function is used to parse the components of the vector into strings, for output, etc. */
template <int dim>
struct rmhd
{

  static const unsigned int n_components  = alphadim;

  static
  std::vector<std::string>
  component_names (int k)
  {

    if(k>=0){
      std::vector<std::string> names;
      std::ostringstream parser;
      parser << "alpha_" << (k);
      names.push_back(parser.str());
      for(unsigned int j=0;j<dim;j++){
	std::ostringstream parser;
	parser << "sigma_" << (k) << "_" << (j) << "";
	names.push_back(parser.str());
      }
      return names;
    } 
    if(k==-1) {
      std::vector<std::string> names;
      std::ostringstream parser;
      parser << "max_alpha";
      names.push_back(parser.str());
      for(unsigned int j=0;j<dim;j++){
	std::ostringstream parser;
	parser << "max_sigma_<" << (j) << ">";
	names.push_back(parser.str());
      }
      return names;
    }
    else{
     std::vector<std::string> names;
      std::ostringstream parser;
      parser << "e_alpha_" << (k);
      names.push_back(parser.str());
      for(unsigned int j=0;j<dim;j++){
	std::ostringstream parser;
	parser << "e_sigma_" << (k) << "_<" << (j) << ">";
	names.push_back(parser.str());
      }
      return names;
    }
  }

  static
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
  component_interpretation ()
  {
    std::vector<DataComponentInterpretation::DataComponentInterpretation> data_component_interpretation;
    data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
    for(unsigned int j=0;j<dim;j++){
      data_component_interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
    }
    return data_component_interpretation;
  }
};
