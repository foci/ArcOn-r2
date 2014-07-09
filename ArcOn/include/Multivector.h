/* Depreciated at present */
template<unsigned int species, typename NUMBER>
class MultiVector {
public:
  MultiVector() {};
  ~MultiVector() {};

  Vector<NUMBER> parts[species];

  Vector<NUMBER>& operator[](unsigned int n) { return parts[n]; }

  const Vector<NUMBER>& operator[](unsigned int n) const { return parts[n]; }

  void reinit(const unsigned int N, bool fast=false) {
    for(unsigned int k=0; k<species; k++) parts[k].reinit(N,fast);
    }

  MultiVector& operator=(const MultiVector& rhs){
    for(unsigned int k=0; k<species; k++) parts[k] = rhs[k];
    return (*this);
  }

  MultiVector& operator=(const NUMBER& rhs){
    for(unsigned int k=0; k<species; k++) parts[k] = rhs;
    return (*this);
  }

  MultiVector& operator-=(const MultiVector& rhs){
    for(unsigned int k=0; k<species; k++) parts[k] -= rhs[k];
    return (*this);
  }

  MultiVector& operator+=(const MultiVector& rhs){
    for(unsigned int k=0; k<species; k++) parts[k] += rhs[k];
    return (*this);
  }

  void add(const NUMBER a, const MultiVector& V){
    for(unsigned int k=0; k<species; k++) parts[k].add(a,V[k]);
  }

  unsigned int size() const {
    if(species>0){
      return parts[0].size();
    } else {
      return 0;
    }
  }
};
