/* Mapping data for primitive basis approach (this is for optimization only) */
template <int dim>
bool arcOn<dim>::pInfoExists(const std::string& name, const int index) const {
  typename arcOn<dim>::pMap::const_iterator it = pmap[index].find(name);
  return( it != pmap[index].end() );
}

template <int dim>
const typename arcOn<dim>::pInfo& arcOn<dim>::pInfoFind(const std::string& name, const int index) const {
  typename arcOn<dim>::pMap::const_iterator it = pmap[index].find(name);
  if( it != pmap[index].end() ){
  } else {
    std::cout << "pInfo for " << name << " NOT FOUND. Current info for: ";
    for ( it=pmap[index].begin() ; it != pmap[index].end(); it++ ){
      std::cout << it->first << "  ";
    }
    std::cout << std::endl;
    Assert( it != pmap[index].end() ,StandardExceptions::ExcNotInitialized() );
  }
  return it->second;
}

template <int dim>
bool arcOn<dim>::MapInfoExists(const std::string& name, const int index) const {
  typename arcOn<dim>::MatMap::const_iterator it = matmap[index].find(name);
  return( it != matmap[index].end() );
}

template <int dim>
const typename arcOn<dim>::MapInfo& arcOn<dim>::MapInfoFind(const std::string& name, const int index) const {
  typename arcOn<dim>::MatMap::const_iterator it = matmap[index].find(name);
  if( it != matmap[index].end() ){
  } else {
    std::cout << "MatInfo for " << name << " NOT FOUND. Current info for: ";
    for ( it=matmap[index].begin() ; it != matmap[index].end(); it++ ){
      std::cout << it->first << "  ";
    }
    std::cout << std::endl;
    Assert( it != matmap[index].end() ,StandardExceptions::ExcNotInitialized() );
  }
  return it->second;
}

template <int dim>
bool arcOn<dim>::PatchInfoExists(const std::string& name, const int index) const {
  typename arcOn<dim>::PatchMap::const_iterator it = patchmap[index].find(name);
  return( it != patchmap[index].end() );
}

template <int dim>
const typename arcOn<dim>::MatPatch& arcOn<dim>::PatchInfoFind(const std::string& name, const int index) const {
  typename arcOn<dim>::PatchMap::const_iterator it = patchmap[index].find(name);
  if( it != patchmap[index].end() ){
  } else {
    std::cout << "PatchInfo for " << name << " NOT FOUND. Current info for: ";
    for ( it=patchmap[index].begin() ; it != patchmap[index].end(); it++ ){
      std::cout << it->first << "  ";
    }
    std::cout << std::endl;
    Assert( it != patchmap[index].end() ,StandardExceptions::ExcNotInitialized() );
  }
  return it->second;
}

