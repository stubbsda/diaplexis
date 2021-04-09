#include "complex.cxx"
#include "complex_dynamics.cxx"
#include "complex_topology.cxx"

namespace DIAPLEXIS {
  template<class kind>
  const int Complex<kind>::ND;
}

template class Complex<SYNARMOSMA::UINT64>;
template class Complex<double>;