#include "spacetime.cxx"
#include "spacetime_disk.cxx"
#include "spacetime_topology.cxx"
#include "spacetime_geometry.cxx"
#include "spacetime_initialization.cxx"
#include "spacetime_explication.cxx"
#include "spacetime_implication.cxx"
#include "spacetime_hyphansis.cxx"
#include "spacetime_optimization.cxx"
#include "spacetime_visualization.cxx"

namespace DIAPLEXIS
{
  template<class kind1,class kind2>
  const int Spacetime<kind1,kind2>::N_EXP;

  template<class kind1,class kind2>
  const int Spacetime<kind1,kind2>::N_IMP;

  template<class kind1,class kind2>
  const std::string Spacetime<kind1,kind2>::EXP_OP[] = {"D","Ux","Ox","R","C","N","A","G","Sg","Sm","Y"};

  template<class kind1,class kind2>
  const std::string Spacetime<kind1,kind2>::IMP_OP[] = {"I","Um","Om","E","F","P","V","Δ"};
}

template class Spacetime<SYNARMOSMA::UINT64,SYNARMOSMA::INT64>;
template class Spacetime<double,double>;
