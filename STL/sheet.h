#include "synarmosma.h"

#ifndef _sheeth
#define _sheeth

namespace DIAPLEXIS {
  class Sheet {
   private:
    int index;
    int parent;
    std::set<int> vx_delta;
    int nstep;
    bool active;
    std::string ops;
    SYNARMOSMA::Homology* H;
    SYNARMOSMA::Homotopy* pi;
    bool pseudomanifold;
    bool boundary;
    bool orientable;  

    void clear();
    void set_topology(const SYNARMOSMA::Homology*,const SYNARMOSMA::Homotopy*,bool,bool,bool);
    void serialize(std::ofstream&) const;
    void deserialize(std::ifstream&);
   public:
    Sheet();
    Sheet(int,SYNARMOSMA::FIELD,SYNARMOSMA::METHOD);
    Sheet(int,int,SYNARMOSMA::FIELD,SYNARMOSMA::METHOD);
    Sheet(const Sheet&);
    Sheet& operator =(const Sheet&);
    ~Sheet();
    friend class Spacetime;
  };
}
#endif
