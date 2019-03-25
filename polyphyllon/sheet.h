#include "synarmosma/homology.h"
#include "synarmosma/homotopy.h"

#ifndef _sheeth
#define _sheeth

namespace DIAPLEXIS {
  class Sheet {
   protected:
    int index = -1;
    int parent = -1;
    int nstep = 0;
    bool active = false;
    bool pseudomanifold = false;
    bool boundary = false;
    bool orientable = false;
    std::string hyphantic_ops = "";
    SYNARMOSMA::Homology* H;
    SYNARMOSMA::Homotopy* pi;

    void clear();
    void set_topology(const SYNARMOSMA::Homology*,const SYNARMOSMA::Homotopy*,bool,bool,bool);
    int serialize(std::ofstream&) const;
    int deserialize(std::ifstream&);
   public:
    Sheet();
    Sheet(int,SYNARMOSMA::Homology::Field,SYNARMOSMA::Homology::Method);
    Sheet(int,int,SYNARMOSMA::Homology::Field,SYNARMOSMA::Homology::Method);
    Sheet(const Sheet&);
    Sheet& operator =(const Sheet&);
    ~Sheet();
    friend class Spacetime;
  };
}
#endif
