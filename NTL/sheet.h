#include "synarmosma.h"

#ifndef _sheeth
#define _sheeth

namespace DIAPLEXIS {
  class Sheet {
   private:
    int index;
    int parent;
    unsigned long colour;
    std::set<int> vx_delta;
    int nstep;
    bool active;
    bool pseudomanifold;
    bool boundary;
    bool orientable;
    std::string ops;
    SYNARMOSMA::Homology* H;
    SYNARMOSMA::Homotopy* pi;

    void clear();
    void serialize(std::ofstream&) const;
    void deserialize(std::ifstream&);
   public:
    Sheet();
    Sheet(int,unsigned long,SYNARMOSMA::FIELD,SYNARMOSMA::METHOD);
    Sheet(int,int,unsigned long,SYNARMOSMA::FIELD,SYNARMOSMA::METHOD);
    Sheet(const Sheet&);
    Sheet& operator =(const Sheet&);
    ~Sheet();
    friend class Spacetime;
  };
}
#endif
