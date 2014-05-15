#include "homology.h"
#include "homotopy.h"

#ifndef _sheeth
#define _sheeth

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
  Homology* H;
  Homotopy* pi;

  void clear();
  void serialize(std::ofstream&) const;
  void deserialize(std::ifstream&);
 public:
  Sheet();
  Sheet(int,unsigned long,FIELD,METHOD);
  Sheet(int,int,unsigned long,FIELD,METHOD);
  Sheet(const Sheet&);
  Sheet& operator =(const Sheet&);
  ~Sheet();
  friend class Spacetime;
};
#endif
