#include "nexus.h"

#ifndef _sheeth
#define _sheeth

class Sheet {
 private:
  int index;
  int parent;
  std::set<int> vx_delta;
  int nstep;
  bool active;
  std::string ops;
  std::vector<Group> HZ;
  Group pi1;
  bool pseudomanifold;
  bool boundary;
  bool orientable;  

  void clear();
  void serialize(std::ofstream&) const;
  void deserialize(std::ifstream&);
 public:
  Sheet();
  Sheet(int);
  Sheet(int,int);
  Sheet(const Sheet&);
  Sheet& operator =(const Sheet&);
  ~Sheet();
  friend class Spacetime;
};
#endif
