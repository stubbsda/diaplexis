#include "nexus.h"

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
  std::vector<Group> HZ;
  Group pi1;

  void clear();
  void serialize(std::ofstream&) const;
  void deserialize(std::ifstream&);
 public:
  Sheet();
  Sheet(int,unsigned long);
  Sheet(int,int,unsigned long);
  Sheet(const Sheet&);
  Sheet& operator =(const Sheet&);
  ~Sheet();
  friend class Spacetime;
};
#endif
