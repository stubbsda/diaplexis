#include "sheet.h"

Sheet::Sheet()
{
  clear();

  active = true;
}

Sheet::Sheet(int n)
{
  clear();

  active = true;
  index = n;
}

Sheet::Sheet(int n,int p)
{
  clear();

  active = true;
  index = n;
  parent = p;
}

Sheet::Sheet(const Sheet& source)
{
  index = source.index;
  parent = source.parent;
  nstep = source.nstep;
  active = source.active;
  ops = source.ops;
  vx_delta = source.vx_delta;
  HZ = source.HZ;
  pi1 = source.pi1;
  pseudomanifold = source.pseudomanifold;
  boundary = source.boundary;
  orientable = source.orientable;
}

Sheet& Sheet::operator =(const Sheet& source)
{
  if (this == &source) return *this;

  index = source.index;
  parent = source.parent;
  active = source.active;
  nstep = source.nstep;
  ops = source.ops;
  vx_delta = source.vx_delta;
  HZ = source.HZ;
  pi1 = source.pi1;
  pseudomanifold = source.pseudomanifold;
  boundary = source.boundary;
  orientable = source.orientable;

  return *this;
}

Sheet::~Sheet()
{

}

void Sheet::clear()
{
  ops = "";
  active = false;
  nstep = 0;  
  parent = -1;
  index = -1;
  vx_delta.clear();
  HZ.clear();
  pi1.clear();
  pseudomanifold = false;
  boundary = false;
  orientable = false;
}

void Sheet::serialize(std::ofstream& s) const
{
  int i,n = (signed) ops.size();

  s.write((char*)(&index),sizeof(int));
  s.write((char*)(&nstep),sizeof(int));
  s.write((char*)(&parent),sizeof(int));
  s.write((char*)(&active),sizeof(bool));
  s.write((char*)(&n),sizeof(int));
  for(i=0; i<n; ++i) {
    s.write((char*)(&ops[i]),sizeof(char));
  }
  // Now the algebraic properties...
  n = (signed) HZ.size();
  s.write((char*)(&n),sizeof(int));
  for(i=0; i<n; ++i) {
    HZ[i].serialize(s);
  }
  pi1.serialize(s);
  s.write((char*)(&pseudomanifold),sizeof(bool));
  s.write((char*)(&boundary),sizeof(bool));
  s.write((char*)(&orientable),sizeof(bool));
}

void Sheet::deserialize(std::ifstream& s) 
{
  int i,n;
  char c;
  Group G;

  clear();

  s.read((char*)(&index),sizeof(int));
  s.read((char*)(&nstep),sizeof(int));
  s.read((char*)(&parent),sizeof(int));
  s.read((char*)(&active),sizeof(bool));
  s.read((char*)(&n),sizeof(int));
  for(i=0; i<n; ++i) {
    s.read((char*)(&c),sizeof(char));
    ops += c;
  }
  // Now the algebraic properties...
  s.read((char*)(&n),sizeof(int));
  for(i=0; i<n; ++i) {
    G.deserialize(s);
    HZ.push_back(G);
  }
  pi1.deserialize(s);
  s.read((char*)(&pseudomanifold),sizeof(bool));
  s.read((char*)(&boundary),sizeof(bool));
  s.read((char*)(&orientable),sizeof(bool));
}


