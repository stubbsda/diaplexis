#include "sheet.h"

using namespace DIAPLEXIS;

Sheet::Sheet()
{
  H = new SYNARMOSMA::Homology(SYNARMOSMA::GF2,SYNARMOSMA::NATIVE);
  pi = new SYNARMOSMA::Homotopy;

  clear();

  active = true;
}

Sheet::Sheet(int n,SYNARMOSMA::FIELD f,SYNARMOSMA::METHOD m)
{
  H = new SYNARMOSMA::Homology(f,m);
  pi = new SYNARMOSMA::Homotopy;

  clear();

  active = true;
  index = n;
}

Sheet::Sheet(int n,int p,SYNARMOSMA::FIELD f,SYNARMOSMA::METHOD m)
{
  H = new SYNARMOSMA::Homology(f,m);
  pi = new SYNARMOSMA::Homotopy;

  clear();

  active = true;
  index = n;
  parent = p;
}

Sheet::Sheet(const Sheet& source)
{
  H = new SYNARMOSMA::Homology(SYNARMOSMA::GF2,SYNARMOSMA::NATIVE);
  pi = new SYNARMOSMA::Homotopy;

  index = source.index;
  parent = source.parent;
  nstep = source.nstep;
  active = source.active;
  ops = source.ops;
  vx_delta = source.vx_delta;
  *H = *source.H;
  *pi = *source.pi;
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
  *H = *source.H;
  *pi = *source.pi;
  pseudomanifold = source.pseudomanifold;
  boundary = source.boundary;
  orientable = source.orientable;

  return *this;
}

Sheet::~Sheet()
{
  delete H;
  delete pi;
}

void Sheet::clear()
{
  ops = "";
  active = false;
  nstep = 0;  
  parent = -1;
  index = -1;
  vx_delta.clear();
  H->clear();
  pi->clear();
  pseudomanifold = false;
  boundary = false;
  orientable = false;
}

void Sheet::set_topology(const SYNARMOSMA::Homology* K,const SYNARMOSMA::Homotopy* p,bool pm,bool bd,bool orient)
{
  delete H;
  H = new SYNARMOSMA::Homology(*K);
  delete pi;
  pi = new SYNARMOSMA::Homotopy(*p); 
  pseudomanifold = pm;
  boundary = bd;
  orientable = orient;
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
  H->serialize(s);
  pi->serialize(s);
  s.write((char*)(&pseudomanifold),sizeof(bool));
  s.write((char*)(&boundary),sizeof(bool));
  s.write((char*)(&orientable),sizeof(bool));
}

void Sheet::deserialize(std::ifstream& s) 
{
  int i,n;
  char c;

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
  H->deserialize(s);
  pi->deserialize(s);
  s.read((char*)(&pseudomanifold),sizeof(bool));
  s.read((char*)(&boundary),sizeof(bool));
  s.read((char*)(&orientable),sizeof(bool));
}


