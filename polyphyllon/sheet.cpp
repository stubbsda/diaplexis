#include "sheet.h"

using namespace DIAPLEXIS;

Sheet::Sheet()
{
  H = new SYNARMOSMA::Homology(SYNARMOSMA::Homology::Field::mod2,SYNARMOSMA::Homology::Method::native);
  pi1 = new SYNARMOSMA::Homotopy;

  active = true;
}

Sheet::Sheet(int n,SYNARMOSMA::Homology::Field f,SYNARMOSMA::Homology::Method m)
{
  H = new SYNARMOSMA::Homology(f,m);
  pi1 = new SYNARMOSMA::Homotopy;

  active = true;
  index = n;
}

Sheet::Sheet(int n,int p,SYNARMOSMA::Homology::Field f,SYNARMOSMA::Homology::Method m)
{
  H = new SYNARMOSMA::Homology(f,m);
  pi1 = new SYNARMOSMA::Homotopy;

  active = true;
  index = n;
  parent = p;
}

Sheet::Sheet(const Sheet& source)
{
  index = source.index;
  parent = source.parent;
  active = source.active;
  hyphantic_ops = source.hyphantic_ops;
  H = new SYNARMOSMA::Homology(*source.H);
  pi1 = new SYNARMOSMA::Homotopy(*source.pi1);
  pseudomanifold = source.pseudomanifold;
  boundary = source.boundary;
  orientable = source.orientable;
}

Sheet& Sheet::operator =(const Sheet& source)
{
  if (this == &source) return *this;

  delete H;
  delete pi1;

  index = source.index;
  parent = source.parent;
  active = source.active;
  hyphantic_ops = source.hyphantic_ops;
  H = new SYNARMOSMA::Homology(*source.H);
  pi1 = new SYNARMOSMA::Homotopy(*source.pi1);
  pseudomanifold = source.pseudomanifold;
  boundary = source.boundary;
  orientable = source.orientable;

  return *this;
}

Sheet::~Sheet()
{
  delete H;
  delete pi1;
}

void Sheet::clear()
{
  hyphantic_ops = "";
  active = false;
  parent = -1;
  index = -1;
  H->clear();
  pi1->clear();
  pseudomanifold = false;
  boundary = false;
  orientable = false;  
}

void Sheet::set_topology(const SYNARMOSMA::Homology* K,const SYNARMOSMA::Homotopy* p,bool pm,bool bd,bool orient)
{
  delete H;
  delete pi1;

  H = new SYNARMOSMA::Homology(*K);
  pi1 = new SYNARMOSMA::Homotopy(*p); 
  pseudomanifold = pm;
  boundary = bd;
  orientable = orient;
}

int Sheet::serialize(std::ofstream& s) const
{
  int i,count = 0,n = (signed) hyphantic_ops.size();

  s.write((char*)(&index),sizeof(int)); count += sizeof(int);
  s.write((char*)(&parent),sizeof(int)); count += sizeof(int);
  s.write((char*)(&active),sizeof(bool)); count += sizeof(bool);
  s.write((char*)(&n),sizeof(int)); count += sizeof(int);
  for(i=0; i<n; ++i) {
    s.write((char*)(&hyphantic_ops[i]),sizeof(char)); count += sizeof(char);
  }
  // Now the algebraic properties...
  count += H->serialize(s);
  count += pi1->serialize(s);
  s.write((char*)(&pseudomanifold),sizeof(bool)); count += sizeof(bool);
  s.write((char*)(&boundary),sizeof(bool)); count += sizeof(bool);
  s.write((char*)(&orientable),sizeof(bool)); count += sizeof(bool);

  return count;
}

int Sheet::deserialize(std::ifstream& s) 
{
  int i,n,count = 0;
  char c;

  clear();

  s.read((char*)(&index),sizeof(int)); count += sizeof(int);
  s.read((char*)(&parent),sizeof(int)); count += sizeof(int);
  s.read((char*)(&active),sizeof(bool)); count += sizeof(bool);
  s.read((char*)(&n),sizeof(int)); count += sizeof(int);
  for(i=0; i<n; ++i) {
    s.read((char*)(&c),sizeof(char)); count += sizeof(char);
    hyphantic_ops += c;
  }
  // Now the algebraic properties...
  count += H->deserialize(s);
  count += pi1->deserialize(s);
  s.read((char*)(&pseudomanifold),sizeof(bool)); count += sizeof(bool);
  s.read((char*)(&boundary),sizeof(bool)); count += sizeof(bool);
  s.read((char*)(&orientable),sizeof(bool)); count += sizeof(bool);

  return count;
}


