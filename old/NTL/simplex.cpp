#include "simplex.h"

using namespace DIAPLEXIS;

Simplex::Simplex()
{
  set_default_values();
}

Simplex::Simplex(int n,unsigned long colour) : Cell(n)
{
  set_default_values();
  ubiquity = colour;
}

Simplex::Simplex(int v1,int v2,unsigned long colour) : Cell(v1,v2)
{
  // A specialized constructor for 0-simplices and 1-simplices
  set_default_values();
  ubiquity = colour;
}

Simplex::Simplex(const std::set<int>& v,unsigned long colour) : Cell(v)
{
  set_default_values();
  ubiquity = colour;
}

Simplex::Simplex(const std::set<int>& v,const NTL::ZZ locus) : Cell(v)
{
  set_default_values();
  ubiquity = locus;
}

Simplex::Simplex(const Simplex& source) : Cell()
{
  vertices = source.vertices;
  ubiquity = source.ubiquity;
  entourage = source.entourage;
  faces = source.faces;
  energy = source.energy;
  volume = source.volume;
  incept = source.incept;
  sq_volume = source.sq_volume;
  orientation = source.orientation;
  modified = source.modified;
}

Simplex& Simplex::operator =(const Simplex& source)
{
  if (this == &source) return *this;

  vertices = source.vertices;
  ubiquity = source.ubiquity;
  entourage = source.entourage;
  faces = source.faces;
  energy = source.energy;
  volume = source.volume;
  incept = source.incept;
  sq_volume = source.sq_volume;
  orientation = source.orientation;
  modified = source.modified;

  return *this;
}

Simplex::~Simplex()
{

}

void Simplex::initialize(int v1,int v2,unsigned long colour)
{
  clear();
  Cell::initialize(v1,v2);
  ubiquity = colour;
}

void Simplex::initialize(int v1,int v2,const NTL::ZZ locus)
{
  clear();
  Cell::initialize(v1,v2);
  ubiquity = locus;
}

void Simplex::initialize(const std::set<int>& vx,unsigned long colour)
{
  clear();
  Cell::initialize(vx);
  ubiquity = colour;
}

void Simplex::initialize(const std::set<int>& vx,const NTL::ZZ locus)
{
  clear();
  Cell::initialize(vx);
  ubiquity = locus;
}

void Simplex::clear()
{
  Cell::clear();
  set_default_values();
}

void Simplex::set_default_values()
{
  ubiquity = 1;
  energy = 0.0;
  volume = 0.0;
  sq_volume = 0.0;
  incept = -1;
  orientation = SYNARMOSMA::DISPARATE;
  modified = true;
}

void Simplex::serialize(std::ofstream& s) const
{
  serialize(s,0);
}

void Simplex::serialize(std::ofstream& s,int nt) const
{
  int i,n;
  long q;
  std::set<int>::const_iterator it;
  NTL::PrimeSeq pseq;
  const int one = 1;
  const int zero = 0; 

  s.write((char*)(&nt),sizeof(int));
  for(i=0; i<nt; ++i) {
    q = pseq.next();
    if (NTL::divide(ubiquity,q) == 1) {
      s.write((char*)(&one),sizeof(int));
    }
    else {
      s.write((char*)(&zero),sizeof(int));
    }
  }
  n = (signed) vertices.size();
  s.write((char*)(&n),sizeof(int));
  for(it=vertices.begin(); it!=vertices.end(); ++it) {
    n = *it;
    s.write((char*)(&n),sizeof(int));
  }
  s.write((char*)(&volume),sizeof(double));
  s.write((char*)(&sq_volume),sizeof(double));
  s.write((char*)(&orientation),sizeof(int));
  s.write((char*)(&energy),sizeof(double));
  s.write((char*)(&incept),sizeof(int));
}

void Simplex::deserialize(std::ifstream& s)
{
  int i,n,m;
  long q;
  NTL::PrimeSeq pseq;

  clear();

  s.read((char*)(&n),sizeof(int));
  ubiquity = 1;
  for(i=0; i<n; ++i) {
    q = pseq.next();
    s.read((char*)(&m),sizeof(int));
    if (m == 1) ubiquity *= q; 
  }
  s.read((char*)(&n),sizeof(int));
  for(i=0; i<n; ++i) {
    s.read((char*)(&m),sizeof(int));
    vertices.insert(m);
  }
  s.read((char*)(&volume),sizeof(double));
  s.read((char*)(&sq_volume),sizeof(double));
  s.read((char*)(&orientation),sizeof(int));
  s.read((char*)(&energy),sizeof(double));
  s.read((char*)(&incept),sizeof(int));
  Cell::calculate_faces();
}

NTL::ZZ Simplex::get_ubiquity() const
{
  return ubiquity;
}

void Simplex::get_faces(std::vector<Simplex>& F) const
{
  // This method grabs the 1+n faces of this n-simplex...
  int i,j,n = (signed) vertices.size();
  std::set<int>::const_iterator it;
  std::set<int> vx;

  F.clear();

  for(i=0; i<n; ++i) {
    j = -1;
    for(it=vertices.begin(); it!=vertices.end(); ++it) {
      j++;
      if (j == i) continue;
      vx.insert(*it);
    }
    F.push_back(Simplex(vx,ubiquity));
    vx.clear();
  }
}

namespace DIAPLEXIS {
  Simplex operator ^(const Simplex& s1,const Simplex& s2)
  {
    std::set<int>::const_iterator it,jt;
    NTL::ZZ locus;
    std::set<int> vx;

    for(it=s1.vertices.begin(); it!=s1.vertices.end(); ++it) {
      jt = std::find(s2.vertices.begin(),s2.vertices.end(),*it);
      if (jt != s2.vertices.end()) vx.insert(*it);
    }
    locus = NTL::GCD(s1.ubiquity,s2.ubiquity);
    Simplex output(vx,locus);

    return output;
  }

  std::ostream& operator <<(std::ostream& s,const Simplex& S)
  {
    s << S.dimension() << std::endl;
    s << SYNARMOSMA::make_key(S.vertices) << std::endl;
    s << S.ubiquity << std::endl;
    return s;
  }
}

