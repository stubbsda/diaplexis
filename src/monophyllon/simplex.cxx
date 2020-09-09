#include "simplex.h"

using namespace DIAPLEXIS;

Simplex::Simplex()
{

}

Simplex::Simplex(int n)
{
  for(int i=0; i<=n; ++i) {
    vertices.insert(i);
  }
  calculate_faces();
}

Simplex::Simplex(int v1,int v2,int d)
{
#ifdef DEBUG
  assert(d == 0 || std::abs(d) == 1);
#endif
  // A specialized constructor for 1-simplices
  vertices.insert(v1);
  faces.push_back(vertices);
  vertices.clear();
  vertices.insert(v2);
  faces.push_back(vertices);
  vertices.insert(v1);
  parity = d;
}

Simplex::Simplex(const std::set<int>& vx,int n)
{
  vertices = vx;
  incept = n;
  calculate_faces();
}

Simplex::Simplex(const Simplex& source)
{
  vertices = source.vertices;
  active = source.active;
  entourage = source.entourage;
  faces = source.faces;
  energy = source.energy;
  volume = source.volume;
  incept = source.incept;
  sq_volume = source.sq_volume;
  parity = source.parity;
  modified = source.modified;
}

Simplex& Simplex::operator =(const Simplex& source)
{
  if (this == &source) return *this;

  vertices = source.vertices;
  active = source.active;
  entourage = source.entourage;
  faces = source.faces;
  energy = source.energy;
  volume = source.volume;
  incept = source.incept;
  sq_volume = source.sq_volume;
  parity = source.parity;
  modified = source.modified;

  return *this;
}

Simplex::~Simplex()
{

}

void Simplex::clear() 
{
  SYNARMOSMA::Cell::clear();

  active = true;
  energy = 0.0;
  volume = 0.0;
  sq_volume = 0.0;
  incept = -1;
  parity = 0;
  modified = true;
}

void Simplex::initialize(int v1,int v2,int d)
{
#ifdef DEBUG
  assert(d == 0 || std::abs(d) == 1);
#endif

  clear();

  SYNARMOSMA::Cell::initialize(v1,v2);
  active = true;
  parity = d;
}

void Simplex::initialize(const std::set<int>& vx)
{
  clear();

  SYNARMOSMA::Cell::initialize(vx);
  active = true;
}

int Simplex::serialize(std::ofstream& s) const
{
  int n,count = 0;
  std::set<int>::const_iterator it;

  s.write((char*)(&active),sizeof(bool)); count += sizeof(bool);
  n = (signed) vertices.size();
  s.write((char*)(&n),sizeof(int)); count += sizeof(int);
  for(it=vertices.begin(); it!=vertices.end(); ++it) {
    n = *it;
    s.write((char*)(&n),sizeof(int)); count += sizeof(int);
  }
  s.write((char*)(&volume),sizeof(double)); count += sizeof(double);
  s.write((char*)(&sq_volume),sizeof(double)); count += sizeof(double);
  s.write((char*)(&parity),sizeof(int)); count += sizeof(int);
  s.write((char*)(&energy),sizeof(double)); count += sizeof(double);
  s.write((char*)(&incept),sizeof(int)); count += sizeof(int);

  return count;
}

int Simplex::deserialize(std::ifstream& s)
{
  int i,n,m,count = 0;

  clear();

  s.read((char*)(&active),sizeof(bool)); count += sizeof(bool);
  s.read((char*)(&n),sizeof(int)); count += sizeof(int);
  for(i=0; i<n; ++i) {
    s.read((char*)(&m),sizeof(int)); count += sizeof(int);
    vertices.insert(m);
  }
  s.read((char*)(&volume),sizeof(double)); count += sizeof(double);
  s.read((char*)(&sq_volume),sizeof(double)); count += sizeof(double);
  s.read((char*)(&parity),sizeof(int)); count += sizeof(int);
  s.read((char*)(&energy),sizeof(double)); count += sizeof(double);
  s.read((char*)(&incept),sizeof(int)); count += sizeof(int);

  SYNARMOSMA::Cell::calculate_faces();

  return count;
}

namespace DIAPLEXIS {
  Simplex operator ^(const Simplex& s1,const Simplex& s2)
  {
    std::set<int>::const_iterator it;
    std::set<int> vx;

    for(it=s1.vertices.begin(); it!=s1.vertices.end(); ++it) {
      if (s2.vertices.count(*it) > 0) vx.insert(*it);
    }
    Simplex output(vx);
    output.active = s1.active && s2.active;
    output.parity = s1.parity * s2.parity;
    return output;
  }

  std::ostream& operator <<(std::ostream& s,const Simplex& S)
  {
    s << S.incept << "  " << S.parity << "  " << S.active << std::endl;
    s << SYNARMOSMA::make_key(S.vertices) << std::endl;
    s << S.energy << "  " << S.sq_volume;
    return s;
  }
}

