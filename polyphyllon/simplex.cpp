#include "simplex.h"

using namespace DIAPLEXIS;

Simplex::Simplex() 
{
  clear();
}

Simplex::Simplex(int n)
{
  clear();

  for(int i=0; i<=n; ++i) {
    vertices.insert(i);
  }
  calculate_faces();  
}

Simplex::Simplex(int v1,int v2,const std::set<int>& locus,int d)
{
  // A specialized constructor for 1-simplices
  clear();

  vertices.insert(v1);
  faces.push_back(vertices);
  vertices.clear();
  vertices.insert(v2);
  faces.push_back(vertices);
  vertices.insert(v1);
  ubiquity = locus;
  if (d == SYNARMOSMA::UNDIRECTED) return;
  orientation = (v1 < v2) ? d : -d;
}

Simplex::Simplex(const std::set<int>& v,const std::set<int>& locus)
{
  clear();

  vertices = v;
  calculate_faces();
  ubiquity = locus;
}

Simplex::Simplex(const Simplex& source)
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

void Simplex::initialize(int v1,int v2,const std::set<int>& locus,int d)
{
  clear();
  SYNARMOSMA::Cell::initialize(v1,v2);
  set_ubiquity(locus);
  if (d == SYNARMOSMA::UNDIRECTED) return;
  orientation = (v1 < v2) ? d : -d;
}

void Simplex::initialize(const std::set<int>& vx,const std::set<int>& locus)
{
  clear();
  SYNARMOSMA::Cell::initialize(vx);
  set_ubiquity(locus);
}

void Simplex::clear()
{
  SYNARMOSMA::Cell::clear();

  ubiquity.clear();
  energy = 0.0;
  volume = 0.0;
  sq_volume = 0.0;
  incept = -1;
  orientation = SYNARMOSMA::UNDIRECTED;
  modified = true;
}

int Simplex::serialize(std::ofstream& s) const
{
  int n,count = 0;
  std::set<int>::const_iterator it;

  n = (signed) ubiquity.size();
  s.write((char*)(&n),sizeof(int)); count += sizeof(int);
  for(it=ubiquity.begin(); it!=ubiquity.end(); ++it) {
    n = *it;
    s.write((char*)(&n),sizeof(int)); count += sizeof(int);
  }
  n = (signed) vertices.size();
  s.write((char*)(&n),sizeof(int)); count += sizeof(int);
  for(it=vertices.begin(); it!=vertices.end(); ++it) {
    n = *it;
    s.write((char*)(&n),sizeof(int)); count += sizeof(int);
  }
  s.write((char*)(&volume),sizeof(double)); count += sizeof(double);
  s.write((char*)(&sq_volume),sizeof(double)); count += sizeof(double);
  s.write((char*)(&orientation),sizeof(int)); count += sizeof(int);
  s.write((char*)(&energy),sizeof(double)); count += sizeof(double);
  s.write((char*)(&incept),sizeof(int)); count += sizeof(int);

  return count;
}

int Simplex::deserialize(std::ifstream& s)
{
  int i,n,m,count = 0;

  clear();

  s.read((char*)(&n),sizeof(int)); count += sizeof(int);
  for(i=0; i<n; ++i) {
    s.read((char*)(&m),sizeof(int)); count += sizeof(int);
    ubiquity.insert(m);
  }
  s.read((char*)(&n),sizeof(int)); count += sizeof(int);
  for(i=0; i<n; ++i) {
    s.read((char*)(&m),sizeof(int)); count += sizeof(int);
    vertices.insert(m);
  }
  s.read((char*)(&volume),sizeof(double)); count += sizeof(double);
  s.read((char*)(&sq_volume),sizeof(double)); count += sizeof(double);
  s.read((char*)(&orientation),sizeof(int)); count += sizeof(int);
  s.read((char*)(&energy),sizeof(double)); count += sizeof(double);
  s.read((char*)(&incept),sizeof(int)); count += sizeof(int);
  SYNARMOSMA::Cell::calculate_faces();

  return count;
}

void Simplex::get_faces(std::vector<Simplex>& F) const
{
  // This method grabs the 1+n faces of this n-simplex...
  int i,j,n = (signed) vertices.size();
  std::set<int>::const_iterator it;
  std::set<int> vx;
  Simplex S;

  F.clear();

  S.ubiquity = ubiquity;

  for(i=0; i<n; ++i) {
    j = -1;
    for(it=vertices.begin(); it!=vertices.end(); ++it) {
      j++;
      if (j == i) continue;
      vx.insert(*it);
    }
    S.vertices = vx;
    S.calculate_faces();
    F.push_back(S);
    vx.clear();
  }
}

namespace DIAPLEXIS {
  Simplex operator ^(const Simplex& s1,const Simplex& s2)
  {
    std::set<int> vx,locus;
    std::set<int>::const_iterator it;

    for(it=s1.vertices.begin(); it!=s1.vertices.end(); ++it) {
      if (s2.vertices.count(*it) > 0) vx.insert(*it);
    }
    
    for(it=s1.ubiquity.begin(); it!=s1.ubiquity.end(); ++it) {
      if (s2.ubiquity.count(*it) > 0) locus.insert(*it);
    }

    Simplex output(vx,locus);
    return output;
  }

  std::ostream& operator <<(std::ostream& s,const Simplex& S)
  {
    std::set<int>::const_iterator it;

    s << S.dimension() << std::endl;
    s << SYNARMOSMA::make_key(S.vertices) << std::endl;
    s << "[ ";
    for(it=S.ubiquity.begin(); it!=S.ubiquity.end(); ++it) {
      s << *it << " ";
    }
    s << "]" << std::endl;
    return s;
  }
}
