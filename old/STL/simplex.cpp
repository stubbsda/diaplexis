#include "simplex.h"

using namespace DIAPLEXIS;

Simplex::Simplex()
{
  set_default_values();
}

Simplex::Simplex(int n,const std::vector<int>& locus) : Cell(n)
{
  set_default_values();
  ubiquity = locus;
}

Simplex::Simplex(int v1,int v2,const std::vector<int>& locus) : Cell(v1,v2)
{
  // A specialized constructor for 0-simplices and 1-simplices
  set_default_values();
  ubiquity = locus;
}

Simplex::Simplex(const std::set<int>& v,const std::vector<int>& locus) : Cell(v)
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

void Simplex::initialize(int v1,int v2,const std::vector<int>& locus)
{
  clear();
  Cell::initialize(v1,v2);
  ubiquity = locus;
}

void Simplex::initialize(const std::set<int>& vx,const std::vector<int>& locus)
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
  ubiquity.clear();
  energy = 0.0;
  volume = 0.0;
  sq_volume = 0.0;
  incept = -1;
  orientation = SYNARMOSMA::DISPARATE;
  modified = true;
}

void Simplex::serialize(std::ofstream& s) const
{
  int i,n;
  std::set<int>::const_iterator it;
 
  n = (signed) ubiquity.size();
  s.write((char*)(&n),sizeof(int));
  for(i=0; i<n; ++i) {
    s.write((char*)(&ubiquity[i]),sizeof(int));
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

  clear();

  s.read((char*)(&n),sizeof(int));
  for(i=0; i<n; ++i) {
    s.read((char*)(&m),sizeof(int));
    ubiquity.push_back(m);
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
    std::vector<int> locus;
    const int n = (signed) s1.ubiquity.size();
    std::set<int> vx;

    for(it=s1.vertices.begin(); it!=s1.vertices.end(); ++it) {
      jt = std::find(s2.vertices.begin(),s2.vertices.end(),*it);
      if (jt != s2.vertices.end()) vx.insert(*it);
    }
    for(int i=0; i<n; ++i) {
      locus.push_back(s1.ubiquity[i] & s2.ubiquity[i]);
    }
    Simplex output(vx,locus);

    return output;
  }

  std::ostream& operator <<(std::ostream& s,const Simplex& S)
  {
    s << S.dimension() << std::endl;
    s << SYNARMOSMA::make_key(S.vertices) << std::endl;
    s << "[";
    for(int i=0; i<(signed) S.ubiquity.size()-1; ++i) {
      s << S.ubiquity[i] << "," << std::endl;
    }
    s << S.ubiquity[S.ubiquity.size()-1] << "]" << std::endl;
    return s;
  }
}

