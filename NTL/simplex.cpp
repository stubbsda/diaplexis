#include "simplex.h"

Simplex::Simplex()
{
  ubiquity = 1;
  energy = 0.0;
  volume = 0.0;
  sq_volume = 0.0;
  incept = -1;
  orientation = SPACELIKE;
  modified = true;
}

Simplex::Simplex(int n,unsigned long colour) : Cell(n)
{
  ubiquity = colour;
  energy = 0.0;
  volume = 0.0;
  sq_volume = 0.0;
  incept = -1;
  orientation = SPACELIKE;
  modified = true;
}

Simplex::Simplex(int v1,int v2,unsigned long colour) : Cell(v1,v2)
{
  // A specialized constructor for 0-simplices and 1-simplices
  ubiquity = colour;
  energy = 0.0;
  volume = 0.0;
  sq_volume = 0.0;
  incept = -1;
  orientation = SPACELIKE;
  modified = true;
}

Simplex::Simplex(const std::set<int>& v,unsigned long colour) : Cell(v)
{
  ubiquity = colour;
  energy = 0.0;
  volume = 0.0;
  sq_volume = 0.0;
  incept = -1;
  orientation = SPACELIKE;
  modified = true;
}

Simplex::Simplex(const std::set<int>& v,const NTL::ZZ locale) : Cell(v)
{
  ubiquity = locale;
  energy = 0.0;
  volume = 0.0;
  sq_volume = 0.0;
  incept = -1;
  orientation = SPACELIKE;
  modified = true;
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

void Simplex::initialize(int v1,int v2,const NTL::ZZ locale)
{
  clear();
  Cell::initialize(v1,v2);
  ubiquity = locale;
}

void Simplex::initialize(const std::set<int>& vx,unsigned long colour)
{
  clear();
  Cell::initialize(vx);
  ubiquity = colour;
}

void Simplex::initialize(const std::set<int>& vx,const NTL::ZZ locale)
{
  clear();
  Cell::initialize(vx);
  ubiquity = locale;
}

void Simplex::clear()
{
  Cell::clear();
  ubiquity = 1;
  energy = 0.0;
  volume = 0.0;
  sq_volume = 0.0;
  incept = -1;
  orientation = SPACELIKE;
  modified = true;
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
  s.write((char*)(&orientation),sizeof(CAUSALITY));
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
  s.read((char*)(&orientation),sizeof(CAUSALITY));
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

inline bool operator ==(const Simplex& s1,const Simplex& s2)
{
  if (s1.vertices == s2.vertices) return true;
  return false;  
}

bool operator !=(const Simplex& s1,const Simplex& s2)
{
  if (s1 == s2) return false;
  return true;  
}

bool operator <=(const Simplex& s1,const Simplex& s2)
{
  if (s2.dimension() < s1.dimension()) return false;
  if (s1 == s2) return true;
  bool output = (s1 < s2) ? true : false;
  return output;
}

bool operator <(const Simplex& s1,const Simplex& s2)
{
  if (s2.dimension() <= s1.dimension()) return false;
  // So s2 is bigger than s1, let's see if s1 is contained 
  // in s2
  return std::includes(s2.vertices.begin(),s2.vertices.end(),s1.vertices.begin(),s1.vertices.end());  
}

Simplex operator ^(const Simplex& s1,const Simplex& s2)
{
  std::set<int>::const_iterator it,jt;
  NTL::ZZ locale;
  std::set<int> vx;

  for(it=s1.vertices.begin(); it!=s1.vertices.end(); ++it) {
    jt = std::find(s2.vertices.begin(),s2.vertices.end(),*it);
    if (jt != s2.vertices.end()) vx.insert(*it);
  }
  locale = NTL::GCD(s1.ubiquity,s2.ubiquity);
  Simplex output(vx,locale);

  return output;
}

std::ostream& operator <<(std::ostream& s,const Simplex& S)
{
  s << S.dimension() << std::endl;
  s << make_key(S.vertices) << std::endl;
  s << S.ubiquity << std::endl;
  return s;
}
