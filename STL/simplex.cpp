#include "simplex.h"

Simplex::Simplex()
{
  energy = 0.0;
  volume = 0.0;
  sq_volume = 0.0;
  incept = -1;
  orientation = SPACELIKE;
  modified = true;
}

Simplex::Simplex(const std::string& vx,const std::vector<int>& locale) : Cell(vx)
{
  // Given a key, create the corresponding simplex...
  energy = 0.0;
  volume = 0.0;
  sq_volume = 0.0;
  incept = -1;
  orientation = SPACELIKE;
  modified = true;
  ubiquity = locale;
}

Simplex::Simplex(int n,const std::vector<int>& locale) : Cell(n)
{
  ubiquity = locale;
  energy = 0.0;
  volume = 0.0;
  sq_volume = 0.0;
  incept = -1;
  orientation = SPACELIKE;
  modified = true;
}

Simplex::Simplex(int v1,int v2,const std::vector<int>& locale) : Cell(v1,v2)
{
  // A specialized constructor for 0-simplices and 1-simplices
  ubiquity = locale;
  energy = 0.0;
  volume = 0.0;
  sq_volume = 0.0;
  incept = -1;
  orientation = SPACELIKE;
  modified = true;
}

Simplex::Simplex(const std::set<int>& v,const std::vector<int>& locale) : Cell(v)
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
  key = source.key;
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
  key = source.key;
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

void Simplex::initialize(int v1,int v2,const std::vector<int>& locale)
{
  std::stringstream s;

  clear();

  if (v1 == -1) {
    vertices.insert(v2);
    s << v2;
  }
  else {
    std::stringstream s1,s2;

    vertices.insert(v1);
    vertices.insert(v2);
    if (v1 < v2) {
      s << v1 << ":" << v2;
    }
    else {
      s << v2 << ":" << v1;
    }
    s1 << v1;
    s2 << v2;
    faces.push_back(s1.str());
    faces.push_back(s2.str());
  }
  key = s.str();
  ubiquity = locale;
}

void Simplex::initialize(const std::string& vx,const std::vector<int>& locale)
{
  clear();
  std::string str;
  boost::char_separator<char> sep(":");
  boost::tokenizer<boost::char_separator<char> > tok(vx,sep);
  for(boost::tokenizer<boost::char_separator<char> >::iterator beg=tok.begin(); beg!=tok.end(); beg++) {
    str = *beg;
    vertices.insert(boost::lexical_cast<int>(str));
  }
  ubiquity = locale;
  string_assembly();
}

void Simplex::initialize(const std::set<int>& vx,const std::vector<int>& locale)
{
  clear();
  vertices = vx;
  ubiquity = locale;
  string_assembly();
}

void Simplex::clear()
{
  vertices.clear();
  entourage.clear();
  faces.clear();
  key = "";
  ubiquity.clear();
  energy = 0.0;
  volume = 0.0;
  sq_volume = 0.0;
  incept = -1;
  orientation = SPACELIKE;
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
  for(it=vertices.begin(); it!=vertices.end(); it++) {
    n = *it;
    s.write((char*)(&n),sizeof(int));
  }
  s.write((char*)(&volume),sizeof(double));
  s.write((char*)(&sq_volume),sizeof(double));
  s.write((char*)(&orientation),sizeof(DIRECTION));
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
  s.read((char*)(&orientation),sizeof(DIRECTION));
  s.read((char*)(&energy),sizeof(double));
  s.read((char*)(&incept),sizeof(int));
  string_assembly();
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
    for(it=vertices.begin(); it!=vertices.end(); it++) {
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
  // in s2, i.e. is s1.key a substring of s2.key?
  std::string::size_type found = s2.key.find(s1.key);
  bool output = (found == std::string::npos) ? false : true;
  return output;
}

Simplex operator ^(const Simplex& s1,const Simplex& s2)
{
  std::set<int>::const_iterator it,jt;
  std::vector<int> locale;
  const int n = (signed) s1.ubiquity.size();
  std::set<int> vx;

  for(it=s1.vertices.begin(); it!=s1.vertices.end(); it++) {
    jt = std::find(s2.vertices.begin(),s2.vertices.end(),*it);
    if (jt != s2.vertices.end()) vx.insert(*it);
  }
  for(int i=0; i<n; ++i) {
    locale.push_back(s1.ubiquity[i] & s2.ubiquity[i]);
  }
  Simplex output(vx,locale);

  return output;
}

std::ostream& operator <<(std::ostream& s,const Simplex& S)
{
  std::set<int>::const_iterator it;
  s << S.dimension() << std::endl;
  s << S.key << std::endl;
  for(it=S.vertices.begin(); it!=S.vertices.end(); it++) {
    s << *it << std::endl;
  }
  return s;
}
