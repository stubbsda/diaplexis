#include "event.h"

using namespace DIAPLEXIS;

Event::Event()
{
  clear();
}

Event::Event(const Event& source)
{
  incept = source.incept;
  topological_dimension = source.topological_dimension;
  deficiency = source.deficiency;
  ubiquity = source.ubiquity;
  entourage = source.entourage;
  energy = source.energy;
  neighbours = source.neighbours;
  anterior = source.anterior;
  posterior = source.posterior;
  theorem = source.theorem;
  entwinement = source.entwinement;
  curvature = source.curvature;
  obliquity = source.obliquity;
  geometric_deficiency = source.geometric_deficiency;
  boundary = source.boundary;
  topology_modified = source.topology_modified;
  geometry_modified = source.geometry_modified;
}

Event& Event::operator =(const Event& source)
{
  if (this == &source) return *this;

  incept = source.incept;
  topological_dimension = source.topological_dimension;
  deficiency = source.deficiency;
  ubiquity = source.ubiquity;
  entourage = source.entourage;
  energy = source.energy;
  neighbours = source.neighbours;
  anterior = source.anterior;
  posterior = source.posterior;
  theorem = source.theorem;
  entwinement = source.entwinement;
  curvature = source.curvature;
  obliquity = source.obliquity;
  geometric_deficiency = source.geometric_deficiency;
  boundary = source.boundary;
  topology_modified = source.topology_modified;
  geometry_modified = source.geometry_modified;

  return *this;
}

Event::~Event()
{

}

void Event::clear()
{
  SYNARMOSMA::Vertex::clear();

  ubiquity.clear();
  boundary = false;
  topology_modified = true;
  geometry_modified = true;
  deficiency = 0.0;
  curvature = 0.0;
  obliquity = 0.0;
  geometric_deficiency = 0.0;
  entwinement.clear();
  theorem.clear();
}

int Event::serialize(std::ofstream& s) const 
{
  int n,count = 0;
  double l;
  std::set<int>::const_iterator it;

  count += SYNARMOSMA::Vertex::serialize(s);

  s.write((char*)(&deficiency),sizeof(double)); count += sizeof(double);
  s.write((char*)(&curvature),sizeof(double)); count += sizeof(double);
  s.write((char*)(&obliquity),sizeof(double)); count += sizeof(double);
  s.write((char*)(&geometric_deficiency),sizeof(double)); count += sizeof(double);
  s.write((char*)(&boundary),sizeof(bool)); count += sizeof(bool);
  n = (signed) ubiquity.size();
  s.write((char*)(&n),sizeof(int)); count += sizeof(int);
  for(it=ubiquity.begin(); it!=ubiquity.end(); ++it) {
    n = *it;
    s.write((char*)(&n),sizeof(int)); count += sizeof(int);
  }
  count += theorem.serialize(s);
  n = (signed) entwinement.size();
  s.write((char*)(&n),sizeof(int)); count += sizeof(int);
  for(int i=0; i<n; ++i) {
    l = entwinement[i];
    s.write((char*)(&l),sizeof(double)); count += sizeof(double);
  }

  return count;
}

int Event::deserialize(std::ifstream& s)
{  
  int i,n,m,count = 0;
  double xc;

  clear();

  count += SYNARMOSMA::Vertex::deserialize(s);

  s.read((char*)(&deficiency),sizeof(double)); count += sizeof(double);
  s.read((char*)(&curvature),sizeof(double)); count += sizeof(double);
  s.read((char*)(&obliquity),sizeof(double)); count += sizeof(double);
  s.read((char*)(&geometric_deficiency),sizeof(double)); count += sizeof(double);
  s.read((char*)(&boundary),sizeof(bool)); count += sizeof(bool);
  s.read((char*)(&n),sizeof(int)); count += sizeof(int);
  for(i=0; i<n; ++i) {
    s.read((char*)(&m),sizeof(int)); count += sizeof(int);
    ubiquity.insert(m);
  }
  count += theorem.deserialize(s);
  s.read((char*)(&n),sizeof(int)); count += sizeof(int);
  for(i=0; i<n; ++i) {
    s.read((char*)(&xc),sizeof(double)); count += sizeof(double);
    entwinement.push_back(xc);
  }

  return count;
}

namespace DIAPLEXIS {

  std::ostream& operator <<(std::ostream& s,const Event& source)
  {
    std::set<int>::const_iterator it;

    s << source.incept << std::endl;
    s << "[  ";
    for(it=source.ubiquity.begin(); it!=source.ubiquity.end(); ++it) {
      s << *it << "  ";
    }
    s << "]" << std::endl;
    s << source.deficiency << "  " << source.energy << "  " << source.curvature << "  " << source.obliquity << "  " << source.geometric_deficiency << std::endl;
    s << source.theorem;

    return s;
  }
}

