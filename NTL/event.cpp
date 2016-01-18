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
  Vertex::clear();
  ubiquity = 1;
  boundary = false;
  topology_modified = true;
  geometry_modified = true;
  deficiency = 0.0;
  curvature = 0.0;
  obliquity = 0.0;
  geometric_deficiency = 0.0;
  entwinement.clear();
}

void Event::write2screen() const
{
  std::cout << incept << std::endl;
  std::cout << deficiency << std::endl;
  std::cout << energy << std::endl;
  std::cout << curvature << std::endl;
  std::cout << obliquity << std::endl;
  std::cout << geometric_deficiency << std::endl;
  std::cout << ubiquity << std::endl;
  std::cout << theorem << std::endl;
}

void Event::serialize(std::ofstream& s) const
{
  serialize(s,0);
}

void Event::serialize(std::ofstream& s,int nt) const
{
  int i,in1;
  long q;
  double l;
  NTL::PrimeSeq pseq;
  const int one = 1;
  const int zero = 0;

  Vertex::serialize(s);

  s.write((char*)(&deficiency),sizeof(double));
  s.write((char*)(&curvature),sizeof(double));
  s.write((char*)(&obliquity),sizeof(double));
  s.write((char*)(&geometric_deficiency),sizeof(double));
  s.write((char*)(&boundary),sizeof(bool)); 

  in1 = (signed) entwinement.size();
  s.write((char*)(&in1),sizeof(int));
  for(i=0; i<in1; ++i) {
    l = entwinement[i];
    s.write((char*)(&l),sizeof(double));
  }

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
}

void Event::deserialize(std::ifstream& s)
{
  int i,j,in1;
  long q;
  double xc;
  NTL::PrimeSeq pseq;

  entwinement.clear();
  topology_modified = true;
  geometry_modified = true;
  ubiquity = 1;
  Vertex::deserialize(s);

  s.read((char*)(&deficiency),sizeof(double));
  s.read((char*)(&curvature),sizeof(double));
  s.read((char*)(&obliquity),sizeof(double));
  s.read((char*)(&geometric_deficiency),sizeof(double));
  s.read((char*)(&boundary),sizeof(bool));

  s.read((char*)(&j),sizeof(int));
  for(i=0; i<j; ++i) {
    s.read((char*)(&xc),sizeof(double));
    entwinement.push_back(xc);
  }

  s.read((char*)(&j),sizeof(int));
  ubiquity = 1;
  for(i=0; i<j; ++i) {
    q = pseq.next();
    s.read((char*)(&in1),sizeof(int));
    if (in1 == 1) ubiquity *= q;    
  }
}

