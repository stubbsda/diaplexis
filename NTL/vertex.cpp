#include "vertex.h"

Vertex::Vertex()
{
  clear();
}

Vertex::Vertex(const Vertex& source)
{
  incept = source.incept;
  global_dimension = source.global_dimension;
  deficiency = source.deficiency;
  ubiquity = source.ubiquity;
  entourage = source.entourage;
  energy = source.energy;
  neighbours = source.neighbours;
  past = source.past;
  future = source.future;
  theorem = source.theorem;
  entwinement = source.entwinement;
  curvature = source.curvature;
  obliquity = source.obliquity;
  geometric_deficiency = source.geometric_deficiency;
  boundary = source.boundary;
  topology_modified = source.topology_modified;
  geometry_modified = source.geometry_modified;
}

Vertex& Vertex::operator =(const Vertex& source)
{
  if (this == &source) return *this;

  incept = source.incept;
  global_dimension = source.global_dimension;
  deficiency = source.deficiency;
  ubiquity = source.ubiquity;
  entourage = source.entourage;
  energy = source.energy;
  neighbours = source.neighbours;
  past = source.past;
  future = source.future;
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

Vertex::~Vertex()
{

}

void Vertex::clear()
{
  ubiquity = 1;
  boundary = false;
  neighbours.clear();
  entourage.clear();
  past.clear();
  future.clear();
  theorem.clear();
  incept = -1;
  global_dimension = -1;
  topology_modified = true;
  geometry_modified = true;
  deficiency = 0.0;
  energy = 0.0;
  curvature = 0.0;
  obliquity = 0.0;
  geometric_deficiency = 0.0;
  entwinement.clear();
}

void Vertex::write2screen() const
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

void Vertex::serialize(std::ofstream& s,int nt) const
{
  int i,in1;
  long q;
  double l;
  std::set<int>::const_iterator it;
  NTL::PrimeSeq pseq;
  const int one = 1;
  const int zero = 0;

  s.write((char*)(&incept),sizeof(int));
  s.write((char*)(&global_dimension),sizeof(int));
  s.write((char*)(&deficiency),sizeof(double));
  s.write((char*)(&energy),sizeof(double));
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

  in1 = (signed) neighbours.size();
  s.write((char*)(&in1),sizeof(int));
  for(it=neighbours.begin(); it!=neighbours.end(); ++it) {
    in1 = *it;
    s.write((char*)(&in1),sizeof(int));
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

  in1 = (signed) past.size();
  s.write((char*)(&in1),sizeof(int));
  for(it=past.begin(); it!=past.end(); ++it) {
    in1 = *it;
    s.write((char*)(&in1),sizeof(int));
  }

  in1 = (signed) future.size();
  s.write((char*)(&in1),sizeof(int));
  for(it=future.begin(); it!=future.end(); ++it) {
    in1 = *it;
    s.write((char*)(&in1),sizeof(int));
  }
  theorem.serialize(s);
}

void Vertex::deserialize(std::ifstream& s)
{
  int i,j,in1;
  long q;
  double xc;
  NTL::PrimeSeq pseq;

  clear();

  s.read((char*)(&incept),sizeof(int));
  s.read((char*)(&global_dimension),sizeof(int));
  s.read((char*)(&deficiency),sizeof(double));
  s.read((char*)(&energy),sizeof(double));
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
  for(i=0; i<j; ++i) {
    s.read((char*)(&in1),sizeof(int));
    neighbours.insert(in1);
  }

  s.read((char*)(&j),sizeof(int));
  ubiquity = 1;
  for(i=0; i<j; ++i) {
    q = pseq.next();
    s.read((char*)(&in1),sizeof(int));
    if (in1 == 1) ubiquity *= q;    
  }

  s.read((char*)(&j),sizeof(int));
  for(i=0; i<j; ++i) {
    s.read((char*)(&in1),sizeof(int));
    past.insert(in1);
  }

  s.read((char*)(&j),sizeof(int));
  for(i=0; i<j; ++i) {
    s.read((char*)(&in1),sizeof(int));
    future.insert(in1);
  }
  theorem.deserialize(s);
}

