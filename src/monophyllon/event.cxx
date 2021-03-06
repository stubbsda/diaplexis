#include "event.h"

using namespace DIAPLEXIS;

template<class kind>
Event<kind>::Event()
{

}

template<class kind>
Event<kind>::Event(const Event<kind>& source)
{
  incept = source.incept;
  topological_dimension = source.topological_dimension;
  deficiency = source.deficiency;
  active = source.active;
  entourage = source.entourage;
  energy = source.energy;
  neighbours = source.neighbours;
  anterior = source.anterior;
  posterior = source.posterior;
  theorem = source.theorem;
  entwinement = source.entwinement;
  obliquity = source.obliquity;
  geometric_deficiency = source.geometric_deficiency;
  boundary = source.boundary;
  topology_modified = source.topology_modified;
  geometry_modified = source.geometry_modified;
}

template<class kind>
Event<kind>& Event<kind>::operator =(const Event<kind>& source)
{
  if (this == &source) return *this;

  incept = source.incept;
  topological_dimension = source.topological_dimension;
  deficiency = source.deficiency;
  active = source.active;
  entourage = source.entourage;
  energy = source.energy;
  neighbours = source.neighbours;
  anterior = source.anterior;
  posterior = source.posterior;
  theorem = source.theorem;
  entwinement = source.entwinement;
  obliquity = source.obliquity;
  geometric_deficiency = source.geometric_deficiency;
  boundary = source.boundary;
  topology_modified = source.topology_modified;
  geometry_modified = source.geometry_modified;

  return *this;
}

template<class kind>
Event<kind>::~Event()
{

}

template<class kind>
void Event<kind>::clear()
{
  SYNARMOSMA::Vertex<kind>::clear();

  active = true;
  boundary = false;
  topology_modified = true;
  geometry_modified = true;
  deficiency = 0.0;
  obliquity = 0.0;
  geometric_deficiency = 0.0;
  entwinement = 0.0;
  theorem.clear();
}

template<class kind>
int Event<kind>::serialize(std::ofstream& s) const
{
  int count = SYNARMOSMA::Vertex<kind>::serialize(s);

  s.write((char*)(&deficiency),sizeof(double)); count += sizeof(double);
  s.write((char*)(&obliquity),sizeof(double)); count += sizeof(double);
  s.write((char*)(&geometric_deficiency),sizeof(double)); count += sizeof(double);
  s.write((char*)(&entwinement),sizeof(double)); count += sizeof(double);
  s.write((char*)(&boundary),sizeof(bool)); count += sizeof(bool);
  s.write((char*)(&active),sizeof(bool)); count += sizeof(bool);
  count += theorem.serialize(s);

  return count;
}

template<class kind>
int Event<kind>::deserialize(std::ifstream& s)
{
  int count = 0;

  clear();

  count += SYNARMOSMA::Vertex<kind>::deserialize(s);

  s.read((char*)(&deficiency),sizeof(double)); count += sizeof(double);
  s.read((char*)(&obliquity),sizeof(double)); count += sizeof(double);
  s.read((char*)(&geometric_deficiency),sizeof(double)); count += sizeof(double);
  s.read((char*)(&entwinement),sizeof(double)); count += sizeof(double);
  s.read((char*)(&boundary),sizeof(bool)); count += sizeof(bool);
  s.read((char*)(&active),sizeof(bool)); count += sizeof(bool);
  count += theorem.deserialize(s);

  return count;
}

namespace DIAPLEXIS {

  template<class kind>
  std::ostream& operator <<(std::ostream& s,const Event<kind>& source)
  {
    s << source.incept << "  " << source.active << std::endl;
    s << source.deficiency << "  " << source.energy << "  " << source.obliquity << "  " << source.geometric_deficiency << std::endl;
    s << source.theorem;

    return s;
  }
}
