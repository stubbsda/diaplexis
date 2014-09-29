#include "synarmosma.h"

#ifndef _vertexh
#define _vertexh

class Vertex {
 // For the event structure, the nodes of spacetime...
 protected:
  int incept;
  int global_dimension;
  bool boundary;
  bool topology_modified;
  bool geometry_modified;
  std::vector<double> entwinement;
  double curvature;
  double obliquity;
  double geometric_deficiency;
  std::vector<int> ubiquity;
  double deficiency;
  double energy;
  std::set<int> neighbours;
  std::set<int> entourage;
  std::set<int> past;
  std::set<int> future;
  Proposition theorem;

 public:
  Vertex();
  Vertex(const Vertex&);
  Vertex& operator =(const Vertex&);
  Vertex(const std::set<int>&);
  virtual ~Vertex();
  void serialize(std::ofstream&) const;
  void deserialize(std::ifstream&);
  void write2screen() const;
  void clear();
  int valence(int) const;
  friend void vertex_difference(int,int,std::vector<double>&);
  friend int operator <(const Vertex&,const Vertex&);
  friend int operator >(const Vertex&,const Vertex&);
  friend class Simplex;
  friend class Sheet;
  friend class Spacetime;
};
#endif

