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
#ifdef UBQ_NTL
  NTL::ZZ ubiquity;
#else
  std::vector<int> ubiquity
#endif
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
  Vertex(const std::set<int>&);
  Vertex& operator =(const Vertex&);
  virtual ~Vertex();
  void serialize(std::ofstream&,int) const;
  void deserialize(std::ifstream&);
  void write2screen() const;
  void clear();
  inline bool active() const; 
  inline bool active(int) const;
  inline void set_active(int);
  inline void set_ubiquity(const std::vector<int>&);
  int valence(int) const;
  friend void vertex_difference(int,int,std::vector<double>&);
  friend class Simplex;
  friend class Sheet;
  friend class Spacetime;
};

inline void Vertex::set_ubiquity(const std::vector<int>& chi)
{
#ifdef UBQ_NTL
  NTL::PrimeSeq pseq;
  long p;
  ubiquity = 1;
  for(int i=0; i<(signed) chi.size(); ++i) {
    p = pseq.next();
    if (chi[i] == 1) ubiquity *= p;
  }
#else
  ubiquity = chi;
#endif
} 

inline void Vertex::set_active(int sheet) 
{
#ifdef UBQ_NTL
  NTL::PrimeSeq pseq;
  long p;
  for(int i=0; i<1+sheet; ++i) {
    p = pseq.next();
  }
  ubiquity = p;
#else
  ubiquity[sheet] = 1;
#endif
}

inline bool Vertex::active(int sheet) const
{
#ifdef UBQ_NTL
  NTL::PrimeSeq pseq;
  long p;
  for(int i=0; i<1+sheet; ++i) {
    p = pseq.next();
  }
  bool output = (NTL::divide(ubiquity,p) == 1) ? true : false;
#else
  bool output = (ubiquity[sheet] == 1) ? true : false;
#endif
  return output;
}

inline bool Vertex::active() const 
{
#ifdef UBQ_NTL
  bool output = (ubiquity == 1) ? false : true;
#else
  bool output = !ghost(ubiquity);
#endif
  return output; 
}
#endif
