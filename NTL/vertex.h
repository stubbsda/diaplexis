#include "synarmosma.h"

#ifndef _vertexh
#define _vertexh

namespace DIAPLEXIS {
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
    NTL::ZZ ubiquity;
    double deficiency;
    double energy;
    std::set<int> neighbours;
    std::set<int> entourage;
    std::set<int> past;
    std::set<int> future;
    SYNARMOSMA::Proposition theorem;

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
    int valence(int) const;
    friend class Spacetime;
  };
}
#endif
