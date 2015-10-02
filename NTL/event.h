#include "synarmosma.h"

#ifndef _eventh
#define _eventh

namespace DIAPLEXIS {
  class Event: public SYNARMOSMA::Vertex {
   // For the event structure, the nodes of spacetime...
   protected:
    bool boundary;
    bool topology_modified;
    bool geometry_modified;
    std::vector<double> entwinement;
    double curvature;
    double obliquity;
    double geometric_deficiency;
    NTL::ZZ ubiquity;
    double deficiency;

   public:
    Event();
    Event(const Event&);
    Event(const std::set<int>&);
    Event& operator =(const Event&);
    virtual ~Event();
    virtual void serialize(std::ofstream&,int) const;
    virtual void serialize(std::ofstream&) const;
    virtual void deserialize(std::ifstream&);
    void write2screen() const;
    virtual void clear();
    int valence(int) const;
    friend class Spacetime;
  };
}
#endif
