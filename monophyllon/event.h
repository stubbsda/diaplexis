#include "synarmosma/vertex.h"
#include "synarmosma/proposition.h"

#ifndef _eventh
#define _eventh

namespace DIAPLEXIS {
  class Event: public SYNARMOSMA::Vertex {
   // For the event structure, the nodes of spacetime...
   protected:
    bool active;
    bool boundary;
    bool topology_modified;
    bool geometry_modified;
    double entwinement;
    double curvature;
    double obliquity;
    double geometric_deficiency;
    double deficiency;
    SYNARMOSMA::Proposition theorem;

    void clear() override;
   public:
    Event();
    Event(const Event&);
    Event(const std::set<int>&);
    Event& operator =(const Event&);
    ~Event() override;
    int serialize(std::ofstream&) const override; 
    int deserialize(std::ifstream&) override;
    int valence(int) const;
    friend std::ostream& operator <<(std::ostream&,const Event&);
    friend class Spacetime;
  };
}
#endif
