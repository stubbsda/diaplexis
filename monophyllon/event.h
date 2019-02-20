#include "synarmosma/vertex.h"
#include "synarmosma/proposition.h"

#ifndef _eventh
#define _eventh

namespace DIAPLEXIS {
  class Event: public SYNARMOSMA::Vertex {
   // For the event structure, the nodes of spacetime...
   protected:
    bool active = true;
    bool boundary = false;
    bool topology_modified = true;
    bool geometry_modified = true;
    double entwinement = 0.0;
    double curvature = 0.0;
    double obliquity = 0.0;
    double geometric_deficiency = 0.0;
    double deficiency = 0.0;
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
    friend class Complex;
    friend class Spacetime;
  };
}
#endif
