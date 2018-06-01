#include "synarmosma/vertex.h"
#include "synarmosma/proposition.h"

#ifndef _eventh
#define _eventh

namespace DIAPLEXIS {
  class Event: public SYNARMOSMA::Vertex {
   // For the event structure, the nodes of spacetime...
   protected:
    bool boundary = false;
    bool topology_modified = true;
    bool geometry_modified = true;
    double curvature = 0.0;
    double obliquity = 0.0;
    double deficiency = 0.0;
    double geometric_deficiency = 0.0;
    std::set<int> ubiquity;
    std::vector<double> entwinement;
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
    inline bool active() const {return !ubiquity.empty();}; 
    inline bool active(int n) const {return (ubiquity.count(n) > 0);};
    inline void set_active(int n) {ubiquity.insert(n);};
    inline void set_inactive(int n) {ubiquity.erase(n);};
    inline void set_ubiquity(const std::set<int>& S) {ubiquity = S;};
    inline void get_ubiquity(std::set<int>& S) const {S = ubiquity;};
    inline int presence() const {return (signed) ubiquity.size();};
    int valence(int) const;
    friend std::ostream& operator <<(std::ostream&,const Event&);    
    friend class Spacetime;
  };
}
#endif
