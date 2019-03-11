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
    inline void activate() {active = true;};
    inline void deactivate() {active = false;};
    inline void clear_entourage() {entourage.clear();};
    inline void set_entourage(const std::set<int>& N) {entourage = N;};
    inline void get_entourage(std::set<int> N) const {N = entourage;};
    inline void get_neighbours(std::set<int>& N) const {N = neighbours;};
    inline void set_neighbours(const std::set<int>& N) {neighbours = N;};
    inline void add_neighbour(int n) {neighbours.insert(n);};
    inline double get_deficiency() const {return deficiency;};
    inline int get_incept() const {return incept;};
    friend std::ostream& operator <<(std::ostream&,const Event&);
    friend class Complex;
  };
}
#endif
