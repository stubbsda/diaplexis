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
    inline void deactivate() {ubiquity.clear(); topology_modified = true;};
    inline void activate(int n) {ubiquity.insert(n); topology_modified = true;};
    inline void deactivate(int n) {ubiquity.erase(n); topology_modified = true;};
    inline void set_ubiquity(const std::set<int>& S) {ubiquity = S; topology_modified = true;};
    inline void get_ubiquity(std::set<int>& S) const {S = ubiquity;};
    inline int presence() const {return (signed) ubiquity.size();};
    inline void clear_entourage() {entourage.clear();};
    inline void set_entourage(const std::set<int>& N) {entourage = N;};
    inline void get_entourage(std::set<int> N) const {N = entourage;};
    inline bool drop_entourage(int);
    inline void get_neighbours(std::set<int>& N) const {N = neighbours;};
    inline void set_neighbours(const std::set<int>& N) {neighbours = N; topology_modified = true;};
    inline bool add_neighbour(int);
    inline bool drop_neighbour(int);
    inline double get_deficiency() const {return deficiency;};
    inline void set_deficiency(double x) {deficiency = x;};
    inline void get_entwinement(std::vector<double>& x) const {x = entwinement;};
    inline void set_entwinement(const std::vector<double>& x) {entwinement = x;};
    inline double get_obliquity() const {return obliquity;};
    inline void set_obliquity(double x) {obliquity = x;};
    inline int get_incept() const {return incept;};
    inline void set_incept(int n) {incept = n;};
    inline bool get_boundary() const {return boundary;};
    inline void set_boundary(bool t) {boundary = t;};
    inline bool get_topology_modified() const {return topology_modified;};
    inline void set_topology_modified(bool t) {topology_modified = t;};
    inline bool get_geometry_modified() const {return geometry_modified;};
    inline void set_geometry_modified(bool t) {geometry_modified = t;};
    inline double get_geometric_deficiency() const {return geometric_deficiency;};
    inline void set_geometric_deficiency(double x) {geometric_deficiency = x;};
    inline int get_topological_dimension() const {return topological_dimension;};
    inline void set_topological_dimension(int n) {topological_dimension = n;};
    friend std::ostream& operator <<(std::ostream&,const Event&);    
    friend class Complex;
  };

  bool Event::add_neighbour(int n)
  {
    if (neighbours.count(n) == 0) {
      neighbours.insert(n);
      topology_modified = true;
      return true;
    }
    return false;
  }

  bool Event::drop_neighbour(int n)
  {
    std::set<int>::const_iterator it = std::find(neighbours.begin(),neighbours.end(),n);
    if (it != neighbours.end()) {
      neighbours.erase(it);
      topology_modified = true;
      return true;
    }
    return false;
  }

  bool Event::drop_entourage(int n)
  {
    std::set<int>::const_iterator it = std::find(entourage.begin(),entourage.end(),n);
    if (it != entourage.end()) {
      entourage.erase(it);
      return true;
    }
    return false;
  }
}
#endif
