#include "synarmosma/cell.h"

#ifndef _simplexh
#define _simplexh

namespace DIAPLEXIS {
  class Simplex: public SYNARMOSMA::Cell {
   private:
    using SYNARMOSMA::Cell::initialize;

   protected:
    double volume;
    double sq_volume;
    double energy;
    bool modified;
    int incept;
    int orientation;
    std::set<int> ubiquity;

    void clear() override;
    void initialize(int,int,const std::set<int>&,int = SYNARMOSMA::UNDIRECTED);
    void initialize(const std::set<int>&,const std::set<int>&);
    inline bool active() const {return !ubiquity.empty();};
    inline bool active(int n) const {return (ubiquity.count(n) > 0);};
    inline void set_inactive(int n) {ubiquity.erase(n);};
    inline void set_active(int n) {ubiquity.insert(n);};
    inline void set_ubiquity(const std::set<int>& S) {ubiquity = S;};
    inline void get_ubiquity(std::set<int>& S) const {S = ubiquity;};
    inline int presence() const {return (signed) ubiquity.size();};
    inline int get_orientation(int,int) const;
    int absolute_embedding() const;
    void get_faces(std::vector<Simplex>&) const;
    int serialize(std::ofstream&) const override;
    int deserialize(std::ifstream&) override;
   public:
    Simplex();
    Simplex(int);
    Simplex(int,int,const std::set<int>&,int = SYNARMOSMA::UNDIRECTED);
    Simplex(const std::set<int>&,const std::set<int>&);
    Simplex(const Simplex&);
    Simplex& operator =(const Simplex&);
    ~Simplex() override;
    friend Simplex operator ^(const Simplex&,const Simplex&);
    friend std::ostream& operator<< (std::ostream&,const Simplex&);
    friend class Spacetime;
  };

  inline int Simplex::get_orientation(int u,int v) const
  {
#ifdef DEBUG
    assert(u != v);
    assert(u >= 0 && v >= 0);
    assert(vertices.size() == 2 && vertices.count(u) == 1 && vertices.count(v) == 1);
#endif
    if (orientation == SYNARMOSMA::UNDIRECTED) return orientation;
    int output = (u < v) ? orientation : -orientation;
    return output;
  }
}
#endif

