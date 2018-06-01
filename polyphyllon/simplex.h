#include "synarmosma/cell.h"

#ifndef _simplexh
#define _simplexh

namespace DIAPLEXIS {
  class Simplex: public SYNARMOSMA::Cell {
   private:
    using SYNARMOSMA::Cell::initialize;

   protected:
    double volume = 0.0;
    double sq_volume = 0.0;
    double energy = 0.0;
    bool modified = true;
    int incept = -1;
    int parity = 0;
    std::set<int> ubiquity;

    void clear() override;
    void initialize(int,int,const std::set<int>&,int = 0);
    void initialize(const std::set<int>&,const std::set<int>&);
    inline bool active() const {return !ubiquity.empty();};
    inline bool active(int n) const {return (ubiquity.count(n) > 0);};
    inline void set_inactive(int n) {ubiquity.erase(n);};
    inline void set_active(int n) {ubiquity.insert(n);};
    inline void set_ubiquity(const std::set<int>& S) {ubiquity = S;};
    inline void get_ubiquity(std::set<int>& S) const {S = ubiquity;};
    inline int presence() const {return (signed) ubiquity.size();};
    inline int get_parity(int,int) const;
    inline bool timelike() const {return (sq_volume < -std::numeric_limits<double>::epsilon());};
    inline bool spacelike() const {return (sq_volume > std::numeric_limits<double>::epsilon());};
    inline bool lightlike() const {return (!timelike() && !spacelike());};
    int absolute_embedding() const;
    void get_faces(std::vector<Simplex>&) const;
    int serialize(std::ofstream&) const override;
    int deserialize(std::ifstream&) override;
   public:
    Simplex();
    Simplex(int);
    Simplex(int,int,const std::set<int>&,int = 0);
    Simplex(const std::set<int>&,const std::set<int>&);
    Simplex(const Simplex&);
    Simplex& operator =(const Simplex&);
    ~Simplex() override;
    friend Simplex operator ^(const Simplex&,const Simplex&);
    friend std::ostream& operator<< (std::ostream&,const Simplex&);
    friend class Spacetime;
  };

  inline int Simplex::get_parity(int u,int v) const
  {
#ifdef DEBUG
    assert(u != v);
    assert(u >= 0 && v >= 0);
    assert(vertices.size() == 2 && vertices.count(u) == 1 && vertices.count(v) == 1);
#endif
    if (parity == 0) return 0;
    int output = (u < v) ? parity : -parity;
    return output;
  }
}
#endif

