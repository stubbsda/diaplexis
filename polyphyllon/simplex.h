#include "synarmosma/cell.h"

#ifndef _simplexh
#define _simplexh

namespace DIAPLEXIS {
  class Simplex: public SYNARMOSMA::Cell {
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
    int absolute_embedding() const;
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
    inline bool active() const {return !ubiquity.empty();};
    inline bool active(int n) const {return (ubiquity.count(n) > 0);};
    inline void deactivate() {ubiquity.clear(); modified = true;};
    inline void deactivate(int n) {ubiquity.erase(n); modified = true;};
    inline void activate(int n) {ubiquity.insert(n); modified = true;};
    inline void set_ubiquity(const std::set<int>& S) {ubiquity = S; modified = true;};
    inline void get_ubiquity(std::set<int>& S) const {S = ubiquity;};
    inline int presence() const {return (signed) ubiquity.size();};
    inline void set_entourage(const std::set<int>& N) {entourage = N;};
    inline void clear_entourage() {entourage.clear();};
    inline int get_parity(int,int) const;
    inline bool get_modified() const {return modified;};
    inline void set_modified(bool t) {modified = t;};
    inline int get_incept() const {return incept;};
    inline void set_incept(int n) {incept = n;};
    inline double get_energy() const {return energy;};
    inline void set_volume(double V) {volume = V;};
    inline double get_volume() const {return volume;};
    inline void set_squared_volume(double V) {sq_volume = V;};
    inline double get_squared_volume() const {return sq_volume;};
    inline bool timelike() const {return (sq_volume < -std::numeric_limits<double>::epsilon());};
    inline bool spacelike() const {return (sq_volume > std::numeric_limits<double>::epsilon());};
    inline bool lightlike() const {return (!timelike() && !spacelike());};
    friend Simplex operator ^(const Simplex&,const Simplex&);
    friend std::ostream& operator<< (std::ostream&,const Simplex&);
    friend class Complex;
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

