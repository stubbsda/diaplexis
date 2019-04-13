#include "synarmosma/cell.h"

#ifndef _simplexh
#define _simplexh

namespace DIAPLEXIS {
  /// A class representing the d-simplices (where \f$d\ge 1\f$) of the spacetime complex.

  /// This class is derived from the Synarmosma library's Cell class, 
  /// whose documentation should therefore also be consulted for more 
  /// details. 
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
    inline int get_parity(int,int) const {return parity;};
    inline bool get_modified() const {return modified;};
    inline void set_modified(bool t) {modified = t;};
    inline int get_incept() const {return incept;};
    inline void set_incept(int n) {incept = n;};
    inline double get_energy() const {return energy;};
    inline double get_volume() const {return volume;};
    inline void set_squared_volume(double V) {sq_volume = V; volume = std::sqrt(std::abs(V));};
    inline double get_squared_volume() const {return sq_volume;};
    inline bool timelike() const {return (sq_volume < -std::numeric_limits<double>::epsilon());};
    inline bool spacelike() const {return (sq_volume > std::numeric_limits<double>::epsilon());};
    inline bool lightlike() const {return (!timelike() && !spacelike());};
    friend Simplex operator ^(const Simplex&,const Simplex&);
    friend std::ostream& operator<< (std::ostream&,const Simplex&);
    friend class Complex;
  };
}
#endif

