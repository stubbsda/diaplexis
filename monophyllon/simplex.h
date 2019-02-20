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
    bool active = true;
    int incept = -1;
    int parity = 0;

    void clear() override;
    void initialize(int,int,int = 0);
    void initialize(const std::set<int>&) override;
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
    Simplex(const Simplex&);
    Simplex(int);
    Simplex(int,int,int = 0);
    Simplex(const std::set<int>&);
    ~Simplex() override;
    Simplex& operator =(const Simplex&);
    friend Simplex operator ^(const Simplex&,const Simplex&);
    friend std::ostream& operator<< (std::ostream&,const Simplex&);
    friend class Complex;
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

