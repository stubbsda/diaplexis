#include "synarmosma/cell.h"

#ifndef _simplexh
#define _simplexh

namespace DIAPLEXIS {
  class Simplex: public SYNARMOSMA::Cell {
   protected:
    double volume;
    double sq_volume;
    double energy;
    bool modified;
    bool active;
    int incept;
    int orientation;

    void clear() override;
    void initialize(int,int,int = SYNARMOSMA::UNDIRECTED);
    void initialize(const std::set<int>&) override;
    inline int get_orientation(int,int) const;
    int absolute_embedding() const;
    void get_faces(std::vector<Simplex>&) const;
    int serialize(std::ofstream&) const override;
    int deserialize(std::ifstream&) override;
   public:
    Simplex();
    Simplex(const Simplex&);
    Simplex(int);
    Simplex(int,int,int = SYNARMOSMA::UNDIRECTED);
    Simplex(const std::set<int>&);
    ~Simplex() override;
    Simplex& operator =(const Simplex&);
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

