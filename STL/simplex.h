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
    int incept;
    std::vector<int> ubiquity;
    SYNARMOSMA::RELATION orientation;

     virtual void set_default_values();
   public:
    Simplex();
    Simplex(const Simplex&);
    Simplex(int,const std::vector<int>&);
    Simplex(int,int,const std::vector<int>&);
    Simplex(const std::set<int>&,const std::vector<int>&);
    virtual ~Simplex();
    Simplex& operator =(const Simplex&);
    void initialize(int,int,const std::vector<int>&);
    void initialize(const std::set<int>&,const std::vector<int>&);
    virtual void clear();
    void get_faces(std::vector<Simplex>&) const;
    virtual void serialize(std::ofstream&) const;
    virtual void deserialize(std::ifstream&);
    friend Simplex operator ^(const Simplex&,const Simplex&);
    friend std::ostream& operator<< (std::ostream&,const Simplex&);
    friend class Spacetime;
  };
}
#endif

