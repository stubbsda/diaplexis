#include "synarmosma.h"

#ifndef _simplexh
#define _simplexh

namespace DIAPLEXIS {
  enum CAUSALITY 
  {
      PAST,
      FUTURE,
      SPACELIKE
  };

  class Simplex: public SYNARMOSMA::Cell {
   protected:
    double volume;
    double sq_volume;
    double energy;
    bool modified;
    int incept;
    NTL::ZZ ubiquity;
    CAUSALITY orientation;

   public:
    Simplex();
    Simplex(const Simplex&);
    Simplex(int,unsigned long);
    Simplex(int,int,unsigned long);
    Simplex(const std::set<int>&,unsigned long);
    Simplex(const std::set<int>&,const NTL::ZZ);
    virtual ~Simplex();
    Simplex& operator =(const Simplex&);
    void initialize(int,int,unsigned long);
    void initialize(int,int,const NTL::ZZ);
    void initialize(const std::set<int>&,unsigned long);
    void initialize(const std::set<int>&,const NTL::ZZ);
    virtual void clear();
    NTL::ZZ get_ubiquity() const;
    void compute_energy();
    int absolute_embedding() const;
    void get_faces(std::vector<Simplex>&) const;
    void serialize(std::ofstream&,int) const;
    virtual void deserialize(std::ifstream&);
    friend Simplex operator ^(const Simplex&,const Simplex&);
    friend std::ostream& operator<< (std::ostream&,const Simplex&);
    friend class Spacetime;
  };
}
#endif

