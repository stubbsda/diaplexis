#include "synarmosma.h"

#ifndef _simplexh
#define _simplexh

enum CAUSALITY 
{
    PAST,
    FUTURE,
    SPACELIKE
};

class Simplex: public Cell {
 protected:
  double volume;
  double sq_volume;
  double energy;
  bool modified;
  int incept;
#ifdef UBQ_NTL
  NTL::ZZ ubiquity;
#else
  std::vector<int> ubiquity;
#endif
  CAUSALITY orientation;

 public:
  Simplex();
  Simplex(const std::string&,unsigned long);
  Simplex(const std::string&,const NTL::ZZ);
  Simplex(const Simplex&);
  Simplex(int,unsigned long);
  Simplex(int,int,unsigned long);
  Simplex(const std::set<int>&,unsigned long);
  Simplex(const std::set<int>&,const NTL::ZZ);
  Simplex(const std::set<int>&,const std::vector<int>&);
  virtual ~Simplex();
  Simplex& operator =(const Simplex&);
  void initialize(const std::string&,unsigned long);
  void initialize(const std::string&,const NTL::ZZ);
  void initialize(int,int,unsigned long);
  void initialize(int,int,const NTL::ZZ);
  void initialize(const std::set<int>&,unsigned long);
  void initialize(const std::set<int>&,const NTL::ZZ);
  void initialize(const std::set<int>&,const std::vector<int>&);
  virtual void clear();
  inline bool active() const;
  inline bool active(int) const;
  inline void set_active(int);
  inline void set_ubiquity(const std::vector<int>&);
  void compute_energy();
  int absolute_embedding() const;
  void get_faces(std::vector<Simplex>&) const;
  void serialize(std::ofstream&,int) const;
  virtual void deserialize(std::ifstream&);
  friend bool operator ==(const Simplex&,const Simplex&);
  friend bool operator !=(const Simplex&,const Simplex&);
  friend bool operator <=(const Simplex&,const Simplex&);
  friend bool operator <(const Simplex&,const Simplex&);
  friend Simplex operator ^(const Simplex&,const Simplex&);
  friend std::ostream& operator<< (std::ostream&,const Simplex&);
  friend class Vertex;
  friend class Sheet;
  friend class Spacetime;
};

inline void Simplex::set_ubiquity(const std::vector<int>& chi)
{
#ifdef UBQ_NTL
  NTL::PrimeSeq pseq;
  long p;
  ubiquity = 1;
  for(int i=0; i<(signed) chi.size(); ++i) {
    p = pseq.next();
    if (chi[i] == 1) ubiquity *= p;
  }
#else
  ubiquity = chi;
#endif
} 

inline void Simplex::set_active(int sheet)
{
#ifdef UBQ_NTL
  NTL::PrimeSeq pseq;
  long p;
  for(int i=0; i<1+sheet; ++i) {
    p = pseq.next();
  }
  ubiquity = p;
#else
  ubiquity[sheet] = 1;
#endif
}

inline bool Simplex::active(int sheet) const
{
#ifdef UBQ_NTL
  NTL::PrimeSeq pseq;
  long p;
  for(int i=0; i<1+sheet; ++i) {
    p = pseq.next();
  }
  bool output = (NTL::divide(ubiquity,p) == 1) ? true : false;
#else
  bool output = (ubiquity[sheet] == 1) ? true : false;
#endif
  return output;
}

inline bool Simplex::active() const
{
#ifdef UBQ_NTL
  bool output = (ubiquity == 1) ? false : true;
#else
  bool output = !ghost(ubiquity);
#endif
  return output;
}
#endif

