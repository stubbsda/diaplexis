#include "cell.h"

class Simplex: public Cell {
 private:
  double volume;
  double sq_volume;
  double energy;
  bool modified;
  int incept;
  std::vector<int> ubiquity;
  DIRECTION orientation;

 public:
  Simplex();
  Simplex(const std::string&,const std::vector<int>&);
  Simplex(const Simplex&);
  Simplex(int,const std::vector<int>&);
  Simplex(int,int,const std::vector<int>&);
  Simplex(const std::set<int>&,const std::vector<int>&);
  virtual ~Simplex();
  Simplex& operator =(const Simplex&);
  void initialize(const std::string&,const std::vector<int>&);
  void initialize(int,int,const std::vector<int>&);
  void initialize(const std::set<int>&,const std::vector<int>&);
  virtual void clear();
  void get_faces(std::vector<Simplex>&) const;
  virtual void serialize(std::ofstream&) const;
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
