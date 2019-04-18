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
    /// This property is the square root of the absolute value of the Simplex::sq_volume 
    /// property and thus non-negative.
    double volume = 0.0;
    /// This property is the square of its volume and may be positive, negative or zero, 
    /// depending on the nature of the spacetime geometry (Euclidean or Lorentzian).
    double sq_volume = 0.0;
    /// This property represents the energy of a simplex and should therefore always be 
    /// non-negative.
    double energy = 0.0;
    /// This property is set to true when the simplex's geometry or topology is modified, 
    /// so that any derived quantities (such as Simplex::sq_volume) can be recomputed.
    bool modified = true;
    /// This property records the relaxation step at which this simplex was originally 
    /// created.
    int incept = -1;
    /// This property stores the parity of this simplex, which can be 0, +1 or -1 and 
    /// is an independent physical property of each 1-simplex. Higher-dimensional simplices 
    /// have a parity value equal to the product of the parity of their edges.
    int parity = 0;
    /// This property is true if this simplex is active and false if the simplex has been 
    /// deleted from the spacetime complex.
    bool active = true;

    /// This method clears all of the instance's extended properties and sets the scalar properties to their default value.
    void clear() override;
    /// This method calls the clear() method and then initializes the instance as a 1-simplex, with the first two arguments the vertices and the optional third argument the value of the Simplex::parity property.
    void initialize(int,int,int = 0);
    /// This method calls the clear() method and then initializes the instance as a simplex whose vertex set is the method's argument.
    void initialize(const std::set<int>&) override;
    /// This method writes the instance properties to a binary disk file and returns the number of bytes written to the file.
    int serialize(std::ofstream&) const override;
    /// This method calls the clear() method on the instance and then reads the properties from a binary disk file and returns the number of bytes read.
    int deserialize(std::ifstream&) override;
   public:
    /// The default constructor which does nothing.
    Simplex();
    /// The standard copy constructor that copies over the properties from the source instance of this class.
    Simplex(const Simplex&);
    /// This constructor that creates an \f$n\f$-simplex with vertices \f$\{0,1,\dots,n\}\f$ where \f$n\f$ is the constructor's argument.
    Simplex(int);
    /// This constructor is specialized for the common case of 1-simplices. The first two arguments are the vertices while the optional third argument is the value of the Simplex::parity property.
    Simplex(int,int,int = 0);
    /// This constructor accepts as its first argument the set of vertices of the simplex while the optional second argument is the value of the Simplex::incept property. 
    Simplex(const std::set<int>&,int = -1);
    /// The destructor which does nothing for this class.
    ~Simplex() override;
    /// The overloaded assignment operator for this class, which behaves exactly like the copy constructor for this class.
    Simplex& operator =(const Simplex&);
    /// This method sets the Simplex::active property to true as well as noting that the topology has been modified.
    inline void activate() {active = true; modified = true;};
    /// This method sets the Simplex::active property to false as well as noting that the topology has been modified.
    inline void deactivate() {active = false; modified = true;};
    /// This method sets the inherited entourage property to the argument.
    inline void set_entourage(const std::set<int>& N) {entourage = N;};
    /// This method clears the inherited entourage property.
    inline void clear_entourage() {entourage.clear();};
    /// This method returns the value of the Simplex::parity property.
    inline int get_parity(int,int) const {return parity;};
    /// This method returns the value of the Simplex::modified property.
    inline bool get_modified() const {return modified;};
    /// This method sets the value of the Simplex::modified property to the argument.
    inline void set_modified(bool t) {modified = t;};
    /// This method returns the value of the Simplex::incept property.
    inline int get_incept() const {return incept;};
    /// This method sets the value of the Simplex::incept property to the argument.
    inline void set_incept(int n) {incept = n;};
    /// This method returns the value of the Simplex::energy property.
    inline double get_energy() const {return energy;};
    /// This method returns the value of the Simplex::volume property.
    inline double get_volume() const {return volume;};
    /// This method sets the value of the Simplex::sq_volume property to the argument and then also recomputes the Simplex::volume property.
    inline void set_squared_volume(double V) {sq_volume = V; volume = std::sqrt(std::abs(V));};
    /// This method returns the value of the Simplex::sq_volume property.
    inline double get_squared_volume() const {return sq_volume;};
    /// This method assumes the geometry is Lorentzian and returns true if the squared volume of the simplex is negative. 
    inline bool timelike() const {return (sq_volume < -std::numeric_limits<double>::epsilon());};
    /// This method assumes the geometry is Lorentzian and returns true if the squared volume of the simplex is positive. 
    inline bool spacelike() const {return (sq_volume > std::numeric_limits<double>::epsilon());};
    /// This method assumes the geometry is Lorentzian and returns true if the simplex is neither timelike nor spacelike, so that its volume is zero.
    inline bool lightlike() const {return (!timelike() && !spacelike());};
    /// This overloaded operator returns a simplex whose vertex set is the intersection of the vertices of its two arguments. The output's Simplex::active property is the conjunction of the same property of its arguments and the Simplex::parity of the output is the product of the same property of its arguments.
    friend Simplex operator ^(const Simplex&,const Simplex&);
    /// This method overrides the ostream operator so as to do a pretty print of an instance of the class.
    friend std::ostream& operator<< (std::ostream&,const Simplex&);
    friend class Complex;
  };
}
#endif

