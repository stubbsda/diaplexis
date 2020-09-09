#include <synarmosma/homology.h>
#include <synarmosma/homotopy.h>

#ifndef _sheeth
#define _sheeth

namespace DIAPLEXIS {
  /// A class representing one sheet of a multi-sheeted spacetime complex, akin to a Riemann surface.
  class Sheet {
   protected:
    /// This property contains the index of this sheet in the Spacetime::codex vector.
    int index = -1;
    /// This property contains the index of this sheet's parent, if it possesses one.
    int parent = -1;
    /// This property indicates whether or not this sheet is currently active or quiescent.
    bool active = false;
    /// This property stores a record of the hyphantic operations that have been successfully executed on this sheet.
    std::string hyphantic_ops = "";
    /// This property is true if this sheet of the spacetime complex satisfies the axioms of a combinatorial pseudomanifold.
    bool pseudomanifold = false;
    /// This property is true if this sheet of the spacetime complex satisfies the axioms of a combinatorial pseudomanifold-with-boundary.
    bool boundary = false;
    /// This property is true if this sheet of the spacetime complex satisfies the axioms of an orientable combinatorial pseudomanifold.
    bool orientable = false;
    /// This property contains the homology of this sheet of the spacetime complex, using a combinatorial presentation for the groups.
    SYNARMOSMA::Homology* H;
    /// This property contains the fundamental group of this sheet of the spacetime complex, using a combinatorial presentation.
    SYNARMOSMA::Homotopy* pi1;

    /// This method calls clear on the Sheet::H and Sheet::pi1 properties and restores the other properties of the class to their default value.
    void clear();
    /// This method sets all of the global topological properties for this sheet: Sheet::H, Sheet::pi1, Sheet::pseudomanifold, Sheet::boundary and Sheet::orientable.
    void set_topology(const SYNARMOSMA::Homology*,const SYNARMOSMA::Homotopy*,bool,bool,bool);
    /// This method writes the instance properties to a binary disk file and returns the number of bytes written to the file.
    int serialize(std::ofstream&) const;
    /// This method calls the clear() method on the instance and then reads the properties from a binary disk file and returns the number of bytes read.
    int deserialize(std::ifstream&);
   public:
    /// The default constructor which allocates the memory for the Sheet::H and Sheet::pi1 properties and sets Sheet::active to true.
    Sheet();
    /// This constructor has as its arguments the Sheet::index property and the base field and computation method for the Sheet::H property. 
    Sheet(int,SYNARMOSMA::Homology::Field,SYNARMOSMA::Homology::Method);
    /// This constructor has as its arguments the Sheet::index and Sheet::parent properties as well as the base field and computation method for the Sheet::H property. 
    Sheet(int,int,SYNARMOSMA::Homology::Field,SYNARMOSMA::Homology::Method);
    /// The standard copy constructor that copies over the properties from the source instance of this class.
    Sheet(const Sheet&);
    /// The overloaded assignment operator for this class, which behaves exactly like the copy constructor for this class.
    Sheet& operator =(const Sheet&);
    /// The destructor which frees the memory from the Sheet::H and Sheet::pi1 properties.
    ~Sheet();
    friend class Spacetime;
  };
}
#endif
