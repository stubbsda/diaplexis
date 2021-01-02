#include <synarmosma/homology.h>
#include <synarmosma/homotopy.h>

#ifndef _sheeth
#define _sheeth

namespace DIAPLEXIS {
  /// A class representing one sheet of a multi-sheeted spacetime complex, akin to a Riemann surface.
  class Sheet {
   protected:
    /// This property contains the index of this sheet in the 
    /// Spacetime::codex vector.
    int index = -1;
    /// This property contains the index of this sheet's parent, 
    /// if it possesses one.
    int parent = -1;
    /// This property is the musical voice (in the polyphonic composition) corresponding to this 
    /// sheet, and is therefore only significant when performing musical hyphansis. 
    int voice = -1;
    /// This property is used to signal whether or not the array Sheet::hyphantic_notes 
    /// has been allocated, which should only be the case when performing musical 
    /// hyphansis. A value of -1 indicates that this array has not been allocated, while 
    /// a positive value indicates the size of the array. 
    int score_allocated = -1;
    /// This property indicates whether or not this sheet is currently active 
    /// or quiescent.
    bool active = false;
    /// This property stores a record of the hyphantic operations that have been successfully 
    /// executed on this sheet.
    std::string hyphantic_ops = "";
    /// This property is true if this sheet of the spacetime complex satisfies the axioms of a 
    /// combinatorial pseudomanifold.
    bool pseudomanifold = false;
    /// This property is true if this sheet of the spacetime complex satisfies the axioms of a 
    /// combinatorial pseudomanifold-with-boundary.
    bool boundary = false;
    /// This property is true if this sheet of the spacetime complex satisfies the axioms of an 
    /// orientable combinatorial pseudomanifold.
    bool orientable = false;
    /// This property contains the homology of this sheet of the spacetime complex, using a combinatorial 
    /// presentation for the groups.
    SYNARMOSMA::Homology* H;
    /// This property contains the fundamental group of this sheet of the spacetime complex, using a 
    /// combinatorial presentation.
    SYNARMOSMA::Homotopy* pi1;
    /// This property stores the set of notes which at each iteration will be 
    /// used to determine the hyphantic operations that will be executed, when 
    /// the Spacetime::weaving property is set to be musical. The score is read 
    /// from the file specified by Spacetime::hyphansis_score and parsed by the 
    /// parse_music_score() method. 
    std::vector<int>* hyphantic_notes;

    /// This method calls clear on the Sheet::H and Sheet::pi1 properties and restores the other properties of the class to their default value.
    void clear();
    /// This method sets all of the global topological properties for this sheet: Sheet::H, Sheet::pi1, Sheet::pseudomanifold, Sheet::boundary and Sheet::orientable.
    void set_topology(const SYNARMOSMA::Homology*,const SYNARMOSMA::Homotopy*,bool,bool,bool);
    /// This method writes the instance properties to a binary disk file and returns the number of bytes written to the file.
    int serialize(std::ofstream&) const;
    /// This method calls the clear() method on the instance and then reads the properties from a binary disk file and returns the number of bytes read.
    int deserialize(std::ifstream&);
    /// This method accepts as its arguments the maximum number of relaxation steps (Spacetime::max_iter) and the filename of the musical score (Spacetime::hyphansis_score). The method parses this file and selects those notes corresponding to its value of Sheet::voice to fill the array Sheet::hyphantic_notes, used by the Spacetime::musical_hyphansis method for the topological weaving of the spacetime complex.
    void parse_music_score(int,std::string&);
    /// This method sets the value of the Sheet::voice property to the method's argument.
    inline void set_voice(int);
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

  void Sheet::set_voice(int n)
  {
    voice = n;
  }
}
#endif
