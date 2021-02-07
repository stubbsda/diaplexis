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
    /// This property indicates whether or not this sheet is currently active, when it is zero, 
    /// or quiescent, when it is greater than zero, indicating how many more relaxation steps 
    /// the sheet will remain quiescent. 
    int sleep = 0;
    /// This floating point property influences the number of child threads that will be created 
    /// by this sheet. It must lie between 0 and 1, with 0 corresponding to sterility. This property 
    /// is reset after each period of quiescence of the sheet, declining as the total number of 
    /// existing sheets grows. 
    double fertility = 0.0;
    /// This floating point property controls whether this sheet will sleep and if so, the number 
    /// of relaxation steps during which it will remain quiescent. It must lie between 0 and 1, with 
    /// 0 corresponding to a sheet which will never sleep. It is reset after each hyphansis stage 
    /// and period of quiescence; in the former case the drowsiness is increased in proportion to 
    /// the number of successful hyphantic operations. 
    double drowsiness = 0.0;
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
    /// This method calls clear on the Sheet::H and Sheet::pi1 properties and restores the other properties of the class to their default value.
    void clear();
    /// This method computes the values of the Sheet::H and (if the second argument is true) Sheet::pi1 properties from the abstract simplicial complex that is its first argument. The method then also computes the Sheet::pseudomanifold, Sheet::boundary and Sheet::orientable properties. 
    void compute_topology(const SYNARMOSMA::Nexus*,bool);
    /// This method sets all of the global topological properties for this sheet: Sheet::H, Sheet::pi1, Sheet::pseudomanifold, Sheet::boundary and Sheet::orientable.
    void set_topology(const SYNARMOSMA::Homology*,const SYNARMOSMA::Homotopy*,bool,bool,bool);
    /// This method writes the instance properties to a binary disk file and returns the number of bytes written to the file.
    int serialize(std::ofstream&) const;
    /// This method calls the clear() method on the instance and then reads the properties from a binary disk file and returns the number of bytes read.
    int deserialize(std::ifstream&);
    /// This method accepts as its arguments the maximum number of relaxation steps (Spacetime::max_iter) and the filename of the musical score (Spacetime::hyphansis_score). The method parses this file and selects those notes corresponding to its value of Sheet::voice to fill the array Sheet::hyphantic_notes, used by the Spacetime::musical_hyphansis method for the topological weaving of the spacetime complex.
    void parse_music_score(int,std::string&);
    /// This method writes to the screen the values of the Sheet::index and Sheet::parent properties; if the argument is true it also writes the value of Sheet::fertility, Sheet::drowsiness and Sheet::sleep, otherwise Sheet::voice.
    void write_parameters(bool) const;
    /// This method determines the number of offspring to be generated at this step, based on a uniform random variate between zero and one (the first argument) and the value of Spacetime::max_children (the second argument); the method returns the number of offspring given the current Sheet::fertility value.
    int reproduction(double,int);
    /// This method writes the output of calling the write method of the Sheet::H property to the method's argument.
    inline void write_homology(std::string&) const;
    /// This method writes the output of calling the write method of the Sheet::pi1 property to the method's argument.
    inline void write_homotopy(std::string&) const;
    /// This method sets the value of the Sheet::hyphantic_ops property to its default (empty) value.
    inline void clear_hyphansis();
    /// This method appends its argument to the Sheet::hyphantic_ops property.
    inline void append_operation(const std::string&);
    /// This method decrements the value of Sheet::sleep by one. 
    inline void hibernate();
    /// This method sets the value of the Sheet::voice property to the method's argument.
    inline void set_voice(int);
    /// This method sets the value of the Sheet::sleep property to the method's argument.
    inline void set_sleep(int);
    /// This method sets the value of the Sheet::drowsiness property to the method's argument.
    inline void set_drowsiness(double);
    /// This method sets the value of the Sheet::fertility property to the method's argument.
    inline void set_fertility(double);
    /// This method returns the current value of the Sheet::voice property. 
    inline int get_voice() const;
    /// This method writes to its first argument the vector of notes in Sheet::hyphantic_notes for the element corresponding to the method's second argument. 
    inline void get_notes(std::vector<int>&,int) const;
    /// This method returns the current value of the Sheet::sleep property. 
    inline int get_sleep() const;
    /// This method returns the current value of the Sheet::pseudomanifold property. 
    inline bool get_pseudomanifold() const;
    /// This method returns the current value of the Sheet::orientable property. 
    inline bool get_orientable() const;
    /// This method returns the current value of the Sheet::boundary property. 
    inline bool get_boundary() const;
    /// This method returns the current value of the Sheet::drowsiness property. 
    inline double get_drowsiness() const;
    /// This method returns the current value of the Sheet::hyphantic_ops property. 
    inline std::string get_hyphantic_operations() const;
  };

  void Sheet::clear_hyphansis()
  {
    hyphantic_ops = "";
  }

  void Sheet::write_homotopy(std::string& output) const
  {
    output = pi1->write();
  }

  void Sheet::write_homology(std::string& output) const
  {
    output = H->write();
  }

  void Sheet::append_operation(const std::string& op)
  {
    hyphantic_ops += op;
  }

  void Sheet::hibernate()
  {
    sleep -= 1;
  }

  void Sheet::set_voice(int n)
  {
    voice = n;
  }

  void Sheet::set_sleep(int n)
  {
    sleep = n;
  }

  void Sheet::set_drowsiness(double x)
  {
    drowsiness = x;
  }

  void Sheet::set_fertility(double x)
  {
    fertility = x;
  }

  bool Sheet::get_pseudomanifold() const
  {
    return pseudomanifold;
  }

  bool Sheet::get_boundary() const
  {
    return boundary;
  }

  bool Sheet::get_orientable() const
  {
    return orientable;
  }

  int Sheet::get_voice() const
  {
    return voice;
  }

  void Sheet::get_notes(std::vector<int>& melody,int n) const
  {
    melody = hyphantic_notes[n];
  }

  int Sheet::get_sleep() const
  {
    return sleep;
  }

  double Sheet::get_drowsiness() const
  {
    return drowsiness;
  }

  std::string Sheet::get_hyphantic_operations() const
  {
    return hyphantic_ops;
  }
}
#endif
