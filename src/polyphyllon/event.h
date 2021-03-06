#include <synarmosma/vertex.h>
#include <synarmosma/proposition.h>

#ifndef _eventh
#define _eventh

namespace DIAPLEXIS {
  template<class kind>
  class Event;

  template<class kind>
  std::ostream& operator <<(std::ostream&,const Event<kind>&);

  template<class kind>
  class Complex;

  /// A class representing the nodes of the spacetime complex, i.e. its event structure.

  /// This class is derived from the Synarmosma library's Vertex class,
  /// whose documentation should therefore also be consulted for more
  /// details.
  template<class kind>
  class Event: public SYNARMOSMA::Vertex<kind> {
   protected:
    /// This property represents a set of logical assertions concerning the event and its
    /// neighbourhood.
    SYNARMOSMA::Proposition theorem;
    /// This property determines whether or not this event lies at the edge of
    /// the spacetime complex.
    bool boundary = false;
    /// This property is set to true when the event's topological structure is modified, so
    /// that its Event::entwinement and Event::deficiency properties need to be recomputed.
    bool topology_modified = true;
    /// This property is set to true when the event's geometry is modified, so that its
    /// Event::obliquity, Event::deficiency and Event::geometric_deficiency properties need
    /// to be recomputed.
    bool geometry_modified = true;
    /// This property represents the extent to which the geometry of this event and
    /// its neighbours diverge from orthogonality, like the Event::entwinement it is
    /// non-negative.
    double obliquity = 0.0;
    /// This property stores the difference between the Event::obliquity and the inherited energy property of
    /// this event, so it can be positive or negative.
    double geometric_deficiency = 0.0;
    /// This property is the difference between the sum of the Event::entwinement and Event::obliquity 
    /// of the event and its inherited energy property.
    double deficiency = 0.0;
    /// This property represents the degree of topological entwinement or tangledness 
    /// of this event - it is a vector since there will be a distinct (non-negative) 
    /// entwinement value for each sheet of the spacetime complex. 
    std::vector<double> entwinement;
    /// This property is a set containing the index number of each of the sheets to which this 
    /// event belongs. 
    std::set<int> ubiquity;
    using SYNARMOSMA::Vertex<kind>::energy;
    using SYNARMOSMA::Vertex<kind>::incept;
    using SYNARMOSMA::Vertex<kind>::topological_dimension;
    using SYNARMOSMA::Vertex<kind>::anterior;
    using SYNARMOSMA::Vertex<kind>::posterior;
    using SYNARMOSMA::Vertex<kind>::neighbours;
    using SYNARMOSMA::Vertex<kind>::entourage;

    /// This method clears all of the instance's extended properties and sets the scalar properties to their default value.
    void clear() override;
   public:
    /// The default constructor which does nothing.
    Event();
    /// The standard copy constructor that copies over the properties from the source instance of this class.
    Event(const Event&);
    /// The overloaded assignment operator for this class, which behaves exactly like the copy constructor for this class.
    Event& operator =(const Event&);
    /// The destructor which does nothing for this class.
    ~Event() override;
    /// This method writes the instance properties to a binary disk file and returns the number of bytes written to the file.
    int serialize(std::ofstream&) const override;
    /// This method calls the clear() method on the instance and then reads the properties from a binary disk file and returns the number of bytes read.
    int deserialize(std::ifstream&) override;
    /// This method returns true if this event belongs to at least one sheet of the spacetime complex and false otherwise.
    bool active() const; 
    /// This method returns true if this event belongs to the sheet whose index is the method's argument and false otherwise.
    bool active(int) const;
    /// This method removes this event from all of the spacetime's sheets and sets the Event::topology_modified property to true.
    void deactivate();
    /// This method adds this event to the sheet whose index is given by the argument and sets the Event::topology_modified property to true.
    void activate(int);
    /// This method removes this event from the sheet whose index is given by the argument and sets the Event::topology_modified property to true.
    void deactivate(int);
    /// This method sets the Event::ubiquity property to the method's argument and sets the Event::topology_modified property to true.
    void set_ubiquity(const std::set<int>&);
    /// This method sets the argument to the Event::ubiquity property.
    void get_ubiquity(std::set<int>&) const;
    /// This method returns the number of sheets of the spacetime complex to which this event belongs.
    int presence() const;
    /// This method clears the inherited posterior property.
    void clear_posterior();
    /// This method sets the argument to the inherited posterior property.
    void get_posterior(std::set<int>&) const;
    /// This method adds the argument to the inherited posterior property and returns true if it was not already a member of this set.
    bool add_posterior(int);
    /// This method clears the inherited anterior property.
    void clear_anterior();
    /// This method sets the argument to the inherited anterior property.
    void get_anterior(std::set<int>&) const;
    /// This method adds the argument to the inherited anterior property and returns true if it was not already a member of this set.
    bool add_anterior(int);
    /// This method clears the inherited entourage property.
    void clear_entourage();
    /// This method sets the argument equal to the inherited entourage property.
    void get_entourage(std::set<int>&) const;
    /// This method sets the inherited entourage property to the argument.
    void set_entourage(const std::set<int>&);
    /// This method removes the argument from the inherited entourage property, returning true if the element was a member of the entourage and false otherwise.
    bool drop_entourage(int);
    /// This method sets the argument to the inherited neighbours property.
    void get_neighbours(std::set<int>&) const;
    /// This method sets the inherited neighbours property to the argument as well as setting Event::topology_modified to true.
    void set_neighbours(const std::set<int>&);
    /// This method returns true if the argument is an element of the inherited neighbours property and false otherwise.
    bool is_neighbour(int) const;
    /// This method adds the argument to the inherited neighbours property and returns true if this is a new neighbour, as well as setting the Event::topology_modified property to true.
    bool add_neighbour(int);
    /// This method deletes the argument from the inherited neighbours property and returns true if this was indeed a neighbour, as well as setting the Event::topology_modified property to true.
    bool drop_neighbour(int);
    /// This method returns the value of the Event::deficiency property.
    double get_deficiency() const;
    /// This method sets the value of the Event::deficiency property to the argument.
    void set_deficiency(double);
    /// This method sets the argument equal to the Event::entwinement property.
    void get_entwinement(std::vector<double>&) const;
    /// This method sets the Event::entwinement property to the argument.
    void set_entwinement(const std::vector<double>&);
    /// This method sets the Event::entwinement property using a C-style array as the method's first argument, along with an integer second argument that is the array's length.
    void set_entwinement(const double*,int);
    /// This method returns the value of the Event::obliquity property.
    double get_obliquity() const;
    /// This method sets the value of the Event::obliquity property to the argument.
    void set_obliquity(double);
    /// This method returns the value of the inherited incept property.
    int get_incept() const;
    /// This method sets the value of the inherited incept property to the argument.
    void set_incept(int);
    /// This method returns the value of the Event::boundary property.
    bool get_boundary() const;
    /// This method sets the value of the Event::boundary property to the argument.
    void set_boundary(bool);
    /// This method returns the value of the Event::topology_modified property.
    bool get_topology_modified() const;
    /// This method sets the value of the Event::topology_modified property to the argument.
    void set_topology_modified(bool);
    /// This method returns the value of the Event::geometry_modified property.
    bool get_geometry_modified() const;
    /// This method sets the value of the Event::geometry_modified property to the argument.
    void set_geometry_modified(bool);
    /// This method returns the value of the Event::geometric_deficiency property.
    double get_geometric_deficiency() const;
    /// This method sets the value of the Event::geometric_deficiency property to the argument.
    void set_geometric_deficiency(double);
    /// This method returns the value of the inherited topological dimension property.
    int get_topological_dimension() const;
    /// This method sets the value of the inherited topological dimension property to the argument.
    void set_topological_dimension(int);
    /// This method overrides the ostream operator so as to do a pretty print of an instance of the class.
    friend std::ostream& operator << <>(std::ostream&,const Event<kind>&);
    friend class Complex<kind>;
  };

  template<class kind>
  inline bool Event<kind>::active() const 
  {
    return !ubiquity.empty();
  } 

  template<class kind>
  inline bool Event<kind>::active(int n) const 
  {
    return (ubiquity.count(n) > 0);
  }

  template<class kind>
  inline void Event<kind>::deactivate() 
  {
    ubiquity.clear(); 
    topology_modified = true;
  }

  template<class kind>
  inline void Event<kind>::activate(int n) 
  {
    ubiquity.insert(n); 
    topology_modified = true;
  }

  template<class kind>
  inline void Event<kind>::deactivate(int n) 
  {
    ubiquity.erase(n); 
    topology_modified = true;
  }

  template<class kind>
  inline void Event<kind>::set_ubiquity(const std::set<int>& S) 
  {
    ubiquity = S; 
    topology_modified = true;
  }

  template<class kind>
  inline void Event<kind>::get_ubiquity(std::set<int>& S) const 
  {
    S = ubiquity;
  }

  template<class kind>
  inline int Event<kind>::presence() const 
  {
    return (signed) ubiquity.size();
  }

  template<class kind>
  inline void Event<kind>::clear_posterior() 
  {
    posterior.clear();
  }

  template<class kind>
  inline void Event<kind>::get_posterior(std::set<int>& S) const 
  {
    S = posterior;
  }

  template<class kind>
  inline bool Event<kind>::add_posterior(int n)
  {
    if (posterior.count(n) == 0) {
      posterior.insert(n);
      return true;
    }
    return false;
  }

  template<class kind>
  inline void Event<kind>::clear_anterior() 
  {
    anterior.clear();
  }

  template<class kind>
  inline void Event<kind>::get_anterior(std::set<int>& S) const 
  {
    S = anterior;
  }

  template<class kind>
  inline bool Event<kind>::add_anterior(int n)
  {
    if (anterior.count(n) == 0) {
      anterior.insert(n);
      return true;
    }
    return false;
  }

  template<class kind>
  inline void Event<kind>::clear_entourage() 
  {
    entourage.clear();
  }

  template<class kind>
  inline void Event<kind>::get_entourage(std::set<int>& S) const 
  {
    S = entourage;
  }

  template<class kind>
  inline void Event<kind>::set_entourage(const std::set<int>& S) 
  {
    entourage = S;
  }

  template<class kind>
  inline bool Event<kind>::drop_entourage(int n)
  {
    std::set<int>::const_iterator it = std::find(entourage.begin(),entourage.end(),n);
    if (it != entourage.end()) {
      entourage.erase(it);
      return true;
    }
    return false;
  }

  template<class kind>
  inline void Event<kind>::get_neighbours(std::set<int>& S) const 
  {
    S = neighbours;
  }

  template<class kind>
  inline void Event<kind>::set_neighbours(const std::set<int>& S) 
  {
    neighbours = S; 
    topology_modified = true;
  }

  template<class kind>
  inline bool Event<kind>::is_neighbour(int n) const 
  {
    return (neighbours.count(n) > 0);
  }

  template<class kind>
  inline bool Event<kind>::add_neighbour(int n)
  {
    if (neighbours.count(n) == 0) {
      neighbours.insert(n);
      topology_modified = true;
      return true;
    }
    return false;
  }

  template<class kind>
  inline bool Event<kind>::drop_neighbour(int n)
  {
    std::set<int>::const_iterator it = std::find(neighbours.begin(),neighbours.end(),n);
    if (it != neighbours.end()) {
      neighbours.erase(it);
      topology_modified = true;
      return true;
    }
    return false;
  }

  template<class kind>
  inline double Event<kind>::get_deficiency() const 
  {
    return deficiency;
  }

  template<class kind>
  inline void Event<kind>::set_deficiency(double x) 
  {
    deficiency = x;
  }

  template<class kind>
  inline void Event<kind>::get_entwinement(std::vector<double>& x) const 
  {
    x = entwinement;
  }

  template<class kind>
  inline void Event<kind>::set_entwinement(const std::vector<double>& x) 
  {
    entwinement = x;
  }

  template<class kind>
  inline void Event<kind>::set_entwinement(const double* x,int n)
  {
    if (n < 1) throw std::invalid_argument("Illegal array length in the Event::set_entwinement method!");

    entwinement.clear(); 
    for(int i=0; i<n; ++i) {
      entwinement.push_back(x[i]);
    }
  }

  template<class kind>
  inline double Event<kind>::get_obliquity() const 
  {
    return obliquity;
  }

  template<class kind>
  inline void Event<kind>::set_obliquity(double x) 
  {
    obliquity = x;
  }

  template<class kind>
  inline int Event<kind>::get_incept() const 
  {
    return incept;
  }

  template<class kind>
  inline void Event<kind>::set_incept(int n) 
  {
    incept = n;
  }

  template<class kind>
  inline bool Event<kind>::get_boundary() const 
  {
    return boundary;
  }

  template<class kind>
  inline void Event<kind>::set_boundary(bool t) 
  {
    boundary = t;
  }

  template<class kind>
  inline bool Event<kind>::get_topology_modified() const 
  {
    return topology_modified;
  }

  template<class kind>
  inline void Event<kind>::set_topology_modified(bool t) 
  {
    topology_modified = t;
  }

  template<class kind>
  inline bool Event<kind>::get_geometry_modified() const 
  {
    return geometry_modified;
  }

  template<class kind>
  inline void Event<kind>::set_geometry_modified(bool t) 
  {
    geometry_modified = t;
  }

  template<class kind>
  inline double Event<kind>::get_geometric_deficiency() const 
  {
    return geometric_deficiency;
  }

  template<class kind>
  inline void Event<kind>::set_geometric_deficiency(double x) 
  {
    geometric_deficiency = x;
  }

  template<class kind>
  inline int Event<kind>::get_topological_dimension() const 
  {
    return topological_dimension;
  }

  template<class kind>
  inline void Event<kind>::set_topological_dimension(int n) 
  {
    topological_dimension = n;
  }
}
#endif
