#include "event.h"
#include "simplex.h"

#include "synarmosma/homology.h"
#include "synarmosma/homotopy.h"

#ifndef _complexh
#define _complexh

namespace DIAPLEXIS {
  /// A class representing an abstract simplicial complex of finite dimension. 
  class Complex {
   private:
    /// This property is true if this complex satisfies the axioms of a combinatorial 
    /// pseudomanifold.
    bool pseudomanifold = false;
    /// This property is true if this complex satisfies the axioms of a 
    /// pseudomanifold-with-boundary.
    bool boundary = false;
    /// This property is true if this complex satisfies the axioms of an 
    /// orientable pseudomanifold.
    bool orientable = false;
    /// This property stores the 0-simplices (events) of the complex as a 
    /// vector of type Event.
    std::vector<Event> events;
    /// This property stores the higher-dimensional simplices, as an array of 
    /// 1 + Complex::ND vectors of type Simplex.
    std::vector<Simplex>* simplices;
    /// This property contains a hash map for each dimension greater than zero 
    /// that allows a rapid lookup of a d-simplex's index in Complex::simplices[d] 
    /// given the simplex's vertex set.
    SYNARMOSMA::hash_map* index_table;
    /// This property stores the homology of the complex, using the 
    /// Homology class of the Synarmosma library.
    SYNARMOSMA::Homology* H;
    /// This property stores the fundamental group of the complex, using the 
    /// Homotopy class of the Synarmosma library.
    SYNARMOSMA::Homotopy* pi1;
    /// The property which handles random number generation for this class and whose default seed 
    /// is set to the current time value (seconds since January 1, 1970).
    SYNARMOSMA::Random* RND;

    /// The maximum dimension of the complex and thus size of the Complex::simplices and 
    /// Complex::index_table arrays.
    static const int ND = 10;
    /// The combinatorial size of the subgraph which is used to compute the Event::entwinement 
    /// property.
    static const int topological_radius = 4;

    /// This method allocates the memory for the properties Complex::simplices, Complex::index_table, Complex::pi1, Complex::H and Complex::RND.
    void allocate();
    void inversion(int);
    /// This method computes the inherited neighbours property of the elements of the Complex::events vector.
    void compute_neighbours();
    void compute_entourages(int);
    bool simplex_addition(int,int,const std::set<int>&,int = -1);
    bool simplex_addition(const std::set<int>&,const std::set<int>&,int = -1);
    void simplex_deletion(int,int,int);
    /// This method calculates which events have had their topology modified based on modified d-simplices (d > 0) and the Complex::topological_radius property, setting the topology_modified property of these events to true.
    void compute_modified_events();
    /// This method accepts a set of modified events and calculates all of the d-simplices (d > 0) in the complex which should therefore also be considered as modified and sets the Simplex::modified property appropriately.
    void compute_dependent_simplices(const std::set<int>&);
    void simplicial_implication(int);
    void simplicial_implication(int,int) const;
    /// This method computes the Event::topological_dimension property of the complex's 0-simplices. 
    void compute_simplicial_dimension();
    /// This method sets the atomic propositions of the Event::theorem property for each active event in the complex. The method's argument is the total number of atomic propositions in the whole complex and the algorithm assigns more atomic propositions to an event's theorem in direct proportion to its energy and topological dimension. The method returns the average number of atomic propositions assigned to each active event.
    double set_logical_atoms(int);
    /// This method computes the conjunction of the Event::theorem property of the argument and the theorem property of each its active neighbour events, calculates the satisfiability of the result and adds this to a running sum. The method then returns this sum divided by the number of neighbouring events.
    double logical_energy(int) const;
    /// This method computes the conjunction of the Event::theorem property of the argument and the theorem property of each of its active neighbour events, finally returning the satisfiability of this conjunction.
    bool logical_conformity(int) const;
    /// This method computes the energy of a d-simplex where the first argument is the dimension and the second is the index of the simplex in Complex::simplices[d]. The energy is defined to be the arithmetic mean of the energy of its vertices.
    void compute_simplex_energy(int,int);
    /// This method computes the parity of a d-simplex where the first argument is the dimension and the second is the index of the simplex in Complex::simplices[d]. The parity is defined to be the product of the parity of its (1+d)*d/2 edges.
    void compute_simplex_parity(int,int);
    bool edge_parity_mutation(int,int);
    bool edge_parity_mutation(int,int,int);
    /// This method accepts the index of a 1-simplex in the complex and then recomputes the parity of all the higher-dimensional simplices which contain this edge.
    void recompute_parity(int);
    /// This method accepts a set of indices of 1-simplices in the complex and recomputes the parity of all the higher-dimensional simplices which contain these edges.
    void recompute_parity(const std::set<int>&);
    /// This method computes the Simplex::parity property for the entire collection of d-simplices (d > 1) in the complex, successively calling the compute_simplex_parity method on each d-simplex whose Simplex::modified property is true. 
    void compute_parity();
    /// This method loops over all the 1-simplices in the simplicial complex - if a 1-simplex is active and at least one of its vertices has a non-zero energy, the corresponding element of the method's argument is given the value 1, otherwise 0. 
    void determine_flexible_edges(std::vector<int>&) const;
    inline bool edge_exists(int,int,int) const;
    int cardinality(int,int) const;
    void compute_degree_distribution(bool,int) const;
    void compute_connectivity_distribution(bool,int) const;
    std::pair<double,double> random_walk(int) const;
    void vertex_degree_statistics(double*,int) const;
    void compute_fvector(std::vector<int>&,std::vector<int>&,int) const;
    void compute_hvector(std::vector<int>&,int) const;
    inline void compute_graph(SYNARMOSMA::Graph*,int) const;
    inline void compute_graph(SYNARMOSMA::Graph*,int,int) const;
    void compute_graph(SYNARMOSMA::Graph*,int,int,int) const;
    void compute_graph(SYNARMOSMA::Graph*,int*,int) const;
    void compute_global_nexus(SYNARMOSMA::Nexus*,int) const;
    void compute_local_nexus(SYNARMOSMA::Nexus*,int,int) const;
    void compute_global_topology(bool);
    /// This method accepts as its first argument an event and then proceeds to compute the number of d-simplices (d > 0) that contain this event, appending this figure to the second argument which is the method's output.
    void simplex_membership(int,std::vector<int>&) const;
    int chromatic_number(int) const;
    int combinatorial_distance(int,int,int) const;
    double entwinement(int) const;
    double cyclic_resistance(int) const;
    /// This method returns the maximum value of the degree - the number of neighbours - over the entire collection of events.
    int max_degree() const;
    int total_dimension(int) const;
    int structural_index(int) const;
    int dimension(int) const;
    int vertex_valence(int,int) const;
    int vertex_dimension(int,int) const;
    int weighted_entourage(int,int) const;
    int cyclicity(int) const;
    double dimensional_stress(int,int) const;
    double dimensional_stress(int,int,int) const;
    double dimensional_frontier(int,int) const;
    double dimensional_uniformity(int,int) const;
    int circuit_rank(int) const;
    int euler_characteristic(int) const;
    int component_analysis(std::vector<int>&,int) const;
    int entourage(int,int) const;
    bool connected(int) const;
    bool consistent(int) const;
    /// This method writes the simplicial inclusion relations for the complex as a directed graph in the DOT file format, where the method's argument is the filename.
    void write_incastrature(const std::string&) const;
    /// This method writes to the console a summary of the topology of the complex, including the simplices that each event belongs to and the number of simplices of a given dimension.
    void write_topology() const;
    /// This method returns false if there is an inactive event whose energy property has a non-zero value.
    bool energy_check() const;
    double total_energy(int) const;
    /// This method carries out the energy diffusion among neighbouring active events in the complex - the method's argument is a parameter that controls the rate of energy transfer.
    void energy_diffusion(double);
    /// This method carries out the energy diffusion among the complex's events using a chip-firing game algorithm, with the method's argument equal to the number of chips to be distributed among the active events.
    void energy_diffusion(int);
    /// This method computes the extent to which a given d-simplex (d is the first argument while the second is the simplex's index in Complex::simplices[d]) can be embedded in Euclidean space by calculating the number of non-zero eigenvalues of a matrix derived from the simplex's edge lengths. The simplex can only be embedded if it is entirely timelike or spacelike, otherwise the method returns -1.
    int simplex_embedding(int,int) const;
    double parity_hamiltonian(double,bool,int) const;
    void write_graph(const std::string&,int) const;
    double distribution_fitness(int*,const std::vector<int>&,int) const;

   public:
    /// The default constructor which calls the allocate() method.
    Complex();
    /// The destructor which frees the memory from the instance's properties (Complex::simplices, Complex::RND etc.).
    ~Complex();
    /// This method attempts to distribute the complex's active events among N processors where N is the argument in such a way as to minimize the output of the distribution_fitness() method.
    void distribute(int) const; 
    /// This method writes the instance properties to a binary disk file and returns the number of bytes written to the file. 
    int serialize(std::ofstream&) const;
    /// This method calls the clear() method on the instance and then reads the properties from a binary disk file and returns the number of bytes read.
    int deserialize(std::ifstream&);
    /// This method clears all of the instance's extended properties and sets the scalar properties to their default value.
    void clear();
    /// This method assembles the vector of neighbour sets from each active event in the complex.
    void get_edge_topology(std::vector<std::set<int> >&) const;
    /// This method writes to the console a set of Event properties for the event whose index is the method's argument.
    void write_event_data(int) const;
    void get_energy_values(std::vector<double>&,int) const;
    void get_deficiency_values(std::vector<double>&,int) const;
    /// This method returns the inherited neighbours property of the event with index equal to the method's argument.
    inline std::set<int> get_neighbours(int n) const {return events[n].neighbours;};
    /// This method returns the homology of the complex as a compact string, the method's argument.
    inline void get_homology(std::string& s) const {s = H->write();};
    /// This method returns a combinatorial presentation of the homotopy group of the complex as a compact string, the method's argument.
    inline void get_homotopy(std::string& s) const {s = pi1->write();};
    /// This method sets the field over which the complex's homology is calculated to the argument.
    inline void set_homology_field(SYNARMOSMA::Homology::Field F) {H->set_field(F);};
    /// This method sets the method by which the complex's homology is calculated to the argument.
    inline void set_homology_method(SYNARMOSMA::Homology::Method M) {H->set_method(M);};
    /// This method returns the field over which the complex's homology is calculated.
    inline SYNARMOSMA::Homology::Field get_homology_field() const {return H->get_field();};
    /// This method returns the method by which the complex's homology is calculated.
    inline SYNARMOSMA::Homology::Method get_homology_method() const {return H->get_method();};
    inline int get_dimension(int sheet) const {return dimension(sheet);};
    inline int get_event_dimension(int n,int sheet) const {return vertex_dimension(n,sheet);};
    inline int get_cardinality(int d,int sheet) const {return cardinality(d,sheet);};
    inline void get_vertex_degree_statistics(double* output,int sheet) const {vertex_degree_statistics(output,sheet);};
    inline void get_fvector(std::vector<int>& f,std::vector<int>& fstar,int sheet) const {compute_fvector(f,fstar,sheet);};
    /// This method returns the vertices of a d-simplex (d > 0) where the arguments are the dimension and index in Complex::simplices[d], while the final argument is an integer array of length 1+d to store the vertices.
    inline void get_simplex_vertices(int d,int n,int* vx) const {simplices[d][n].get_vertices(vx);};
    /// This method returns a string of a d-simplex's vertices (where d > 0) separated by a colon in ascending order.  
    inline std::string get_simplex_key(int d,int n) const {return SYNARMOSMA::make_key(simplices[d][n].vertices);};
    /// This method returns the index of the 1-simplex whose vertex set is the method's argument; if the edge does not exist the method returns -1. 
    inline int get_edge_index(const std::set<int>&) const;
    inline int get_cyclicity(int sheet) const {return cyclicity(sheet);};
    inline int get_circuit_rank(int sheet) const {return circuit_rank(sheet);};
    /// This method returns whether or not the simplicial complex is orientable, assuming it is a combinatorial pseudomanifold.
    inline bool is_orientable() const {return orientable;};
    inline int get_euler_characteristic(int sheet) const {return euler_characteristic(sheet);};
    /// This method returns the Event::obliquity property of the event with index equal to the method's argument.
    inline double get_event_obliquity(int n) const {return events[n].obliquity;};
    /// This method returns the inherited energy property of the event with index equal to the method's argument.
    inline double get_event_energy(int n) const {return events[n].get_energy();};
    /// This method returns the Event::deficiency property of the event with index equal to the method's argument.
    inline double get_event_deficiency(int n) const {return events[n].deficiency;};
    inline double get_event_entwinement(int n,int sheet) const {return events[n].entwinement[sheet];};
    /// This method returns the number of elements in the Complex::events vector.
    inline int get_events() const {return (signed) events.size();};
    /// This method returns the output of the total_energy() method.
    inline double get_total_energy(int sheet) const {return total_energy(sheet);};
    /// This method returns the cardinality of the entourage property of the event (if the first argument is zero) or simplex (if it is greater than zero) whose index is specified by the second argument. 
    inline int get_entourage_cardinality(int d,int n) const {int output = (d == 0) ? (signed) events[n].neighbours.size() : (signed) simplices[d][n].entourage.size(); return output;}; 
    /// This method returns the incept property of an Event or Simplex - depending on the dimension (the first argument) - of the element of Complex::events or Complex::simplices[d] with index value equal to the second argument.
    inline int get_incept(int d,int n) const {int output = (d == 0) ? events[n].incept : simplices[d][n].incept; return output;};
    /// This method returns the value of the Complex::pseudomanifold property and sets the argument to the value of the Complex::boundary property.
    inline bool is_pseudomanifold(bool*) const;
    friend class Spacetime;
  };

  bool Complex::edge_exists(int u,int v,int sheet) const
  {
    std::set<int> vx;
    vx.insert(u); vx.insert(v);
    SYNARMOSMA::hash_map::const_iterator qt = index_table[1].find(vx);  
    if (qt == index_table[1].end()) return false;
    if (sheet == -1) {
      if (!simplices[1][qt->second].active()) return false;
    }
    else {
      if (!simplices[1][qt->second].active(sheet)) return false;
    }
    return true;
  }

  int Complex::get_edge_index(const std::set<int>& vx) const 
  {
#ifdef DEBUG
    assert(vx.size() == 2);
#endif    
    SYNARMOSMA::hash_map::const_iterator qt = index_table[1].find(vx); 
    if (qt == index_table[1].end()) return -1;
    return qt->second;
  }

  bool Complex::is_pseudomanifold(bool* bdry) const 
  {
    *bdry = boundary;
    return pseudomanifold;
  }

  void Complex::compute_graph(SYNARMOSMA::Graph* G,int sheet) const
  {
    if (events.empty()) return;

    int offset[events.size()];
    compute_graph(G,offset,sheet);
  }

  void Complex::compute_graph(SYNARMOSMA::Graph* G,int base,int sheet) const
  {
    compute_graph(G,base,Complex::topological_radius,sheet);
  }
}
#endif
