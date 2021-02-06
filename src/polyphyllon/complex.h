#include <synarmosma/homology.h>
#include <synarmosma/homotopy.h>

#include "event.h"
#include "simplex.h"

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
    /// The combinatorial size of the subgraph which is used to compute the Event::entwinement 
    /// property.
    int topological_radius = 4;
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

    /// This method allocates the memory for the properties Complex::simplices, Complex::index_table, Complex::pi1, Complex::H and Complex::RND.
    void allocate();
    /// This method transforms the complex by converting it to a graph (one-dimensional complex) whose topology is the inverse of the current topology (using a different random sheet for each newly created 1-simplex), calling the compute_entourages() method with the argument -1. The method's argument is the number of sheets in the complex. 
    void inversion(int);
    /// This method computes the inherited neighbours property of the elements of the Complex::events vector.
    void compute_neighbours();
    /// This method calculates the value of the inherited entourage property for the elements of the Complex::events and Complex::simplices arrays and a given sheet, the method's unique argument.
    void compute_entourages(int);
    /// This method adds a 1-simplex to the complex, with the first two arguments specifying the events connected by this new edge, the third argument the ubiquity of the edge and the optional final argument the incept property of the Simplex class. The method returns false if the edge already exists in the complex and has the same ubiquity, true otherwise.
    bool simplex_addition(int,int,const std::set<int>&,int = -1);
    /// This method adds a d-simplex to the complex, with the first argument specifying the vertex set of the new simplex, the second argument its ubiquity and the optional final argument the incept property of the simplex. The method returns false if this simplex already exists in the complex and has the same ubiquity, true otherwise.
    bool simplex_addition(const std::set<int>&,const std::set<int>&,int = -1);
    /// This method deletes a d-simplex by removing the sheet index (the method's third argument) from the Simplex::ubiquity property, where d is the first argument and the second argument is the index of the simplex in Complex::simplices[d]. The method then recursively deletes all the higher-dimensional simplices which depend on this simplex from the sheet. The method returns false if the simplex doesn't belong to the sheet in question and true otherwise. 
    bool simplex_deletion(int,int,int);
    /// This method calculates which events have had their topology modified based on modified d-simplices (d > 0) and the Complex::topological_radius property, setting the topology_modified property of these events to true.
    void compute_modified_events();
    /// This method accepts a set of modified events and calculates all of the d-simplices (d > 0) in the complex which should therefore also be considered as modified and sets the Simplex::modified property appropriately.
    void compute_dependent_simplices(const std::set<int>&);
    /// This method ensures that the entailment property among simplices is respected, so that if a d-simplex is belongs to the sheet indicated by the argument, then so do all of its (d-1)-dimensional faces and so forth.
    void simplicial_implication(int);
    /// This method will calculate all of the n-simplices (n > 1) belonging to the sheet indicated by the second argument that are entailed by the base event (the method's first argument) and its neighbours (via their mutual edges) and list them, as well as checking to see if they already exist in the complex.
    void simplicial_implication(int,int) const;
    /// This method computes the inherited topological_dimension property of the Event class for the complex's 0-simplices. 
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
    /// This method seeks to change the parity of the edge joining together the two events specified by the method's first two arguments, existing on the sheet indicated by its third argument. Assuming the edge exists, the method alters its parity (changing it to +1 or -1 with equal probability if it is zero and multiplying it by -1 otherwise), calls recompute_parity() and returns true.
    bool edge_parity_mutation(int,int,int);
    /// This method accepts the index of a 1-simplex in the complex and then recomputes the parity of all the higher-dimensional simplices which contain this edge.
    void recompute_parity(int);
    /// This method accepts a set of indices of 1-simplices in the complex and recomputes the parity of all the higher-dimensional simplices which contain these edges.
    void recompute_parity(const std::set<int>&);
    /// This method computes the Simplex::parity property for the entire collection of d-simplices (d > 1) in the complex, successively calling the compute_simplex_parity method on each d-simplex whose Simplex::modified property is true. 
    void compute_parity();
    /// This method loops over all the 1-simplices in the simplicial complex - if a 1-simplex is active and at least one of its vertices has a non-zero energy, the corresponding element of the method's argument is given the value 1, otherwise 0. 
    void determine_flexible_edges(std::vector<int>&) const;
    /// This method returns true if there is an active edge belonging to the sheet indicated in the third argument which connects the two vertices that are the method's first two arguments.
    bool edge_exists(int,int,int) const;
    /// This method returns the number of active d-simplices belonging to the sheet indicated by the final argument, where d is the method's first argument.
    int cardinality(int,int) const;
    /// This method computes the 1-skeleton of the sheet indicated by the second argument, calls the method Graph::degree_distribution from the Synarmosma library and prints out the vertex degree histogram to the console. If the method's first argument is true the histogram is logarithmic and linear otherwise. 
    void compute_degree_distribution(bool,int) const;
    /// This method computes and writes to the console the histogram of combinatorial distances among the complex's events belonging to the sheet indicated by the second argument. The first argument controls the algorithm used; if true it assumes memory is abundant relative to the number of events and uses a direct method for the calculation.
    void compute_connectivity_distribution(bool,int) const;
    /// This method computes the 1-skeleton of the sheet indicated by the unique argument and using the compute_graph() method and then returns the output of the random_walk method of the Graph class in the Synarmosma library.
    std::pair<double,double> random_walk(int) const;
    /// This method accepts as its argument an array of three double precision numbers and computes the 1-skeleton of the sheet indicated by the second argument and then uses Graph class methods to assign the maximum degree, minimum degree and average degree to the elements of the method's first argument.
    void vertex_degree_statistics(double*,int) const;
    /// This method computes the \f$f\f$-vector and \f$f^*\f$-vector of a given sheet (the final argument) of the complex, the first of which is the number of d-simplices for each d (d >= 0) and the second the dual of the first, storing the output in the method's first two arguments.
    void compute_fvector(std::vector<int>&,std::vector<int>&,int) const;
    /// This method computes the \f$h\f$-vector of a given sheet (the final argument) of the simplicial complex, after first computing the \f$f\f$-vector and \f$f^*\f$-vector upon which it depends; the output is written to the method's first argument.
    void compute_hvector(std::vector<int>&,int) const;
    /// This method calls the compute_graph method after allocating the memory for the offset array and results in the calculation of the 1-skeleton corresponding to a sheet of the complex, indicated by the method's final argument.
    void compute_graph(SYNARMOSMA::Graph*,int) const;
    /// This method calls the compute_graph method and stores the offset array used to deal with the problem of inactive events, so that the 1-skeleton of a sheet (the last argument) the complex is stored in the first argument.
    void compute_graph(SYNARMOSMA::Graph*,int*,int) const;
    /// This method computes the 1-skeleton of a sheet (the third argument) of the simplicial complex centred on a base vertex (the second argument) and with a combinatorial distance equal to Complex::topological_radius.
    void compute_graph(SYNARMOSMA::Graph*,int,int) const;
    /// This method computes the 1-skeleton of a sheet (the fourth argument) of the simplicial complex centred on a base vertex (the second argument) and with a combinatorial distance equal to the third argument.
    void compute_graph(SYNARMOSMA::Graph*,int,int,int) const;
    /// This method computes the instance of the Synarmosma library's Nexus class corresponding to a given sheet (the second argument) of the complex.
    void compute_global_nexus(SYNARMOSMA::Nexus*,int) const;
    /// This method computes the instance of the Synarmosma library's Nexus class corresponding to those simplices belonging to a given sheet (the final argument) and which contain the event whose index is the second argument of this method.
    void compute_local_nexus(SYNARMOSMA::Nexus*,int,int) const;
    /// This method computes the Nexus class from the Synarmosma library corresponding to all active d-simplices in the complex and then computes the homology which is stored in the Complex::H property. If the method's argument is true, it also computes the fundamental group of the Nexus and stores it in Complex::pi1. Finally this method computes the Complex::pseudomanifold, Complex::orientable and Complex::boundary properties.
    void compute_global_topology(bool);
    /// This method accepts as its first argument an event and then proceeds to compute the number of d-simplices (d > 0) that contain this event, appending this figure to the second argument which is the method's output.
    void simplex_membership(int,std::vector<int>&) const;
    /// This method computes the 1-skeleton of the sheet indicated by the method's argument as an instance of the Synarmosma library's Graph class and returns the chromatic_number method of this class.
    int chromatic_number(int) const;
    /// This method calculates the combinatorial distance - the minimal number of event to event hops - between the two events whose index in the Complex::events vector is given by the first two arguments, on the sheet indicated by the final argument.
    int combinatorial_distance(int,int,int) const;
    /// This method calls the compute_graph method for the sheet indicated by the argument and returns the output of the entwinement method of the Synarmosma library's Graph class.     
    double entwinement(int) const;
    /// This method computes the 1-skeleton of the sheet indicated by the method's argument as an instance of the Synarmosma library's Graph class and returns the cyclic_resistance method of this class.
    double cyclic_resistance(int) const;
    /// This method returns the maximum value of the degree - the number of neighbours - over the entire collection of events.
    int max_degree() const;
    /// This method returns the sum of the vertex_dimension() for each active event in the sheet indicated by the method's argument.
    int total_dimension(int) const;
    /// This method returns the sum of the number of d-simplices (d >= 0) belonging to the sheet indicated by the argument, weighted by the number of vertices in each simplex.
    int structural_index(int) const;
    /// This method returns the highest dimension D of this complex which has a D-simplex belonging to the sheet indicated by the argument. If there are no d-simplices belonging to the sheet for any non-negative d, the method returns -1.
    int dimension(int) const;
    /// This method returns the number of 1-simplices belonging to the sheet indicated by the second argument that are in the entourage of the event given by the first argument.
    int vertex_valence(int,int) const;
    /// This method returns the dimension of the highest-dimensional simplex belonging to the sheet indicated by the second argument and which contains the event specified by the method's first argument; if the event is inactive the method returns -1.
    int vertex_dimension(int,int) const;
    /// This method accepts two events as its arguments and returns the sum of active d-simplices that contain both events, weighted by the dimension.
    int weighted_entourage(int,int) const;
    /// This method computes the 1-skeleton of the sheet indicated by the method's argument as an instance of the Synarmosma library's Graph class and if the graph is connected (otherwise it returns zero), it returns the number of edges less the number of bridges in the graph.
    int cyclicity(int) const;
    /// This method returns the sum of the dimensional stress associated with each d-simplex (d > 0) that belongs to the sheet indicated by the method's second argument and contains the event whose index is the method's first argument.
    double dimensional_stress(int,int) const;
    /// This method measures the standard deviation of the simplicial dimensions on a given sheet (the third argument) of the vertices of a given d-simplex, with the first argument the dimension d and the second the index in the array Complex::simplices[d].
    double dimensional_stress(int,int,int) const;
    /// This method returns the percentage of 1-simplices belonging to the sheet indicated by the second argument whose two vertices have different topological dimensions, assuming they are greater than the method's first argument.
    double dimensional_frontier(int,int) const;
    /// This method computes the sum of the differences between the topological dimension of each event belonging to the sheet (indicated by the method's second argument) and the method's first argument (or the dimension of the complex if it is greater). The method then returns this sum divided by the number of active events. 
    double dimensional_uniformity(int,int) const;
    /// This method returns the circuit rank of a given sheet of the complex, i.e. \f$T = c - v + e\f$ where \f$c\f$ is the number of connected components, \f$v\f$ the number of events and \f$e\f$ the number of 1-simplices. The method's argument is the index of the sheet.
    int circuit_rank(int) const;
    /// This method returns the Euler characteristic of a given sheet (indicated by the argument) of the complex, i.e. the alternating sum of the elements of the \f$f\f$-vector, \f$\chi = \sum_{n=0}^D (-1)^n f_n\f$.
    int euler_characteristic(int) const;
    /// This method computes the number of connected components (the method's return value) of the sheet indicated by the final argument as well as the number of events in each component, stored in the method's first argument.
    int component_analysis(std::vector<int>&,int) const;
    /// This method computes the number of active d-simplices belonging to the sheet indicated by the second argument and which contain the event whose index is the first argument, multiplied by the dimension d. If the event is inactive the method returns -1. 
    int entourage(int,int) const;
    /// This method returns true if the 1-skeleton of events belonging to the sheet indicated by its argument is connected and false otherwise.
    bool connected(int) const;
    /// This method carries out a series of tests of the internal consistency of various Event and Simplex properties for those elements which belong to the sheet indicated by the method's argument. The method returns false if any problems are discovered.
    bool consistent(int) const;
    /// This method writes the simplicial inclusion relations for the complex as a directed graph in the DOT file format, where the method's first argument is the filename and the second the sheet to write out.
    void write_incastrature(const std::string&,int) const;
    /// This method writes to the console a summary of the topology of a given sheet (the method's argument) of the complex, including the simplices that each event belongs to and the number of simplices of a given dimension.
    void write_topology(int) const;
    /// This method returns false if there is an inactive event whose energy property has a non-zero value.
    bool energy_check() const;
    /// This method returns the sum of the energy property for every event belonging to the sheet whose index is the method's argument.
    double total_energy(int) const;
    /// This method carries out the energy diffusion among neighbouring active events in the complex - the method's first argument is a parameter that controls the rate of energy transfer, the second is a threshold for distinguishing energy conservation during diffusion.
    void energy_diffusion(double,double);
    /// This method carries out the energy diffusion among the complex's events using a chip-firing game algorithm, with the method's argument equal to the number of chips to be distributed among the active events.
    void energy_diffusion(int);
    /// This method computes the extent to which a given d-simplex (d is the first argument while the second is the simplex's index in Complex::simplices[d]) can be embedded in Euclidean space by calculating the number of non-zero eigenvalues of a matrix derived from the simplex's edge lengths. The simplex can only be embedded if it is entirely timelike or spacelike, otherwise the method returns -1.
    int simplex_embedding(int,int) const;
    /// This method returns the value of a simple Ising-like Hamiltonian based on the Simplex::parity property of the complex's d-simplices (d > 0) belonging to a given sheet (whose index is the method's final argument). The first argument is the strength of the magnetic field and the second determines whether it is ferromagnetic or not.
    double parity_hamiltonian(double,bool,int) const;
    /// This method writes the complex's 1-skeleton as a binary file whose name is the first argument with the sheet index the second argument. The file format is a list of the energy of all events belonging to the sheet followed by a pair of integers for each edge belonging to the sheet.
    void write_graph(const std::string&,int) const;
    /// This method returns the "fitness" of a distribution of the complex's events among a set of N processors, where N is the method's final argument. The first argument is the number of events assigned to each processor while the second is the assignment of each event to a processor. The method attempts to ensure that neighbouring events are assigned to the same processor while also have an equal distribution of events among the processors.
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
    /// This method writes the energy property of all the events belonging to the sheet indicated in the second argument to the method's first argument.
    void get_energy_values(std::vector<double>&,int) const;
    /// This method writes the deficiency property of all the active events belonging to the sheet indicated in the second argument to the method's first argument.
    void get_deficiency_values(std::vector<double>&,int) const;
    /// This method returns the inherited neighbours property of the event with index equal to the method's argument.
    void get_neighbours(int,std::set<int>&) const;
    /// This method returns the homology of the complex as a compact string, the method's argument.
    void get_homology(std::string&) const;
    /// This method returns a combinatorial presentation of the homotopy group of the complex as a compact string, the method's argument.
    void get_homotopy(std::string&) const;
    /// This method sets the field over which the complex's homology is calculated to the argument.
    void set_homology_field(SYNARMOSMA::Homology::Field);
    /// This method sets the method by which the complex's homology is calculated to the argument.
    void set_homology_method(SYNARMOSMA::Homology::Method);
    /// This method returns the field over which the complex's homology is calculated.
    SYNARMOSMA::Homology::Field get_homology_field() const;
    /// This method returns the method by which the complex's homology is calculated.
    SYNARMOSMA::Homology::Method get_homology_method() const;
    /// This method returns true if the element of Complex::events given by the method's argument is active, false otherwise.
    bool active_event(int,int) const;
    /// This method returns true if the element of Complex::events (for d = 0) or Complex::simplices (for d > 0) given by the method's second argument is active.
    bool active_simplex(int,int,int) const;
    /// This method returns the output of the dimension() method for a given sheet.
    int get_dimension(int) const;
    /// This method returns the inherited topological_dimension of the event with index equal to the method's first argument and sheet equal to the second argument.
    int get_event_dimension(int,int) const;
    /// This method returns the number of d-simplices belonging to the sheet indicated by the second argument.
    int get_cardinality(int,int) const;
    /// This method sets the first argument, an array of three double precision numbers, to the maximum, minimum and average degree of the events belonging to the sheet indicated by the second argument.    
    void get_vertex_degree_statistics(double*,int) const;
    /// This method calls the compute_fvector() method using its three arguments.     
    void get_fvector(std::vector<int>&,std::vector<int>&,int) const;
    /// This method returns the vertices of a d-simplex (d > 0) where the arguments are the dimension and index in Complex::simplices[d], while the final argument is an integer array of length 1+d to store the vertices.
    void get_simplex_vertices(int,int,int*) const;
    /// This method returns a string of a d-simplex's vertices (where d > 0) separated by a colon in ascending order.  
    std::string get_simplex_key(int,int) const;
    /// This method returns the index of the 1-simplex whose vertex set is the method's argument; if the edge does not exist the method returns -1. 
    int get_edge_index(const std::set<int>&) const;
    /// This method returns the output of the cyclicity() method for a given sheet.
    int get_cyclicity(int) const;
    /// This method returns the output of the circuit_rank() method for a given sheet.
    int get_circuit_rank(int) const;
    /// This method returns whether or not the simplicial complex is orientable, assuming it is a combinatorial pseudomanifold.
    bool is_orientable() const;
    /// This method returns the Euler characteristic of a given sheet of the complex.
    int get_euler_characteristic(int) const;
    /// This method returns the Event::obliquity property of the event with index equal to the method's argument.
    double get_event_obliquity(int) const;
    /// This method returns the inherited energy property of the event with index equal to the method's argument.
    double get_event_energy(int) const;
    /// This method returns the Event::deficiency property of the event with index equal to the method's argument.
    double get_event_deficiency(int) const;
    /// This method returns the Event::entwinement property for a particular sheet (the second argument) of the event with index equal to the method's first argument.
    double get_event_entwinement(int,int) const;
    /// This method returns the number of elements in the Complex::events vector.
    int get_events() const;
    /// This method returns the output of the total_energy() method for the given sheet.
    double get_total_energy(int) const;
    /// This method returns the cardinality of the entourage property of the event (if the first argument is zero) or simplex (if it is greater than zero) whose index is specified by the second argument. 
    int get_entourage_cardinality(int,int) const; 
    /// This method returns the incept property of an Event or Simplex - depending on the dimension (the first argument) - of the element of Complex::events or Complex::simplices[d] with index value equal to the second argument.
    int get_incept(int,int) const;
    /// This method sets the first element of the argument to the Complex::pseudomanifold property and second element to the value of the Complex::boundary property.
    void is_pseudomanifold(std::pair<bool,bool>&) const;
    friend class Spacetime;
  };

  inline void Complex::get_neighbours(int n,std::set<int>& S) const
  {
    S = events[n].neighbours;
  }

  inline void Complex::get_homology(std::string& s) const 
  {
    s = H->write();
  }

  inline void Complex::get_homotopy(std::string& s) const 
  {
    s = pi1->write();
  }

  inline void Complex::set_homology_field(SYNARMOSMA::Homology::Field F) 
  {
    H->set_field(F);
  }

  inline void Complex::set_homology_method(SYNARMOSMA::Homology::Method M) 
  {
    H->set_method(M);
  }

  inline SYNARMOSMA::Homology::Field Complex::get_homology_field() const 
  {
    return H->get_field();
  }

  inline SYNARMOSMA::Homology::Method Complex::get_homology_method() const 
  {
    return H->get_method();
  }

  inline bool Complex::active_event(int n,int sheet) const 
  {
    if (sheet == -1) return (!events[n].ubiquity.empty());

    return (events[n].ubiquity.count(sheet) > 0);
  }

  inline bool Complex::active_simplex(int d,int n,int sheet) const 
  {
    if (d == 0) return active_event(n,sheet);
    if (sheet == -1) return (!simplices[d][n].ubiquity.empty());
    return (simplices[d][n].ubiquity.count(sheet) > 0); 
  }

  inline int Complex::get_dimension(int sheet) const
  {
    return dimension(sheet);
  }

  inline int Complex::get_event_dimension(int n,int sheet) const
  {
    return vertex_dimension(n,sheet);
  }

  inline int Complex::get_cardinality(int d,int sheet) const 
  {
    return cardinality(d,sheet);
  }

  inline void Complex::get_vertex_degree_statistics(double* output,int sheet) const 
  {
    vertex_degree_statistics(output,sheet);
  }

  inline void Complex::get_fvector(std::vector<int>& f,std::vector<int>& fstar,int sheet) const
  {
    compute_fvector(f,fstar,sheet);
  }

  inline void Complex::get_simplex_vertices(int d,int n,int* vx) const 
  {
    simplices[d][n].get_vertices(vx);
  }

  inline std::string Complex::get_simplex_key(int d,int n) const 
  {
    return SYNARMOSMA::make_key(simplices[d][n].vertices);
  }

  inline int Complex::get_cyclicity(int sheet) const
  {
    return cyclicity(sheet);
  }

  inline int Complex::get_circuit_rank(int sheet) const
  {
    return circuit_rank(sheet);
  }

  inline bool Complex::is_orientable() const 
  {
    return orientable;
  }

  inline int Complex::get_euler_characteristic(int sheet) const
  {
    return euler_characteristic(sheet);
  }

  inline double Complex::get_event_obliquity(int n) const 
  {
    return events[n].obliquity;
  }

  inline double Complex::get_event_energy(int n) const 
  {
    return events[n].get_energy();
  }

  inline double Complex::get_event_deficiency(int n) const 
  {
    return events[n].deficiency;
  }

  inline double Complex::get_event_entwinement(int n,int sheet) const 
  {
    return events[n].entwinement[sheet];
  }

  inline int Complex::get_events() const
  {
    return (signed) events.size();
  }

  inline double Complex::get_total_energy(int sheet) const 
  {
    return total_energy(sheet);
  }

  inline int Complex::get_entourage_cardinality(int d,int n) const
  {
    int output = (d == 0) ? (signed) events[n].neighbours.size() : (signed) simplices[d][n].entourage.size(); 
    return output;
  }

  inline int Complex::get_incept(int d,int n) const 
  {
    int output = (d == 0) ? events[n].incept : simplices[d][n].incept; 
    return output;
  }

  inline bool Complex::edge_exists(int u,int v,int sheet) const
  {
    if (u == v || u < 0 || v < 0) throw std::invalid_argument("Illegal vertex values in the Complex::edge_exists method!");

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

  inline int Complex::get_edge_index(const std::set<int>& vx) const 
  {
    if (vx.size() != 2) throw std::invalid_argument("Illegal set cardinality in the Complex::get_edge_index method!");
    
    SYNARMOSMA::hash_map::const_iterator qt = index_table[1].find(vx); 
    if (qt == index_table[1].end()) return -1;
    return qt->second;
  }

  inline void Complex::is_pseudomanifold(std::pair<bool,bool>& output) const 
  {
    output.first = pseudomanifold;
    output.second = boundary;
  }

  inline void Complex::compute_graph(SYNARMOSMA::Graph* G,int sheet) const
  {
    if (events.empty()) return;

    int offset[events.size()];
    compute_graph(G,offset,sheet);
  }

  inline void Complex::compute_graph(SYNARMOSMA::Graph* G,int base,int sheet) const
  {
    compute_graph(G,base,Complex::topological_radius,sheet);
  }
}
#endif
