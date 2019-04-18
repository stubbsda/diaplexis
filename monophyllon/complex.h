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
    /// This method transforms the complex by converting it to a graph (one-dimensional complex) whose topology is the inverse of the current topology, calling the Complex::compute_entourages method. 
    void inversion();
    /// This method computes the inherited neighbours property of the elements of the Complex::events vector.
    void compute_neighbours();
    /// This method calculates the value of the inherited entourage property for the elements of the Complex::events and Complex::simplices arrays.
    void compute_entourages();
    bool simplex_addition(int,int,int);
    bool simplex_addition(const std::set<int>&,int = -1);
    bool simplex_addition(const std::set<int>&,std::set<int>&);
    /// This method deletes a d-simplex by setting its active property to false, where d is the first argument and the second argument is the index of the simplex in Complex::simplices[d]. The method then recursively deletes all the higher-dimensional simplices which depend on this simplex. 
    void simplex_deletion(int,int);
    void compute_modified_vertices();
    void compute_dependent_simplices(const std::set<int>&);
    void simplicial_implication();
    void simplicial_implication(int);
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
    bool edge_parity_mutation(int);
    bool edge_parity_mutation(int,int);
    void recompute_parity(int);
    void recompute_parity(const std::set<int>&);
    /// This method computes the Simplex::parity property for the entire collection of d-simplices (d > 1) in the complex, successively calling the compute_simplex_parity method on each d-simplex whose Simplex::modified property is true. 
    void compute_parity();
    /// This method loops over all the 1-simplices in the simplicial complex - if a 1-simplex is active and at least one of its vertices has a non-zero energy, the corresponding element of the method's argument is given the value 1, otherwise 0. 
    void determine_flexible_edges(std::vector<int>&);
    /// This method returns true if there is an active edge connecting the two vertices that are the method's arguments.
    bool edge_exists(int,int) const;
    /// This method returns the number of active d-simplices in the complex, where d is the method's argument.
    int cardinality(int) const;
    void compute_degree_distribution(bool) const;
    void compute_connectivity_distribution(bool) const;
    /// This method computes the 1-skeleton of the complex using the compute_graph method and then returns the output of the random_walk method of the Graph class in the Synarmosma library.
    std::pair<double,double> random_walk() const;
    void vertex_degree_statistics(double*) const;
    void compute_fvector(std::vector<int>&,std::vector<int>&) const;
    void compute_hvector(std::vector<int>&) const;
    /// This method calls the compute_graph method after allocating the memory for the offset array and results in the calculation of the 1-skeleton corresponding to the complex.
    void compute_graph(SYNARMOSMA::Graph*) const;
    /// This method calls the compute_graph method and stores the offset array used to deal with the problem of inactive events, so that the 1-skeleton of the complex is stored in the first argument.
    void compute_graph(SYNARMOSMA::Graph*,int*) const;
    /// This method computes the 1-skeleton of the simplicial complex centred on a base vertex (the second argument) and with a combinatorial distance equal to Complex::topological_radius.
    void compute_graph(SYNARMOSMA::Graph* G,int) const
    /// This method computes the 1-skeleton of the simplicial complex centred on a base vertex (the second argument) and with a combinatorial distance equal to the third argument.
    void compute_graph(SYNARMOSMA::Graph*,int,int) const;
    /// This method computes the instance of the Synarmosma library's Nexus class corresponding to this complex.
    void compute_global_nexus(SYNARMOSMA::Nexus*) const;
    /// This method computes the instance of the Synarmosma library's Nexus class corresponding to those simplices which contain the event whose index is the second argument of this method.
    void compute_local_nexus(SYNARMOSMA::Nexus*,int) const;
    /// This method allocates an instance of the Synarmosma library's Nexus class, calls the compute_global_nexus method and then computes the Complex::H and (if the method's argument is true) Complex::pi1 properties, along with the related Complex::pseudomanifold, Complex::boundary and Complex::orientable properties. 
    void compute_global_topology(bool);
    /// This method accepts as its first argument an event and then proceeds to compute the number of d-simplices (d > 0) that contain this event, appending this figure to the second argument which is the method's output.
    void simplex_membership(int,std::vector<int>&) const;
    /// This method computes the 1-skeleton of the complex as an instance of the Synarmosma library's Graph class and returns the chromatic_number method of this class.
    int chromatic_number() const;
    /// This method calculates the combinatorial distance - the minimal number of event to event hops - between the two events whose index in the Complex::events vector is given by the arguments
    int combinatorial_distance(int,int) const;
    /// This method calls the compute_graph method and returns the output of the entwinement method of the Synarmosma library's Graph class. 
    double entwinement() const;
    double cyclic_resistance() const;
    /// This method returns the maximum value of the degree - the number of neighbours - over the entire collection of events.
    int max_degree() const;
    /// This method returns the sum of the vertex_dimension() for each active event in the complex.
    int total_dimension() const;
    /// This method returns the sum of the number of active d-simplices (d >= 0), weighted by the number of vertices in each simplex.
    int structural_index() const;
    /// This method returns the highest dimension D of this complex which has an active D-simplex. If there are no active d-simplices for any non-negative d, the method returns -1.
    int dimension() const;
    /// This method returns the number of active 1-simplices that are in the entourage of the event given by the argument.
    int vertex_valence(int) const;
    /// This method returns the dimension of the highest-dimensional simplex in the complex which contains the event specified by the method's argument; if the event is inactive the method returns -1.
    int vertex_dimension(int) const;
    int weighted_entourage(int,int) const;
    int cyclicity() const;
    /// This method returns the sum of the dimensional stress associated with each d-simplex (d > 0) that is active and contains the event whose index is the method's argument.
    double dimensional_stress(int) const;
    /// This method measures the standard deviation of the simplicial dimensions of the vertices of a given d-simplex, with the first argument the dimension d and the second the index in the array Complex::simplices[d].
    double dimensional_stress(int,int) const;
    /// This method returns the percentage of active 1-simplices in the complex whose two vertices have different topological dimensions, assuming they are greater than the method's argument.
    double dimensional_frontier(int) const;
    /// This method computes the sum of the differences between the topological dimension of each active event and the method's argument (or the dimension of the complex if it is greater). The method then returns this sum divided by the number of active events. 
    double dimensional_uniformity(int) const;
    /// This method returns the circuit rank of the complex, i.e. \f$T = c - v + e\f$ where \f$c\f$ is the number of connected components, \f$v\f$ the number of events and \f$e\f$ the number of 1-simplices.
    int circuit_rank() const;
    /// This method returns the Euler characteristic of the complex, i.e. the alternating sum of the elements of the \f$f\f$-vector, \f$\chi = \sum_{n=0}^D (-1)^n f_n\f$.
    int euler_characteristic() const;
    /// This method computes the number of connected components of the complex (the method's return value) as well as the number of events in each component, stored in the method's argument.
    int component_analysis(std::vector<int>&) const;
    /// This method computes the number of active d-simplices which contain the event whose index is the argument, multiplied by the dimension d. If the event is inactive the method returns -1. 
    int entourage(int) const;
    /// This method returns true if the complex is connected and false otherwise. 
    bool connected() const;
    /// This method carries out a series of tests of the internal consistency of various Event and Simplex properties and returns false if any problems are discovered.
    bool consistent() const;
    /// This method writes the simplicial inclusion relations for the complex as a directed graph in the DOT file format, where the method's argument is the filename.
    void write_incastrature(const std::string&) const;
    /// This method writes to the screen a summary of the topology of the complex, including the simplices that each event belongs to and the number of simplices of a given dimension.
    void write_topology() const;
    /// This method returns false if there is an inactive event whose energy property has a non-zero value.
    bool energy_check() const;
    /// This method returns the sum of the energy property for every active event in the complex.
    double total_energy() const;
    void energy_diffusion(double);
    void energy_diffusion(int);
    int simplex_embedding(int,int) const;
    double parity_hamiltonian(double,bool) const;
    /// This method writes the complex's 1-skeleton as a binary file whose name is the argument. The file format is a list of the energy of all active events followed by a pair of integers for each active edge.
    void write_graph(const std::string&) const;
    inline double distribution_fitness(int*,const std::vector<int>&,int) const;

   public:
    /// The default constructor which calls the allocate() method.
    Complex();
    /// The destructor which frees the memory from the instance's properties (Complex::simplices, Complex::RND etc.).
    ~Complex();
    void distribute(int) const; 
    /// This method writes the instance properties to a binary disk file and returns the number of bytes written to the file.    
    int serialize(std::ofstream&) const;
    /// This method calls the clear() method on the instance and then reads the properties from a binary disk file and returns the number of bytes read.
    int deserialize(std::ifstream&);
    /// This method clears all of the instance's extended properties and sets the scalar properties to their default value.
    void clear();
    /// This method assembles the vector of neighbour sets from each active event in the complex.
    void get_edge_topology(std::vector<std::set<int> >&) const;
    void write_vertex_data(int) const;
    void get_energy_values(std::vector<double>&) const;
    void get_deficiency_values(std::vector<double>&) const;
    inline std::set<int> get_neighbours(int n) const {return events[n].neighbours;};
    inline void get_homotopy(std::string& s) const {s = pi1->write();};
    inline void get_homology(std::string& s) const {s = H->write();};
    inline void set_homology_field(SYNARMOSMA::Homology::Field F) {H->set_field(F);};
    inline void set_homology_method(SYNARMOSMA::Homology::Method M) {H->set_method(M);};
    inline SYNARMOSMA::Homology::Field get_homology_field() const {return H->get_field();};
    inline SYNARMOSMA::Homology::Method get_homology_method() const {return H->get_method();};
    /// This method returns true if the element of Complex::events given by the method's argument is active, false otherwise.
    inline bool active_event(int n) const {return events[n].active;};
    /// This method returns true if the element of Complex::events (for d = 0) or Complex::simplices (for d > 0) given by the method's second argument is active.
    inline bool active_simplex(int d,int n) const {bool output = (d == 0) ? events[n].active : simplices[d][n].active; return output;};
    /// This method returns the value of the dimension() method.
    inline int get_dimension() const {return dimension();};
    inline int get_event_dimension(int n) const {return vertex_dimension(n);};
    inline int get_cardinality(int d) const {return cardinality(d);};
    inline void get_vertex_degree_statistics(double* output) const {vertex_degree_statistics(output);};
    inline void get_fvector(std::vector<int>& f,std::vector<int>& fstar) const {compute_fvector(f,fstar);};
    inline void get_simplex_vertices(int d,int n,int* vx) const {simplices[d][n].get_vertices(vx);};
    inline std::string get_simplex_key(int d,int n) const {return SYNARMOSMA::make_key(simplices[d][n].vertices);};
    inline int get_edge_index(const std::set<int>& vx) const {SYNARMOSMA::hash_map::const_iterator qt = index_table[1].find(vx); return qt->second;};
    inline int get_cyclicity() const {return cyclicity();};
    inline int get_circuit_rank() const {return circuit_rank();};
    inline bool is_orientable() const {return orientable;};
    inline int get_euler_characteristic() const {return euler_characteristic();};
    inline double get_event_obliquity(int n) const {return events[n].obliquity;};
    inline double get_event_energy(int n) const {return events[n].get_energy();};
    inline double get_event_deficiency(int n) const {return events[n].deficiency;};
    inline double get_event_entwinement(int n) const {return events[n].entwinement;};
    inline int get_events() const {return (signed) events.size();};
    inline double get_total_energy() const {return total_energy();}; 
    inline int get_entourage_cardinality(int d,int n) const {int output = (d == 0) ? (signed) events[n].neighbours.size() : (signed) simplices[d][n].entourage.size(); return output;}; 
    inline int get_incept(int d,int n) const {int output = (d == 0) ? events[n].incept : simplices[d][n].incept; return output;};
    bool is_pseudomanifold(bool*) const;
    friend class Spacetime;
  };

  inline bool Complex::edge_exists(int u,int v) const
  {
    std::set<int> vx;
    vx.insert(u); vx.insert(v);
    SYNARMOSMA::hash_map::const_iterator qt = index_table[1].find(vx);
    if (qt == index_table[1].end()) return false;
    if (!simplices[1][qt->second].active) return false;
    return true;
  }

  inline bool Complex::is_pseudomanifold(bool* bdry) const 
  {
    *bdry = boundary;
    return pseudomanifold;
  }

  inline void Complex::compute_graph(SYNARMOSMA::Graph* G,int base) const
  {
    compute_graph(G,base,Complex::topological_radius);
  }

  inline double Complex::distribution_fitness(int* volume,const std::vector<int>& affinity,int nprocs) const
  {
    int i,sum = 0,bcount = 0;
    double mu,sigma = 0.0;
    std::set<int>::const_iterator it;
    const int nv = (signed) events.size();

    for(i=0; i<nprocs; ++i) {
      sum += volume[i];
    }
    mu = double(sum)/double(nprocs);
    for(i=0; i<nprocs; ++i) {
      sigma += (volume[i] - mu)*(volume[i] - mu);
    }
    sigma = std::sqrt(sigma/double(nprocs));
    for(i=0; i<nv; ++i) {
      if (!events[i].active) continue;
      for(it=events[i].neighbours.begin(); it!=events[i].neighbours.end(); ++it) {
        if (*it < i) continue;
        if (affinity[*it] != affinity[i]) bcount++;
      }
    }
    return sigma + double(bcount);
  }
}
#endif
