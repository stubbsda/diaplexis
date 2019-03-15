#include "event.h"
#include "simplex.h"

#include "synarmosma/homology.h"
#include "synarmosma/homotopy.h"

#ifndef _complexh
#define _complexh

namespace DIAPLEXIS {
  class Complex {
   private:
    bool pseudomanifold = false;
    bool boundary = false;
    bool orientable = false;
    std::vector<Event> events;
    std::vector<Simplex>* simplices;
    SYNARMOSMA::hash_map* index_table;
    SYNARMOSMA::Homology* H;
    SYNARMOSMA::Homotopy* pi1;
    SYNARMOSMA::Random* RND;

    // Maximum combinatorial dimension of the spacetime
    static const int ND = 10;
    // The combinatorial size of the subgraph used to compute
    // the topological entwinement at a vertex
    static const int topological_radius = 4;

    void allocate();
    void inversion();
    void compute_neighbours();
    void compute_entourages();
    bool simplex_addition(int,int,int);
    bool simplex_addition(const std::set<int>&,int = -1);
    bool simplex_addition(const std::set<int>&,std::set<int>&);
    void simplex_deletion(int,int);
    void compute_geometric_dependency(const std::set<int>&);
    void compute_topological_dependency(const std::set<int>&);
    void simplicial_implication();
    void simplicial_implication(int);
    void compute_simplicial_dimension();
    double set_logical_atoms(int);
    double logical_energy(int) const;
    bool logical_conformity(int) const;
    void compute_simplex_energy(int,int);
    void compute_simplex_parity(int,int);
    bool edge_parity_mutation(int);
    bool edge_parity_mutation(int,int);
    void recompute_parity(int);
    void recompute_parity(const std::set<int>&);
    void compute_parity();
    void determine_flexible_edges(std::vector<int>&);
    bool edge_exists(int,int) const;
    int cardinality(int) const;
    int cardinality_safe(int) const;
    void compute_graph(SYNARMOSMA::Graph*,int) const;
    void compute_degree_distribution(bool) const;
    void compute_connectivity_distribution(bool) const;
    std::pair<double,double> random_walk() const;
    void vertex_degree_statistics(double*) const;
    void compute_fvector(std::vector<int>&,std::vector<int>&) const;
    void compute_hvector(std::vector<int>&) const;
    void compute_graph(SYNARMOSMA::Graph*,int,int) const;
    void compute_graph(SYNARMOSMA::Graph*) const;
    void compute_graph(SYNARMOSMA::Graph*,int*) const;
    void compute_global_nexus(SYNARMOSMA::Nexus*) const;
    void compute_local_nexus(SYNARMOSMA::Nexus*,int) const;
    void compute_global_topology(bool);
    void simplex_membership(int,std::vector<int>&) const;
    int chromatic_number() const;
    double dimensional_stress(int,int) const;
    int combinatorial_distance(int,int) const;
    double entwinement() const;
    double cyclic_resistance() const;
    int max_degree() const;
    int total_dimension() const;
    int structural_index() const;
    int dimension() const;
    int vertex_valence(int) const;
    int vertex_dimension(int) const;
    int weighted_entourage(int,int) const;
    int cyclicity() const;
    double dimensional_frontier(int) const;
    double dimensional_uniformity(int) const;
    int circuit_rank() const;
    int euler_characteristic() const;
    int component_analysis(std::vector<int>&) const;
    int entourage(int) const;
    bool connected() const;
    bool consistent() const;
    void write_incastrature(const std::string&) const;
    void write_topology() const;
    bool energy_check() const;
    double total_energy() const;
    void energy_diffusion(double);
    void energy_diffusion(int);
    void simplicial_implication(int) const;
    int simplex_embedding(int,int) const;
    double dimensional_stress(int) const;
    void compute_delta(std::set<int>&);

    double parity_hamiltonian(double,bool) const;
    void write_graph(const std::string&) const;
    inline double distribution_fitness(int*,const std::vector<int>&,int) const;

   public:
    Complex();
    Complex(const Complex&);
    ~Complex();
    void distribute(int) const; 
    int serialize(std::ofstream&) const;
    int deserialize(std::ifstream&);
    void clear();
    void get_edge_topology(std::vector<std::set<int> >&) const;
    bool active_element(int,int) const;
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
    inline bool active_event(int n) const {return events[n].active;};
    inline bool active_simplex(int d,int n) const {return simplices[d][n].active;};
    inline int get_dimension() const {return dimension();};
    inline int get_event_dimension(int n) const {return vertex_dimension(n);};
    inline int get_cardinality(int d) const {return cardinality_safe(d);};
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

  inline int Complex::cardinality(int d) const
  {
    int i,n = 0;
    if (d == 0) {
      const int M = (signed) events.size();
      for(i=0; i<M; ++i) {
        if (!events[i].active) continue;
        n++;
      }
    }
    else {
      const int M = (signed) simplices[d].size();
      for(i=0; i<M; ++i) {
        if (!simplices[d][i].active) continue;
        n++;
      }
    }
    return n;
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
