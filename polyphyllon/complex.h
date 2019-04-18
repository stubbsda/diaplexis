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
    void inversion(int);
    void compute_neighbours();
    void compute_entourages(int);
    bool simplex_addition(int,int,const std::set<int>&,int = -1);
    bool simplex_addition(const std::set<int>&,const std::set<int>&,int = -1);
    void simplex_deletion(int,int,int);
    void compute_modified_vertices();
    void compute_dependent_simplices(const std::set<int>&);
    void simplicial_implication(int);
    void simplicial_implication(int,int);
    void compute_simplicial_dimension();
    double set_logical_atoms(int);
    double logical_energy(int) const;
    bool logical_conformity(int) const;
    void compute_simplex_energy(int,int);
    void compute_simplex_parity(int,int);
    bool edge_parity_mutation(int,int);
    bool edge_parity_mutation(int,int,int);
    void recompute_parity(int);
    void recompute_parity(const std::set<int>&);
    void compute_parity();
    void determine_flexible_edges(std::vector<int>&) const;
    void write_topology(int) const;
    void write_incastrature(const std::string&,int) const;
    bool edge_exists(int,int,int) const;
    int cardinality(int,int) const;
    void compute_degree_distribution(bool,int) const;
    void compute_connectivity_distribution(bool,int) const;
    std::pair<double,double> random_walk(int) const;
    void vertex_degree_statistics(double*,int) const;
    void compute_fvector(std::vector<int>&,std::vector<int>&,int) const;
    void compute_hvector(std::vector<int>&,int) const;
    void compute_graph(SYNARMOSMA::Graph*,int) const;
    void compute_graph(SYNARMOSMA::Graph*,int,int) const;
    void compute_graph(SYNARMOSMA::Graph*,int,int,int) const;
    void compute_graph(SYNARMOSMA::Graph*,int*,int) const;
    void compute_global_nexus(SYNARMOSMA::Nexus*,int) const;
    void compute_local_nexus(SYNARMOSMA::Nexus*,int,int) const;
    void compute_global_topology(bool);
    void simplex_membership(int,std::vector<int>&) const;
    int chromatic_number(int) const;
    double dimensional_stress(int,int) const;
    int combinatorial_distance(int,int,int) const;
    double entwinement(int) const;
    double cyclic_resistance(int) const;
    int max_degree() const;
    int total_dimension(int) const;
    int structural_index(int) const;
    int dimension(int) const;
    int vertex_valence(int,int) const;
    int vertex_dimension(int,int) const;
    int weighted_entourage(int,int) const;
    int cyclicity(int) const;
    double dimensional_frontier(int,int) const;
    double dimensional_uniformity(int,int) const;
    bool active_simplex(int,int,int) const;
    int circuit_rank(int) const;
    int euler_characteristic(int) const;
    int component_analysis(std::vector<int>&,int) const;
    int entourage(int,int) const;
    bool connected(int) const;
    bool consistent(int) const;
    void write_incastrature(const std::string&) const;
    void write_topology() const;
    bool energy_check() const;
    double total_energy(int) const;
    void energy_diffusion(double);
    void energy_diffusion(int);
    void simplicial_implication(int,int) const;
    int simplex_embedding(int,int) const;
    double dimensional_stress(int,int,int) const;
    void compute_delta(std::set<int>&);
    double parity_hamiltonian(double,bool,int) const;
    void write_graph(const std::string&,int) const;
    double distribution_fitness(int*,const std::vector<int>&,int) const;

   public:
    Complex();
    ~Complex();
    void distribute(int) const; 
    int serialize(std::ofstream&) const;
    int deserialize(std::ifstream&);
    void clear();
    void get_edge_topology(std::vector<std::set<int> >&) const;
    bool active_element(int,int) const;
    void write_vertex_data(int) const;
    void get_energy_values(std::vector<double>&,int) const;
    void get_deficiency_values(std::vector<double>&,int) const;
    inline std::set<int> get_neighbours(int n) const {return events[n].neighbours;};
    inline void get_homotopy(std::string& s) const {s = pi1->write();};
    inline void get_homology(std::string& s) const {s = H->write();};
    inline void set_homology_field(SYNARMOSMA::Homology::Field F) {H->set_field(F);};
    inline void set_homology_method(SYNARMOSMA::Homology::Method M) {H->set_method(M);};
    inline SYNARMOSMA::Homology::Field get_homology_field() const {return H->get_field();};
    inline SYNARMOSMA::Homology::Method get_homology_method() const {return H->get_method();};
    inline int get_dimension(int sheet) const {return dimension(sheet);};
    inline int get_cardinality(int d,int sheet) const {return cardinality(d,sheet);};
    inline int get_event_dimension(int n,int sheet) const {return vertex_dimension(n,sheet);};
    inline double get_event_obliquity(int n) const {return events[n].obliquity;};
    inline double get_event_energy(int n) const {return events[n].get_energy();};
    inline double get_event_deficiency(int n) const {return events[n].deficiency;};
    inline double get_event_entwinement(int n,int sheet) const {return events[n].entwinement[sheet];};
    inline std::string get_simplex_key(int d,int n) const {return SYNARMOSMA::make_key(simplices[d][n].vertices);};
    inline void get_simplex_vertices(int d,int n,int* vx) const {simplices[d][n].get_vertices(vx);};
    inline void get_vertex_degree_statistics(double* output,int sheet) const {vertex_degree_statistics(output,sheet);};
    inline void get_fvector(std::vector<int>& f,std::vector<int>& fstar,int sheet) const {compute_fvector(f,fstar,sheet);};
    inline int get_edge_index(const std::set<int>& vx) const {SYNARMOSMA::hash_map::const_iterator qt = index_table[1].find(vx); return qt->second;};
    inline int get_cyclicity(int sheet) const {return cyclicity(sheet);};
    inline int get_circuit_rank(int sheet) const {return circuit_rank(sheet);};
    inline bool is_orientable() const {return orientable;};
    inline double get_total_energy(int sheet) const {return total_energy(sheet);};
    inline int get_euler_characteristic(int sheet) const {return euler_characteristic(sheet);};
    inline int get_events() const {return (signed) events.size();};
    inline int get_entourage_cardinality(int d,int n) const {int output = (d == 0) ? (signed) events[n].neighbours.size() : (signed) simplices[d][n].entourage.size(); return output;}; 
    inline int get_incept(int d,int n) const {int output = (d == 0) ? events[n].incept : simplices[d][n].incept; return output;};
    inline bool is_pseudomanifold(bool*) const;
    friend class Spacetime;
  };

  inline bool Complex::edge_exists(int u,int v,int sheet) const
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

  inline bool Complex::is_pseudomanifold(bool* bdry) const 
  {
    *bdry = boundary;
    return pseudomanifold;
  }

  inline void Complex::compute_graph(SYNARMOSMA::Graph* G,int base,int sheet) const
  {
    compute_graph(G,base,Complex::topological_radius,sheet);
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
      if (!events[i].active()) continue;
      for(it=events[i].neighbours.begin(); it!=events[i].neighbours.end(); ++it) {
        if (*it < i) continue;
        if (affinity[*it] != affinity[i]) bcount++;
      }
    }
    return sigma + double(bcount);
  }
}
#endif
