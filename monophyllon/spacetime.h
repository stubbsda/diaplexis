// For XML processing...
#include <pugixml.hpp>

#include <boost/timer/timer.hpp>
#include <boost/filesystem.hpp>

#include "synarmosma/geometry.h"
#include "synarmosma/homology.h"
#include "synarmosma/homotopy.h"
#include "synarmosma/directed_graph.h"

#include "event.h"
#include "simplex.h"

#ifndef _spacetimeh
#define _spacetimeh

namespace DIAPLEXIS {
  typedef struct {
    std::vector<Event> events;
    std::vector<Simplex>* simplices;
    SYNARMOSMA::hash_map* index_table;
    SYNARMOSMA::Geometry geometry;
  } Plegma;

  class Spacetime {
   protected:
    enum Initial_Topology
    {
        random,
        monoplex,
        cartesian,
        singleton,
        diskfile
    };

    enum Geometry_Solver
    {
        minimal,
        mechanical,
        evolutionary,
        annealing,
        simplex
    };

    enum class Integrator
    {
        euler,
        rk4
    };

    enum class Hyphansis
    {
        dynamic,
        musical
    };

    // The main (variable) properties of the Spacetime class
    int iterations;
    int system_size;
    int nactive;
    double error;
    double global_deficiency;
    double topology_delta;
    double geometry_delta;
    double energy_delta;
    Plegma anterior;
    std::vector<Simplex>* simplices;
    SYNARMOSMA::hash_map* index_table;
    std::vector<Event> events;
    std::string ops;
    std::set<int> vx_delta;
    SYNARMOSMA::Random* RND;
    SYNARMOSMA::Geometry* geometry;
    SYNARMOSMA::Homology* H;
    SYNARMOSMA::Homotopy* pi;
    bool pseudomanifold;
    bool boundary;
    bool orientable;
    std::vector<int> flexible_edge;

    // The global parameters
    Initial_Topology initial_state;
    Initial_Topology original_state;
    Hyphansis weaving;
    int initial_size;
    int max_iter;
    int initial_dim;
    boost::posix_time::ptime start_time;
    std::string date_string;
    std::string pid_string;
    std::string state_file;
    std::string log_file;
    std::string input_file;
    double edge_probability;
    double edge_reorientability;
    bool diskless;
    bool perturb_topology;
    bool perturb_geometry;
    bool perturb_energy;
    bool superposable;
    bool compressible;
    bool converged;
    bool high_memory;
    bool instrument_convergence;
    int checkpoint_frequency;
    std::string hyphansis_file;
    std::string hyphansis_score;

    // Now the parameters associated with the
    // geometry solver
    Geometry_Solver solver;
    Integrator engine;
    double geometry_tolerance;
    double thermalization;
    double thermal_variance;
    double simplex_alpha;
    double simplex_gamma;
    double simplex_rho;
    double simplex_sigma;
    double spring_constant;
    double repulsion_constant;
    double damping_constant;
    double step_size;
    double edge_flexibility_threshold;
    bool cgradient_refinement;
    int solver_its;
    int ngenerations;
    int pool_size;
    int njousts;
    int max_int_steps;
    int max_CG_steps;
    int max_LS_steps;
    int thermal_sweep;
    int annealing_steps;

    // Stuff for the implicative/explicative operators:
    static const int N_EXP = 11;
    static const int N_IMP = 8;
    static const std::string EXP_OP[N_EXP];
    static const std::string IMP_OP[N_IMP];
    // The combinatorial size of the subgraph used to compute
    // the topological entwinement at a vertex
    static const int topological_radius = 4;
    // Maximum combinatorial dimension of the spacetime
    static const int ND = 10;

    // Error tolerance for convergence
    static const double convergence_threshold;
    // The initial temperature for cosmic annealing
    static const double T_zero;
    // The rate of thermal decay in the annealing schedule
    static const double kappa;
    // The coupling constant between the topological-geometric torsion 
    // and the energy
    static const double Lambda;

    bool edge_exists(int,int) const;
    int cardinality(int) const;
    int cardinality_safe(int) const;
    void compute_graph(SYNARMOSMA::Graph*,int) const;
    void compute_degree_distribution(bool) const;
    void compute_connectivity_distribution() const;
    void random_walk(double*,double*) const;
    void arclength_statistics(double*) const;
    void vertex_degree_statistics(double*) const;
    void compute_fvector(std::vector<int>&,std::vector<int>&) const;
    void compute_hvector(std::vector<int>&) const;
    void compute_graph(SYNARMOSMA::Graph*,int,int) const;
    void compute_graph(SYNARMOSMA::Graph*) const;
    void compute_graph(SYNARMOSMA::Graph*,int*) const;
    void compute_causal_graph(SYNARMOSMA::Directed_Graph*,int) const;
    void compute_global_nexus(SYNARMOSMA::Nexus*) const;
    void compute_local_nexus(SYNARMOSMA::Nexus*,int) const;
    void simplex_membership(int,std::vector<int>&) const;
    int chromatic_number() const;
    double dimensional_stress(int,int) const;
    int combinatorial_distance(int,int) const;
    double entwinement() const;
    double cyclic_resistance() const;
    int max_degree() const;
    bool delaunay() const;
    int total_dimension() const;
    int structural_index() const;
    int dimension() const;
    int vertex_valence(int) const;
    int vertex_dimension(int) const;
    int weighted_entourage(int,int) const;
    int cyclicity() const;
    double dimensional_frontier(int) const;
    double dimensional_uniformity() const;
    bool active_simplex(int,int) const;
    int circuit_rank() const;
    int euler_characteristic() const;
    int component_analysis(std::vector<int>&) const;
    int entourage(int) const;
    bool connected() const;
    bool consistent() const;
    bool correctness();
    void clean() const;
    bool energy_check() const;
    double total_energy() const;
    void simplicial_implication(int) const;
    int simplex_embedding(int,int) const;
    double dimensional_stress(int) const;
    void write_graph(const std::string&) const;

    void compute_geometric_gradient(std::vector<double>&,bool);
    void compute_geometric_dependency(const std::set<int>&);
    void compute_topological_dependency(const std::set<int>&);
    void simplicial_implication();
    void reciprocate();
    void compute_simplex_energy(int,int);
    void compute_simplex_orientation(int,int);
    // The various methods needed for the hyphantic operators
    void hyphansis();
    void dynamic_hyphansis(const std::vector<std::pair<int,double> >&);
    void musical_hyphansis(const std::vector<std::pair<int,double> >&);
    std::string implicative_scale(int,std::vector<double>&) const;
    std::string explicative_scale(int,std::vector<double>&) const;
    int select_vertex(const std::vector<int>&,double) const;
    int vertex_addition(const std::vector<double>&);
    int vertex_addition(const std::set<int>&);
    int vertex_addition(int);
    bool vertex_deletion(int);
    bool simplex_addition(const std::set<int>&);
    void simplex_deletion(int,int);
    bool interplication(int,double,int);
    bool circumvolution();
    bool circumvolution(int);
    int compression(double,std::set<int>&);
    bool unravel(int);
    bool reduction(int);
    bool contraction(int,double);
    bool compensation_m(int);
    bool compensation_g(int);
    bool expansion(int);
    bool expansion(int,double);
    bool foliation_m(int);
    bool foliation_x(int);
    bool amputation(int,double);
    bool fusion_x(int,double);
    bool fusion_m(int);
    bool fission(int,double);
    bool inflation(int,double);
    bool deflation(int);
    bool perforation(int,int);
    bool correction(int);
    bool germination(int);
    bool vertex_twist();
    bool stellar_addition(int);
    bool stellar_deletion(int);
    bool edge_reorientation(int);
    bool edge_reorientation(int,int);
    void vertex_fusion(int,int);
    void recompute_orientation(int);
    void recompute_orientation(const std::set<int>&);
    void compute_orientation();
    void compute_neighbours();
    void compute_entourages();
    void superposition_fusion(std::set<int>&);
    void superposition_fission(std::set<int>&);
    void assemble_indices();
    void regularization(bool);
    void entourage(std::vector<int>&) const;

    void inversion();
    bool realizable(int,int) const;
    void write_topology() const;
    void write_incastrature(const std::string&) const;
    void determine_flexible_edges();
    void compute_simplicial_dimension();
    void compute_volume();
    void compute_lengths();
    void compute_curvature();
    void compute_obliquity();
    double compute_abnormality() const;
    double compute_abnormality(const std::vector<double>&) const;
    void mechanical_force(const std::vector<int>&,const std::vector<double>&,double*) const;
    double minimize_lengths(const std::vector<int>&,const std::vector<int>&,int*) const;
    void structural_deficiency();
    void compute_total_lightcone(int,std::set<int>&,std::set<int>&) const;
    void compute_lightcones();
    double compute_temporal_vorticity(int) const;
    double compute_temporal_nonlinearity() const;
    void set_logical_atoms(int);
    double logical_energy(int) const;
    bool logical_conformity(int) const;
    double representational_energy(bool) const;

    void implication(std::string&) const;
    void explication(std::string&) const;
    void compute_delta();
    bool global_operations();
    void write_distribution(const std::vector<int>&) const;
    void compute_colours(std::vector<unsigned char>&,bool) const;
    void compute_global_topology();
    void build_initial_state();
    void write_log() const;
    void read_complex(std::ifstream&);
    void read_state();
    void write_complex(std::ofstream&) const;
    void write_state() const;
    void read_parameters(const char*);
    void set_default_values();
    void update_viewer();
    void energy_diffusion();
    void energy_diffusion(int);
    void optimize();
    bool adjust_dimension();
    void analyze_convergence();
    void test_harness(int,int);
    void condense();
    void initialize();
    void clear();
    void write(Spacetime&) const;
    void read(const Spacetime&);
    inline double distribution_fitness(int*,const std::vector<int>&,int) const;

   public:
    Spacetime();
    Spacetime(bool);
    Spacetime(const char*);
    Spacetime(const char*,bool);
    ~Spacetime();
    void read_state(const std::string&);
    bool step_forwards();
    void evolve();
    void chorogenesis(int);
    void distribute(int) const; 
    void restart(const char*,bool);
    void export_visual_data(std::vector<float>&,std::vector<float>&,std::vector<int>&,int*,bool) const;
    void export_visual_data(std::vector<float>&,std::vector<int>&,int*) const;
    void set_checkpoint_frequency(int);
    bool active_element(int,int) const;
    void get_edge_topology(std::vector<std::set<int> >&) const;
    void get_energy_values(std::vector<double>&) const;
    void get_deficiency_values(std::vector<double>&) const;
    void get_coordinates(int,std::vector<double>&) const;
    void get_coordinates(std::vector<double>&) const;
    double get_geometric_distance(int,int,bool) const;
    bool is_pseudomanifold(bool*) const;
    void get_energy_extrema(double*) const;
    void get_deficiency_extrema(double*) const;
    void write_vertex_data(int) const;
    inline int get_background_dimension() const {return geometry->dimension();};
    inline int get_cardinality(int d) const {return cardinality_safe(d);};
    inline int get_event_dimension(int n) const {return vertex_dimension(n);};
    inline double get_event_obliquity(int n) const {return events[n].obliquity;};
    inline double get_event_energy(int n) const {return events[n].get_energy();};
    inline double get_event_curvature(int n) const {return events[n].curvature;};
    inline double get_event_deficiency(int n) const {return events[n].deficiency;};
    inline double get_event_entwinement(int n) const {return events[n].entwinement;};
    inline std::string get_state_file() const {return state_file;};
    inline std::string get_simplex_key(int d,int n) const {return SYNARMOSMA::make_key(simplices[d][n].vertices);};
    inline std::string get_ops() const {return ops;};
    inline void get_simplex_vertices(int d,int n,int* vx) const {simplices[d][n].get_vertices(vx);};
    inline void get_arclength_statistics(double* output) const {arclength_statistics(output);};
    inline void get_vertex_degree_statistics(double* output) const {vertex_degree_statistics(output);};
    inline void get_fvector(std::vector<int>& f,std::vector<int>& fstar) const {compute_fvector(f,fstar);};
    inline void get_delta(double* output) const {output[0] = topology_delta; output[1] = geometry_delta; output[2] = energy_delta;};
    inline int get_edge_index(const std::set<int>& vx) const {SYNARMOSMA::hash_map::const_iterator qt = index_table[1].find(vx); return qt->second;};
    inline int get_iterations() const {return iterations;};
    inline int get_dimension() const {return dimension();};
    inline int get_cyclicity() const {return cyclicity();};
    inline int get_circuit_rank() const {return circuit_rank();};
    inline double get_error() const {return error;};
    inline double get_total_energy() const {return total_energy();}; 
    inline int get_maximum_iterations() const {return max_iter;};
    inline bool is_converged() const {return converged;};
    inline bool is_orientable() const {return orientable;};
    inline int get_euler_characteristic() const {return euler_characteristic();};
    inline int get_events() const {return (signed) events.size();};
    inline int get_entourage_cardinality(int d,int n) const {int output = (d == 0) ? (signed) events[n].neighbours.size() : (signed) simplices[d][n].entourage.size(); return output;}; 
    inline int get_incept(int d,int n) const {int output = (d == 0) ? events[n].incept : simplices[d][n].incept; return output;};
  };

  inline bool Spacetime::edge_exists(int u,int v) const
  {
    std::set<int> vx;
    vx.insert(u); vx.insert(v);
    SYNARMOSMA::hash_map::const_iterator qt = index_table[1].find(vx);
    if (qt == index_table[1].end()) return false;
    if (!simplices[1][qt->second].active) return false;
    return true;
  }

  inline bool Spacetime::is_pseudomanifold(bool* bdry) const 
  {
    *bdry = boundary;
    return pseudomanifold;
  }

  inline double Spacetime::distribution_fitness(int* volume,const std::vector<int>& affinity,int nprocs) const
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

  inline void Spacetime::compute_graph(SYNARMOSMA::Graph* G,int base) const
  {
    compute_graph(G,base,Spacetime::topological_radius);
  }

  inline int Spacetime::cardinality(int d) const
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
}
#endif

