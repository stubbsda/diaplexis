// For XML processing...
#include <pugixml.hpp>

#include <boost/timer/timer.hpp>
#include <boost/filesystem.hpp>

#include "synarmosma/geometry.h"
#include "synarmosma/directed_graph.h"

#include "sheet.h"
#include "complex.h"

#ifndef _spacetimeh
#define _spacetimeh

namespace DIAPLEXIS {
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
    int iterations = 0;
    int system_size = 0;
    int nactive = 1;
    bool reversible = false;
    double error = 0.0;
    double global_deficiency = 0.0;
    std::set<int> modified_vertices;
    Complex* skeleton;
    SYNARMOSMA::Geometry* geometry;
    std::vector<Sheet> codex;

    // The global parameters
    Initial_Topology initial_state = Initial_Topology::random;
    Initial_Topology original_state = Initial_Topology::random;
    Hyphansis weaving = Hyphansis::dynamic;
    int initial_size = 10;
    int max_iter = 50;
    int nt_initial = 1;
    unsigned int initial_dim = 4;
    boost::posix_time::ptime start_time;
    std::string date_string = "";
    std::string pid_string = "";
    std::string state_file = "data/spacetime";
    std::string log_file = "data/spacetime";
    std::string input_file = "";
    double edge_probability = std::log(500.0)/250.0;
    double parity_mutation = 0.01;
    bool diskless = false;
    bool perturb_topology = false;
    bool perturb_geometry = false;
    bool perturb_energy = true;
    bool superposable = false;
    bool compressible = false;
    bool permutable = false;
    bool converged = false;
    bool foliodynamics = false;
    bool high_memory = true;
    int checkpoint_frequency = 50;
    std::string hyphansis_file = "data/hyphansis";
    std::string hyphansis_score = "";

    // Now the parameters associated with the
    // geometry solver
    Geometry_Solver solver = Geometry_Solver::minimal;
    Integrator engine = Integrator::rk4;
    double geometry_tolerance = 0.0001;
    // Minimal
    int solver_its = 50;
    // Evolutionary
    int ngenerations = 1000;
    int njousts = 20;
    int pool_size = 100;
    // Annealing
    double thermal_variance = 0.5;
    int thermal_sweep = 1000;
    int annealing_steps = 500;
    double thermalization = 0.001;
    // Mechanical
    int max_int_steps = 10000;
    double step_size = 0.05;
    double spring_constant = -1.5;
    double repulsion_constant = 1.0;
    double damping_constant = 0.85;
    double edge_flexibility_threshold = 2.0;
    bool cgradient_refinement = true;
    // Conjugate gradient
    int max_CG_steps = 10;
    int max_LS_steps = 20;
    // Simplex
    double simplex_alpha = 1.0;
    double simplex_gamma = 2.0;
    double simplex_rho = 0.5;
    double simplex_sigma = 0.5;

    // Stuff for the implicative/explicative operators:
    static const int N_EXP = 11;
    static const int N_IMP = 8;
    static const std::string EXP_OP[N_EXP];
    static const std::string IMP_OP[N_IMP];

    // The intensity of branching in the polycosmos, 0 < x < 1
    static const double ramosity;
    // Error tolerance for convergence
    static const double convergence_threshold;
    // The initial temperature for cosmic annealing
    static const double T_zero;
    // The rate of thermal decay in the annealing schedule
    static const double kappa;
    // The coupling constant between the topological-geometric torsion 
    // and the energy
    static const double Lambda;

    std::string sheet_activity() const;
    void compute_geometric_gradient(std::vector<double>&,bool);
    // The various methods needed for the hyphantic operators
    void hyphansis(int);
    void dynamic_hyphansis(const std::vector<std::pair<int,double> >&,int);
    void musical_hyphansis(const std::vector<std::pair<int,double> >&,int);
    std::string implicative_scale(int,std::vector<double>&) const;
    std::string explicative_scale(int,std::vector<double>&) const;
    int select_vertex(const std::vector<int>&,double,int) const;
    int vertex_addition(const std::vector<double>&,int);
    int vertex_addition(const std::set<int>&,int);
    int vertex_addition(int,int);
    bool vertex_deletion(int,int);
    bool simplex_addition(const std::set<int>&,int);
    void simplex_deletion(int,int,int);
    bool interplication(int,double,int,int);
    bool circumvolution(int);
    bool circumvolution(int,int);
    int compression(double,std::set<int>&);
    bool unravel(int,int);
    bool reduction(int,int);
    bool contraction(int,double,int);
    bool compensation_m(int,int);
    bool compensation_g(int,int);
    bool expansion(int,int);
    bool expansion(int,double,int);
    bool foliation_m(int,int);
    bool foliation_x(int,int);
    bool amputation(int,double,int);
    bool fusion_x(int,double,int);
    bool fusion_m(int,int);
    bool fission(int,double,int);
    bool inflation(int,double,int);
    bool deflation(int,int);
    bool perforation(int,int,int);
    bool correction(int,int);
    bool germination(int,int);
    bool vertex_twist(int);
    bool stellar_addition(int,int);
    bool stellar_deletion(int,int);
    bool edge_parity_mutation(int,int);
    bool edge_parity_mutation(int,int,int);
    void vertex_fusion(int,int,int);
    void recompute_parity(int);
    void recompute_parity(const std::set<int>&);
    void compute_parity();
    void compute_neighbours();
    void compute_entourages(int);
    void superposition_fusion(std::set<int>&);
    void superposition_fission(std::set<int>&);
    void assemble_indices();
    void regularization(bool,int);
    void entourage(std::vector<int>&) const;

    void inversion();
    bool realizable(int,int) const;
    void compute_simplicial_dimension();
    void compute_volume();
    void compute_lengths();
    void compute_obliquity();
    double compute_abnormality() const;
    double compute_abnormality(const std::vector<double>&) const;
    void mechanical_force(const std::vector<int>&,const std::vector<double>&,double*) const;
    double minimize_lengths(const std::vector<int>&,const std::vector<int>&,int*) const;
    void structural_deficiency();
    void compute_lightcones();
    void compute_total_lightcone(int,std::set<int>&,std::set<int>&) const;
    double compute_temporal_vorticity(int,int) const;
    double compute_temporal_nonlinearity() const;
    double representational_energy(bool) const;

    void implication(std::string&) const;
    void explication(std::string&) const;
    void compute_delta();
    int sheet_fission(int);
    int sheet_dynamics();
    bool global_operations();
    void compute_colours(std::vector<unsigned char>&,bool,bool) const;
    void compute_global_topology(int);
    void build_initial_state(const std::set<int>&);
    void write_log() const;
    void read_complex(std::ifstream&);
    void read_state();
    void write_complex(std::ofstream&) const;
    void write_state() const;
    void read_parameters(const std::string&);
    int ubiquity_permutation(double,std::set<int>&);
    void update_viewer();
    void mechanical_solver();
    void evolutionary_solver();
    void annealing_solver();
    void simplex_solver();
    void optimize();
    bool adjust_dimension();
    void condense();
    void initialize();
    void allocate();
    void clear();
    void write(Spacetime&) const;
    void read(const Spacetime&);

   public:
    Spacetime();
    Spacetime(bool);
    Spacetime(const std::string&);
    Spacetime(const std::string&,bool);
    ~Spacetime();
    void read_state(const std::string&);
    void fallback();
    bool advance();
    void evolve();
    void chorogenesis(int);
    void restart(const std::string&,bool);
    void export_visual_data(std::vector<float>&,std::vector<float>&,std::vector<int>&,int*,bool) const;
    void export_visual_data(std::vector<float>&,std::vector<int>&,int*,int) const;
    bool active_element(int,int) const;
    void get_ubiquity(int,int,std::string&) const;
    void get_edge_topology(std::vector<std::set<int> >&) const;
    void get_ubiquity_vector(std::vector<int>&) const;
    void get_energy_values(std::vector<double>&,int) const;
    void get_deficiency_values(std::vector<double>&,int) const;
    void get_coordinates(int,std::vector<double>&) const;
    void get_coordinates(std::vector<double>&) const;
    bool is_pseudomanifold(bool*,int) const;
    void get_energy_extrema(double*) const;
    void get_deficiency_extrema(double*) const;
    void write_vertex_data(int) const;

    inline void set_checkpoint_frequency(int n) {checkpoint_frequency = n;};
    inline double get_geometric_distance(int,int) const {return geometry->get_squared_distance(n,m,false);};
    inline int get_background_dimension() const {return geometry->dimension();};
    inline int get_cardinality(int d,int sheet) const {return cardinality_safe(d,sheet);};
    inline int get_event_dimension(int n,int sheet) const {return vertex_dimension(n,sheet);};
    inline double get_event_obliquity(int n) const {return events[n].obliquity;};
    inline double get_event_energy(int n) const {return events[n].get_energy();};
    inline double get_event_deficiency(int n) const {return events[n].deficiency;};
    inline double get_event_entwinement(int n,int sheet) const {return events[n].entwinement[sheet];};
    inline std::string get_state_file() const {return state_file;};
    inline std::string get_simplex_key(int d,int n) const {return SYNARMOSMA::make_key(simplices[d][n].vertices);};
    inline std::string get_sheet_ops(int n) const {return codex[n].ops;};
    inline std::string get_sheet_activity() const {return sheet_activity();};
    inline void get_simplex_vertices(int d,int n,int* vx) const {simplices[d][n].get_vertices(vx);};
    inline void get_arclength_statistics(double* output,int sheet) const {arclength_statistics(output,sheet);};
    inline void get_vertex_degree_statistics(double* output,int sheet) const {vertex_degree_statistics(output,sheet);};
    inline void get_fvector(std::vector<int>& f,std::vector<int>& fstar,int sheet) const {compute_fvector(f,fstar,sheet);};
    inline void get_delta(double* output) const {output[0] = topology_delta; output[1] = geometry_delta; output[2] = energy_delta;};
    inline int get_codex_size() const {return (signed) codex.size();};
    inline int get_edge_index(const std::set<int>& vx) const {SYNARMOSMA::hash_map::const_iterator qt = index_table[1].find(vx); return qt->second;};
    inline int get_iterations() const {return iterations;};
    inline int get_dimension(int sheet) const {return dimension(sheet);};
    inline int get_cyclicity(int sheet) const {return cyclicity(sheet);};
    inline int get_circuit_rank(int sheet) const {return circuit_rank(sheet);};
    inline double get_error() const {return error;};
    inline double get_total_energy(int sheet) const {return total_energy(sheet);};
    inline int get_maximum_iterations() const {return max_iter;};
    inline bool is_converged() const {return converged;};
    inline bool is_orientable(int sheet) const {bool output = (sheet == -1) ? orientable : codex[sheet].orientable; return output;};
    inline int get_euler_characteristic(int sheet) const {return euler_characteristic(sheet);};
    inline int get_events() const {return (signed) events.size();};
    inline int get_entourage_cardinality(int d,int n) const {int output = (d == 0) ? (signed) events[n].neighbours.size() : (signed) simplices[d][n].entourage.size(); return output;}; 
    inline int get_incept(int d,int n) const {int output = (d == 0) ? events[n].incept : simplices[d][n].incept; return output;};
  };

  inline std::string Spacetime::sheet_activity() const
  {
    std::string out = "(";
    for(int i=0; i<(signed) codex.size()-1; ++i) {
      out += (boost::lexical_cast<std::string>(codex[i].active) + ",");
    }
    out += boost::lexical_cast<std::string>(codex[codex.size()-1].active);
    out += ")";
    return out;
  }

  inline bool Spacetime::is_pseudomanifold(bool* bdry,int sheet) const 
  {
    if (sheet == -1) {
      return skeleton->is_pseudomanifold(bdry);
    }
    else {
      *bdry = codex[sheet].boundary;
      return codex[sheet].pseudomanifold;
    }
  }
}
#endif

