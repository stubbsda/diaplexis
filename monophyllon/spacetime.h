// For XML processing...
#include <pugixml.hpp>

#include <boost/timer/timer.hpp>
#include <boost/filesystem.hpp>

#include "synarmosma/geometry.h"
#include "synarmosma/directed_graph.h"

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
    double topology_delta = 0.0;
    double geometry_delta = 0.0;
    double energy_delta = 0.0;
    std::string hyphantic_ops = "";
    std::set<int> vx_delta;
    std::vector<int> flexible_edge;
    Complex* skeleton;
    SYNARMOSMA::Geometry* geometry;

    // The global parameters
    Initial_Topology initial_state = Initial_Topology::random;
    Initial_Topology original_state = Initial_Topology::random;
    Hyphansis weaving = Hyphansis::dynamic;
    int initial_size = 10;
    int max_iter = 50;
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
    bool converged = false;
    bool high_memory = true;
    bool instrument_convergence = true;
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

    // Error tolerance for convergence
    static const double convergence_threshold;
    // The initial temperature for cosmic annealing
    static const double T_zero;
    // The rate of thermal decay in the annealing schedule
    static const double kappa;
    // The coupling constant between the topological-geometric torsion 
    // and the energy
    static const double Lambda;

    // The various methods needed for the hyphantic operators
    void hyphansis();
    void dynamic_hyphansis(const std::vector<std::pair<int,double> >&);
    void musical_hyphansis(const std::vector<std::pair<int,double> >&);
    std::string implicative_scale(int,std::vector<double>&) const;
    std::string explicative_scale(int,std::vector<double>&) const;
    void implication(std::string&) const;
    void explication(std::string&) const;
    int select_vertex(const std::vector<int>&,double) const;
    int vertex_addition(const std::vector<double>&);
    int vertex_addition(const std::set<int>&);
    int vertex_addition(int);
    bool vertex_deletion(int);
    bool interplication(int,double,int);
    bool circumvolution();
    bool circumvolution(int);
    int compression(double,std::set<int>&);
    bool unravel(int);
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
    void vertex_fusion(int,int);
    void superposition_fusion(std::set<int>&);
    void superposition_fission(std::set<int>&);
    void assemble_indices();
    void regularization(bool);
    void entourage(std::vector<int>&) const;

    bool realizable(int,int) const;
    void compute_volume();
    void compute_lengths();
    void compute_obliquity();
    double compute_abnormality() const;
    double compute_abnormality(const std::vector<double>&) const;
    void compute_geometric_gradient(std::vector<double>&,bool);
    void mechanical_force(const std::vector<int>&,const std::vector<double>&,double*) const;
    void mechanical_solver();
    void annealing_solver();
    void simplex_solver();
    void evolutionary_solver();
    void optimize();
    double minimize_lengths(const std::vector<int>&,const std::vector<int>&,int*) const;
    void structural_deficiency();
    void compute_total_lightcone(int,std::set<int>&,std::set<int>&) const;
    void compute_lightcones();
    void compute_causal_graph(SYNARMOSMA::Directed_Graph*,int) const;
    double compute_temporal_vorticity(int) const;
    double compute_temporal_nonlinearity() const;
    double representational_energy(bool) const;

    bool global_operations();
    void write_distribution(const std::vector<int>&) const;
    void compute_colours(std::vector<unsigned char>&,bool) const;
    void build_initial_state();
    void write_log() const;
    void read_state();
    void write_state() const;
    void read_parameters(const std::string&);
    void set_default_values();
    void update_viewer();
    bool adjust_dimension();
    void analyze_convergence();
    void test_harness(int,int);
    void arclength_statistics(double*) const;
    void condense();
    void initialize();
    void allocate();
    void clear();
    void clean() const;
    bool correctness();
    void write(Spacetime&) const;
    void read(const Spacetime&);

   public:
    Spacetime();
    Spacetime(bool);
    Spacetime(const std::string&);
    Spacetime(const std::string&,bool);
    ~Spacetime();
    void read_state(const std::string&);
    bool advance();
    void fallback();
    void evolve();
    void chorogenesis(int);
    void distribute(int) const; 
    void restart(const std::string&,bool);
    void export_visual_data(std::vector<float>&,std::vector<float>&,std::vector<int>&,int*,bool) const;
    void export_visual_data(std::vector<float>&,std::vector<int>&,int*) const;
    void set_checkpoint_frequency(int);
    void get_coordinates(int,std::vector<double>&) const;
    void get_coordinates(std::vector<double>&) const;
    void get_energy_extrema(double*) const;
    void get_deficiency_extrema(double*) const;
    inline double get_geometric_distance(int n,int m) const {return geometry->get_squared_distance(n,m,false);};
    inline int get_background_dimension() const {return geometry->dimension();};
    inline std::string get_state_file() const {return state_file;};
    inline std::string get_hyphantic_operations() const {return hyphantic_ops;};
    inline void get_arclength_statistics(double* output) const {arclength_statistics(output);};
    inline void get_delta(double* output) const {output[0] = topology_delta; output[1] = geometry_delta; output[2] = energy_delta;};
    inline int get_iterations() const {return iterations;};
    inline double get_error() const {return error;};
    inline int get_maximum_iterations() const {return max_iter;};
    inline bool is_converged() const {return converged;};
  };
}
#endif

