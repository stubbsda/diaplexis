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
  /// A class representing the spacetime, a combination of dynamic geometric and topological properties. 
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
    /// This property represents the number of relaxation steps 
    /// that have been taken thus far in the simulation.
    int iterations = 0;
    /// This property is the size of the spacetime complex, adjusted 
    /// for dimensionality - the sum of the active events, each event 
    /// multiplied by its geometric dimensionality. 
    int system_size = 0;
    /// This property stores the number of active events in the 
    /// spacetime complex.
    int nactive = 1;
    /// This property is true if the current relaxation step can be 
    /// reversed and false otherwise. The simulation always stores a 
    /// snapshot of the preceding relaxation step and can thus reverse 
    /// itself by a single step (but only a single step). 
    bool reversible = false;
    /// This property stores the total error of the spacetime: the sum 
    /// of the square of the Event::deficiency property for all active 
    /// events and the Spacetime::global_deficiency property, all divided 
    /// by the number of active events.
    double error = 0.0;
    /// This property stores the global component of the spacetime error: 
    /// the sum of the output of representational_energy() and \f$2\pi\f$
    /// times the output of Complex::euler_characteristic, less the sum 
    /// of the energy property of all active events. 
    double global_deficiency = 0.0;
    /// This property contains a list of all of the hyphantic operations 
    /// successfully carried out on the spacetime complex since the start 
    /// of the simulation. 
    std::string hyphantic_ops = "";
    /// This property contains all of the spacetime's topology, stored in a 
    /// pointer to an instance of the Complex class.
    Complex* skeleton;
    /// This property contains the spacetime's geometry, stored in a pointer to 
    /// an instance of the Synarmosma library's Geometry class.
    SYNARMOSMA::Geometry* geometry;

    // The global parameters
    /// This property controls the sort of initial state from which the simulation 
    /// begins, including eventually from a checkpoint file. 
    Initial_Topology initial_state = Initial_Topology::random;
    /// If the initial state of the simulation is Initial_Topology::diskfile, then 
    /// this property records that the original state was indeed "diskfile" while the 
    /// checkpoint file alters the value of Spacetime::initial_state.
    Initial_Topology original_state = Initial_Topology::random;
    /// This property determines the way in which the hyphantic operators which 
    /// are applied at each relaxation step are chosen. If this property is set to 
    /// Hyphansis::musical then a musical score is used to control the operators, 
    /// otherwise a dynamic algorithm is used based on Event::deficiency value 
    /// of the spacetime's active events.  
    Hyphansis weaving = Hyphansis::dynamic;
    /// This property stores the initial number of events in the 
    /// spacetime complex.
    int initial_size = 10;
    /// This property stores the maximum number of relaxation steps 
    /// that will be carried out in the simulation.
    int max_iter = 50;
    /// This property is the initial topological dimension of the 
    /// spacetime complex, i.e. the output of the Complex::dimension() 
    /// method.
    int initial_dim = 4;
    /// This property stores the time (as measured in seconds since 
    /// the beginning of the Unix era on January 1, 1970) at which the 
    /// simulation started. 
    boost::posix_time::ptime start_time;
    /// This property stores the date (in the format YYYYMMDD) of the 
    /// simulation's beginning, which is used in the name of the output 
    /// files for this particular simulation. 
    std::string date_string = "";
    /// This property stores the Unix process ID of the simulation, 
    /// which is used in the name of the output files for this particular 
    /// simulation. 
    std::string pid_string = "";
    /// The initial part of the name of the checkpoint file(s) created 
    /// during the simulation - these are placed in the "data" directory 
    /// with the name spacetime_date_pid_step.dat
    std::string state_file = "data/spacetime";
    /// The initial part of the name of the log file created during the 
    /// simulation - this file is placed in the "data" directory with the 
    /// name spacetime_date_pid.log
    std::string log_file = "data/spacetime";
    /// The name of the input file from which the simulation will read its 
    /// state assuming the value of Spacetime::initial_state is set to 
    /// Initial_Topology::diskfile. 
    std::string input_file = "";
    /// If the spacetime's initial topology is random, then this property 
    /// controls the probability that there is an edge connecting any pair 
    /// of distinct events. 
    double edge_probability = std::log(500.0)/250.0;
    /// This property, a percentage between zero and one, controls the frequency 
    /// with which the spacetime complex's d-simplices (d > 0) have their Simplex::parity
    /// property mutated at each relaxation step. 
    double parity_mutation = 0.01;
    /// If this property is true the simulation does not carry out any disk 
    /// I/O operations - there is no log file, no checkpoint files are created 
    /// and so forth. 
    bool diskless = false;
    /// This property determines whether or not the spacetime complex's 
    /// initial topology is perturbed at the start of the simulation. 
    bool perturb_topology = false;
    /// This property determines whether or not the spacetime complex's 
    /// initial geometry is perturbed at the start of the simulation. 
    bool perturb_geometry = false;
    /// This property determines whether or not the spacetime complex's 
    /// initial energy distribution is perturbed at the start of the 
    /// simulation. 
    bool perturb_energy = true;
    /// If this property is set to true then the superposition_fusion() 
    /// and superposition_fission() methods are called at each relaxation 
    /// step. 
    bool superposable = false;
    /// If this property is set to true then the compression() method is 
    /// called at each relaxation step.
    bool compressible = false;
    /// If the Spacetime::error property falls below the value of Spacetime::convergence_threshold,
    /// then this property is set to true.
    bool converged = false;
    /// If the amount of physical memory available on the computer is plentiful 
    /// relative to the size of the spacetime complex, then this property should 
    /// be set to true to improve the program's performance (at the expense of 
    /// significantly increasing its memory footprint).
    bool high_memory = true;
    /// This property controls how frequently the the method write_state() 
    /// is invoked during the simulation, in terms of how many relaxation 
    /// steps elapse before a new checkpoint is made. The value must be greater 
    /// than zero and if it is greater than Spacetime::max_iter then no 
    /// checkpoint files will be created.
    int checkpoint_frequency = 50;
    /// The initial part of the name of the hyphansis log file created during 
    /// the simulation - this file is placed in the "data" directory with the 
    /// name hyphansis_date_pid.log
    std::string hyphansis_file = "data/hyphansis";
    /// If the Spacetime::weaving property is set to Hyphansis::musical, then 
    /// this property stores the filename of the score file that is used to 
    /// control the selection of hyphantic operators that will be applied at 
    /// each relaxation step. 
    std::string hyphansis_score = "";

    // Now the parameters associated with the
    // geometry solver
    /// This property controls which optimization algorithm is used to try and 
    /// ensure the geometry meshes well with the topology and energy distribution 
    /// at a given relaxation step. 
    Geometry_Solver solver = Geometry_Solver::minimal;
    /// This is the convergence threshold used for determining if a geometry solver 
    /// has succeeded and the algorithm can exit. 
    double geometry_tolerance = 0.0001;
    /// This property is only meaningful for the Geometry_Solver::minimal case and 
    /// controls the number of random event coordinate changes attempted before exiting 
    /// the geometry solver. This property should be positive. 
    int solver_its = 50;
    /// This property is only meaningful for the Geometry_Solver::evolutionary case and 
    /// controls the number of generations of artificial evolution that will be attempted 
    /// before exiting the geometry solver. This property should be positive.
    int ngenerations = 1000;
    /// This property is only meaningful for the Geometry_Solver::evolutionary case and 
    /// controls the number of one-on-one contests carried out at each generation to determine 
    /// which individuals will survive to the next generation. This property should be 
    /// positive and no greater than Spacetime::pool_size.  
    int njousts = 20;
    /// This property is only meaningful for the Geometry_Solver::evolutionary case and 
    /// controls the number of individuals in a generation; a higher number means more 
    /// genetic diversity but at the price of consuming more memory. This property should 
    /// be positive.
    int pool_size = 100;
    // Annealing
    /// This positive property controls the variance used in mutating the geometry of an 
    /// individual event in the simulated annealing geometry solver and so is only relevant 
    /// when Spacetime::solver is set to Geometry_Solver::annealing.
    double thermal_variance = 0.5;
    /// This positive property determines the number of iterations of the thermal sweep loop, 
    /// used in testing whether the spacetime geometry has reached thermal equilibrium at a 
    /// given temperature. It is thus only relevant when Spacetime::solver is set to 
    /// Geometry_Solver::annealing.
    int thermal_sweep = 1000;
    /// This positive property is the number of iterations performed by the main annealing loop, 
    /// which gradually lowers the temperature, ensuring the spacetime geometry reaches thermal 
    /// equilibrium at each temperature. It is thus only relevant when Spacetime::solver is set 
    /// to Geometry_Solver::annealing.
    int annealing_steps = 500;
    /// This positive property that should normally be less than unity is what determines whether 
    /// or not the spacetime geometry has reached thermal equilibrium, by calculating the variance 
    /// of energies over the ensemble of thermal sweeps. If the variance is less than this threshold 
    /// then the spacetime geometry has reached thermal equilibrium. This property is naturally 
    /// meaningful only when Spacetime::solver is set to Geometry_Solver::annealing.
    double thermalization = 0.001;
    // Mechanical
    /// This property controls which ODE integration algorithm is used (forward Euler 
    /// or fourth order Runge-Kutta) if the Spacetime::solver property is set to 
    /// Geometry_Solver::mechanical. 
    Integrator engine = Integrator::rk4;
    /// This property sets the maximum number of integration steps that will be carried 
    /// out for the Geometry_Solver::mechanical solver model, which uses a velocity norm 
    /// convergence criterion. This property should be positive.
    int max_int_steps = 10000;
    /// This property is the step size for the integration algorithm used to numerically 
    /// solve the ODE system in the Geometry_Solver::mechanical solver and it should normally 
    /// be greater than zero and (significantly) less than one. 
    double step_size = 0.05;
    /// This property determines in part the damped spring model which is used to compute 
    /// event geometry in the Geometry_Solver::mechanical engine. This property should be negative 
    /// and expresses the spring's stiffness. 
    double spring_constant = -1.5;
    /// This property determines in part the damped spring model which is used to compute 
    /// event geometry in the Geometry_Solver::mechanical engine. This property should be positive 
    /// and expresses the strength of the repulsive force between events. 
    double repulsion_constant = 1.0;
    /// This property determines in part the damped spring model which is used to compute 
    /// event geometry in the Geometry_Solver::mechanical engine. This property should be positive 
    /// and expresses the strength of the friction in the system, causing the model to gradually 
    /// lose its kinetic energy.  
    double damping_constant = 0.85;
    /// This property is used to compute the degree of "abnormality" in the spacetime geometry; 
    /// if the distance between two connected events is less than this value, then the distance is 
    /// ignored in the abnormality calculation. This property should be positive and is only 
    /// meaningful if Spacetime::solver is Geometry_Solver::mechanical and Spacetime::cgradient_refinement 
    /// is true.
    double edge_flexibility_threshold = 2.0;
    /// This property controls whether or not a the Geometry_Solver::mechanical engine will attempt 
    /// to perform a conjugate gradient optimization designed to normalize as much as possible the 
    /// inter-event distance so that it is close to unity. 
    bool cgradient_refinement = true;
    /// This property is only meaningful if Spacetime::solver is set to Geometry_Solver::mechanical 
    /// and Spacetime::cgradient_refinement is true, in which case it controls the maximum number of 
    /// conjugate gradient steps. 
    int max_CG_steps = 10;
    /// This property is only meaningful if Spacetime::solver is set to Geometry_Solver::mechanical 
    /// and Spacetime::cgradient_refinement is true, in which case it controls the maximum number of 
    /// line search steps. 
    int max_LS_steps = 20;
    // Simplex
    /// This positive property controls the reflection transformation used when the Spacetime::solver 
    /// property is set to Geometry_Solver::simplex.
    double simplex_alpha = 1.0;
    /// This property controls the expansion transformation used when the Spacetime::solver property
    /// is set to Geometry_Solver::simplex; it must be greater than the maximum of 1.0 and Spacetime::simplex_alpha. 
    double simplex_gamma = 2.0;
    /// This positive property controls the external contraction transformation used when the Spacetime::solver 
    /// property is set to Geometry_Solver::simplex. This property must be less than unity.
    double simplex_rho = 0.5;
    /// This positive property controls the internal contraction transformation used when the Spacetime::solver 
    /// property is set to Geometry_Solver::simplex. This property must be less than unity.
    double simplex_sigma = 0.5;

    // Stuff for the implicative/explicative operators:
    /// The total number of available explicative operators 
    /// for the hyphansis phase of each relaxation step.
    static const int N_EXP = 11;
    /// The total number of available implicative operators 
    /// for the hyphansis phase of each relaxation step.
    static const int N_IMP = 8;
    /// This array of strings stores the abbreviations that are 
    /// used to represent the explicative operators in the 
    /// log file for hyphansis as well as the Spacetime::hyphantic_ops 
    /// property.
    static const std::string EXP_OP[N_EXP];
    /// This array of strings stores the abbreviations that are 
    /// used to represent the implicative operators in the 
    /// log file for hyphansis as well as the Spacetime::hyphantic_ops 
    /// property.
    static const std::string IMP_OP[N_IMP];

    /// This property stores the error tolerance for the 
    /// convergence of the relaxation process.
    static const double convergence_threshold;
    /// The coupling constant between the topological-geometric torsion 
    /// and the energy in the structure equation.
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
    void vertex_fusion(int,int);
    bool vertex_twist();

    bool circumvolution();
    bool circumvolution(int);
    bool contraction(int,double);
    bool compensation_m(int);
    bool compensation_g(int);
    bool expansion(int);
    bool expansion(int,double);
    bool foliation_m(int);
    bool foliation_x(int);
    bool reduction(int);
    bool amputation(int,double);
    bool fusion_x(int,double);
    bool fusion_m(int);
    bool fission(int,double);
    bool inflation(int,double);
    bool deflation(int);
    bool perforation(int,int);
    bool correction(int);
    bool germination(int);
    bool stellar_addition(int);
    bool stellar_deletion(int);

    bool interplication(int,double,int);
    int compression(double);
    void superposition_fusion(double);
    void superposition_fission(int);
    void regularization(bool);
    void entourage(std::vector<int>&) const;

    bool realizable(int,int) const;
    void compute_volume();
    void compute_lengths();
    void compute_obliquity();
    double compute_abnormality(const std::vector<int>&) const;
    double compute_abnormality(const std::vector<double>&,const std::vector<int>&) const;
    void compute_geometric_gradient(std::vector<double>&,bool,const std::vector<int>&);
    void mechanical_force(const std::vector<int>&,const std::vector<int>&,const std::vector<double>&,const std::vector<double>&,double*) const;
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
    void compute_colours(std::vector<unsigned char>&,bool) const;
    void build_initial_state();
    void write_log() const;
    void write(Spacetime&) const;
    void read(const Spacetime&);
    void read_state(const std::string&);
    void write_state(const std::string& = "") const;
    void read_parameters(const std::string&);
    bool adjust_dimension();
    void arclength_statistics(double*) const;
    void condense();
    void initialize();
    void allocate();
    void clear();
    void clean() const;
    bool correctness();

   public:
    Spacetime();
    Spacetime(bool);
    Spacetime(const std::string&);
    Spacetime(const std::string&,bool);
    ~Spacetime();
    bool advance();
    void fallback();
    void evolve();
    void chorogenesis(int);
    void restart(const std::string&,bool);
    void export_visual_data(std::vector<float>&,std::vector<float>&,std::vector<int>&,int*,bool) const;
    void export_visual_data(std::vector<float>&,std::vector<int>&,int*) const;
    inline void set_checkpoint_frequency(int n) {checkpoint_frequency = n;};
    void get_coordinates(int,std::vector<double>&) const;
    void get_coordinates(std::vector<double>&) const;
    void get_energy_extrema(double*) const;
    void get_deficiency_extrema(double*) const;
    inline double get_geometric_distance(int n,int m) const {return geometry->get_squared_distance(n,m,false);};
    inline int get_background_dimension() const {return geometry->dimension();};
    inline std::string get_state_file() const {return state_file;};
    inline std::string get_hyphantic_operations() const {return hyphantic_ops;};
    inline void get_arclength_statistics(double* output) const {arclength_statistics(output);};
    inline int get_iterations() const {return iterations;};
    inline double get_error() const {return error;};
    inline int get_maximum_iterations() const {return max_iter;};
    inline bool is_converged() const {return converged;};
  };
}
#endif

