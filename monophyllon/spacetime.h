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
    /// This enumerated class lists the five initial states for the spacetime: random (the Spacetime::edge_probability 
    /// property is used to randomly distribute edges among the events), monoplex (a d-simplex where d is the value of 
    /// the Spacetime::initial_dim property), cartesian (the initial state is an orthonormal lattice of events), singleton 
    /// (a single event and thus a sort of "Big Bang" model) and diskfile, where the initial state is read from a binary 
    /// disk file for checkpoint-restart. 
    enum Initial_Topology
    {
        random,
        monoplex,
        cartesian,
        singleton,
        diskfile
    };

    /// This enumerated class lists the five different geometry solvers available for fitting the spacetime geometry to 
    /// its current topology: the minimal solver (random event geometry changes that are accepted if the error is reduced),
    /// the mechanical solver (using a damped spring model and an ODE integrator for assigning event geometry along with 
    /// a conjugate gradient algorithm to reduce 1-simplex abnormality), the evolutionary solver (using an evolutionary 
    /// model with a population where each individual has a different geometry), the annealing solver (using a simulated 
    /// annealing algorithm to minimize the Spacetime::error value) and the simplex solver (using the Nelder-Mead downhill 
    /// simplex optimization algorithm).  
    enum Geometry_Solver
    {
        minimal,
        mechanical,
        evolutionary,
        annealing,
        simplex
    };

    /// This enumerated class lists the two different algorithms available for solving the ODE system when using the 
    /// mechanical geometry solver: either a forward Euler method (fast but rather inaccurate and unstable) or the fourth-order 
    /// Runge-Kutta method which is more costly but also more reliable. 
    enum class Integrator
    {
        euler,
        rk4
    };

    /// This enumerated class lists the two different algorithms according to which the hyphansis step in the spacetime 
    /// simulation is carried out: either a dynamic method based on the magnitude and sign of the Event::deficiency property 
    /// for each active event or else one which is based on a music score. 
    enum class Hyphansis
    {
        dynamic,
        musical
    };

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
    /// steps elapse before a new checkpoint is made. If the value is zero
    /// (or greater than Spacetime::max_iter) then no checkpoint files will 
    /// be written to disk.
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

    /// This method is called by the advance() method to carry out the hyphansis step of topological change.
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

    /// This method tests where the d-simplex specified by the method's two arguments (the dimension d and index of the simplex in Complex::simplices[d]) can be realized given the geometry of the simplex's events; the method returns true if the d-simplex has a real volume and false otherwise. 
    bool realizable(int,int) const;
    /// This method computes the volume of all of the spacetime's active d-simplices (d > 0) using the current geometry; to handle the 1-simplices, it calls compute_lengths().
    void compute_volume();
    /// This method computes the length of all of the spacetime's active 1-simplices (edges) using the current geometry.
    void compute_lengths();
    /// This method computes the Event::obliquity property for all of the spacetime's active events using the current geometry.
    void compute_obliquity();
    double compute_abnormality(const std::vector<int>&) const;
    double compute_abnormality(const std::vector<double>&,const std::vector<int>&) const;
    void compute_geometric_gradient(std::vector<double>&,bool,const std::vector<int>&);
    void mechanical_force(const std::vector<int>&,const std::vector<int>&,const std::vector<double>&,const std::vector<double>&,double*) const;
    void mechanical_solver();
    /// This method optimizes the spacetime geometry using a simulated annealing algorithm.
    void annealing_solver();
    /// This method optimizes the spacetime geometry using the Nelder-Mead downhill simplex algorithm.
    void simplex_solver();
    /// This method optimizes the spacetime geometry using an evolutionary algorithm, with a population of geometries undergoing selection based on the value of the Spacetime::error property.
    void evolutionary_solver();
    /// This method carries out the geometry optimization, calling other methods (e.g. mechanical_solver(), simplex_solver() etc.) as necessary.
    void optimize();
    double minimize_lengths(const std::vector<int>&,const std::vector<int>&,int*) const;
    /// This method computes the current value of the terms in the structure equation for the spacetime, including setting the value of the Spacetime::error property. 
    void structural_deficiency();
    void compute_total_lightcone(int,std::set<int>&,std::set<int>&) const;
    void compute_lightcones();
    void compute_causal_graph(SYNARMOSMA::Directed_Graph*,int) const;
    double compute_temporal_vorticity(int) const;
    double compute_temporal_nonlinearity() const;
    double representational_energy(bool) const;

    /// This method carries out the various global operations that must be performed at each relaxation step.
    bool global_operations();
    /// This method computes the colours to be used for the events when visualizing the spacetime; a set of three RGB float values (lying between zero and unity) is used for each event - the inactive events are made invisible - and the criterion used for the colouring is either the energy (second argument true) or deficiency property (false) of the event. 
    void compute_colours(std::vector<unsigned char>&,bool) const;
    /// This method computes the initial configuration of the spacetime based on its parameters.
    void build_initial_state();
    /// This method writes out information on the current configuration of the Spacetime instance to the log file, whose name is stored in the Spacetime::log_file property.
    void write_log() const;
    void write(Spacetime&) const;
    void read(const Spacetime&);
    /// This method calls clear() and then reads an instance of the Spacetime class from the binary file whose name is the method's unique argument.
    void read_state(const std::string&);
    /// This method writes the Spacetime instance to a binary file; the file's name is the method's argument if it has one, otherwise it uses the Spacetime::state_file property.
    void write_state(const std::string& = "") const;
    /// This method reads the run-time parameters of the Spacetime class from an XML file whose name is the method's unique argument.
    void read_parameters(const std::string&);
    /// This method adjusts the dimension of each spacetime event according to its simplicial membership (topological dimension) and updates the value of Spacetime::system_size accordingly. It returns the output from calling the Geometry::adjust_dimension method of the Synarmosma library.
    bool adjust_dimension();
    /// This method computes the maximum, minimum and arithmetic mean of the absolute value of the length of the complex's active 1-simplices and writes it to the method's argument, an array with three elements. 
    void arclength_statistics(double*) const;
    /// This method reduces memory pressure in the simulation by eliminating inactive events and 1-simplices from the spacetime, if the ratio of active to inactive is less than half for both categories. The method returns the maximum of the current ratios.  
    double condense();
    /// This method attributes a value to a variety of Spacetime properties and the calls the build_initial_state() method.  
    void initialize();
    /// This method allocates the resources for the Spacetime::geometry and Spacetime::skeleton properties.
    void allocate();
    /// This method calls the clear method on the extended properties of the class and sets Spacetime::iterations to zero and Spacetime::hyphantic_ops to the empty string. 
    void clear();
    /// This method returns true if all of the spacetime's d-simplices (d >= 0) have a modified property equal to false; it will otherwise return false. 
    bool clean() const;
    /// This method computes the value of Spacetime::error in the spacetime's current state, then sets all of the d-simplices (d >= 0) as modified and re-computes the Spacetime::error property; if there is any difference in the two values it returns false. 
    bool correctness();
    bool advance();
    void fallback();
    void restart(const std::string&,bool);

   public:
    /// This constructor begins by calling the allocate() method; if the optional argument is true, it sets Spacetime::diskless to true and Spacetime::checkpoint_frequency to zero, then calls initialize(). 
    Spacetime(bool = false);
    /// This constructor accepts a string argument which contains the XML parameter file which is parsed using the read_parameters() method, after calling the allocate() method. If the optional second argument is true, it sets Spacetime::diskless to true and Spacetime::checkpoint_frequency to zero, then calls initialize(). 
    Spacetime(const std::string&,bool = false);
    /// The destructor frees the resources associated with the Spacetime::geometry and Spacetime::skeleton properties.
    ~Spacetime();
    /// This is the main public method of this class, which evolves the spacetime over a set of relaxation steps until Spacetime::converged is true or the number of steps exceeds Spacetime::max_iter. The method returns the value of Spacetime::error at completion.
    double evolve();
    /// This method is similar to evolve(), in that it carries out a simulation of the spacetime but this time with the goal of adapting the topology to the geometric constraints (e.g. dimensionality). The spacetime is assumed to start in a state with many edges, some of which must be deleted to match the geometry and degree distribution for a regular space. The method returns the value of Spacetime::error at completion. 
    double chorogenesis();
    /// This method writes (using an offset to handle inactive events and edges) the coordinates of the spacetime's active events in the first argument, an edge table (listing the two event indices) in the second argument and the third argument is a pair of integers containing the number of active events and active edges.
    void export_visual_data(std::vector<float>&,std::vector<int>&,std::pair<int,int>&) const;
    /// This method writes the coordinates (second argument), edge table (third argument) and the total number of events and edges (the pair of integers that is the fourth argument). The first argument is a vector of colours for each event with inactive events and edges set to be invisible; the final argument determines the basis for the colouring scheme, using energy or deficiency values.
    void export_visual_data(std::vector<float>&,std::vector<float>&,std::vector<int>&,std::pair<int,int>&,bool) const;
    /// This method sets the value of Spacetime::checkpoint_frequency to the method's argument.
    inline void set_checkpoint_frequency(int n) {checkpoint_frequency = n;};
    /// This method sets the second argument to the coordinates of the event whose index is the first argument.
    inline void get_coordinates(int n,std::vector<double>& x) const {geometry->get_coordinates(n,x);};
    /// This method collects the coordinates of all active events in the spacetime, stored in the argument as a vector of a length equal to the product of the background dimension and the number of active events.
    void get_coordinates(std::vector<double>&) const;
    /// This method computes the maximum and minimum energy over the set of active events and sets the argument's two elements to them. 
    void get_energy_extrema(std::pair<double,double>&) const;
    /// This method computes the maximum and minimum deficiency over the set of active events and sets the argument's two elements to them. 
    void get_deficiency_extrema(std::pair<double,double>&) const;
    /// This method returns the square of the distance between the two events that are the method's arguments. 
    inline double get_geometric_distance(int n,int m) const {return geometry->get_squared_distance(n,m,false);};
    /// This method returns the background dimension of the Spacetime::geometry property.
    inline int get_background_dimension() const {return geometry->dimension();};
    /// This method returns the value of Spacetime::state_file.
    inline std::string get_state_file() const {return state_file;};
    /// This method returns the value of Spacetime::hyphantic_ops.
    inline std::string get_hyphantic_operations() const {return hyphantic_ops;};
    /// This method calls the arclength_statistics() method with the same argument.
    inline void get_arclength_statistics(double* output) const {arclength_statistics(output);};
    /// This method returns the value of Spacetime::iterations.
    inline int get_iterations() const {return iterations;};
    /// This method returns the value of Spacetime::error.
    inline double get_error() const {return error;};
    /// This method returns the value of Spacetime::max_iter.
    inline int get_maximum_iterations() const {return max_iter;};
    /// This method returns the value of Spacetime::converged.
    inline bool is_converged() const {return converged;};
  };
}
#endif

