// For XML processing...
#include <pugixml.hpp>

#include <synarmosma/geometry.h>
#include <synarmosma/directed_graph.h>

#include "complex.h"

#ifndef _spacetimeh
#define _spacetimeh

namespace DIAPLEXIS {
  /// A class representing the spacetime, a combination of dynamic geometric and topological properties.
  template<class kind1,class kind2>
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
    /// the sum of the output of representational_energy() and compute_temporal_nonlinearity() 
    /// and \f$2\pi\f$ times the output of Complex::euler_characteristic, less 
    /// the sum of the energy property of all active events.
    double global_deficiency = 0.0;
    /// This property contains all of the spacetime's topology, stored in a 
    /// pointer to an instance of the Complex class.
    Complex<kind1>* skeleton;
    /// This property contains the spacetime's geometry, stored in a pointer to 
    /// an instance of the Synarmosma library's Geometry class.
    SYNARMOSMA::Geometry<kind2>* geometry;
    /// This property contains a list of all of the hyphantic operations 
    /// successfully carried out on the spacetime complex since the start 
    /// of the simulation. 
    std::string hyphantic_ops = "";
    /// This property stores the set of notes which at each iteration will be 
    /// used to determine the hyphantic operations that will be executed, when 
    /// the Spacetime::weaving property is set to be musical. The score is read 
    /// from the file specified by Spacetime::hyphansis_score and parsed by the 
    /// read_parameters() method. 
    std::vector<int>* hyphantic_notes;
    /// This property stores the indices of the hyphantic music scale's 25 notes 
    /// by their key value, stretching from the lowest pitch (most explicative) 
    /// value of 21 (0) to the highest pitch (most implicative) value of 59 (24).
    /// It is used by the musical_hyphansis() method to determine from which bin 
    /// of events to select a base event for a given hyphantic operation, the bins 
    /// having been populated by pitch_mapping().
    int key_mapping[60];

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
    /// This property stores the error tolerance for the 
    /// convergence of the relaxation process and geometry 
    /// solvers.
    double convergence_threshold = 0.00001;
    /// The coupling constant between the topological-geometric torsion 
    /// and the energy in the structure equation.
    double coupling_constant = 0.2;

    /// This property is the initial topological dimension of the 
    /// spacetime complex, i.e. the output of the Complex::dimension() 
    /// method.
    int initial_dim = 4;
    /// This property stores the time (as measured in seconds since 
    /// the beginning of the Unix era on January 1, 1970) at which the 
    /// simulation started. 
    time_t start_time;
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
    /// This positive floating point property controls the degree of tolerance 
    /// for determining the normality of inter-event distances, i.e. the extent to 
    /// which two events have a separation whose absolute value is equal to unity.  
    double abnormality_threshold = 0.1;
    /// This positive floating point property is the criterion for determining if 
    /// two events will be fused together during superposition: to be a candidate 
    /// for such a fusion, the pair of events must have a distance such that the 
    /// absolute value of its square is less than this cut-off.
    double superposition_threshold = 0.05;
    /// This property controls the number of attempts made at event fission by the 
    /// event_fission() method, as a percentage of the total number of active events 
    /// in the spacetime complex and for that reason it must lie between zero and 
    /// and one. 
    double event_fission_percentage = 0.1;
    /// If this property is true the simulation does not carry out any disk 
    /// I/O operations - there is no log file, no checkpoint files are created 
    /// and so forth. 
    bool diskless = false;
    /// This property is used to signal whether or not the array Spacetime::hyphantic_notes 
    /// has been allocated, which should only be the case when performing musical 
    /// hyphansis.  
    bool score_allocated = false;
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
    // Minimal Method
    /// This property is only meaningful for the Geometry_Solver::minimal case and 
    /// controls the number of random event coordinate changes attempted before exiting 
    /// the geometry solver. This property should be positive. 
    int solver_its = 50;
    // Evolutionary Method
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
    // Annealing Method
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
    // Mechanical Method
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
    // Simplex Method
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

    /// This method is called by the advance() method to carry out the hyphansis step of topological change; it assembles a list of active events and their deficiency and then calls either dynamic_hyphansis() or musical_hyphansis().
    void hyphansis();
    /// This method carries out the hyphansis step according to a purely dynamic scheme, based on the magnitude and sign of a event's deficiency. The method returns the number of successful hyphantic operations performed.
    int dynamic_hyphansis();
    /// This method carries out the hyphansis step according to a scheme based on a musical composition (in the Spacetime::hyphansis_score file), using a mapping between the notes and the hyphantic operators; the operators are chosen based on the magnitude and sign of a event's deficiency. The method returns the number of successful hyphantic operations performed.
    int musical_hyphansis();
    /// This method puts active events in an array of 25 bins (the method's argument) according to the deficiency value of event, with the width of the bins set by the pitch intervals of the hyphantic musical scale. 
    void pitch_mapping(std::set<int>*) const;
    /// This method is used in the musical_hyphansis() method and converts a (higher-pitched) key - a musical note - into an implicative hyphantic operator (the string output) along with the parameter value for its use, if necessary (the second argument).  
    std::string implicative_scale(int,std::vector<double>&) const;
    /// This method is used in the musical_hyphansis() method and converts a (lower-pitched) key - a musical note - into an explicative hyphantic operator (the string output) along with the parameter value for its use, if necessary (the second argument).  
    std::string explicative_scale(int,std::vector<double>&) const;
    /// This method is used by the dynamic_hyphansis() method and assigns an implicative hyphantic operator to the method's unique argument, based in part on the current value of Spacetime::iterations as well as hasard, through pseudo-random numbers.
    void implication(std::string&) const;
    /// This method is used by the dynamic_hyphansis() method and assigns an explicative hyphantic operator to the method's unique argument, based in part on the current value of Spacetime::iterations as well as hasard, through pseudo-random numbers.
    void explication(std::string&) const;
    /// This method adds a new event to the spacetime, where the argument is the coordinate vector for the new event; the method returns the index of the new event. 
    int event_addition(const std::vector<double>&);
    /// This method adds a new event to the spacetime, where the argument is a set of parent events for the new event; the method returns the index of the new event. 
    int event_addition(const std::set<int>&);
    /// This method adds a new event to the spacetime, where the argument is the parent event for the new event; the method returns the index of the new event.     
    int event_addition(int);
    /// This method deletes the event whose index is the argument, by setting its Event::active property to false. The method returns false if the event was already inactive, otherwise true. 
    bool event_deletion(int);
    /// This method fuses the second argument into the first, both interpreted as indices of events; it returns true if the fusion succeeded and false otherwise.
    bool event_fusion(int,int);
    /// This method carries out an event fusion that is designed to twist the topology of the complex and create non-orientability; it returns true if the event fusion succeeded, false otherwise.
    bool event_twist();
    /// This method attempts a circumvolution using active boundary edges; it is normally only called when the global energy reaches a critical threshold. The method returns true if it is successful.
    bool circumvolution();
    /// This method seeks to fuse together two d-simplices (d > 0), one of which contains the event whose index is the method's argument; it returns true if successful.
    bool circumvolution(int);
    /// This method looks for active 1-simplices connected to a particular event (the method's first argument) and whose length exceeds a certain threshold (the second argument); if the method finds such edges, it attempts to delete one and returns true if successful.
    bool contraction(int,double);
    /// This method attempts to remove a 1-simplex from the event whose index is the argument, favouring edges whose length is much greater than unity and whose connecting event also has an excessive degree. The method returns true if it succeeds in deleting an edge.
    bool compensation_m(int);
    /// This method adds or removes a 1-simplex from the event whose index is the argument, depending on the event's degree which has a goal of being twice the background dimension. When choosing an edge the method favours shorter over longer edges and ensuring that the partner event has a similar degree deficiency or excess. If an edge is added or removed the method returns true. 
    bool compensation_g(int);
    /// This method creates a new d-simplex (d > 0) containing the event whose index is the method's argument and returning true if successful. This method creates new events for the simplex's other d vertices.
    bool expansion(int);
    /// This method creates a new d-simplex (d > 0) containing the event whose index is the method's first argument and returning true if successful. The second argument is the percentage of the new simplex's vertex-events which should be newly created. 
    bool expansion(int,double);
    /// This method examines all the active neighbours with non-zero deficiency that are connected to the event whose index is the method's argument. Two such neighbours are chosen and along with the base event a 2-simplex is created; if the simplex creation is successful, the method returns true. 
    bool foliation_m(int);
    /// This method examines all the active neighbours with non-zero deficiency that are connected to the event whose index is the method's argument. One such neighbour is chosen and the connecting 1-simplex is deleted; if the edge deletion is successful, the method returns true. 
    bool foliation_x(int);
    /// This method calculates the edges connecting the event whose index is the method's argument and selects one of those which is used in a d-simplex (d > 1) so that it can be deleted. If the edge deletion is successful the method returns true.
    bool reduction(int);
    /// This method deletes an event drawn from the base event whose index is the method's first argument and its neighbours, assuming their deficiency property exceeds the value of the second argument. The method also deletes all the higher-dimensional simplices depending on this event; it returns true if the deletion is successful.
    bool amputation(int,double);
    /// This method fuses an active event with non-zero deficiency and which lies within a distance L (the second argument) of the base event (the first argument), to the base event. It returns true if the fusion is successful. 
    bool fusion_x(int,double);
    /// This method accepts as its unique argument the index of an event and chooses at random one of the active events among its neighbours and then fuses this neighbour to the original event. It returns true if this fusion is successful.
    bool fusion_m(int);
    /// This method causes an event to undergo fission into one or more events, duplicating the original event's neighbour connections depending on the value of the second argument which should lie between zero and unity. If the first argument is non-negative, it is the index of the event which undergoes fission otherwise a random active event is chosen. The method returns true if the fission is successful.
    bool fission(int,double);
    /// This method inflates a d-simplex into an n-simplex, where n > d. If the first argument is non-negative, it is assumed to be the index of an event which must belong to the d-simplex, otherwise an active d-simplex is chosen at random. The second argument argument is a creativity index controlling whether or not new events are created to supply the other vertices of the n-simplex. The method returns true if it succeeds in inflating a d-simplex.
    bool inflation(int,double);
    /// This method accepts the index of an active event as its argument and if the topological dimension of this event is greater than unity, it deletes a d-simplex (d > 1) containing this event, thereby reducing the event's topological dimension to d-1. 
    bool deflation(int);
    /// This method attempts to create a hole or perforation in a simplex. If the first argument is non-negative the method tries to create this hole in a d-simplex (d > 1) containing the base event (the first argument), where d is the second argument. If the first argument is negative, the method is global and the second argument is ignored; it tries to find a random d-simplex (d > 1) in which to create a hole. The method returns true if it is successful in creating such a perforation in the spacetime complex. 
    bool perforation(int,int);
    /// This method needs to loop over all active events and then find those which are capable of adding another event at a distance of (roughly) unity from the base event (the argument) and which is orthogonal to the base event's current set of edges.
    bool correction(int);
    /// This method constructs new neighbour events w_i for the base event whose index is the method's argument and which are unit distance from the base event and orthogonal to the base event's existing edges, if possible. It returns true if it is successful and false otherwise.
    bool germination(int);
    /// This method accepts as its argument the index of an event and, if this event is a member of a 2-simplex, carries out a Δ => Y transformation, returning true if it is successful.
    bool stellar_addition(int);
    /// This method accepts as its argument the index of an event and, if this event has a degree of at least three and is not a member of a d-simplex (d > 1), carries out a Y => Δ transformation, returning true if it is successful.
    bool stellar_deletion(int);

    /// This method first deletes events in the centre of the spacetime and then inserts a combinatorial black hole (i.e. a compact, highly-entwined knot) in the location of these inactive events. The first argument is the event on which the knot should be centred, the second argument the approximate radius of the knot and the final argument its dimensionality. The method returns true if it succeeds.
    bool interplication(int,double,int);
    /// This method deletes 1-simplices whose length exceeds the method's argument, if it satisfies a Boltzmann criterion. If the spacetime complex is disconnected it then adds the fewest and shortest possible edges to re-connect it. The method returns the number of 1-simplices that were originally deleted.
    int compression(double);
    /// This method fuses together pairs of events when the absolute value of their squared distance is less than the square of Spacetime::superposition_threshold, using the event_fusion() method; the method returns the number of pairs of events fused together.
    int superposition_fusion();
    /// This method does the opposite of superposition_fusion() - it randomly selects N active events, where N is equal to product of Spacetime::event_fission_percentage and the total number of active events. These events must then satisfy a Boltzmann criterion with a temperature equal to 100, after which they undergo fission using the fission() method.
    int event_fission();
    /// This method ensures that the spacetime complex is connected, consistent and that it satisfies the entailment axiom of simplicial complexes. If the argument is false, the Complex::simplicial_implication() method is called at the beginning of the method.
    void regularization(bool);

    /// This method tests whether the d-simplex specified by the method's two arguments (the dimension d and index of the simplex in Complex::simplices[d]) can be realized given the geometry of the simplex's events; the method returns true if the d-simplex has a positive squared volume and false otherwise. 
    bool realizable(int,int) const;
    /// This method computes the volume of all of the spacetime's active d-simplices (d > 0) using the current geometry; to handle the 1-simplices, it calls compute_lengths().
    void compute_volume();
    /// This method computes the length of all of the spacetime's active 1-simplices (edges) using the current geometry.
    void compute_lengths();
    /// This method computes the Event::obliquity property for all of the spacetime's active events using the current geometry.
    void compute_obliquity();
    /// This method returns the abormality of the lengths of the spacetime's 1-simplices, while ignoring those 1-simplices deemed "flexible", which are listed in the method's unique argument.
    double compute_abnormality(const std::vector<int>&) const;
    /// This method returns the abormality of the lengths of the spacetime's 1-simplices, while ignoring those 1-simplices deemed "flexible", which are listed in the method's second argument; the first argument is the vector of current event coordinates. 
    double compute_abnormality(const std::vector<double>&,const std::vector<int>&) const;
    /// This method computes the gradient of the edge length abnormality with respect to the event geometry; the gradient is written to the method's first argument, while the third argument is the vector listing which edges are deemed "flexible". If the second argument is true this method computes the negative gradient.   
    void compute_geometric_gradient(std::vector<double>&,bool,const std::vector<int>&);
    /// This method computes the force on each event in the damped spring mechanical model used by the mechanical_solver() method. The method's arguments are offset mappings for distinguishing active and inactive events (the first two arguments), a vector of energy values for active events, the current event positions and finally the calculated forces on these events. 
    void mechanical_force(const std::vector<int>&,const std::vector<int>&,const std::vector<double>&,const std::vector<double>&,double*) const;
    /// This method optimizes the spacetime geometry using a damped spring mechanical model for the 1-skeleton.
    void mechanical_solver();
    /// This method optimizes the spacetime geometry using a simulated annealing algorithm.
    void annealing_solver();
    /// This method optimizes the spacetime geometry using the Nelder-Mead downhill simplex algorithm.
    void simplex_solver();
    /// This method optimizes the spacetime geometry using an evolutionary algorithm, with a population of geometries undergoing selection based on the value of the Spacetime::error property.
    void evolutionary_solver();
    /// This method carries out the geometry optimization, calling other methods (e.g. mechanical_solver(), simplex_solver() etc.) as necessary.
    void optimize();
    /// This method calculates and returns the minimum of the absolute value of the distances separating the two vectors of events that are its first arguments. The final argument will contain the index of the two events that are the closest. 
    double minimize_lengths(const std::vector<int>&,const std::vector<int>&,int*) const;
    /// This method computes the current value of the terms in the structure equation for the spacetime, including setting the value of the Spacetime::error property. 
    void structural_deficiency();
    /// This method accepts as its first argument the index of an event and then calculates the "total" (i.e. all events that lie in its past and future, including those which are not direct neighbours of the original event) past and future light cones of this event, which are output as the method's second and third arguments.
    void compute_total_lightcone(int,std::set<int>&,std::set<int>&) const;
    /// This method computes the anterior and posterior properties of the active spacetime events.
    void compute_lightcones();
    /// This method computes a partially directed graph in which directed edges reflect timelike separation between events, with the graph based on a particular event (the method's second argument). 
    void compute_causal_graph(SYNARMOSMA::Directed_Graph*,int) const;
    /// This method computes and returns the extent to which the chrono-geometry matches the topology in the vicinity of a given event (the method's argument); a value of zero indicates that the geometry and topology match well, while a positive value indicates more entwinement is needed while a negative value less entwinement.  
    double compute_temporal_vorticity(int) const;
    /// This method calls compute_total_lightcone() for each active event and then checks for temporal loops, sinks and sources as well as the cyclicity of the causal graph of this event if the event is a source or sink. The method's return value is an index of such exotic causal behaviour in the spacetime as a whole.
    double compute_temporal_nonlinearity() const;
    /// This method computes a measure of how well the spacetime's geometry corresponds to its topological 1-skeleton and returns this value. If the argument is true the spatial embedding is weighted by the absolute value of the length of the 1-simplices.  
    double representational_energy(bool) const;

    /// This method carries out the various global operations that must be performed at each relaxation step.
    bool global_operations();
    /// This method computes the colours to be used for the events when visualizing the spacetime; a set of three RGB 8 bit unsigned integers (lying between 0 and 255) is used for each event - the inactive events are made invisible - and the criterion used for the colouring is either the energy (second argument true) or deficiency property (false) of the event. 
    void compute_colours(std::vector<unsigned char>&,bool) const;
    /// This method computes the initial configuration of the spacetime based on its parameters.
    void build_initial_state();
    /// This method writes out information on the current configuration of the Spacetime instance to the log file, whose name is stored in the Spacetime::log_file property.
    void write_log() const;
    /// This method calls clear() and then reads an instance of the Spacetime class from the binary file whose name is the method's unique argument.
    void read_state(const std::string&);
    /// This method writes the Spacetime instance to a binary file; the file's name is the method's argument if it has one, otherwise it uses the Spacetime::state_file property.
    void write_state(const std::string& = "") const;
    /// This method reads the run-time parameters of the Spacetime class from an XML file whose name is the method's unique argument.
    void read_parameters(const std::string&);
    /// This method adjusts the dimension of each spacetime event according to its simplicial membership (topological dimension) and updates the value of Spacetime::system_size accordingly. It returns the output from calling the Geometry::adjust_dimension method of the Synarmosma library.
    bool adjust_dimension();
    /// This method computes the maximum, minimum and arithmetic mean of the absolute value of the length of the spacetime's active 1-simplices and writes it to the method's argument, an array with three elements. 
    void arclength_statistics(double*) const;
    /// This method reduces memory pressure in the simulation by eliminating inactive events and 1-simplices from the spacetime, if the ratio of active to inactive is less than half for both categories. The method returns the maximum of the current ratios.  
    double condense();
    /// This method attributes a value to a variety of Spacetime properties and the calls the build_initial_state() method.  
    void initialize();
    /// This method allocates the resources for the Spacetime::geometry and Spacetime::skeleton properties.
    void allocate();
    /// This method first calls the Complex::consistent method, returning false if it does, and then checks that all the floating point properties of the Event class are sensible, returning false otherwise.
    bool consistent() const;
    /// This method calls the clear method on the extended properties of the class and sets Spacetime::iterations to zero and Spacetime::hyphantic_ops to the empty string. 
    void clear();
    /// This method returns true if all of the spacetime's d-simplices (d >= 0) have a modified property equal to false; it will otherwise return false. 
    bool clean() const;
    /// This method computes the value of Spacetime::error in the spacetime's current state, then sets all of the d-simplices (d >= 0) as modified and re-computes the Spacetime::error property; if there is any difference in the two values it returns false. 
    bool correctness();
    /// This method steps forward in the relaxation process, successively calling the hyphansis(), regularization(), condense() and global_operations() methods. The output of the latter is the output of this method.
    bool advance();
    /// If Spacetime::reversible is true, this method will call read_state() using the previous step's snapshot on disk to restore the simulation. 
    void fallback();
    /// This method calls clear(), re-reads the parameters from its first argument (an XML file) and then calls the initialize() method. The second argument determines whether or not it reuses the existing random number seed.
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
    /// This method is similar to evolve(), in that it carries out a simulation of the spacetime but in this case with the goal of adapting the topology to the geometric constraints (e.g. dimensionality). The spacetime is assumed to start in a state with many edges, some of which must be deleted to match the geometry and degree distribution for a regular space. The method returns the value of Spacetime::error at completion. 
    double chorogenesis();
    /// This method writes (using an offset to handle inactive events and edges) the coordinates of the spacetime's active events in the first argument, an edge table (listing the two event indices) in the second argument and the third argument is a pair of integers containing the number of active events and active edges.
    void export_visual_data(std::vector<float>&,std::vector<int>&,std::pair<int,int>&) const;
    /// This method writes the coordinates (second argument), edge table (third argument) and the total number of events and edges (the pair of integers that is the fourth argument). The first argument is a vector of colours for each event with inactive events and edges set to be invisible; the final argument determines the basis for the colouring scheme, using energy or deficiency values.
    void export_visual_data(std::vector<float>&,std::vector<float>&,std::vector<int>&,std::pair<int,int>&,bool) const;
    /// This method sets the value of Spacetime::checkpoint_frequency to the method's argument.
    void set_checkpoint_frequency(int);
    /// This method sets the second argument to the coordinates of the event whose index is the first argument.
    void get_coordinates(int,std::vector<double>&) const;
    /// This method collects the coordinates of all active events in the spacetime, stored in the argument as a vector of a length equal to the product of the background dimension and the number of active events.
    void get_coordinates(std::vector<double>&) const;
    /// This method computes the maximum and minimum energy over the set of active events and sets the argument's two elements to them. 
    void get_energy_extrema(std::pair<double,double>&) const;
    /// This method computes the maximum and minimum deficiency over the set of active events and sets the argument's two elements to them. 
    void get_deficiency_extrema(std::pair<double,double>&) const;
    /// This method returns the square of the distance between the two events that are the method's arguments. 
    double get_geometric_distance(int,int) const;
    /// This method returns the background dimension of the Spacetime::geometry property.
    int get_background_dimension() const;
    /// This method returns the value of Spacetime::state_file.
    std::string get_state_file() const;
    /// This method returns the value of Spacetime::hyphantic_ops.
    std::string get_hyphantic_operations() const;
    /// This method furnishes a public version of the arclength_statistics() method, using the same argument.
    void get_arclength_statistics(double* output) const;
    /// This method returns the value of Spacetime::iterations.
    int get_iterations() const;
    /// This method returns the value of Spacetime::error.
    double get_error() const;
    /// This method returns the value of Spacetime::max_iter.
    int get_maximum_iterations() const;
    /// This method returns the value of Spacetime::converged.
    bool finished() const;
  };

  template<class kind1,class kind2>
  inline void Spacetime<kind1,kind2>::set_checkpoint_frequency(int n) 
  {
    checkpoint_frequency = n;
  }

  template<class kind1,class kind2>
  inline void Spacetime<kind1,kind2>::get_coordinates(int n,std::vector<double>& x) const 
  {
    geometry->get_coordinates(n,x);
  }

  template<class kind1,class kind2>
  inline double Spacetime<kind1,kind2>::get_geometric_distance(int n,int m) const
  {
    return geometry->get_squared_distance(n,m,false);
  }

  template<class kind1,class kind2>
  inline int Spacetime<kind1,kind2>::get_background_dimension() const
  {
    return geometry->dimension();
  }

  template<class kind1,class kind2>
  inline std::string Spacetime<kind1,kind2>::get_state_file() const 
  {
    return state_file;
  }

  template<class kind1,class kind2>
  inline std::string Spacetime<kind1,kind2>::get_hyphantic_operations() const
  {
    return hyphantic_ops;
  }

  template<class kind1,class kind2>
  inline void Spacetime<kind1,kind2>::get_arclength_statistics(double* output) const 
  {
    arclength_statistics(output);
  }

  template<class kind1,class kind2>
  inline int Spacetime<kind1,kind2>::get_iterations() const
  {
    return iterations;
  }

  template<class kind1,class kind2>
  inline double Spacetime<kind1,kind2>::get_error() const
  {
    return error;
  }

  template<class kind1,class kind2>
  inline int Spacetime<kind1,kind2>::get_maximum_iterations() const 
  {
    return max_iter;
  }

  template<class kind1,class kind2>
  inline bool Spacetime<kind1,kind2>::finished() const 
  {
    return converged;
  }
}
#endif

