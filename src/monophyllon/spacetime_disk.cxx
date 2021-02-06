#include "spacetime.h"

using namespace DIAPLEXIS;

void Spacetime::read_parameters(const std::string& filename)
{
  int q,tvalue,n_is = 0,n_so = 0,D = 0;
  unsigned long rs;
  bool euclidean = false,relational = false,uniform = false;
  std::string name,value;
  pugi::xml_document pfile;
  pugi::xml_node global,gsolver;

  // Open the file
  if (!(pfile.load_file(filename.c_str()))) throw std::invalid_argument("Unable to parse parameter file!");

  global = pfile.child("Parameters").child("Global");
  for(pugi::xml_node params = global.first_child(); params; params = params.next_sibling()) {
    name = params.name();
    value = params.first_child().value();

    if (name == "InitialEvents") {
      initial_size = std::stoi(value);
    }
    else if (name == "InitialState") {
      std::transform(value.begin(),value.end(),value.begin(),::toupper);
      if (value == "CARTESIAN") {
        initial_state = Initial_Topology::cartesian;
      }
      else if (value == "RANDOM") {
        initial_state = Initial_Topology::random;
      }
      else if (value == "SINGLETON") {
        initial_state = Initial_Topology::singleton;
      }
      else if (value == "MONOPLEX") {
        initial_state = Initial_Topology::monoplex;
      }
      else if (value == "DISKFILE") {
        initial_state = Initial_Topology::diskfile;
      }
      else {
        throw std::runtime_error("Unrecognized initial state " + value + "!");
      }
      n_is++;
    }
    else if (name == "InputFile") {
      input_file = value;
    }
    else if (name == "EuclideanGeometry") {
      tvalue = std::stoi(value);
      assert(tvalue == 0 || tvalue == 1);
      euclidean = (tvalue == 1) ? true : false;
    }
    else if (name == "RelationalGeometry") {
      tvalue = std::stoi(value);
      assert(tvalue == 0 || tvalue == 1);
      relational = (tvalue == 1) ? true : false;
    }
    else if (name == "DimensionalUniformity") {
      tvalue = std::stoi(value);
      assert(tvalue == 0 || tvalue == 1);
      uniform = (tvalue == 1) ? true : false;
    }
    else if (name == "RandomSeed") {
      rs = (unsigned) std::stol(value);
      if (rs == 0) rs = (unsigned) std::time(nullptr);
      skeleton->RND->set_seed(rs);
    }
    else if (name == "ConvergenceThreshold") {
      convergence_threshold = std::stod(value);
    }
    else if (name == "CouplingConstant") {
      coupling_constant = std::stod(value);
    }
    else if (name == "TopologicalRadius") {
      skeleton->topological_radius = std::stoi(value);
    }
    else if (name == "BackgroundDimension") {
      D = std::stoi(value);
    }
    else if (name == "EdgeProbability") {
      edge_probability = std::stod(value);
    }
    else if (name == "AbnormalityThreshold") {
      abnormality_threshold = std::stod(value);
    }
    else if (name == "SuperpositionThreshold") {
      superposition_threshold = std::stod(value);
    }
    else if (name == "MaximumIterations") {
      max_iter = std::stoi(value);
    }
    else if (name == "InitialDimension") {
      initial_dim = std::stoi(value);
    }
    else if (name == "PerturbTopology") {
      tvalue = std::stoi(value);
      assert(tvalue == 0 || tvalue == 1);
      perturb_topology = (tvalue == 1) ? true : false;
    }
    else if (name == "PerturbGeometry") {
      tvalue = std::stoi(value);
      assert(tvalue == 0 || tvalue == 1);
      perturb_geometry = (tvalue == 1) ? true : false;
    }
    else if (name == "PerturbEnergy") {
      tvalue = std::stoi(value);
      assert(tvalue == 0 || tvalue == 1);
      perturb_energy = (tvalue == 1) ? true : false;
    }
    else if (name == "MemoryFootprint") {
      std::transform(value.begin(),value.end(),value.begin(),::toupper);
      if (value == "HIGH") {
        high_memory = true;
      }
      else if (value == "LOW") {
        high_memory = false;
      }
      else {
        throw std::runtime_error("Unrecognized memory footprint " + value + "!");
      }
    }
    else if (name == "Hyphansis") {
      std::transform(value.begin(),value.end(),value.begin(),::toupper);
      if (value == "MUSICAL") {
        weaving = Hyphansis::musical;
      }
      else if (value == "DYNAMIC") {
        weaving = Hyphansis::dynamic;
      }
      else {
        throw std::runtime_error("Unrecognized hyphansis method " + value + "!");
      }
    }
    else if (name == "HyphansisScore") {
      hyphansis_score = value;
    }
    else if (name == "ParityMutation") {
      parity_mutation = std::stod(value);
    }
    else if (name == "HomologyMethod") {
      std::transform(value.begin(),value.end(),value.begin(),::toupper);
      if (value == "GAP") {
        skeleton->set_homology_method(SYNARMOSMA::Homology::Method::gap);
      }
      else if (value == "NATIVE") {
        skeleton->set_homology_method(SYNARMOSMA::Homology::Method::native);
      }
      else {
        throw std::runtime_error("Unrecognized homology method " + value + "!");
      }
    }
    else if (name == "HomologyBase") {
      std::transform(value.begin(),value.end(),value.begin(),::toupper);
      if (value == "INT") {
        skeleton->set_homology_field(SYNARMOSMA::Homology::Field::int32);
      }
      else if (value == "NTL::ZZ") {
        skeleton->set_homology_field(SYNARMOSMA::Homology::Field::multiprecision);
      }
      else if (value == "GF2") {
        skeleton->set_homology_field(SYNARMOSMA::Homology::Field::mod2);
      }
      else {
        throw std::runtime_error("Unrecognized homology base field " + value + "!");
      }
    }
    else if (name == "Compressible") {
      tvalue = std::stoi(value);
      assert(tvalue == 0 || tvalue == 1);
      compressible = (tvalue == 1) ? true : false;
    }
    else if (name == "Superposable") {
      tvalue = std::stoi(value);
      assert(tvalue == 0 || tvalue == 1);
      superposable = (tvalue == 1) ? true : false;
    }
    else if (name == "CheckpointFrequency") {
      checkpoint_frequency = std::stoi(value);
    }
  }

  gsolver = pfile.child("Parameters").child("GeometrySolver");
  for(pugi::xml_node params = gsolver.first_child(); params; params = params.next_sibling()) {
    name = params.name();
    value = params.first_child().value();

    if (name == "SolverType") {
      std::transform(value.begin(),value.end(),value.begin(),::toupper);
      if (value == "MINIMAL") {
        solver = Geometry_Solver::minimal;
      }
      else if (value == "MECHANICAL") {
        solver = Geometry_Solver::mechanical;
      }
      else if (value == "EVOLUTIONARY") {
        solver = Geometry_Solver::evolutionary;
      }
      else if (value == "ANNEALING") {
        solver = Geometry_Solver::annealing;
      }
      else if (value == "SIMPLEX") {
        solver = Geometry_Solver::simplex;
      }
      else {
        throw std::runtime_error("Unrecognized geometry solver " + value + "!");
      }
      n_so++;
    }
    else if (name == "SolverIterations") {
      solver_its = std::stoi(value);
    }
    else if (name == "MaximumGenerations") {
      ngenerations = std::stoi(value);
    }
    else if (name == "MaximumJousts") {
      njousts = std::stoi(value);
    }
    else if (name == "PoolSize") {
      pool_size = std::stoi(value);
    }
    else if (name == "ThermalSweeps") {
      thermal_sweep = std::stoi(value);
    }
    else if (name == "AnnealingSteps") {
      annealing_steps = std::stoi(value);
    }
    else if (name == "ThermalVariance") {
      thermal_variance = std::stod(value);
    }
    else if (name == "ThermalizationCriterion") {
      thermalization = std::stod(value);
    }
    else if (name == "IntegrationEngine") {
      std::transform(value.begin(),value.end(),value.begin(),::toupper);
      if (value == "EULER") {
        engine = Integrator::euler;
      }
      else if (value == "RK4") {
        engine = Integrator::rk4;
      }
      else {
        throw std::runtime_error("Unrecognized integration engine " + value + "!");
      }
    }
    else if (name == "StepSize") {
      step_size = std::stod(value);
    }
    else if (name == "MaximumIntegrationSteps") {
      max_int_steps = std::stoi(value);
    }
    else if (name == "DampingConstant") {
      damping_constant = std::stod(value);
    }
    else if (name == "SpringConstant") {
      spring_constant = std::stod(value);
    }
    else if (name == "RepulsionConstant") {
      repulsion_constant = std::stod(value);
    } 
    else if (name == "ConjugateGradientRefinement") {
      tvalue = std::stoi(value);
      assert(tvalue == 0 || tvalue == 1);
      cgradient_refinement = (tvalue == 1) ? true : false;
    }
    else if (name == "MaximumConjugateGradientSteps") {
      max_CG_steps = std::stoi(value);
    }
    else if (name == "MaximumLineSolverSteps") {
      max_LS_steps = std::stoi(value);
    }
    else if (name == "EdgeFlexibilityThreshold") {
      edge_flexibility_threshold = std::stod(value);
    }
    else if (name == "ReflectionCoefficient") {
      simplex_alpha = std::stod(value);
    }
    else if (name == "ExpansionCoefficient") {
      simplex_gamma = std::stod(value);
    }
    else if (name == "ContractionCoefficient") {
      simplex_rho = std::stod(value);
    }
    else if (name == "ShrinkageCoefficient") {
      simplex_sigma = std::stod(value);
    }
  }

  // Now a series of tests to make sure that the parameters
  // aren't entirely crazy...
  assert(convergence_threshold > std::numeric_limits<double>::epsilon());
  assert(coupling_constant > std::numeric_limits<double>::epsilon());
  assert(abnormality_threshold > std::numeric_limits<double>::epsilon());
  assert(max_iter >= 0);
  assert(skeleton->topological_radius > 0);
  assert(checkpoint_frequency > 0);
  // One and only one initial state should be chosen...
  assert(n_is == 1);
  // And similarly for the solver type...
  assert(n_so == 1);
  if (superposable) assert(superposition_threshold > std::numeric_limits<double>::epsilon());
  // Finally for the background dimension...
  assert(D > 0);

  if (weaving == Hyphansis::dynamic) {
    assert(parity_mutation > -std::numeric_limits<double>::epsilon());
    assert((parity_mutation - 1.0) < std::numeric_limits<double>::epsilon());
  }
  else {
    // Make sure the score file exists and has the right structure...
    int i,v,n,its,voice = -1,mv = -1,nsilent = 0;
    bool prolonged = false;
    std::string line;
    std::vector<std::string> elements;
    std::set<int> voices,legal_notes;
    std::ifstream mscore;
    const char delimiter = '/';

    score_allocated = true;
    hyphantic_notes = new std::vector<int>[max_iter+1];

    // Implicative notes...
    legal_notes.insert(45); legal_notes.insert(47); legal_notes.insert(48); legal_notes.insert(49);
    legal_notes.insert(50); legal_notes.insert(52); legal_notes.insert(53); legal_notes.insert(54);
    legal_notes.insert(56); legal_notes.insert(57); legal_notes.insert(58); legal_notes.insert(59);

    // Explicative notes...
    legal_notes.insert(21); legal_notes.insert(23); legal_notes.insert(25); legal_notes.insert(26);
    legal_notes.insert(28); legal_notes.insert(29); legal_notes.insert(30); legal_notes.insert(32);
    legal_notes.insert(33); legal_notes.insert(34); legal_notes.insert(35); legal_notes.insert(37);

    // The unique neutral note...
    legal_notes.insert(40);

    mscore.exceptions(std::ifstream::badbit);
    // Open the file containing the hyphantic score 
    try {
      mscore.open(hyphansis_score);

      // Now read the measure that corresponds to this iteration and sheet...
      while(mscore.good()) {
        getline(mscore,line);
        // If the line is empty or doesn't contain a forward slash, ignore it...
        if (line.empty()) continue;
        if (line.find(delimiter) == std::string::npos) continue;
        // Tokenize the line at the forward slash...
        SYNARMOSMA::split(line,delimiter,elements);
        assert(elements.size() == 3);
        its = std::stoi(elements[0]) - 1;
        if (its > max_iter) {
          prolonged = true;
          if (its > mv) mv = its;
          continue;
        }
        // Now read the voice...
        v = std::stoi(elements[1]);
        if (voices.count(v) == 0) voices.insert(v);
        // Initialize the voice variable to be the first voice we find in the score file; 
        // note that we assume the voice indices to be non-negative...
        if (voice < 0) voice = v;
        // There should only be one voice in the case of the Monophyllon version of this library...
        if (v != voice) continue;
        // Check the legality of the note...
        n = std::stoi(elements[2]);
        assert(legal_notes.count(n) > 0);
        hyphantic_notes[its].push_back(n);
      }
    }
    catch (const std::ifstream::failure& e) {
      std::cout << "Error in opening or reading the " << hyphansis_score << " file!" << std::endl;
    }
    // Close the score file
    mscore.close();

    if (voices.size() > 1) std::cout << "Warning: The musical score is polyphonic, with " << voices.size() << " voices, but only voice " << voice << " will be used!" << std::endl;
    if (prolonged) std::cout << "Warning: The length of this musical score (" << mv << ") exceeds the maximum number of iterations (" << max_iter << ")!" << std::endl;
    // Check for iterations during which there is no hyphansis...
    for(i=0; i<=max_iter; ++i) {
      if (hyphantic_notes[i].empty()) nsilent++;
    }
    if (nsilent > 0) std::cout << "Warning: There are " << nsilent << " relaxation steps during which there will be no hyphantic operations!" << std::endl;
  }

  geometry->initialize(euclidean,relational,uniform,high_memory,D);

  if (initial_state == Initial_Topology::random) {
    assert(edge_probability > std::numeric_limits<double>::epsilon() && (edge_probability - 1.0) < -std::numeric_limits<double>::epsilon());
    assert(initial_size > 1);
  }
  else if (initial_state == Initial_Topology::cartesian) {
    assert(initial_size > 1);
    // If initial_size != n^dimension, n \in Z+, we have a problem
    std::vector<std::pair<long,int> > factors;
    SYNARMOSMA::factorize(initial_size,factors);
    for(q=0; q<(signed) factors.size(); ++q) {
      if (factors[q].second % D != 0) throw std::invalid_argument("Grid size inconsistent with the background dimension!");
    }
  }
  else if (initial_state == Initial_Topology::monoplex) {
    assert(initial_dim > 0);
    assert(initial_dim <= Complex::ND);
    initial_size = 1 + initial_dim;
  }
  else if (initial_state == Initial_Topology::singleton) {
    initial_size = 1;
  }

  if (solver == Geometry_Solver::minimal) {
    assert(solver_its > 0);
  }
  else if (solver == Geometry_Solver::evolutionary) {
    assert(pool_size > 0);
    assert(ngenerations > 0);
    assert(njousts > 0);
  }
  else if (solver == Geometry_Solver::annealing) {
    assert(annealing_steps > 0);
    assert(thermalization > std::numeric_limits<double>::epsilon());
    assert(thermal_variance > std::numeric_limits<double>::epsilon());
    assert(thermal_sweep > 0);
  }
  else if (solver == Geometry_Solver::mechanical) {
    if (!relational && !uniform) throw std::invalid_argument("Mechanical solver requires dimensional uniformity when geometry is absolute!");
    assert(step_size > std::numeric_limits<double>::epsilon() && (step_size - 1.0) < -std::numeric_limits<double>::epsilon());
    assert(spring_constant < -std::numeric_limits<double>::epsilon());
    assert(repulsion_constant > std::numeric_limits<double>::epsilon());
    assert(damping_constant > std::numeric_limits<double>::epsilon());
    assert(max_int_steps > 0);
    if (cgradient_refinement) {
      assert(max_CG_steps > 0);
      assert(max_LS_steps > 0);
      assert(edge_flexibility_threshold > std::numeric_limits<double>::epsilon());
    }
  }
  else if (solver == Geometry_Solver::simplex) {
    assert(simplex_alpha > std::numeric_limits<double>::epsilon());
    assert((simplex_gamma - 1.0) > std::numeric_limits<double>::epsilon());
    assert(simplex_gamma > simplex_alpha);
    assert(simplex_rho > std::numeric_limits<double>::epsilon() && (simplex_rho - 1.0) < -std::numeric_limits<double>::epsilon());
    assert(simplex_sigma > std::numeric_limits<double>::epsilon() && (simplex_sigma - 1.0) < -std::numeric_limits<double>::epsilon());
  }
  original_state = initial_state;
}

void Spacetime::write_log() const
{
  int i,v,nn,ne,ntime,nspace,nf,np,nat,vd_avg,vd_min,vd_max;
  int nsource,nsink,nnull,nch,in1,sum,mx,mn;
  double w,avg_ven,max_ven,min_ven,vdata[3];
  bool atemporal;
  static bool fcall = true;
  std::string nvalue;
  pugi::xml_document logfile;
  pugi::xml_node rstep,sheet,atom;
  SYNARMOSMA::hash_map::const_iterator qt;
  std::set<int> S,edge;
  std::set<int>::const_iterator it;
  std::string group_string,bvalue[] = {"False","True"};
  const int Nv = (signed) skeleton->events.size();
  const int Ne = (signed) skeleton->simplices[1].size();
  int vdimension[Nv];

  if (fcall) {
    char hostname[80];
    std::string sname[] = {"RANDOM","MONOPLEX","CARTESIAN","SINGLETON","DISKFILE"};
    std::string hname[] = {"MINIMAL","MECHANICAL","EVOLUTIONARY","ANNEALING","SIMPLEX"};

    gethostname(hostname,80);

    std:: string cdate = std::string(std::ctime(&start_time));
    cdate.erase(std::remove(cdate.begin(),cdate.end(),'\n'),cdate.end());

    std::ofstream s(log_file,std::ios::out | std::ios::trunc);
    s << "<?xml version=\"1.0\" ?>" << std::endl;
    s << "<LogFile>" << std::endl;
    s << "<Library>Diaplexis</Library>" << std::endl;
    s << "<Hostname>" << hostname << "</Hostname>" << std::endl;
    s << "<StartTime>" << cdate << "</StartTime>" << std::endl;
    s << "<CompileTimeParameters>" << std::endl;
    s << "<MaximumSimplicialDimension>" << Complex::ND << "</MaximumSimplicialDimension>" << std::endl;
    s << "<AtomicPropositions>" << SYNARMOSMA::Proposition::get_clause_size() << "</AtomicPropositions>" << std::endl;
    s << "</CompileTimeParameters>" << std::endl;
    s << "<RunTimeParameters>" << std::endl;
    s << "<Global>" << std::endl;
    if (initial_state != Initial_Topology::singleton && initial_state != Initial_Topology::monoplex) s << "<InitialEvents>" << initial_size << "</InitialEvents>" << std::endl;
    s << "<MaximumIterations>" << max_iter << "</MaximumIterations>" << std::endl;
    s << "<EuclideanGeometry>" << bvalue[geometry->get_euclidean()] << "</EuclideanGeometry>" << std::endl;
    s << "<RelationalGeometry>" << bvalue[geometry->get_relational()] << "</RelationalGeometry>" << std::endl;
    s << "<DimensionalUniformity>" << bvalue[geometry->get_uniform()] << "</DimensionalUniformity>" << std::endl;
    s << "<BackgroundDimension>" << geometry->dimension() << "</BackgroundDimension>" << std::endl;
    s << "<RandomSeed>" << skeleton->RND->get_seed() << "</RandomSeed>" << std::endl;
    s << "<ConvergenceThreshold>" << convergence_threshold << "</ConvergenceThreshold>" << std::endl;
    s << "<CouplingConstant>" << coupling_constant << "</CouplingConstant>" << std::endl;
    s << "<TopologicalRadius>" << skeleton->topological_radius << "</TopologicalRadius>" << std::endl;
    s << "<AbnormalityThreshold>" << abnormality_threshold << "</AbnormalityThreshold>" << std::endl;
    s << "<Superposable>" << bvalue[superposable] << "</Superposable>" << std::endl;
    if (superposable) s << "<SuperpositionThreshold>" << superposition_threshold << "</SuperpositionThreshold>" << std::endl;
    s << "<Compressible>" << bvalue[compressible] << "</Compressible>" << std::endl;
    if (initial_state == Initial_Topology::diskfile) {
      s << "<InitialState>" << sname[initial_state] << " (" << sname[original_state] << ")</InitialState>" << std::endl;
      s << "<InputFile>" << input_file << "</InputFile>" << std::endl;
    }
    else {
      s << "<InitialState>" << sname[initial_state] << "</InitialState>" << std::endl;
    }
    if (initial_state == Initial_Topology::random) s << "<EdgeProbability>" << edge_probability << "</EdgeProbability>" << std::endl;
    if (initial_state == Initial_Topology::monoplex) s << "<InitialDimension>" << initial_dim << "</InitialDimension>" << std::endl;
    s << "<CheckpointFrequency>" << checkpoint_frequency << "</CheckpointFrequency>" << std::endl;
    if (initial_state == Initial_Topology::cartesian) {
      s << "<PerturbTopology>" << bvalue[perturb_topology] << "</PerturbTopology>" << std::endl;
      s << "<PerturbGeometry>" << bvalue[perturb_geometry] << "</PerturbGeometry>" << std::endl;
    }
    if (initial_state != Initial_Topology::random) s << "<PerturbEnergy>" << bvalue[perturb_energy] << "</PerturbEnergy>" << std::endl;
    if (skeleton->get_homology_method() == SYNARMOSMA::Homology::Method::gap) {
      s << "<HomologyMethod>GAP</HomologyMethod>" << std::endl;
    }
    else if (skeleton->get_homology_method() == SYNARMOSMA::Homology::Method::native) {
      s << "<HomologyMethod>Native</HomologyMethod>" << std::endl;
    }
    if (skeleton->get_homology_field() == SYNARMOSMA::Homology::Field::int32) {
      s << "<HomologyBase>INT</HomologyBase>" << std::endl;
    }
    else if (skeleton->get_homology_field() == SYNARMOSMA::Homology::Field::multiprecision) {
      s << "<HomologyBase>NTL::ZZ</HomologyBase>" << std::endl;
    }
    else if (skeleton->get_homology_field() == SYNARMOSMA::Homology::Field::mod2) {
      s << "<HomologyBase>GF2</HomologyBase>" << std::endl;
    }
    if (weaving == Hyphansis::dynamic) {
      s << "<Hyphansis>DYNAMIC</Hyphansis>" << std::endl;
      s << "<ParityMutation>" << parity_mutation << "</ParityMutation>" << std::endl; 
    }
    else {
      s << "<Hyphansis>MUSICAL</Hyphansis>" << std::endl;
      s << "<HyphansisScore>" << hyphansis_score << "</HyphansisScore>" << std::endl;
    }
    s << "</Global>" << std::endl;
    s << "<Geometry>" << std::endl;
    s << "<SolverType>" << hname[solver] << "</SolverType>" << std::endl;
    if (solver == Geometry_Solver::minimal) {
      s << "<SolverIterations>" << solver_its << "</SolverIterations>" << std::endl;
    }
    else if (solver == Geometry_Solver::evolutionary) {
      s << "<MaximumGenerations>" << ngenerations << "</MaximumGenerations>" << std::endl;
      s << "<MaximumJousts>" << njousts << "</MaximumJousts>" << std::endl;
      s << "<PoolSize>" << pool_size << "<</PoolSize>" << std::endl;
    }
    else if (solver == Geometry_Solver::annealing) {
      s << "<ThermalVariance>" << thermal_variance << "</ThermalVariance>" << std::endl;
      s << "<ThermalSweeps>" << thermal_sweep << "</ThermalSweeps>" << std::endl;
      s << "<AnnealingSteps>" << annealing_steps << "</AnnealingSteps>" << std::endl;
      s << "<ThermalizationCriterion>" << thermalization << "</ThermalizationCriterion>" << std::endl;
    }
    else if (solver == Geometry_Solver::mechanical) {
      if (engine == Integrator::euler) {
        s << "<IntegrationEngine>EULER</IntegrationEngine>" << std::endl;
      }
      else {
        s << "<IntegrationEngine>RK4</IntegrationEngine>" << std::endl;
      }
      s << "<StepSize>" << step_size << "</StepSize>" << std::endl;
      s << "<MaximumIntegrationSteps>" << max_int_steps << "</MaximumIntegrationSteps>" << std::endl;
      s << "<DampingConstant>" << damping_constant << "</DampingConstant>" << std::endl;
      s << "<SpringConstant>" << spring_constant << "</SpringConstant>" << std::endl;
      s << "<RepulsionConstant>" << repulsion_constant << "</RepulsionConstant>" << std::endl;
      s << "<ConjugateGradientRefinement>" << bvalue[cgradient_refinement] << "</ConjugateGradientRefinement>" << std::endl;
      if (cgradient_refinement) {
        s << "<MaximumConjugateGradientSteps>" << max_CG_steps << "</MaximumConjugateGradientSteps>" << std::endl;
        s << "<MaximumLineSolverSteps>" << max_LS_steps << "</MaximumLineSolverSteps>" << std::endl;
        s << "<EdgeFlexibilityThreshold>" << edge_flexibility_threshold << "</EdgeFlexibilityThreshold>" << std::endl;
      }
    }
    else if (solver == Geometry_Solver::simplex) {
      s << "<ReflectionCoefficient>" << simplex_alpha << "</ReflectionCoefficient>" << std::endl;
      s << "<ExpansionCoefficient>" << simplex_gamma << "</ExpansionCoefficient>" << std::endl;
      s << "<ContractionCoefficient>" << simplex_rho << "</ContractionCoefficient>" << std::endl;
      s << "<ShrinkageCoefficient>" << simplex_sigma << "</ShrinkageCoefficient>" << std::endl;
    }
    s << "</Geometry>" << std::endl;
    s << "</RunTimeParameters>" << std::endl;
    s << "</LogFile>" << std::endl;
    s.close();
    fcall = false;
  }

  nn = skeleton->cardinality(0);
  avg_ven = 0.0;
  max_ven = 0.0;
  min_ven = 10000.0;
  for(i=0; i<Nv; ++i) {
    if (!skeleton->active_event(i)) continue;
    w = skeleton->events[i].get_energy();
    if (w > max_ven) max_ven = w;
    if (w < min_ven) min_ven = w;
    avg_ven += w;
  }
  avg_ven /= double(nn);

  nsink = 0;
  nsource = 0;
  nch = 0;
  nat = 0;
  for(i=0; i<Nv; ++i) {
    // Is this event atemporal (no timelike edges) or chiral
    // (number of FUTURE edges != number of PAST edges)?
    if (!skeleton->active_event(i)) continue;
    nf = 0;
    np = 0;
    atemporal = true;
    skeleton->events[i].get_neighbours(S);
    for(it=S.begin(); it!=S.end(); ++it) {
      in1 = *it;
      edge.clear();
      edge.insert(i); edge.insert(in1);
      qt = skeleton->index_table[1].find(edge);
      if (!skeleton->simplices[1][qt->second].timelike()) continue;
      atemporal = false;
      if (geometry->get_temporal_order(i,in1) == SYNARMOSMA::Relation::before) {
        nf++;
      }
      else {
        np++;
      }
    }
    if (atemporal) {
      nat++;
    }
    else {
      if (nf != np) nch++;
      if (nf == 0) nsink++;
      if (np == 0) nsource++;
    }
  }
  nspace = 0;
  ntime = 0;
  nnull = 0;
  for(i=0; i<Ne; ++i) {
    if (!skeleton->active_simplex(1,i)) continue;
    if (skeleton->simplices[1][i].timelike()) {
      ntime++;
    }
    else if (skeleton->simplices[1][i].spacelike()) {
      nspace++;
    }
    else {
      nnull++;
    }
  }

  logfile.load_file(log_file.c_str());
  rstep = logfile.child("LogFile").append_child("RelaxationStep");

  atom = rstep.append_child("Iteration");
  nvalue = std::to_string(iterations);
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  atom = rstep.append_child("StructuralDeficiency");
  nvalue = std::to_string(error);
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  atom = rstep.append_child("Pseudomanifold");
  nvalue = (skeleton->pseudomanifold) ? "True" : "False";
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  if (skeleton->pseudomanifold) {
    atom = rstep.append_child("Boundary");
    nvalue = (skeleton->boundary) ? "True" : "False";
    atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

    atom = rstep.append_child("Orientable");
    nvalue = (skeleton->orientable) ? "True" : "False";
    atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());
  }
  else {
    atom = rstep.append_child("Boundary");
    atom.append_child(pugi::node_pcdata).set_value("NULL");

    atom = rstep.append_child("Orientable");
    atom.append_child(pugi::node_pcdata).set_value("NULL");
  }

  atom = rstep.append_child("Homology");
  skeleton->get_homology(nvalue);
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  if (high_memory) {
    atom = rstep.append_child("Homotopy");
    skeleton->get_homotopy(nvalue);
    atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());
  }

  ne = skeleton->cardinality(1);
  atom = rstep.append_child("CircuitRank");
  nvalue = std::to_string(ne-nn+1);
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  atom = rstep.append_child("ActiveEvents");
  nvalue = std::to_string(nn);
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  atom = rstep.append_child("SourceEvents");
  nvalue = std::to_string(nsource);
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  atom = rstep.append_child("SinkEvents");
  nvalue = std::to_string(nsink);
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  atom = rstep.append_child("ActiveEdges");
  nvalue = std::to_string(ne);
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  atom = rstep.append_child("SpacelikeEdges");
  nvalue = std::to_string(nspace);
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  atom = rstep.append_child("TimelikeEdges");
  nvalue = std::to_string(ntime);
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  atom = rstep.append_child("NullEdges");
  nvalue = std::to_string(nnull);
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  atom = rstep.append_child("CyclicEdges");
  nvalue = std::to_string(skeleton->cyclicity());
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  arclength_statistics(vdata);
  atom = rstep.append_child("MinimumArcLength");
  nvalue = std::to_string(vdata[1]);
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  atom = rstep.append_child("MeanArcLength");
  nvalue = std::to_string(vdata[2]);
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  atom = rstep.append_child("MaximumArcLength");
  nvalue = std::to_string(vdata[0]);
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  skeleton->vertex_degree_statistics(vdata);
  atom = rstep.append_child("MinimumEventDegree");
  nvalue = std::to_string(vdata[1]);
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  atom = rstep.append_child("MeanEventDegree");
  nvalue = std::to_string(vdata[2]);
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  atom = rstep.append_child("MaximumEventDegree");
  nvalue = std::to_string(vdata[0]);
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  atom = rstep.append_child("TotalEnergy");
  nvalue = std::to_string(skeleton->total_energy());
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  atom = rstep.append_child("MinimumEventEnergy");
  nvalue = std::to_string(min_ven);
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  atom = rstep.append_child("MeanEventEnergy");
  nvalue = std::to_string(avg_ven);
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  atom = rstep.append_child("MaximumEventEnergy");
  nvalue = std::to_string(max_ven);
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  atom = rstep.append_child("AtemporalEvents");
  nvalue = std::to_string(nat);
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  atom = rstep.append_child("ChiralEvents");
  nvalue = std::to_string(nch);
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  atom = rstep.append_child("EulerCharacteristic");
  nvalue = std::to_string(skeleton->euler_characteristic());
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  atom = rstep.append_child("HyphanticSequence");
  if (hyphantic_ops == "") {
    atom.append_child(pugi::node_pcdata).set_value("NULL");
  }
  else {
    atom.append_child(pugi::node_pcdata).set_value(hyphantic_ops.c_str());
  }

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i)
#endif
  for(i=0; i<Nv; ++i) {
    vdimension[i] = (!skeleton->active_event(i)) ? -1 : skeleton->vertex_dimension(i);    
  }

  sum = 0;
  mx = 0;
  mn = Complex::ND;
  for(i=0; i<Nv; ++i) {
    v = vdimension[i];
    if (v == -1) continue;
    sum += v;
    if (v > mx) mx = v;
    if (v < mn) mn = v;
  }
  vd_avg = double(sum)/double(skeleton->cardinality(0));
  vd_min = mn;
  vd_max = mx;

  atom = rstep.append_child("MinimumEventDimension");
  nvalue = std::to_string(vd_min);
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  atom = rstep.append_child("MeanEventDimension");
  nvalue = std::to_string(vd_avg);
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  atom = rstep.append_child("MaximumEventDimension");
  nvalue = std::to_string(vd_max);
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  logfile.save_file(log_file.c_str());
}

void Spacetime::read_state(const std::string& filename)
{
  int i,n;
  char c;
  std::string cmodel,fmodel;

  std::ifstream s(filename,std::ios::in | std::ios::binary);
  if (!s.is_open()) throw std::invalid_argument("The file " + filename + " cannot be opened!");

  clear();

  s.seekg(0);

  s.read((char*)(&n),sizeof(int));
  if (n != 1) {
    s.close();
    throw std::runtime_error("Reading in a state file for non-monophyllon version of the Diaplexis library!"); 
  }

  s.read((char*)(&n),sizeof(int));
  if (n != Complex::ND) {
    s.close();
    throw std::runtime_error("The compiled binary's maximum simplicial dimension " + std::to_string(Complex::ND) + " does not match that (" + std::to_string(n) + ") of the data file!");
  }

  s.read((char*)(&n),sizeof(int));
  if (n != (signed) SYNARMOSMA::Proposition::get_clause_size()) {
    s.close();
    throw std::runtime_error("The compiled binary's propositional clause size " + std::to_string(SYNARMOSMA::Proposition::get_clause_size()) + " does not match that (" + std::to_string(n) + ") of the data file!");
  }

  s.read((char*)(&initial_size),sizeof(int));
  s.read((char*)(&initial_dim),sizeof(int));
  s.read((char*)(&max_iter),sizeof(int));
  s.read((char*)(&checkpoint_frequency),sizeof(int));
  s.read((char*)(&n),sizeof(int));
  skeleton->topological_radius = n;
  // Skip changing the initial_state as it should remain DISKFILE...
  s.read((char*)(&original_state),sizeof(Initial_Topology));
  s.read((char*)(&high_memory),sizeof(bool));
  s.read((char*)(&superposable),sizeof(bool));
  s.read((char*)(&compressible),sizeof(bool));
  s.read((char*)(&perturb_topology),sizeof(bool));
  s.read((char*)(&perturb_geometry),sizeof(bool));
  s.read((char*)(&perturb_energy),sizeof(bool));
  s.read((char*)(&edge_probability),sizeof(double));
  s.read((char*)(&abnormality_threshold),sizeof(double));
  s.read((char*)(&superposition_threshold),sizeof(double));
  s.read((char*)(&parity_mutation),sizeof(double));
  s.read((char*)(&weaving),sizeof(Hyphansis));
  if (weaving == Hyphansis::musical) {
    int j,q;
    std::vector<int> bar;

    if (score_allocated) delete[] hyphantic_notes;
    hyphantic_notes = new std::vector<int>[max_iter+1];
    // Get the score file...
    s.read((char*)(&n),sizeof(int));
    hyphansis_score.resize(n);
    s.read((char*)(&hyphansis_score[0]),n);
    // And now the actual notes...
    for(i=0; i<=max_iter; ++i) {
      s.read((char*)(&n),sizeof(int));
      for(j=0; j<n; ++j) {
        s.read((char*)(&q),sizeof(int));
        bar.push_back(q);
      }
      hyphantic_notes[i] = bar;
      bar.clear();
    }
  }

  s.read((char*)(&solver),sizeof(Geometry_Solver));
  s.read((char*)(&engine),sizeof(Integrator));
  s.read((char*)(&convergence_threshold),sizeof(double));
  s.read((char*)(&coupling_constant),sizeof(double));
  s.read((char*)(&solver_its),sizeof(int));
  s.read((char*)(&ngenerations),sizeof(int));
  s.read((char*)(&pool_size),sizeof(int));
  s.read((char*)(&njousts),sizeof(int));
  s.read((char*)(&annealing_steps),sizeof(int));
  s.read((char*)(&thermal_sweep),sizeof(int));
  s.read((char*)(&thermal_variance),sizeof(double));
  s.read((char*)(&thermalization),sizeof(double));
  s.read((char*)(&max_CG_steps),sizeof(int));
  s.read((char*)(&max_LS_steps),sizeof(int));
  s.read((char*)(&spring_constant),sizeof(double));
  s.read((char*)(&repulsion_constant),sizeof(double));
  s.read((char*)(&damping_constant),sizeof(double));
  s.read((char*)(&step_size),sizeof(double));
  s.read((char*)(&edge_flexibility_threshold),sizeof(double));
  s.read((char*)(&cgradient_refinement),sizeof(bool));
  s.read((char*)(&simplex_alpha),sizeof(double));
  s.read((char*)(&simplex_gamma),sizeof(double));
  s.read((char*)(&simplex_rho),sizeof(double));
  s.read((char*)(&simplex_sigma),sizeof(double));

  s.read((char*)(&system_size),sizeof(int));
  s.read((char*)(&iterations),sizeof(int));
  s.read((char*)(&global_deficiency),sizeof(double));
  s.read((char*)(&error),sizeof(double));
  s.read((char*)(&converged),sizeof(bool));

  s.read((char*)(&n),sizeof(int));
  for(i=0; i<n; ++i) {
    s.read((char*)(&c),sizeof(char));
    hyphantic_ops += c;
  }

  geometry->deserialize(s);
  skeleton->deserialize(s);

  s.close();
}

void Spacetime::write_state(const std::string& filename) const
{
  int i,n = SYNARMOSMA::Proposition::get_clause_size();
  char c;
  std::string datafile;

  // This is a monophyllon file, so the first value is 1...
  const int ftype = 1;

  if (filename == "") {
    datafile = state_file + "_" + std::to_string(iterations) + ".dat";
  }
  else {
    datafile = filename;
  }
  std::ofstream s(datafile,std::ios::out | std::ios::trunc | std::ios::binary);

  // First the global parameters...
  s.write((char*)(&ftype),sizeof(int));
  s.write((char*)(&Complex::ND),sizeof(int));
  s.write((char*)(&n),sizeof(int));

  // Now the global runtime parameters specific to the
  // single Spacetime instance...
  s.write((char*)(&initial_size),sizeof(int));
  s.write((char*)(&initial_dim),sizeof(int));
  s.write((char*)(&max_iter),sizeof(int));
  s.write((char*)(&checkpoint_frequency),sizeof(int));
  n = skeleton->topological_radius;
  s.write((char*)(&n),sizeof(int));
  s.write((char*)(&initial_state),sizeof(Initial_Topology));
  s.write((char*)(&high_memory),sizeof(bool));
  s.write((char*)(&superposable),sizeof(bool));
  s.write((char*)(&compressible),sizeof(bool));
  s.write((char*)(&perturb_topology),sizeof(bool));
  s.write((char*)(&perturb_geometry),sizeof(bool));
  s.write((char*)(&perturb_energy),sizeof(bool));
  s.write((char*)(&edge_probability),sizeof(double));
  s.write((char*)(&abnormality_threshold),sizeof(double));
  s.write((char*)(&superposition_threshold),sizeof(double));
  s.write((char*)(&parity_mutation),sizeof(double));
  s.write((char*)(&weaving),sizeof(Hyphansis));
  if (weaving == Hyphansis::musical) {
    int j,m = (signed) hyphansis_score.size();
    s.write((char*)(&m),sizeof(int));
    s.write((char*)(&hyphansis_score[0]),m);
    for(i=0; i<=max_iter; ++i) {
      n = (signed) hyphantic_notes[i].size();
      s.write((char*)(&n),sizeof(int));
      for(j=0; j<n; ++j) {
        m = hyphantic_notes[i][j];
        s.write((char*)(&m),sizeof(int));
      }
    }
  }

  // Geometric runtime constants...
  s.write((char*)(&solver),sizeof(Geometry_Solver));
  s.write((char*)(&engine),sizeof(Integrator));
  s.write((char*)(&convergence_threshold),sizeof(double));
  s.write((char*)(&coupling_constant),sizeof(double));
  s.write((char*)(&solver_its),sizeof(int));
  s.write((char*)(&ngenerations),sizeof(int));
  s.write((char*)(&pool_size),sizeof(int));
  s.write((char*)(&njousts),sizeof(int));
  s.write((char*)(&annealing_steps),sizeof(int));
  s.write((char*)(&thermal_sweep),sizeof(int));
  s.write((char*)(&thermal_variance),sizeof(double));
  s.write((char*)(&thermalization),sizeof(double));
  s.write((char*)(&max_CG_steps),sizeof(int));
  s.write((char*)(&max_LS_steps),sizeof(int));
  s.write((char*)(&spring_constant),sizeof(double));
  s.write((char*)(&repulsion_constant),sizeof(double));
  s.write((char*)(&damping_constant),sizeof(double));
  s.write((char*)(&step_size),sizeof(double));
  s.write((char*)(&edge_flexibility_threshold),sizeof(double));
  s.write((char*)(&cgradient_refinement),sizeof(bool));
  s.write((char*)(&simplex_alpha),sizeof(double));
  s.write((char*)(&simplex_gamma),sizeof(double));
  s.write((char*)(&simplex_rho),sizeof(double));
  s.write((char*)(&simplex_sigma),sizeof(double));

  // Variables indicating the current state of the
  // Spacetime instance...
  s.write((char*)(&system_size),sizeof(int));
  s.write((char*)(&iterations),sizeof(int));
  s.write((char*)(&global_deficiency),sizeof(double));
  s.write((char*)(&error),sizeof(double));
  s.write((char*)(&converged),sizeof(bool));

  // Now write out the hyphantic operations that have been performed...
  n = hyphantic_ops.length();
  s.write((char*)(&n),sizeof(int));
  for(i=0; i<n; ++i) {
    c = hyphantic_ops[i];
    s.write((char*)(&c),sizeof(char));
  }

  // The principal body of the file...
  geometry->serialize(s);
  skeleton->serialize(s);

  // Close the file and return...
  s.close();
}
