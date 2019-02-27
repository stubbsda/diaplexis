#include "spacetime.h"

using namespace DIAPLEXIS;

void Spacetime::read_parameters(const std::string& filename)
{
  pugi::xml_document pfile;
  pugi::xml_node global,gsolver;
  std::string name,value;
  unsigned int rs;
  int q,n_is = 0,n_so = 0,D = 0;
  bool euclidean = false,relational = false,uniform = false;

  // Open the file
  if (!(pfile.load_file(filename.c_str()))) throw std::invalid_argument("Unable to parse parameter file!");

  global = pfile.child("Parameters").child("Global");
  for(pugi::xml_node params = global.first_child(); params; params = params.next_sibling()) {
    name = params.name();
    value = params.first_child().value();

    if (name == "InitialEvents") {
      initial_size = boost::lexical_cast<int>(value);
    }
    else if (name == "InitialState") {
      boost::to_upper(value);
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
      n_is++;
    }
    else if (name == "InputFile") {
      input_file = value;
    }
    else if (name == "EuclideanGeometry") {
      boost::to_upper(value);
      euclidean = (value == "YES") ? true : false;
    }
    else if (name == "RelationalGeometry") {
      boost::to_upper(value);
      relational = (value == "YES") ? true : false;
    }
    else if (name == "DimensionalUniformity") {
      boost::to_upper(value);
      uniform = (value == "YES") ? true : false;
    }
    else if (name == "RandomSeed") {
      rs = boost::lexical_cast<unsigned int>(value);
      if (rs == 0) rs = (unsigned) std::time(nullptr);
      skeleton->RND->set_seed(rs);
    }
    else if (name == "BackgroundDimension") {
      D = boost::lexical_cast<int>(value);
    }
    else if (name == "EdgeProbability") {
      edge_probability = boost::lexical_cast<double>(value);
    }
    else if (name == "MaximumIterations") {
      max_iter = boost::lexical_cast<int>(value);
    }
    else if (name == "InitialDimension") {
      initial_dim = boost::lexical_cast<int>(value);
    }
    else if (name == "PerturbTopology") {
      boost::to_upper(value);
      perturb_topology = (value == "YES") ? true : false;
    }
    else if (name == "PerturbGeometry") {
      boost::to_upper(value);
      perturb_geometry = (value == "YES") ? true : false;
    }
    else if (name == "PerturbEnergy") {
      boost::to_upper(value);
      perturb_energy = (value == "YES") ? true : false;
    }
    else if (name == "MemoryFootprint") {
      boost::to_upper(value);
      high_memory = (value == "HIGH") ? true : false;
    }
    else if (name == "InstrumentConvergence") {
      boost::to_upper(value);
      instrument_convergence = (value == "YES") ? true : false;
    }
    else if (name == "Hyphansis") {
      boost::to_upper(value);
      if (value == "MUSICAL") {
        weaving = Hyphansis::musical;
      }
      else if (value == "DYNAMIC") {
        weaving = Hyphansis::dynamic;
      }
    }
    else if (name == "HyphansisScore") {
      hyphansis_score = value;
    }
    else if (name == "ParityMutation") {
      parity_mutation = boost::lexical_cast<double>(value);
    }
    else if (name == "HomologyMethod") {
      boost::to_upper(value);
      if (value == "GAP") {
        skeleton->set_homology_method(SYNARMOSMA::Homology::Method::gap); 
      }
      else if (value == "NATIVE") {
        skeleton->set_homology_method(SYNARMOSMA::Homology::Method::native);
      }
    }
    else if (name == "HomologyBase") {
      boost::to_upper(value);
      if (value == "INT") {
        skeleton->set_homology_field(SYNARMOSMA::Homology::Field::int32);
      }
      else if (value == "NTL::ZZ") {
        skeleton->set_homology_field(SYNARMOSMA::Homology::Field::multiprecision);
      }
      else if (value == "GF2") {
        skeleton->set_homology_field(SYNARMOSMA::Homology::Field::mod2);
      }
    }
    else if (name == "Compressible") {
      boost::to_upper(value);
      compressible = (value == "YES") ? true : false;
    }
    else if (name == "Superposable") {
      boost::to_upper(value);
      superposable = (value == "YES") ? true : false;
    }
    else if (name == "CheckpointFrequency") {
      checkpoint_frequency = boost::lexical_cast<int>(value);
    }
  }

  gsolver = pfile.child("Parameters").child("GeometrySolver");
  for(pugi::xml_node params = gsolver.first_child(); params; params = params.next_sibling()) {
    name = params.name();
    value = params.first_child().value();

    if (name == "SolverType") {
      boost::to_upper(value);
      if (value == "MECHANICAL") {
        solver = Geometry_Solver::mechanical;
      }
      else if (value == "EVOLUTIONARY") {
        solver = Geometry_Solver::evolutionary;
      }
      else if (value == "ANNEALING") {
        solver = Geometry_Solver::annealing;
      }
      else if (value == "MINIMAL") {
        solver = Geometry_Solver::minimal;
      }
      else if (value == "SIMPLEX") {
        solver = Geometry_Solver::simplex;
      }
      n_so++;
    }
    else if (name == "SolverIterations") {
      solver_its = boost::lexical_cast<int>(value);
    }
    else if (name == "MaximumGenerations") {
      ngenerations = boost::lexical_cast<int>(value);
    }
    else if (name == "MaximumJousts") {
      njousts = boost::lexical_cast<int>(value);
    }
    else if (name == "PoolSize") {
      pool_size = boost::lexical_cast<int>(value);
    }
    else if (name == "ThermalSweeps") {
      thermal_sweep = boost::lexical_cast<int>(value);
    }
    else if (name == "AnnealingSteps") {
      annealing_steps = boost::lexical_cast<int>(value);
    }
    else if (name == "ThermalVariance") {
      thermal_variance = boost::lexical_cast<double>(value);
    }
    else if (name == "ThermalizationCriterion") {
      thermalization = boost::lexical_cast<double>(value);
    }
    else if (name == "IntegrationEngine") {
      boost::to_upper(value);
      if (value == "EULER") {
        engine = Integrator::euler;
      }
      else if (value == "RK4") {
        engine = Integrator::rk4;
      }
    }
    else if (name == "StepSize") {
      step_size = boost::lexical_cast<double>(value);
    }
    else if (name == "MaximumIntegrationSteps") {
      max_int_steps = boost::lexical_cast<int>(value);
    }
    else if (name == "DampingConstant") {
      damping_constant = boost::lexical_cast<double>(value);
    }
    else if (name == "SpringConstant") {
      spring_constant = boost::lexical_cast<double>(value);
    }
    else if (name == "RepulsionConstant") {
      repulsion_constant = boost::lexical_cast<double>(value);
    } 
    else if (name == "ConjugateGradientRefinement") {
      boost::to_upper(value);
      cgradient_refinement = (value == "YES") ? true : false;
    }
    else if (name == "MaximumConjugateGradientSteps") {
      max_CG_steps = boost::lexical_cast<int>(value);
    }
    else if (name == "MaximumLineSolverSteps") {
      max_LS_steps = boost::lexical_cast<int>(value);
    }
    else if (name == "EdgeFlexibilityThreshold") {
      edge_flexibility_threshold = boost::lexical_cast<double>(value);
    }
    else if (name == "ReflectionCoefficient") {
      simplex_alpha = boost::lexical_cast<double>(value);
    }
    else if (name == "ExpansionCoefficient") {
      simplex_gamma = boost::lexical_cast<double>(value);
    }
    else if (name == "ContractionCoefficient") {
      simplex_rho = boost::lexical_cast<double>(value);
    }
    else if (name == "ShrinkageCoefficient") {
      simplex_sigma = boost::lexical_cast<double>(value);
    }
    else if (name == "GeometryTolerance") {
      geometry_tolerance = boost::lexical_cast<double>(value);
    }
  }

  // Now a series of tests to make sure that the parameters
  // aren't entirely crazy...
  assert(geometry_tolerance > std::numeric_limits<double>::epsilon());
  assert(max_iter >= 0);
  assert(initial_size > 0);
  // One and only one initial state should be chosen...
  assert(n_is == 1);
  // And similarly for the solver type...
  assert(n_so == 1);
  // Finally for the background dimension...
  assert(D > 0);

  if (weaving == Hyphansis::dynamic) {
    assert(parity_mutation >= 0.0);
    assert(parity_mutation <= 1.0);
  }
  else {
    // Make sure the score file exists
    assert(boost::filesystem::exists(hyphansis_score));
  }

  geometry->initialize(euclidean,relational,uniform,high_memory,D);

  if (initial_state == Initial_Topology::random) {
    assert(edge_probability > std::numeric_limits<double>::epsilon() && (edge_probability - 1.0) < -std::numeric_limits<double>::epsilon());
    assert(initial_size > 1);
  }
  else if (initial_state == Initial_Topology::cartesian) {
    // If initial_size != n^dimension, n \in Z+, we have a problem
    std::vector<std::pair<long,int> > factors;
    SYNARMOSMA::factorize(initial_size,factors);
    for(q=0; q<(signed) factors.size(); ++q) {
      if (factors[q].second % D != 0) throw std::invalid_argument("Grid size inconsistent with the background dimension!");
    }
  }
  else if (initial_state == Initial_Topology::monoplex) {
    assert(initial_dim <= Complex::ND);
  }
  else if (initial_state == Initial_Topology::singleton) {
    assert(initial_size == 1);
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
  std::set<int> S;
  std::set<int>::const_iterator it;
  std::string group_string,bvalue[] = {"False","True"};
  const int Nv = (signed) skeleton->events.size();
  const int Ne = (signed) skeleton->simplices[1].size();
  int vdimension[Nv];

  if (fcall) {
    char hostname[80];
    std::string sname[] = {"RANDOM","MONOPLEX","CARTESIAN","SINGLETON","DISKFILE"};
    std::string hname[] = {"Minimal","Conjugate Gradient","Evolutionary","Simulated Annealing","Simplex"};

    gethostname(hostname,80);

    std::ofstream s(log_file,std::ios::out | std::ios::trunc);
    s << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << std::endl;
    s << "<LogFile>" << std::endl;
    s << "<Library>Diaplexis</Library>" << std::endl;
    s << "<Hostname>" << hostname << "</Hostname>" << std::endl;
    s << "<StartDate>" << start_time.date() << "</StartDate>" << std::endl;
    s << "<StartTime>" << start_time.time_of_day() << "</StartTime>" << std::endl;
    s << "<CompileTimeParameters>" << std::endl;
    s << "<MaximumDimension>" << Complex::ND << "</MaximumDimension>" << std::endl;
    s << "<AtomicPropositions>" << SYNARMOSMA::Proposition::get_clause_size() << "</AtomicPropositions>" << std::endl;
    s << "<TopologicalRadius>" << Complex::topological_radius << "</TopologicalRadius>" << std::endl;
    s << "<ConvergenceThreshold>" << Spacetime::convergence_threshold << "</ConvergenceThreshold>" << std::endl;
    s << "<InitialAnnealingTemperature>" << Spacetime::T_zero << "</InitialAnnealingTemperature>" << std::endl;
    s << "<ThermalDecayRate>" << Spacetime::kappa << "</ThermalDecayRate>" << std::endl;
    s << "<EnergyCouplingConstant>" << Spacetime::Lambda << "</EnergyCouplingConstant>" << std::endl;
    s << "</CompileTimeParameters>" << std::endl;
    s << "<RunTimeParameters>" << std::endl;
    s << "<Global>" << std::endl;
    if (initial_state != Initial_Topology::singleton) s << "<InitialEvents>" << initial_size << "</InitialEvents>" << std::endl;
    s << "<MaximumIterations>" << max_iter << "</MaximumIterations>" << std::endl;
    s << "<EuclideanGeometry>" << bvalue[geometry->get_euclidean()] << "</EuclideanGeometry>" << std::endl;
    s << "<RelationalGeometry>" << bvalue[geometry->get_relational()] << "</RelationalGeometry>" << std::endl;
    s << "<DimensionalUniformity>" << bvalue[geometry->get_uniform()] << "</DimensionalUniformity>" << std::endl;
    s << "<BackgroundDimension>" << geometry->dimension() << "</BackgroundDimension>" << std::endl;
    s << "<RandomSeed>" << skeleton->RND->get_seed() << "</RandomSeed>" << std::endl;
    s << "<Superposable>" << bvalue[superposable] << "</Superposable>" << std::endl;
    s << "<Compressible>" << bvalue[compressible] << "</Compressible>" << std::endl;
    if (initial_state == Initial_Topology::diskfile) {
      s << "<InitialState>" << sname[initial_state] << " (" << sname[original_state] << ")</InitialState>" << std::endl;
      s << "<InputFile>" << input_file << "</InputFile>" << std::endl;
    }
    else {
      s << "<InitialState>" << sname[initial_state] << "</InitialState>" << std::endl;
    }
    if (initial_state == Initial_Topology::random) s << "<EdgeProbability>" << edge_probability << "</EdgeProbability>" << std::endl;
    s << "<CheckpointFrequency>" << checkpoint_frequency << "</CheckpointFrequency>" << std::endl;
    if (initial_state == Initial_Topology::cartesian) {
      s << "<PerturbTopology>" << bvalue[perturb_topology] << "</PerturbTopology>" << std::endl;
      s << "<PerturbGeometry>" << bvalue[perturb_geometry] << "</PerturbGeometry>" << std::endl;
      s << "<PerturbEnergy>" << bvalue[perturb_energy] << "</PerturbEnergy>" << std::endl;
    }
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
    }
    else {
      s << "<Hyphansis>MUSICAL</Hyphansis>" << std::endl;
      s << "<HyphansisScore>" << hyphansis_score << "</HyphansisScore>" << std::endl;
    }
    s << "</Global>" << std::endl;
    s << "<Geometry>" << std::endl;
    s << "<GeometryTolerance>" << geometry_tolerance << "</GeometryTolerance>" << std::endl;
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
    // Is this vertex atemporal (no timelike edges) or chiral
    // (number of FUTURE edges != number of PAST edges)?
    if (!skeleton->active_event(i)) continue;
    nf = 0;
    np = 0;
    atemporal = true;
    for(it=skeleton->events[i].neighbours.begin(); it!=skeleton->events[i].neighbours.end(); ++it) {
      in1 = *it;
      S.clear();
      S.insert(i); S.insert(in1);
      qt = skeleton->index_table[1].find(S);
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
  nvalue = boost::lexical_cast<std::string>(iterations);
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  atom = rstep.append_child("StructuralDeficiency");
  nvalue = boost::lexical_cast<std::string>(error);
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  if (instrument_convergence) {
    atom = rstep.append_child("TopologyDelta");
    nvalue = boost::lexical_cast<std::string>(topology_delta);
    atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

    atom = rstep.append_child("GeometryDelta");
    nvalue = boost::lexical_cast<std::string>(geometry_delta);
    atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

    atom = rstep.append_child("EnergyDelta");
    nvalue = boost::lexical_cast<std::string>(energy_delta);
    atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());
  }

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
  nvalue = boost::lexical_cast<std::string>(ne-nn+1);
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  atom = rstep.append_child("ActiveEvents");
  nvalue = boost::lexical_cast<std::string>(nn);
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  atom = rstep.append_child("SourceEvents");
  nvalue = boost::lexical_cast<std::string>(nsource);
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  atom = rstep.append_child("SinkEvents");
  nvalue = boost::lexical_cast<std::string>(nsink);
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  atom = rstep.append_child("ActiveEdges");
  nvalue = boost::lexical_cast<std::string>(ne);
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  atom = rstep.append_child("SpacelikeEdges");
  nvalue = boost::lexical_cast<std::string>(nspace);
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  atom = rstep.append_child("TimelikeEdges");
  nvalue = boost::lexical_cast<std::string>(ntime);
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  atom = rstep.append_child("NullEdges");
  nvalue = boost::lexical_cast<std::string>(nnull);
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  atom = rstep.append_child("CyclicEdges");
  nvalue = boost::lexical_cast<std::string>(skeleton->cyclicity());
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  arclength_statistics(vdata);
  atom = rstep.append_child("MinimumArcLength");
  nvalue = boost::lexical_cast<std::string>(vdata[1]);
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  atom = rstep.append_child("MeanArcLength");
  nvalue = boost::lexical_cast<std::string>(vdata[2]);
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  atom = rstep.append_child("MaximumArcLength");
  nvalue = boost::lexical_cast<std::string>(vdata[0]);
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  skeleton->vertex_degree_statistics(vdata);
  atom = rstep.append_child("MinimumEventDegree");
  nvalue = boost::lexical_cast<std::string>(vdata[1]);
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  atom = rstep.append_child("MeanEventDegree");
  nvalue = boost::lexical_cast<std::string>(vdata[2]);
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  atom = rstep.append_child("MaximumEventDegree");
  nvalue = boost::lexical_cast<std::string>(vdata[0]);
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  atom = rstep.append_child("TotalEnergy");
  nvalue = boost::lexical_cast<std::string>(skeleton->total_energy());
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  atom = rstep.append_child("MinimumEventEnergy");
  nvalue = boost::lexical_cast<std::string>(min_ven);
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  atom = rstep.append_child("MeanEventEnergy");
  nvalue = boost::lexical_cast<std::string>(avg_ven);
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  atom = rstep.append_child("MaximumEventEnergy");
  nvalue = boost::lexical_cast<std::string>(max_ven);
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  atom = rstep.append_child("AtemporalEvents");
  nvalue = boost::lexical_cast<std::string>(nat);
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  atom = rstep.append_child("ChiralEvents");
  nvalue = boost::lexical_cast<std::string>(nch);
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  atom = rstep.append_child("EulerCharacteristic");
  nvalue = boost::lexical_cast<std::string>(skeleton->euler_characteristic());
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
  nvalue = boost::lexical_cast<std::string>(vd_min);
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  atom = rstep.append_child("MeanEventDimension");
  nvalue = boost::lexical_cast<std::string>(vd_avg);
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  atom = rstep.append_child("MaximumEventDimension");
  nvalue = boost::lexical_cast<std::string>(vd_max);
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  logfile.save_file(log_file.c_str());
}

void Spacetime::write(Spacetime& state) const
{
  state.system_size = system_size;
  state.error = error;
  state.global_deficiency = global_deficiency;
  state.skeleton->events = skeleton->events;
  for(int i=0; i<=Complex::ND; ++i) {
    state.skeleton->simplices[i] = skeleton->simplices[i];
    state.skeleton->index_table[i] = skeleton->index_table[i];
  }
}

void Spacetime::read(const Spacetime& source)
{
  system_size = source.system_size;
  error = source.error;
  global_deficiency = source.global_deficiency;
  skeleton->events = source.skeleton->events;
  for(int i=0; i<=Complex::ND; ++i) {
    skeleton->simplices[i] = source.skeleton->simplices[i];
    skeleton->index_table[i] = source.skeleton->index_table[i];
  }
}

void Spacetime::read_state(const std::string& filename)
{
  int i,n;
  char c;
  double x;
  std::string cmodel,fmodel;

  std::ifstream s(filename,std::ios::in | std::ios::binary);
  if (!s.is_open()) {
    // File doesn't exist, print an error message and die
    std::cerr << "The file " << filename << " cannot be found." << std::endl;
    std::cerr << "Exiting..." << std::endl;
    std::exit(1);
  }

  clear();

  s.seekg(0);

  s.read((char*)(&n),sizeof(int));
  if (n != 1) {
    s.close();
    std::cerr << "Reading in a state file for non-monophyllon version of the Diaplexis library." << std::endl;
    std::cerr << "Exiting..." << std::endl;
    std::exit(1);
  }

  s.read((char*)(&n),sizeof(int));
  if (n != Complex::ND) {
    s.close();
    std::cerr << "The compiled binary's maximum simplicial dimension " << Complex::ND << " does not match that (" << n << ") of the data file." << std::endl;
    std::cerr << "Exiting..." << std::endl;
    std::exit(1);
  }

  s.read((char*)(&n),sizeof(int));
  if (n != (signed) SYNARMOSMA::Proposition::get_clause_size()) {
    s.close();
    std::cerr << "The compiled binary's atomic clause number " << SYNARMOSMA::Proposition::get_clause_size() << " does not match that (" << n << ") of the data file." << std::endl;
    std::cerr << "Exiting..." << std::endl;
    std::exit(1);
  }

  s.read((char*)(&n),sizeof(int));
  if (n != Complex::topological_radius) {
    s.close();
    std::cerr << "The compiled binary's topological radius " << Complex::topological_radius << " does not match that (" << n << ") of the data file." << std::endl;
    std::cerr << "Exiting..." << std::endl;
    std::exit(1);
  }

  s.read((char*)(&x),sizeof(double));
  if (!SYNARMOSMA::double_equality(x,Spacetime::convergence_threshold)) {
    s.close();
    std::cerr << "The compiled binary's convergence threshold " << Spacetime::convergence_threshold << " does not match that (" << x << ") of the data file." << std::endl;
    std::cerr << "Exiting..." << std::endl;
    std::exit(1);
  }

  s.read((char*)(&x),sizeof(double));
  if (!SYNARMOSMA::double_equality(x,Spacetime::T_zero)) {
    s.close();
    std::cerr << "The compiled binary's T_zero " << Spacetime::T_zero << " does not match that (" << x << ") of the data file." << std::endl;
    std::cerr << "Exiting..." << std::endl;
    std::exit(1);
  }

  s.read((char*)(&x),sizeof(double));
  if (!SYNARMOSMA::double_equality(x,Spacetime::kappa)) {
    s.close();
    std::cerr << "The compiled binary's kappa " << Spacetime::kappa << " does not match that (" << x << ") of the data file." << std::endl;
    std::cerr << "Exiting..." << std::endl;
    std::exit(1);
  }

  s.read((char*)(&x),sizeof(double));
  if (!SYNARMOSMA::double_equality(x,Spacetime::Lambda)) {
    s.close();
    std::cerr << "The compiled binary's Lambda " << Spacetime::Lambda << " does not match that (" << x << ") of the data file." << std::endl;
    std::cerr << "Exiting..." << std::endl;
    std::exit(1);
  }

  s.read((char*)(&initial_size),sizeof(int));
  s.read((char*)(&initial_dim),sizeof(int));
  s.read((char*)(&max_iter),sizeof(int));
  s.read((char*)(&checkpoint_frequency),sizeof(int));
  // Skip changing the initial_state as it should remain DISKFILE...
  s.read((char*)(&original_state),sizeof(Initial_Topology));
  s.read((char*)(&instrument_convergence),sizeof(bool));
  s.read((char*)(&high_memory),sizeof(bool));
  s.read((char*)(&superposable),sizeof(bool));
  s.read((char*)(&compressible),sizeof(bool));
  s.read((char*)(&perturb_topology),sizeof(bool));
  s.read((char*)(&perturb_geometry),sizeof(bool));
  s.read((char*)(&perturb_energy),sizeof(bool));
  s.read((char*)(&edge_probability),sizeof(double));
  s.read((char*)(&parity_mutation),sizeof(double));
  s.read((char*)(&weaving),sizeof(Hyphansis));
  if (weaving == Hyphansis::musical) {
    s.read((char*)(&n),sizeof(int));
    hyphansis_score.resize(n);
    s.read((char*)(&hyphansis_score[0]),n);
  }

  s.read((char*)(&solver),sizeof(Geometry_Solver));
  s.read((char*)(&engine),sizeof(Integrator));
  s.read((char*)(&geometry_tolerance),sizeof(double));
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

void Spacetime::write_state() const
{
  int i,n = SYNARMOSMA::Proposition::get_clause_size();
  char c;

  // This is a monophyllon file, so the first value is 1...
  const int ftype = 1;

  std::stringstream sstream;
  sstream << iterations;
  std::string filename = state_file + "_" + sstream.str() + ".dat";

  std::ofstream s(filename,std::ios::out | std::ios::trunc | std::ios::binary);

  // First the global parameters...
  s.write((char*)(&ftype),sizeof(int));
  s.write((char*)(&Complex::ND),sizeof(int));
  s.write((char*)(&n),sizeof(int));
  s.write((char*)(&Complex::topological_radius),sizeof(int));
  s.write((char*)(&Spacetime::convergence_threshold),sizeof(double));
  s.write((char*)(&Spacetime::T_zero),sizeof(double));
  s.write((char*)(&Spacetime::kappa),sizeof(double));
  s.write((char*)(&Spacetime::Lambda),sizeof(double));

  // Now the global runtime parameters specific to the
  // single Spacetime instance...
  s.write((char*)(&initial_size),sizeof(int));
  s.write((char*)(&initial_dim),sizeof(int));
  s.write((char*)(&max_iter),sizeof(int));
  s.write((char*)(&checkpoint_frequency),sizeof(int));
  s.write((char*)(&initial_state),sizeof(Initial_Topology));
  s.write((char*)(&instrument_convergence),sizeof(bool));
  s.write((char*)(&high_memory),sizeof(bool));
  s.write((char*)(&superposable),sizeof(bool));
  s.write((char*)(&compressible),sizeof(bool));
  s.write((char*)(&perturb_topology),sizeof(bool));
  s.write((char*)(&perturb_geometry),sizeof(bool));
  s.write((char*)(&perturb_energy),sizeof(bool));
  s.write((char*)(&edge_probability),sizeof(double));
  s.write((char*)(&parity_mutation),sizeof(double));
  s.write((char*)(&weaving),sizeof(Hyphansis));
  if (weaving == Hyphansis::musical) {
    int m = (signed) hyphansis_score.size();
    s.write((char*)(&m),sizeof(int));
    s.write((char*)(&hyphansis_score[0]),m);
  }

  // Geometric runtime constants...
  s.write((char*)(&solver),sizeof(Geometry_Solver));
  s.write((char*)(&engine),sizeof(Integrator));
  s.write((char*)(&geometry_tolerance),sizeof(double));
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
