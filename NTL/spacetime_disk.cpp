#include "spacetime.h"

extern SYNARMOSMA::Random RND;

using namespace DIAPLEXIS;

void Spacetime::read_parameters(const char* filename)
{
  pugi::xml_document pfile;
  pugi::xml_node global,gsolver;
  std::string name,value;
  unsigned int rs;
  int q,n_is = 0,n_so = 0,D = 0;
  bool euclidean = false,relational = false,uniform = false;

  // Open the file
  if (!(pfile.load_file(filename))) {
    std::cerr << "The file " << filename << " either does not exist or could not be loaded correctly." << std::endl;
    std::cerr << "Exiting..." << std::endl;
    std::exit(1);
  }

  global = pfile.child("Parameters").child("Global");
  for(pugi::xml_node params = global.first_child(); params; params = params.next_sibling()) {
    name = params.name();
    value = params.first_child().value();

    if (name == "InitialVertices") {
      initial_size = boost::lexical_cast<int>(value);
    }
    else if (name == "InitialState") {
      if (value == "CARTESIAN") {
        initial_state = CARTESIAN;
      }
      else if (value == "RANDOM") {
        initial_state = RANDOM;
      }
      else if (value == "SINGLETON") {
        initial_state = SINGLETON;
      }
      else if (value == "MONOPLEX") {
        initial_state = MONOPLEX;
      }
      else if (value == "DISKFILE") {
        initial_state = DISKFILE;
      }
      n_is++;
    }
    else if (name == "InputFile") {
      input_file = value;
    }
    else if (name == "EuclideanGeometry") {
      if (value == "True") euclidean = true;
    }
    else if (name == "RelationalGeometry") {
      if (value == "True") relational = true;
    }
    else if (name == "DimensionalUniformity") {
      if (value == "True") uniform = true;
    }
    else if (name == "RandomSeed") {
      rs = boost::lexical_cast<unsigned int>(value);
      if (rs == 0) {
        RND.set_seed(static_cast<unsigned int>(std::time(0)));
      }
      else {
        RND.set_seed(rs);
      }
    }
    else if (name == "BackgroundDimension") {
      D = boost::lexical_cast<int>(value);
    }
    else if (name == "EdgeProbability") {
      edge_probability = boost::lexical_cast<double>(value);
    }
    else if (name == "MaxIterations") {
      max_iter = boost::lexical_cast<int>(value);
    }
    else if (name == "InitialSheets") {
      nt_initial = boost::lexical_cast<int>(value);
    }
    else if (name == "InitialDimension") {
      initial_dim = boost::lexical_cast<int>(value);
    }
    else if (name == "PerturbTopology") {
      perturb_topology = (value == "True") ? true : false;
    }
    else if (name == "PerturbGeometry") {
      perturb_geometry = (value == "True") ? true : false;
    }
    else if (name == "PerturbEnergy") {
      perturb_energy = (value == "True") ? true : false;
    }
    else if (name == "MemoryFootprint") {
      high_memory = (value == "High") ? true : false;
    }
    else if (name == "InstrumentConvergence") {
      instrument_convergence = (value == "True") ? true : false;
    }
    else if (name == "Hyphansis") {
      if (value == "MUSICAL") weaving = MUSICAL;
    }
    else if (name == "HyphansisScore") {
      hyphansis_score = value;
    }
    else if (name == "HomologyMethod") {
      if (value == "GAP") {
        H->set_method(SYNARMOSMA::GAP); 
      }
      else if (value == "Native") {
        H->set_method(SYNARMOSMA::NATIVE);
      }
    }
    else if (name == "HomologyBase") {
      if (value == "INT") {
        H->set_field(SYNARMOSMA::INT);
      }
      else if (value == "NTL::ZZ") {
        H->set_field(SYNARMOSMA::ZZ);
      }
      else {
        // This should be the Galois field GF(2)
        H->set_field(SYNARMOSMA::GF2);
      }
    }
    else if (name == "Compressible") {
      compressible = (value == "True") ? true : false;
    }
    else if (name == "Permutable") {
      permutable = (value == "True") ? true : false;
    }
    else if (name == "Superposable") {
      superposable = (value == "True") ? true : false;
    }
    else if (name == "SheetDynamics") {
      foliodynamics = (value == "True") ? true : false;
    }
    else if (name == "Checkpoint") {
      checkpoint = (value == "True") ? true : false;
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
      if (value == "MECHANICAL") {
        solver = MECHANICAL;
      }
      else if (value == "EVOLUTIONARY") {
        solver = EVOLUTIONARY;
      }
      else if (value == "ANNEALING") {
        solver = ANNEALING;
      }
      else if (value == "MINIMAL") {
        solver = MINIMAL;
      }
      else if (value == "SIMPLEX") {
        solver = SIMPLEX;
      }
      n_so++;
    }
    else if (name == "SolverIterations") {
      solver_its = boost::lexical_cast<int>(value);
    }
    else if (name == "MaxGenerations") {
      ngenerations = boost::lexical_cast<int>(value);
    }
    else if (name == "MaxJousts") {
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
      int_engine = value;
    }
    else if (name == "StepSize") {
      step_size = boost::lexical_cast<double>(value);
    }
    else if (name == "MaxIntegrationSteps") {
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
      cgradient_refinement = (value == "True") ? true : false;
    }
    else if (name == "MaxConjugateGradientSteps") {
      max_CG_steps = boost::lexical_cast<int>(value);
    }
    else if (name == "MaxLineSolverSteps") {
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
    else if (name == "GeometryCutoff") {
      geometry_cutoff = boost::lexical_cast<double>(value);
    }
  }

  // Now a series of tests to make sure that the parameters
  // aren't entirely crazy...
  assert(geometry_cutoff > 0.0);
  assert(max_iter >= 0);
  assert(nt_initial > 0);
  assert(initial_size > 0);
  // One and only one initial state should be chosen...
  assert(n_is == 1);
  // And similarly for the solver type...
  assert(n_so == 1);
  // Finally for the background dimension...
  assert(D > 0);

  geometry->initialize(euclidean,relational,uniform,high_memory,D);
  if (instrument_convergence) anterior.geometry.initialize(euclidean,relational,uniform,high_memory,D);

  if (initial_state == RANDOM) {
    assert(edge_probability > 0.0 && edge_probability < 1.0);
    assert(initial_size > 1);
  }
  else if (initial_state == CARTESIAN) {
    // If initial_size != n^dimension, n \in Z+, we have a problem
    std::vector<std::pair<long,int> > factors;
    SYNARMOSMA::factorize(initial_size,factors);
    for(q=0; q<(signed) factors.size(); ++q) {
      if (factors[q].second % D != 0) {
        std::cerr << "The grid size " << initial_size << " is not consistent with the background dimension " << D << "." << std::endl;
        std::cerr << "Exiting..." << std::endl;
        std::exit(1);
      }
    }
  }
  else if (initial_state == MONOPLEX) {
    assert(initial_dim <= Spacetime::ND);
  }
  else if (initial_state == SINGLETON) {
    assert(initial_size == 1);
  }

  if (solver == MINIMAL) {
    assert(solver_its > 0);
  }
  else if (solver == EVOLUTIONARY) {
    assert(pool_size > 0);
    assert(ngenerations > 0);
    assert(njousts > 0);
  }
  else if (solver == ANNEALING) {
    assert(annealing_steps > 0);
    assert(thermalization > 0.0);
    assert(thermal_variance > 0.0);
    assert(thermal_sweep > 0);
  }
  else if (solver == MECHANICAL) {
    if (!relational && !uniform) {
      std::cerr << "If the geometry model is absolute it must be dimensionally uniform in order to use the MECHANICAL geometry solver!" << std::endl;
      std::exit(1);
    }
    assert(int_engine == "EULER" || int_engine == "RK4");
    assert(step_size > 0.0 && step_size < 1.0);
    assert(spring_constant < 0.0);
    assert(repulsion_constant > 0.0);
    assert(damping_constant > 0.0);
    assert(max_int_steps > 0);
    if (cgradient_refinement) {
      assert(max_CG_steps > 0);
      assert(max_LS_steps > 0);
      assert(edge_flexibility_threshold > 0.0);
    }
  }
  else if (solver == SIMPLEX) {
    assert(simplex_alpha > 0.0);
    assert(simplex_gamma > 1.0);
    assert(simplex_gamma > simplex_alpha);
    assert(simplex_rho > 0.0 && simplex_rho < 1.0);
    assert(simplex_sigma > 0.0 && simplex_sigma < 1.0);
  }
  original_state = initial_state;
}

void Spacetime::write_log() const
{
  int i,j,l,v,nn,ne,max_val,min_val,ntime,nspace,ncyclic,nf,np,nat;
  int nsource,nsink,nnull,nch,in1,sum,mx,mn,vx[2];
  double w,wm,avg_val,avg_ven,max_ven,min_ven,vdata[3];
  double max_length,min_length,avg_length;
  bool atemporal;
  static bool fcall = true;
  SYNARMOSMA::RELATION rho;
  std::string nvalue;
  pugi::xml_document logfile;
  pugi::xml_node rstep,sheet,atom;
  SYNARMOSMA::hash_map::const_iterator qt;
  std::set<int>::const_iterator it;
  std::string group_string,bvalue[] = {"False","True"};
  const int Nv = (signed) events.size();
  const int Ne = (signed) simplices[1].size();
  const int Nt = (signed) codex.size();
  double vd_avg[Nt];
  int vd_max[Nt],vd_min[Nt];
  int vdimension[Nv*Nt];

  if (fcall) {
    char hostname[80];
    std::string sname[] = {"RANDOM","MONOPLEX","CARTESIAN","SINGLETON","DISKFILE"};
    std::string hname[] = {"Minimal","Conjugate Gradient","Evolutionary","Simulated Annealing","Simplex"};

    gethostname(hostname,80);

    std::ofstream s(log_file.c_str(),std::ios::trunc);
    s << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << std::endl;
    s << "<LogFile>" << std::endl;
    s << "<Library>Diaplexis (NTL)</Library>" << std::endl;
    s << "<Hostname>" << hostname << "</Hostname>" << std::endl;
    s << "<StartDate>" << start_time.date() << "</StartDate>" << std::endl;
    s << "<StartTime>" << start_time.time_of_day() << "</StartTime>" << std::endl;
    s << "<CompileTimeParameters>" << std::endl;
    s << "<MaximumDimension>" << Spacetime::ND << "</MaximumDimension>" << std::endl;
    s << "<AtomicPropositions>" << SYNARMOSMA::Proposition::get_clause_size() << "</AtomicPropositions>" << std::endl;
    s << "<TopologicalRadius>" << Spacetime::topological_radius << "</TopologicalRadius>" << std::endl;
    s << "<PolycosmicRamosity>" << Spacetime::ramosity << "</PolycosmicRamosity>" << std::endl;
    s << "<MachineEpsilon>" << Spacetime::epsilon << "</MachineEpsilon>" << std::endl;
    s << "<InitialAnnealingTemperature>" << Spacetime::T_zero << "</InitialAnnealingTemperature>" << std::endl;
    s << "<ThermalDecayRate>" << Spacetime::kappa << "</ThermalDecayRate>" << std::endl;
    s << "<EnergyCouplingConstant>" << Spacetime::Lambda << "</EnergyCouplingConstant>" << std::endl;
    s << "</CompileTimeParameters>" << std::endl;
    s << "<RunTimeParameters>" << std::endl;
    s << "<Global>" << std::endl;
    if (initial_state != SINGLETON) s << "<InitialVertices>" << initial_size << "</InitialVertices>" << std::endl;
    s << "<InitialSheets>" << nt_initial << "</InitialSheets>" << std::endl;
    s << "<MaxIterations>" << max_iter << "</MaxIterations>" << std::endl;
    s << "<EuclideanGeometry>" << bvalue[geometry->get_euclidean()] << "</EuclideanGeometry>" << std::endl;
    s << "<RelationalGeometry>" << bvalue[geometry->get_relational()] << "</RelationalGeometry>" << std::endl;
    s << "<DimensionalUniformity>" << bvalue[geometry->get_uniform()] << "</DimensionalUniformity>" << std::endl;
    s << "<BackgroundDimension>" << geometry->dimension() << "</BackgroundDimension>" << std::endl;
    s << "<SheetDynamics>" << bvalue[foliodynamics] << "</SheetDynamics>" << std::endl;
    s << "<RandomSeed>" << RND.get_seed() << "</RandomSeed>" << std::endl;
    s << "<Superposable>" << bvalue[superposable] << "</Superposable>" << std::endl;
    s << "<Compressible>" << bvalue[compressible] << "</Compressible>" << std::endl;
    s << "<Permutable>" << bvalue[permutable] << "</Permutable>" << std::endl;
    if (initial_state == DISKFILE) {
      s << "<InitialState>" << sname[initial_state] << " (" << sname[original_state] << ")</InitialState>" << std::endl;
    }
    else {
      s << "<InitialState>" << sname[initial_state] << "</InitialState>" << std::endl;
    }
    if (initial_state == RANDOM) s << "<EdgeProbability>" << edge_probability << "</EdgeProbability>" << std::endl;
    s << "<Checkpoint>" << bvalue[checkpoint] << "</Checkpoint>" << std::endl;
    if (checkpoint) s << "<CheckpointFrequency>" << checkpoint_frequency << "</CheckpointFrequency>" << std::endl;
    if (initial_state == DISKFILE) s << "<InputFile>" << input_file << "</InputFile>" << std::endl;
    if (initial_state == CARTESIAN) {
      s << "<PerturbTopology>" << bvalue[perturb_topology] << "</PerturbTopology>" << std::endl;
      s << "<PerturbGeometry>" << bvalue[perturb_geometry] << "</PerturbGeometry>" << std::endl;
      s << "<PerturbEnergy>" << bvalue[perturb_energy] << "</PerturbEnergy>" << std::endl;
    }
    if (H->get_method() == SYNARMOSMA::GAP) {
      s << "<HomologyMethod>GAP</HomologyMethod>" << std::endl;
    }
    else {
      s << "<HomologyMethod>Native</HomologyMethod>" << std::endl;
    }
    if (H->get_field() == SYNARMOSMA::INT) {
      s << "<HomologyBase>INT</HomologyBase>" << std::endl;
    }
    else if (H->get_field() == SYNARMOSMA::ZZ) {
      s << "<HomologyBase>NTL::ZZ</HomologyBase>" << std::endl;
    }
    else if (H->get_field() == SYNARMOSMA::GF2) {
      s << "<HomologyBase>GF2</HomologyBase>" << std::endl;
    }
    if (weaving == DYNAMIC) {
      s << "<Hyphansis>DYNAMIC</Hyphansis>" << std::endl; 
    }
    else {
      s << "<Hyphansis>MUSICAL</Hyphansis>" << std::endl;
      s << "<HyphansisScore>" << hyphansis_score << "</HyphansisScore>" << std::endl;
    }
    s << "</Global>" << std::endl;
    s << "<Geometry>" << std::endl;
    s << "<GeometryCutoff>" << geometry_cutoff << "</GeometryCutoff>" << std::endl;
    s << "<SolverType>" << hname[solver] << "</SolverType>" << std::endl;
    if (solver == MINIMAL) {
      s << "<SolverIterations>" << solver_its << "</SolverIterations>" << std::endl;
    }
    else if (solver == EVOLUTIONARY) {
      s << "<MaxGenerations>" << ngenerations << "</MaxGenerations>" << std::endl;
      s << "<MaxJousts>" << njousts << "</MaxJousts>" << std::endl;
      s << "<PoolSize>" << pool_size << "<</PoolSize>" << std::endl;
    }
    else if (solver == ANNEALING) {
      s << "<ThermalVariance>" << thermal_variance << "</ThermalVariance>" << std::endl;
      s << "<ThermalSweeps>" << thermal_sweep << "</ThermalSweeps>" << std::endl;
      s << "<AnnealingSteps>" << annealing_steps << "</AnnealingSteps>" << std::endl;
      s << "<ThermalizationCriterion>" << thermalization << "</ThermalizationCriterion>" << std::endl;
    }
    else if (solver == MECHANICAL) {
      s << "<IntegrationEngine>" << int_engine << "</IntegrationEngine>" << std::endl;
      s << "<StepSize>" << step_size << "</StepSize>" << std::endl;
      s << "<MaxIntegrationSteps>" << max_int_steps << "</MaxIntegrationSteps>" << std::endl;
      s << "<DampingConstant>" << damping_constant << "</DampingConstant>" << std::endl;
      s << "<SpringConstant>" << spring_constant << "</SpringConstant>" << std::endl;
      s << "<RepulsionConstant>" << repulsion_constant << "</RepulsionConstant>" << std::endl;
      s << "<ConjugateGradientRefinement>" << bvalue[cgradient_refinement] << "</ConjugateGradientRefinement>" << std::endl;
      if (cgradient_refinement) {
        s << "<MaxConjugateGradientSteps>" << max_CG_steps << "</MaxConjugateGradientSteps>" << std::endl;
        s << "<MaxLineSolverSteps>" << max_LS_steps << "</MaxLineSolverSteps>" << std::endl;
        s << "<EdgeFlexibilityThreshold>" << edge_flexibility_threshold << "</EdgeFlexibilityThreshold>" << std::endl;
      }
    }
    else if (solver == SIMPLEX) {
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

  nn = cardinality(0,-1);
  avg_ven = 0.0;
  max_ven = 0.0;
  min_ven = 10000.0;
  for(i=0; i<Nv; ++i) {
    if (events[i].ubiquity == 1) continue;
    w = events[i].energy;
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
    // (number of SYNARMOSMA::AFTER edges != number of SYNARMOSMA::BEFORE edges)?
    if (events[1].ubiquity == 1) continue;
    nf = 0;
    np = 0;
    atemporal = true;
    for(it=events[i].entourage.begin(); it!=events[i].entourage.end(); ++it) {
      in1 = *it;
      if (simplices[1][in1].orientation == SYNARMOSMA::DISPARATE) continue;
      rho = simplices[1][in1].orientation;
      simplices[1][in1].get_vertices(vx);
      if (i == vx[1]) rho = (rho == SYNARMOSMA::AFTER) ? SYNARMOSMA::BEFORE : SYNARMOSMA::AFTER;
      atemporal = false;
      if (rho == SYNARMOSMA::AFTER) {
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
    if (simplices[1][i].ubiquity == 1) continue;
    wm = simplices[1][i].volume;
    if (wm < Spacetime::epsilon) {
      nnull++;
    }
    else if (simplices[1][i].orientation == SYNARMOSMA::DISPARATE) {
      nspace++;
    }
    else {
      ntime++;
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
  nvalue = (pseudomanifold) ? "True" : "False";
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  if (pseudomanifold) {
    atom = rstep.append_child("Boundary");
    nvalue = (boundary) ? "True" : "False";
    atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

    atom = rstep.append_child("Orientable");
    nvalue = (orientable) ? "True" : "False";
    atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());
  }
  else {
    atom = rstep.append_child("Boundary");
    atom.append_child(pugi::node_pcdata).set_value("NULL");

    atom = rstep.append_child("Orientable");
    atom.append_child(pugi::node_pcdata).set_value("NULL");
  }

  atom = rstep.append_child("Homology");
  nvalue = H->write();
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  if (high_memory) {
    atom = rstep.append_child("Homotopy");
    nvalue = pi->write();
    atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());
  }

  ne = cardinality(1,-1);
  atom = rstep.append_child("Cyclomaticity");
  nvalue = boost::lexical_cast<std::string>(ne-nn+1);
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  atom = rstep.append_child("ActiveVertices");
  nvalue = boost::lexical_cast<std::string>(nn);
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  atom = rstep.append_child("SourceVertices");
  nvalue = boost::lexical_cast<std::string>(nsource);
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  atom = rstep.append_child("SinkVertices");
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
  nvalue = boost::lexical_cast<std::string>(cyclicity(-1));
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  arclength_statistics(vdata,-1);
  atom = rstep.append_child("MinArcLength");
  nvalue = boost::lexical_cast<std::string>(vdata[1]);
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  atom = rstep.append_child("AvgArcLength");
  nvalue = boost::lexical_cast<std::string>(vdata[2]);
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  atom = rstep.append_child("MaxArcLength");
  nvalue = boost::lexical_cast<std::string>(vdata[0]);
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  vertex_degree_statistics(vdata,-1);
  atom = rstep.append_child("MinVertexDegree");
  nvalue = boost::lexical_cast<std::string>(vdata[1]);
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  atom = rstep.append_child("AvgVertexDegree");
  nvalue = boost::lexical_cast<std::string>(vdata[2]);
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  atom = rstep.append_child("MaxVertexDegree");
  nvalue = boost::lexical_cast<std::string>(vdata[0]);
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  atom = rstep.append_child("TotalEnergy");
  nvalue = boost::lexical_cast<std::string>(total_energy(-1));
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  atom = rstep.append_child("MinVertexEnergy");
  nvalue = boost::lexical_cast<std::string>(min_ven);
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  atom = rstep.append_child("AvgVertexEnergy");
  nvalue = boost::lexical_cast<std::string>(avg_ven);
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  atom = rstep.append_child("MaxVertexEnergy");
  nvalue = boost::lexical_cast<std::string>(max_ven);
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  atom = rstep.append_child("AtemporalVertices");
  nvalue = boost::lexical_cast<std::string>(nat);
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  atom = rstep.append_child("ChiralVertices");
  nvalue = boost::lexical_cast<std::string>(nch);
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

  atom = rstep.append_child("EulerCharacteristic");
  nvalue = boost::lexical_cast<std::string>(euler_characteristic(-1));
  atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i,j)
#endif
  for(i=0; i<Nv; ++i) {
    for(j=0; j<Nt; ++j) {
      vdimension[Nt*i+j] = (NTL::divide(events[i].ubiquity,codex[j].colour) == 0) ? -1 : vertex_dimension(i,j);
    }
  }

  for(i=0; i<Nt; ++i) {
    sum = 0;
    mx = 0;
    mn = Spacetime::ND;
    for(j=0; j<Nv; ++j) {
      v = vdimension[Nt*j+i];
      if (v == -1) continue;
      sum += v;
      if (v > mx) mx = v;
      if (v < mn) mn = v;
    }
    vd_avg[i] = double(sum)/double(cardinality(0,i));
    vd_min[i] = mn;
    vd_max[i] = mx;
  }

  for(i=0; i<Nt; ++i) {
    sheet = rstep.append_child("Sheet");

    atom = sheet.append_child("Index");
    nvalue = boost::lexical_cast<std::string>(i+1);
    atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

    atom = sheet.append_child("State");
    atom.append_child(pugi::node_pcdata).set_value(bvalue[codex[i].active].c_str());

    nn = cardinality(0,i);
    ne = cardinality(1,i);
    ncyclic = cyclicity(i);

    avg_ven = 0.0;
    avg_val = 0.0;

    max_ven = 0.0;
    max_val = 0;

    min_ven = 10000.0;
    min_val = 1000;
    for(j=0; j<Nv; ++j) {
      if (NTL::divide(events[j].ubiquity,codex[i].colour) == 0) continue;
      min_ven = events[j].energy;
      min_val = vertex_valence(j,i);
      break;
    }
    for(j=0; j<Nv; ++j) {
      if (NTL::divide(events[j].ubiquity,codex[i].colour) == 0) continue;
      w = events[j].energy;
      if (w > max_ven) max_ven = w;
      if (w < min_ven) min_ven = w;
      avg_ven += w;
      l = vertex_valence(j,i);
      if (l > max_val) max_val = l;
      if (l < min_val) min_val = l;
      avg_val += double(l);
    }
    avg_val /= double(nn);
    avg_ven /= double(nn);

    nsink = 0;
    nsource = 0;
    nch = 0;
    nat = 0;
    for(j=0; j<Nv; ++j) {
      // Is this vertex atemporal (no timelike edges) or chiral
      // (number of SYNARMOSMA::AFTER edges != number of SYNARMOSMA::BEFORE edges)?
      if (NTL::divide(events[j].ubiquity,codex[i].colour) == 0) continue;
      nf = 0;
      np = 0;
      atemporal = true;
      for(it=events[j].entourage.begin(); it!=events[j].entourage.end(); ++it) {
        in1 = *it;
        if (simplices[1][in1].orientation == SYNARMOSMA::DISPARATE) continue;
        rho = simplices[1][in1].orientation;
        simplices[1][in1].get_vertices(vx);
        if (j == vx[1]) rho = (rho == SYNARMOSMA::AFTER) ? SYNARMOSMA::BEFORE : SYNARMOSMA::AFTER;
        atemporal = false;
        if (rho == SYNARMOSMA::AFTER) {
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
    max_length = 0.0;
    min_length = 1000.0;
    avg_length = 0.0;
    if (ne > 0) {
      for(j=0; j<Ne; ++j) {
        if (NTL::divide(simplices[1][j].ubiquity,codex[i].colour) == 0) continue;
        wm = simplices[1][j].volume;
        avg_length += wm;
        if (wm < Spacetime::epsilon) {
          nnull++;
        }
        else if (simplices[1][j].orientation == SYNARMOSMA::DISPARATE) {
          nspace++;
        }
        else {
          ntime++;
        }
      }
      avg_length /= double(ne);
      for(j=0; j<Ne; ++j) {
        if (NTL::divide(simplices[1][j].ubiquity,codex[i].colour) == 0) continue;
        wm = simplices[1][j].volume;
        if (wm > max_length) max_length = wm;
        if (wm < min_length) min_length = wm;
      }
    }

    atom = sheet.append_child("Pseudomanifold");
    nvalue = (codex[i].pseudomanifold) ? "True" : "False";
    atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

    if (codex[i].pseudomanifold) {
      atom = sheet.append_child("Boundary");
      nvalue = (codex[i].boundary) ? "True" : "False";
      atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

      atom = sheet.append_child("Orientable");
      nvalue = (codex[i].orientable) ? "True" : "False";
      atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());
    }
    else {
      atom = sheet.append_child("Boundary");
      atom.append_child(pugi::node_pcdata).set_value("NULL");

      atom = sheet.append_child("Orientable");
      atom.append_child(pugi::node_pcdata).set_value("NULL");
    }

    atom = sheet.append_child("Homology");
    nvalue = codex[i].H->write();
    atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

    if (high_memory) {
      atom = sheet.append_child("Homotopy");
      nvalue = codex[i].pi->write();
      atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());
    }

    atom = sheet.append_child("Cyclomaticity");
    nvalue = boost::lexical_cast<std::string>(ne-nn+1);
    atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

    atom = sheet.append_child("ActiveVertices");
    nvalue = boost::lexical_cast<std::string>(nn);
    atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

    atom = sheet.append_child("SourceVertices");
    nvalue = boost::lexical_cast<std::string>(nsource);
    atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

    atom = sheet.append_child("SinkVertices");
    nvalue = boost::lexical_cast<std::string>(nsink);
    atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

    atom = sheet.append_child("ActiveEdges");
    nvalue = boost::lexical_cast<std::string>(ne);
    atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

    atom = sheet.append_child("SpacelikeEdges");
    nvalue = boost::lexical_cast<std::string>(nspace);
    atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

    atom = sheet.append_child("TimelikeEdges");
    nvalue = boost::lexical_cast<std::string>(ntime);
    atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

    atom = sheet.append_child("NullEdges");
    nvalue = boost::lexical_cast<std::string>(nnull);
    atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

    atom = sheet.append_child("CyclicEdges");
    nvalue = boost::lexical_cast<std::string>(ncyclic);
    atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

    atom = sheet.append_child("MinArcLength");
    nvalue = boost::lexical_cast<std::string>(min_length);
    atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

    atom = sheet.append_child("AvgArcLength");
    nvalue = boost::lexical_cast<std::string>(avg_length);
    atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

    atom = sheet.append_child("MaxArcLength");
    nvalue = boost::lexical_cast<std::string>(max_length);
    atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

    atom = sheet.append_child("MinVertexDegree");
    nvalue = boost::lexical_cast<std::string>(min_val);
    atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

    atom = sheet.append_child("AvgVertexDegree");
    nvalue = boost::lexical_cast<std::string>(avg_val);
    atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

    atom = sheet.append_child("MaxVertexDegree");
    nvalue = boost::lexical_cast<std::string>(max_val);
    atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

    atom = sheet.append_child("MinVertexDimension");
    nvalue = boost::lexical_cast<std::string>(vd_min[i]);
    atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

    atom = sheet.append_child("AvgVertexDimension");
    nvalue = boost::lexical_cast<std::string>(vd_avg[i]);
    atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

    atom = sheet.append_child("MaxVertexDimension");
    nvalue = boost::lexical_cast<std::string>(vd_max[i]);
    atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

    atom = sheet.append_child("TotalEnergy");
    nvalue = boost::lexical_cast<std::string>(total_energy(i));
    atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

    atom = sheet.append_child("MinVertexEnergy");
    nvalue = boost::lexical_cast<std::string>(min_ven);
    atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

    atom = sheet.append_child("AvgVertexEnergy");
    nvalue = boost::lexical_cast<std::string>(avg_ven);
    atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

    atom = sheet.append_child("MaxVertexEnergy");
    nvalue = boost::lexical_cast<std::string>(max_ven);
    atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

    atom = sheet.append_child("AtemporalVertices");
    nvalue = boost::lexical_cast<std::string>(nat);
    atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

    atom = sheet.append_child("ChiralVertices");
    nvalue = boost::lexical_cast<std::string>(nch);
    atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

    atom = sheet.append_child("EulerCharacteristic");
    nvalue = boost::lexical_cast<std::string>(euler_characteristic(i));
    atom.append_child(pugi::node_pcdata).set_value(nvalue.c_str());

    atom = sheet.append_child("HyphanticSequence");
    if (codex[i].ops == "") {
      atom.append_child(pugi::node_pcdata).set_value("NULL");
    }
    else {
      atom.append_child(pugi::node_pcdata).set_value(codex[i].ops.c_str());
    }
  }
  logfile.save_file(log_file.c_str());
}

void Spacetime::read_complex(std::ifstream& s)
{
  int i,j,n;
  Event v;
  Simplex S;
  Sheet t;

  s.read((char*)(&n),sizeof(int));
  for(i=0; i<n; ++i) {
    v.deserialize(s);
    events.push_back(v);
  }

  for(i=1; i<=Spacetime::ND; ++i) {
    s.read((char*)(&n),sizeof(int));
    for(j=0; j<n; ++j) {
      S.deserialize(s);
      simplices[i].push_back(S);
    }
  }

  for(i=1; i<=Spacetime::ND; ++i) {
    for(j=0; j<(signed) simplices[i].size(); ++j) {
      index_table[i][simplices[i][j].vertices] = j;
    }
  }
  compute_entourages(-1);

  // Now the algebraic properties...
  H->deserialize(s);
  pi->deserialize(s);

  s.read((char*)(&j),sizeof(int));
  for(i=0; i<j; ++i) {
    codex.push_back(t);
  }
  nactive = 0;
  for(i=0; i<j; ++i) {
    psequence.next();
    codex[i].deserialize(s);
    if (codex[i].active) nactive++;
  }
}

void Spacetime::read_state(const std::string& filename)
{
  int i,j,n;
  std::string cmodel,fmodel;
  unsigned int q;
  double x;

  std::ifstream s(filename.c_str(),std::ios::binary | std::ios::in);
  if (!s.is_open()) {
    // File doesn't exist, print an error message and die
    std::cerr << "The file " << filename << " cannot be found." << std::endl;
    std::cerr << "Exiting..." << std::endl;
    std::exit(1);
  }

  clear();

  s.seekg(0,std::ios_base::beg);

  s.read((char*)(&n),sizeof(int));
  if (n != Spacetime::ND) {
    s.close();
    std::cerr << "The compiled binary's maximum simplicial dimension " << Spacetime::ND << " does not match that (" << n << ") of the data file." << std::endl;
    std::cerr << "Exiting..." << std::endl;
    std::exit(1);
  }

  s.read((char*)(&n),sizeof(int));
  if (n != SYNARMOSMA::Proposition::get_clause_size()) {
    s.close();
    std::cerr << "The compiled binary's atomic clause number " << SYNARMOSMA::Proposition::get_clause_size() << " does not match that (" << n << ") of the data file." << std::endl;
    std::cerr << "Exiting..." << std::endl;
    std::exit(1);
  }

  s.read((char*)(&n),sizeof(int));
  if (n != Spacetime::topological_radius) {
    s.close();
    std::cerr << "The compiled binary's topological radius " << Spacetime::topological_radius << " does not match that (" << n << ") of the data file." << std::endl;
    std::cerr << "Exiting..." << std::endl;
    std::exit(1);
  }

  s.read((char*)(&x),sizeof(double));
  if (!SYNARMOSMA::double_equality(x,Spacetime::ramosity)) {
    s.close();
    std::cerr << "The compiled binary's ramosity " << Spacetime::ramosity << " does not match that (" << x << ") of the data file." << std::endl;
    std::cerr << "Exiting..." << std::endl;
    std::exit(1);
  }

  s.read((char*)(&x),sizeof(double));
  if (!SYNARMOSMA::double_equality(x,Spacetime::epsilon)) {
    s.close();
    std::cerr << "The compiled binary's epsilon " << Spacetime::epsilon << " does not match that (" << x << ") of the data file." << std::endl;
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
  s.read((char*)(&nt_initial),sizeof(int));
  s.read((char*)(&initial_dim),sizeof(int));
  s.read((char*)(&max_iter),sizeof(int));
  s.read((char*)(&q),sizeof(int));
  s.read((char*)(&checkpoint),sizeof(bool));
  s.read((char*)(&checkpoint_frequency),sizeof(int));
  // Skip changing the initial_state as it should remain DISKFILE...
  s.read((char*)(&original_state),sizeof(TOPOLOGY));
  s.read((char*)(&foliodynamics),sizeof(bool));
  s.read((char*)(&instrument_convergence),sizeof(bool));
  s.read((char*)(&high_memory),sizeof(bool));
  s.read((char*)(&superposable),sizeof(bool));
  s.read((char*)(&compressible),sizeof(bool));
  s.read((char*)(&permutable),sizeof(bool));
  s.read((char*)(&perturb_topology),sizeof(bool));
  s.read((char*)(&perturb_geometry),sizeof(bool));
  s.read((char*)(&perturb_energy),sizeof(bool));
  s.read((char*)(&edge_probability),sizeof(double));

  s.read((char*)(&solver),sizeof(ALGORITHM));
  s.read((char*)(&geometry_cutoff),sizeof(double));
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
  s.read((char*)(&simplex_alpha),sizeof(double));
  s.read((char*)(&simplex_gamma),sizeof(double));
  s.read((char*)(&simplex_rho),sizeof(double));
  s.read((char*)(&simplex_sigma),sizeof(double));

  s.read((char*)(&system_size),sizeof(int));
  s.read((char*)(&iterations),sizeof(int));
  s.read((char*)(&global_deficiency),sizeof(double));
  s.read((char*)(&error),sizeof(double));
  s.read((char*)(&converged),sizeof(bool));
  s.read((char*)(&pseudomanifold),sizeof(bool));
  s.read((char*)(&boundary),sizeof(bool));
  s.read((char*)(&orientable),sizeof(bool));

  geometry->deserialize(s);

  read_complex(s);

  // Now read the anterior spacetime state and calculate the deltas...
  if (instrument_convergence) {
    Event v;
    Simplex S;

    s.read((char*)(&topology_delta),sizeof(double));
    s.read((char*)(&geometry_delta),sizeof(double));
    s.read((char*)(&energy_delta),sizeof(double));

    s.read((char*)(&n),sizeof(int));
    for(i=0; i<n; ++i) {
      v.deserialize(s);
      anterior.events.push_back(v);
    }

    for(i=1; i<=Spacetime::ND; ++i) {
      s.read((char*)(&n),sizeof(int));
      for(j=0; j<n; ++j) {
        S.deserialize(s);
        anterior.simplices[i].push_back(S);
      }
    }

    // Regenerate the anterior index table...
    for(i=1; i<=Spacetime::ND; ++i) {
      for(j=0; j<(signed) anterior.simplices[i].size(); ++j) {
        anterior.index_table[i][anterior.simplices[i][j].vertices] = j;
      }
    }
  }
  s.close();
  RND.set_seed(q);
}

void Spacetime::write_graph(const std::string& filename,int sheet) const
{
  int i,j,vx[2],N1 = 0;
  double E;
  std::vector<int> offset;
  const int nv = (signed) events.size();
  const int ne = (signed) simplices[1].size();

  for(i=0; i<nv; ++i) {
    offset.push_back(-1);
  }

  std::ofstream s(filename.c_str(),std::ios::binary | std::ios::trunc);
  i = cardinality(0,sheet);
  s.write((char*)(&i),sizeof(int));
  i = cardinality(1,sheet);
  s.write((char*)(&i),sizeof(int));
  if (sheet == -1) {
    // First calculate the appropriate vertex offset and write out the energy...
    for(i=0; i<nv; ++i) {
      if (events[i].ubiquity == 1) continue;
      offset[i] = N1;
      N1++;
      E = events[i].energy;
      s.write((char*)(&E),sizeof(double));
    }
    // Finally the edges...
    for(i=0; i<ne; ++i) {
      if (simplices[1][i].ubiquity == 1) continue;
      simplices[1][i].get_vertices(vx);
      j = offset[vx[0]];
      s.write((char*)(&j),sizeof(int));
      j = offset[vx[1]];
      s.write((char*)(&j),sizeof(int));
    }
  }
  else {
    // First calculate the appropriate vertex offset and write out the energy...
    for(i=0; i<nv; ++i) {
      if (NTL::divide(events[i].ubiquity,codex[sheet].colour) == 0) continue;
      offset[i] = N1;
      N1++;
      E = events[i].energy;
      s.write((char*)(&E),sizeof(double));
    }
    // Finally the edges...
    for(i=0; i<ne; ++i) {
      if (NTL::divide(simplices[1][i].ubiquity,codex[sheet].colour) == 0) continue;
      simplices[1][i].get_vertices(vx);
      j = offset[vx[0]];
      s.write((char*)(&j),sizeof(int));
      j = offset[vx[1]];
      s.write((char*)(&j),sizeof(int));
    }
  }
  s.close();
}

void Spacetime::write_complex(std::ofstream& s) const
{
  int i,j,n;
  const int nt = (signed) codex.size();

  n = (signed) events.size();
  s.write((char*)(&n),sizeof(int));
  for(i=0; i<n; ++i) {
    events[i].serialize(s,nt);
  }

  for(i=1; i<=Spacetime::ND; ++i) {
    n = (signed) simplices[i].size();
    s.write((char*)(&n),sizeof(int));
    for(j=0; j<n; ++j) {
      simplices[i][j].serialize(s,nt);
    }
  }
  // Now the algebraic properties...
  H->serialize(s);
  pi->serialize(s);

  // Finally the data for each Sheet instance...
  j = (signed) codex.size();
  s.write((char*)(&j),sizeof(int));
  for(i=0; i<j; ++i) {
    codex[i].serialize(s);
  }
}

void Spacetime::write_state() const
{
  int i,j,n = SYNARMOSMA::Proposition::get_clause_size();
  const int nt = (signed) codex.size();

  std::stringstream sstream;
  sstream << iterations;
  std::string filename = state_file + "_" + sstream.str() + ".dat";

  std::ofstream s(filename.c_str(),std::ios::binary | std::ios::trunc);
  s.seekp(0,std::ios_base::beg);

  // First the global parameters...
  s.write((char*)(&Spacetime::ND),sizeof(int));
  s.write((char*)(&n),sizeof(int));
  s.write((char*)(&Spacetime::topological_radius),sizeof(int));
  s.write((char*)(&Spacetime::ramosity),sizeof(double));
  s.write((char*)(&Spacetime::epsilon),sizeof(double));
  s.write((char*)(&Spacetime::T_zero),sizeof(double));
  s.write((char*)(&Spacetime::kappa),sizeof(double));
  s.write((char*)(&Spacetime::Lambda),sizeof(double));

  // Now the global runtime parameters specific to the
  // single Spacetime instance...
  s.write((char*)(&initial_size),sizeof(int));
  s.write((char*)(&nt_initial),sizeof(int));
  s.write((char*)(&initial_dim),sizeof(int));
  s.write((char*)(&max_iter),sizeof(int));
  // Write the initial random number seed to disk...
  unsigned int q = RND.get_seed();
  s.write((char*)(&q),sizeof(int));
  s.write((char*)(&checkpoint),sizeof(bool));
  s.write((char*)(&checkpoint_frequency),sizeof(int));
  s.write((char*)(&initial_state),sizeof(TOPOLOGY));
  s.write((char*)(&foliodynamics),sizeof(bool));
  s.write((char*)(&instrument_convergence),sizeof(bool));
  s.write((char*)(&high_memory),sizeof(bool));
  s.write((char*)(&superposable),sizeof(bool));
  s.write((char*)(&compressible),sizeof(bool));
  s.write((char*)(&permutable),sizeof(bool));
  s.write((char*)(&perturb_topology),sizeof(bool));
  s.write((char*)(&perturb_geometry),sizeof(bool));
  s.write((char*)(&perturb_energy),sizeof(bool));
  s.write((char*)(&edge_probability),sizeof(double));

  // Geometric runtime constants...
  s.write((char*)(&solver),sizeof(ALGORITHM));
  s.write((char*)(&geometry_cutoff),sizeof(double));
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
  s.write((char*)(&pseudomanifold),sizeof(bool));
  s.write((char*)(&boundary),sizeof(bool));
  s.write((char*)(&orientable),sizeof(bool));

  // The principal body of the file...
  geometry->serialize(s);

  write_complex(s);

  // Now write the deltas and the anterior spacetime state...
  if (instrument_convergence) {
    s.write((char*)(&topology_delta),sizeof(double));
    s.write((char*)(&geometry_delta),sizeof(double));
    s.write((char*)(&energy_delta),sizeof(double));

    n = (signed) anterior.events.size();
    s.write((char*)(&n),sizeof(int));
    for(i=0; i<n; ++i) {
      anterior.events[i].serialize(s,nt);
    }

    for(i=1; i<=Spacetime::ND; ++i) {
      n = (signed) anterior.simplices[i].size();
      s.write((char*)(&n),sizeof(int));
      for(j=0; j<n; ++j) {
        anterior.simplices[i][j].serialize(s,nt);
      }
    }
  }

  // Close the file and return...
  s.close();
}
