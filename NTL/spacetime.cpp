#include "spacetime.h"

Random RND;

const double Spacetime::ramosity = 0.15;
const double Spacetime::epsilon = 0.00001;
const double Spacetime::T_zero = 500.0;
const double Spacetime::kappa = 1.35;
const int Spacetime::N_EXP;
const int Spacetime::N_IMP;
const char Spacetime::EXP_OP[] = {'D','U','S','R','C','N','A','G'};
const char Spacetime::IMP_OP[] = {'F','U','O','E','I','P','V'};
const int Spacetime::topological_radius;

const double seqn_weights[] = {1.0,0.0,0.2,0.2};

Spacetime::Spacetime()
{
  set_default_values();
  initialize();
}

Spacetime::Spacetime(bool no_disk)
{
  set_default_values();
  diskless = no_disk;
  if (diskless) checkpoint = false;
  initialize();
}

Spacetime::Spacetime(const char* filename)
{
  set_default_values();
  read_parameters(filename);
  initialize();
}

Spacetime::Spacetime(const char* filename,bool no_disk)
{
  set_default_values();
  read_parameters(filename);
  diskless = no_disk;
  if (diskless) checkpoint = false;
  initialize();
}

Spacetime::~Spacetime()
{
  events.clear();
  codex.clear();
  HZ.clear();
  pi1.clear();
  for(int i=1; i<=ND; ++i) {
    simplices[i].clear();
    index_table[i].clear();
  }
  anterior.events.clear();
  for(int i=1; i<=ND; ++i) {
    anterior.simplices[i].clear();
    anterior.index_table[i].clear();
  }
  delete geometry;
}

void Spacetime::set_default_values()
{
  perturb_geometry = false;
  perturb_topology = false;
  perturb_energy = true;
  ngenerations = 1000;
  solver_its = 50;
  njousts = 20;
  pool_size = 100;
  thermal_variance = 0.5;
  thermal_sweep = 1000;
  annealing_steps = 500;
  thermalization = 0.001;
  int_engine = "RK4";
  max_int_steps = 10000;
  step_size = 0.05;
  spring_constant = -1.5;
  repulsion_constant = 1.0;
  damping_constant = 0.85;
  cgradient_refinement = true;
  max_CG_steps = 10;
  max_LS_steps = 20;
  edge_flexibility_threshold = 2.0;
  simplex_alpha = 1.0;
  simplex_gamma = 2.0;
  simplex_rho = 0.5;
  simplex_sigma = 0.5;
  system_size = 0;
  error = 0.0;
  global_deficiency = 0.0;
  topology_delta = 0.0;
  geometry_delta = 0.0;
  energy_delta = 0.0;
  state_file = std::string("data/spacetime");
  log_file = std::string("data/spacetime");
  input_file = std::string("");
  initial_size = 500;
  initial_dim = 4;
  initial_state = RANDOM;
  solver = MINIMAL;
  converged = false;
  permutable = false;
  compressible = false;
  superposable = false;
  foliodynamics = false;
  checkpoint = true;
  checkpoint_frequency = 50;
  iterations = 0;
  geometry_cutoff = 0.0001;
  max_iter = 50;
  edge_probability = std::log(500.0)/250.0;
  nactive = 1;
  nt_initial = 1;
  homology_method = NATIVE;
  homology_base = GF2;
  pseudomanifold = false;
  boundary = false;
  orientable = false;
  diskless = false;
  original_state = RANDOM;
  geometry = new Geometry(false);
}

void Spacetime::set_checkpoint_frequency(int a)
{
  checkpoint = true;
  checkpoint_frequency = a;
}

void Spacetime::restart(const char* filename,bool save_seed)
{
  if (save_seed) {
    unsigned int n = RND.get_seed();
    clear();
    read_parameters(filename);
    RND.set_seed(n);
  }
  else {
    clear();
    read_parameters(filename);
  }
  initialize();
}

void Spacetime::clear()
{
  events.clear();
  codex.clear();
  HZ.clear();
  pi1.clear();
  geometry->clear();
  for(int i=1; i<=ND; ++i) {
    simplices[i].clear();
    index_table[i].clear();
  }
  anterior.events.clear();
  for(int i=1; i<=ND; ++i) {
    anterior.simplices[i].clear();
    anterior.index_table[i].clear();
  }
  psequence.reset(1);
  iterations = 0;
}

void Spacetime::distribute(int nprocs) const
{
  assert(nprocs > 0);
  const int nv = (signed) events.size();
  int i,j,k,n,p,p_old,ecount,bcount,its = 0,affinity[nv],volume[nprocs],max_dim = 0,current = -1,cproc = 0,nreal = 0;
  bool done,good,bdry;
  double cost,current_cost;
  std::set<int> vx,neg,next;
  std::set<int>::const_iterator it;
  const int D = dimension(-1);

  // The algorithm should begin by making current equal to the 
  // vertex with the highest simplicial dimension...
  for(i=0; i<nv; ++i) {
    affinity[i] = -1;
    if (events[i].ubiquity == 1) continue;
    nreal++;
    if (events[i].global_dimension > max_dim) {
      current = i;
      max_dim = events[i].global_dimension;
    }
  }
  for(i=0; i<nprocs; ++i) {
    volume[i] = 0;
  }
  std::cout << "Distributing " << nreal << " vertices across " << nprocs << " processor elements, with maximum simplicial dimension of " << max_dim << "." << std::endl;
  // First do the higher-dimensional vertices...
  do {
    // So take all of the n-simplices (n > 1) that are active and which 
    // contain the current vertex and assign all of their vertices to the 
    // current processor...
    if (events[current].global_dimension > 1) {
      affinity[current] = cproc;
      volume[cproc] += 1;
      // Need to loop over all d-simplices, d > 1
      for(i=2; i<=events[current].global_dimension; ++i) {
        for(j=0; j<(signed) simplices[i].size(); ++j) {
          if (simplices[i][j].ubiquity == 1) continue;
          if (simplices[i][j].contains(current)) {
            simplices[i][j].get_vertices(vx);
            for(it=vx.begin(); it!=vx.end(); ++it) {
              if (affinity[*it] == -1) {
                affinity[*it] = cproc;
                volume[cproc] += 1;
              }
              else {
                assert(affinity[*it] == cproc);
              }
            }
          }
        }
      }
      // Now handle the neighbouring d-simplices (d > 1)...
      do {
        for(i=2; i<=D; ++i) {
          for(j=0; j<(signed) simplices[i].size(); ++j) {
            if (simplices[i][j].ubiquity == 1) continue;
            simplices[i][j].get_vertices(vx);
            bdry = false;
            for(it=vx.begin(); it!=vx.end(); ++it) {
              if (affinity[*it] == -1) {
                neg.insert(*it);
              }
              else {
                assert(affinity[*it] == cproc);
                bdry = true;
              }
            }
            if (!neg.empty() && bdry) {
              for(it=neg.begin(); it!=neg.end(); ++it) {
                next.insert(*it);
              }
            }
            neg.clear(); 
          }
        }
        if (next.empty()) break;
        for(it=next.begin(); it!=next.end(); ++it) {
          affinity[*it] = cproc;
        }  
        next.clear();
      } while(true);
    }    
    done = true;
    for(i=0; i<nv; ++i) {
      if (events[i].ubiquity == 1) continue;
      if (affinity[i] == -1 && events[i].global_dimension > 1) {
        done = false;
        current = i;
        break;
      }
    }
    cproc = (cproc + 1)%nprocs;
  } while(!done);
  // Now everything else, adding vertices to equilibrate the population counts 
  // and minimize the number of boundary edges...
  for(i=0; i<nprocs; ++i) {
    if (volume[i] == 0) {
      // Find an initial vertex, ideally far from any existing vertices that have 
      // been assigned to a processor...
      do {
        n = RND.irandom(nv);
        if (events[n].ubiquity == 1) continue;
        if (affinity[n] > -1) continue;
        good = true;
        for(it=events[n].neighbours.begin(); it!=events[n].neighbours.end(); ++it) {
          if (affinity[*it] > -1) {
            good = false;
            break;
          }
        }
        if (good) break;
      } while(true);
      affinity[n] = i;
      volume[i] += 1;
    }
    done = false;
    do {
      next.clear();
      for(j=0; j<nv; ++j) {
        if (events[j].ubiquity == 1) continue;
        if (affinity[j] > -1) continue;
        for(it=events[j].neighbours.begin(); it!=events[j].neighbours.end(); ++it) {
          if (affinity[*it] == i) next.insert(j);
        }
      }
      if (next.empty()) break;
      for(it=next.begin(); it!=next.end(); ++it) {
        assert(affinity[*it] == -1);
        affinity[*it] = i;
        volume[i] += 1;
        if (volume[i] >= nreal/nprocs) {
          done = true;
          break;
        }
      }
      if (done) break;  
    } while(true);
  }
  for(i=0; i<nv; ++i) {
    if (events[i].ubiquity == 1) continue;
    if (affinity[i] == -1) {
      for(it=events[i].neighbours.begin(); it!=events[i].neighbours.end(); ++it) {
        p = affinity[*it];
        if (p > -1) {
          affinity[i] = p;
          volume[p] += 1;
          break;
        }
      }
    }
  }
  current_cost = distribution_fitness(volume,affinity,nprocs);
  std::cout << "At iteration 0 the cost is " << current_cost << std::endl;
  do {
    do {
      n = RND.irandom(nv);
      if (events[n].ubiquity == 1) continue;
      if (events[n].global_dimension > 1) continue;
      break;
    } while(true);
    p_old = affinity[n];
    do {
      p = RND.irandom(nprocs);
      if (p != p_old) break;
    } while(true);
    affinity[n] = p;
    volume[p] += 1;
    volume[p_old] -= 1;
    cost = distribution_fitness(volume,affinity,nprocs);
    if (cost < current_cost) {
      current_cost = cost;
    }
    else {
      affinity[n] = p_old;
      volume[p] -= 1;
      volume[p_old] += 1;
    }
    its++;
    if (its % 250 == 0) std::cout << "At iteration " << its << " the cost is " << current_cost << std::endl; 
  } while(its < 10000);
  // Some basic sanity checks: every vertex has a processor...
  n = 0;
  for(i=0; i<nv; ++i) {
    if (events[i].ubiquity == 1) continue;
    if (affinity[i] == -1) n++;
  }
  assert(n == 0);
  // And the processors haven't overcounted the vertices.
  n = 0;
  for(i=0; i<nprocs; ++i) {
    n += volume[i];
  }
  assert(n == nreal);
  // Analysis of the distribution of vertices and edges among the processors...
  for(i=0; i<nprocs; ++i) {
    ecount = 0; bcount = 0;
    for(j=0; j<nv; ++j) {
      if (events[j].ubiquity == 1) continue;
      if (affinity[j] != i) continue;
      for(it=events[j].neighbours.begin(); it!=events[j].neighbours.end(); ++it) {
        k = *it;
        if (affinity[k] == i) {
          ecount++;
        }
        else {
          bcount++;
        }
      }
    } 
    std::cout << "Processor " << i << " owns " << volume[i] << " vertices and " << ecount << " edges, with " << bcount << " boundary edges." << std::endl;
  }
}

void Spacetime::get_ubiquity_vector(std::vector<int>& output) const
{
  int i,j,n;
  const int nv = (signed) events.size();
  const int ne = (signed) simplices[1].size();
  const int nt = (signed) codex.size();

  output.clear();

  for(i=0; i<nv; ++i) {
    if (events[i].ubiquity == 1) continue;
    for(j=0; j<nt; ++j) {
      n = (NTL::divide(events[i].ubiquity,codex[j].colour) == 1) ? 1 : 0;
      output.push_back(n);
    }
  }
  for(i=0; i<ne; ++i) {
    if (simplices[1][i].ubiquity == 1) continue;
    for(j=0; j<nt; ++j) {
      n = (NTL::divide(simplices[1][i].ubiquity,codex[j].colour) == 1) ? 1 : 0;
      output.push_back(n);
    }
  }
}

void Spacetime::get_energy_values(std::vector<double>& output,int sheet) const
{
  int i;
  const int nv = (signed) events.size();

  output.clear();

  if (sheet == -1) {
    for(i=0; i<nv; ++i) {
      if (events[i].ubiquity == 1) continue;
      output.push_back(events[i].energy);
    }
  } 
  else {
    for(i=0; i<nv; ++i) {
      if (NTL::divide(events[i].ubiquity,codex[sheet].colour) == 0) continue;
      output.push_back(events[i].energy);
    }
  }
}

void Spacetime::get_deficiency_values(std::vector<double>& output,int sheet) const
{
  int i;
  const int nv = (signed) events.size();

  output.clear();

  if (sheet == -1) {
    for(i=0; i<nv; ++i) {
      if (events[i].ubiquity == 1) continue;
      output.push_back(events[i].deficiency);
    }
  }
  else {
    for(i=0; i<nv; ++i) {
      if (NTL::divide(events[i].ubiquity,codex[sheet].colour) == 0) continue;
      output.push_back(events[i].deficiency);
    }
  }
}

bool Spacetime::edge_exists(int u,int v,int sheet) const
{
  hash_map::const_iterator qt = index_table[1].find(make_key(u,v));
  if (qt == index_table[1].end()) return false;
  if (sheet == -1) {
    if (simplices[1][qt->second].ubiquity == 1) return false;
  }
  else {
    if (NTL::divide(simplices[1][qt->second].ubiquity,codex[sheet].colour) == 0) return false;
  }
  return true;
}

std::string Spacetime::sheet_activity() const
{
  std::string out = "(";
  for(int i=0; i<(signed) codex.size()-1; ++i) {
    out += (boost::lexical_cast<std::string>(codex[i].active) + ",");
  }
  out += boost::lexical_cast<std::string>(codex[codex.size()-1].active);
  out += ")";
  return out;
}

char Spacetime::implication() const
{
  // Should return one of {F,U,O,E,I,P,V}
  double alpha;
  if (iterations < 50) {
    alpha = RND.drandom();
    if (alpha < 0.3) {
      return 'F';
    }
    else if (alpha < 0.6) {
      if (RND.drandom() < 0.5) {
        return 'E';
      }
      else {
        return 'I';
      }
    }
    else if (alpha < 0.75) {
      return 'O';
    }
    else if (alpha < 0.9) {
      return 'U';
    }
    else {
      if (RND.drandom() < 0.33) {
        return 'P';
      }
      else {
        return 'V';
      }
    }
  }
  else {
    if (RND.drandom() < 0.5) {
      if (RND.drandom() < 0.5) {
        return 'P';
      }
      else {
        return 'V';
      }
    }
    else {
      alpha = RND.drandom();
      if (alpha < 0.4) {
        return 'O';
      }
      else if (alpha < 0.6) {
        return 'U';
      }
      else if (alpha < 0.8) {
        return 'F';
      }
      else {
        if (RND.drandom() < 0.67) {
          return 'E';
        }
        else {
          return 'I';
        }
      }
    }
  }
}

char Spacetime::explication() const
{
  // Should return one of {U,D,S,R,G,A,C,N}
  //if (RND.drandom() < 0.1) return 'G';
  double alpha;
  if (iterations < 50) {
    if (RND.drandom() < (0.35 + 1.0/(1+iterations/2))) {
      return 'S';
    }
    else {
      if (iterations <= 10) return 'C';
      if (RND.drandom() < 0.5) {
        return 'C';
      }
      else {
        return 'A';
      }
    }
  }
  else {
    alpha = RND.drandom();
    if (alpha < 0.25) {
      return 'C';
    }
    else if (alpha < 0.5) {
      return 'S';
    }
    else if (alpha < 0.75) {
      return 'N';
    }
    else {
      return 'U';
    }
  }
}

void Spacetime::clean() const
{
  int i,j;
  for(i=0; i<(signed) events.size(); ++i) {
    assert(!events[i].topology_modified);
  }
  for(i=1; i<=ND; ++i) {
    for(j=0; j<(signed) simplices[i].size(); ++j) {
      assert(!simplices[i][j].modified);
    }
  }
}

bool Spacetime::step_forwards()
{
  int i,n = (signed) codex.size();
  std::vector<int> order;
  timeval Z;
  double t1,t2;
  static double htime = 0.0;
  static double gtime = 0.0;

  gettimeofday(&Z,NULL);
  t1 = Z.tv_sec + (Z.tv_usec/1000000.0);
  RND.shuffle(order,n);
  for(i=0; i<n; ++i) {
    if (!codex[order[i]].active) continue;
    hyphansis(order[i]);
  }
  regularization(false,-1);
  geometry->compute_distances();
  gettimeofday(&Z,NULL);
  t2 = Z.tv_sec + (Z.tv_usec/1000000.0);
  htime += (t2 - t1);
#ifdef VERBOSE
  std::cout << "Calling global operations..." << std::endl;
#endif
  gettimeofday(&Z,NULL);
  t1 = Z.tv_sec + (Z.tv_usec/1000000.0);
  bool out = global_operations();
  gettimeofday(&Z,NULL);
  t2 = Z.tv_sec + (Z.tv_usec/1000000.0);
  gtime += (t2 - t1);
  if (out) {
    std::cout << boost::format("%-45s %8.3f %-10s\n") % "Time required for topological hyphansis was" % htime % "seconds.";
    std::cout << boost::format("%-45s %8.3f %-10s\n") % "Time required for global operations was" % gtime % "seconds.";
  }
  return out;
}

void Spacetime::test_harness(int type,int n)
{
  int i,j,d,nv;
  double alpha;
  hash_map::const_iterator qt;
  std::set<int> vx;
  std::vector<double> xc;
  unsigned long sheet = 2;
  Simplex S;

  for(i=1; i<=ND; ++i) {
    simplices[i].clear();
    index_table[i].clear();
  }

  if (type == 0) {
    nv = (signed) events.size();

    for(i=0; i<nv; ++i) {
      events[i].entourage.clear();
      events[i].topology_modified = true;
      events[i].geometry_modified = true;
    }

    // Randomize the 2D coordinates of the initial vertices...
    for(i=0; i<Geometry::background_dimension; ++i) {
      xc.push_back(0.0);
    }
    for(i=0; i<nv; ++i) {
      for(j=0; j<Geometry::background_dimension; ++j) {
        xc[j] = RND.nrandom(0.0,2.5);
      }
      geometry->set_coordinates(i,xc);
    }
  }
  else {
    Vertex vt;

    nv = n;
    events.clear();

    vt.ubiquity = sheet;
    for(i=0; i<n; ++i) {
      events.push_back(vt);
      geometry->add_vertex(-1);
    }
  }

  S.ubiquity = sheet;
  do {
    alpha = RND.drandom();
    if (alpha < 0.5) {
      d = 2;
    }
    else if (alpha < 0.75) {
      d = 3;
    }
    else if (alpha < 0.9) {
      d = 4;
    }
    else if (alpha < 0.98) {
      d = 5;
    }
    else {
      d = 6;
    }
    do {
      vx.insert(RND.irandom(nv));
    } while((signed) vx.size() < (1+d));
    S.initialize(vx,sheet);
    qt = index_table[d].find(S.key);
    if (qt == index_table[d].end()) {
      simplices[d].push_back(S);
      index_table[d][S.key] = (signed) simplices[d].size() - 1;
      simplicial_implication(0);
      simplicial_implication(-1);
      compute_neighbours();
    }
    vx.clear();
  } while (!connected(-1));
  regularization(false,0);
  regularization(false,-1);
  assert(consistent(-1));
  adjust_dimension();

  compute_volume();
  compute_curvature();
  compute_obliquity();
  structural_deficiency();
#ifdef VERBOSE
  std::cout << "At relaxation step " << iterations << " the global error is " << error << std::endl;
  std::cout << "Sheet activity " << sheet_activity() << std::endl;
#endif
}

void Spacetime::evolve()
{
  do {
    if (step_forwards()) break;
  } while(true);

  if (diskless) return;

  // Write the spacetime state to disk
  write_state();

  // A final message for the logfile: the completion time...
  std::string nvalue;
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();

  pugi::xml_document logfile;
  pugi::xml_node rstep;

  logfile.load_file(log_file.c_str());
  rstep = logfile.child("LogFile").insert_child_after("FinishDate",logfile.child("LogFile").child("StartTime"));
  nvalue = boost::lexical_cast<std::string>(now.date());
  rstep.append_child(pugi::node_pcdata).set_value(nvalue.c_str());
  logfile.save_file(log_file.c_str());

  logfile.load_file(log_file.c_str());
  rstep = logfile.child("LogFile").insert_child_after("FinishTime",logfile.child("LogFile").child("FinishDate"));
  nvalue = boost::lexical_cast<std::string>(now.time_of_day());
  rstep.append_child(pugi::node_pcdata).set_value(nvalue.c_str());
  logfile.save_file(log_file.c_str());

  logfile.load_file(log_file.c_str());
  rstep = logfile.child("LogFile").insert_child_after("Converged",logfile.child("LogFile").child("FinishTime"));
  nvalue = (converged) ? "True" : "False";
  rstep.append_child(pugi::node_pcdata).set_value(nvalue.c_str());
  logfile.save_file(log_file.c_str());
}

bool Spacetime::correctness()
{
  int i,j;
  const int nv = (signed) events.size();

  compute_volume();
  compute_curvature();
  compute_obliquity();
  structural_deficiency();
  double anterior_error = error;

  for(i=0; i<nv; ++i) {
    events[i].topology_modified = true;
    events[i].geometry_modified = true;
  }
  for(i=1; i<=ND; ++i) {
    for(j=0; j<(signed) simplices[i].size(); ++j) {
      simplices[i][j].modified = true;
    }
  }

  compute_volume();
  compute_curvature();
  compute_obliquity();
  structural_deficiency();
  if (std::abs(anterior_error - error) < Spacetime::epsilon) return true;
  std::cout << "Error difference is " << anterior_error << "  " << error << "  " << std::abs(anterior_error - error) << std::endl;
  return false;
}

void Spacetime::write_vertex_data(int v) const
{
  if (v < 0 || v >= (signed) events.size()) return;
  std::cout << "For vertex " << v << " we have:" << std::endl;
  std::cout << "    Incept = " << events[v].incept << std::endl;
  std::cout << "    Structural deficiency = " << events[v].deficiency << std::endl;
  std::cout << "    Energy = " << events[v].energy << std::endl;
  std::cout << "    Orthogonality = " << events[v].obliquity << std::endl;
  std::cout << std::endl;
}

void Spacetime::write_topology(int sheet) const
{
  int i,j,n;
  std::vector<int> vx;
  UINT64 q;
  const int ns = cardinality(0,sheet);
  const int ulimit = dimension(sheet);

  if (sheet == -1) {
    for(i=0; i<(signed) events.size(); ++i) {
      if (events[i].ubiquity == 1) continue;
      simplex_membership(i,vx);
      std::cout << i+1 << ": [";
      for(j=1; j<ulimit; ++j) {
        std::cout << vx[j-1] << ",";
      }
      std::cout << vx[ulimit-1] << "]" << std::endl;
    }
    for(i=ND; i>=1; i--) {
      n = cardinality(i,sheet);
      if (n > 0) {
        if (i > 0) {
          q = 1;
          for(j=ns; j>=ns-i; j--) {
            q *= (UINT64) j;
          }
          q = q/factorial(i+1);
        }
        else {
          q = (UINT64) ns;
        }
        std::cout << "There are " << n << " (" << q << ") " << i << "-simplices in this complex." << std::endl;
      }
    }
  }
  else {
    for(i=0; i<(signed) events.size(); ++i) {
      if (NTL::divide(events[i].ubiquity,codex[sheet].colour) == 0) continue;
      simplex_membership(i,vx);
      std::cout << i+1 << ": [";
      for(j=1; j<ulimit; ++j) {
        std::cout << vx[j-1] << ",";
      }
      std::cout << vx[ulimit-1] << "]" << std::endl;
    }
    for(i=ND; i>=1; i--) {
      n = cardinality(i,sheet);
      if (n > 0) {
        if (i > 0) {
          q = 1;
          for(j=ns; j>=ns-i; j--) {
            q *= (UINT64) j;
          }
          q = q/factorial(i+1);
        }
        else {
          q = (UINT64) ns;
        }
        std::cout << "There are " << n << " (" << q << ") " << i << "-simplices in this complex." << std::endl;
      }
    }
  }
  std::cout << "There are " << ns << " (" << ns << ") 0-simplices in this complex." << std::endl;
  std::cout << "The Euler characteristic is " << euler_characteristic(sheet) << std::endl;
}

void Spacetime::write_incastrature(const std::string& filename,int sheet) const
{
  // A method that generates the Hasse diagram corresponding to the
  // complex's simplicial structure, with the diagram stored in the
  // PDF "hasse.pdf"; the method assumes that the Graphviz library
  // has been installed on the system.
  int i,j,k;
  std::string nm;

  std::ofstream s(filename.c_str(),std::ios::trunc);

  s << "digraph G {" << std::endl;
  if (sheet == -1) {
    for(i=ND; i>0; i--) {
      for(j=0; j<(signed) simplices[i].size(); ++j) {
        if (simplices[i][j].ubiquity == 1) continue;
        nm = simplices[i][j].key;
        for(k=0; k<1+i; ++k) {
          s << "  \"" << nm << "\" -> \"" << simplices[i][j].faces[k] << "\";" << std::endl;
        }
      }
    }
    for(i=0; i<(signed) events.size(); ++i) {
      if (events[i].ubiquity == 1) continue;
      s << "  \"" << i << "\" -> \"NULL\";" << std::endl;
    }
  }
  else {
    for(i=ND; i>0; i--) {
      for(j=0; j<(signed) simplices[i].size(); ++j) {
        if (NTL::divide(simplices[i][j].ubiquity,codex[sheet].colour) == 0) continue;
        nm = simplices[i][j].key;
        for(k=0; k<1+i; ++k) {
          s << "  \"" << nm << "\" -> \"" << simplices[i][j].faces[k] << "\";" << std::endl;
        }
      }
    }
    for(i=0; i<(signed) events.size(); ++i) {
      if (NTL::divide(events[i].ubiquity,codex[sheet].colour) == 0) continue;
      s << "  \"" << i << "\" -> \"NULL\";" << std::endl;
    }
  }
  s << "}" << std::endl;
  s.close();
}

bool Spacetime::energy_check() const
{
  bool output = true;
  const int nvertex = (signed) events.size();

  for(int i=0; i<nvertex; ++i) {
    if (events[i].energy > 0.0) {
      if (events[i].ubiquity == 1) {
        std::cout << "Potential problem here: " << i << "  " << events[i].energy << std::endl;
        output = false;
      }
    }
  }
  return output;
}

void Spacetime::structural_deficiency()
{
  /*
  \begin{equation}
  \Lamba_G(v) = \frac {1}{NT} \sum_{i=1}^NT \lamba_i(v)
  \end{equation}
  where $\lambda_i(v)$ is the entwinement of the graph centred on the vertex v
  for the sheet $i$ and $NT$ is the total number of sheets.

  The topological quantities must be calculated over the whole ensemble
  of sheets, to be used in conjunction with the global quantities: the
  geometry (vertex coordinates) and energy.
  */
  int i,j,k,c,n;
  double sum,sum1,sum2,l,l_inv,d1,d2,delta,E_G,E_total = 0.0;
  bool found;
  std::set<int>::const_iterator it;
  hash_map::const_iterator qt;
  Graph* G = new Graph;
  const double lambda = 1.0;
  const double pfactor = 5.5;
  const double na = double(cardinality(0,-1));
  const int nv = (signed) events.size();
  const int nt = (signed) codex.size();
  double R[nv],gvalue[nv],Lambda[nv],rho[nv];

  for(i=0; i<nv; ++i) {
    events[i].deficiency = 0.0;
    events[i].geometric_deficiency = 0.0;
    Lambda[i] = 0.0;
    gvalue[i] = 0.0;
    R[i] = 0.0;
    rho[i] = 0.0;
  }

  for(i=0; i<nv; ++i) {
    if (events[i].ubiquity == 1) {
      n = (signed) events[i].entwinement.size();
      for(j=n; j<nt; ++j) {
        events[i].entwinement.push_back(0.0);
      }
      continue;
    }
    if (!events[i].topology_modified) continue;
    events[i].entwinement.clear();
    for(j=0; j<nt; ++j) {
      if (NTL::divide(events[i].ubiquity,codex[j].colour) == 0) {
        events[i].entwinement.push_back(0.0);
        continue;
      }
      //compute_graph(G,i,j);
      //sum = seqn_weights[0]*G->completeness();
      //compute_graph(G,i,j);
      //sum += seqn_weights[1]*G->entwinement()/double(G->nvertex - 1);
      sum = 0.5*double(vertex_dimension(i,j) - 1); // *G->entwinement();
      events[i].entwinement.push_back(sum);
    }
    events[i].global_dimension = vertex_dimension(i,-1);
    events[i].topology_modified = false;
  }

  for(i=0; i<nv; ++i) {
    sum = 0.0;
    for(j=0; j<nt; ++j) {
      sum += events[i].entwinement[j];
    }
    gvalue[i] = sum/double(nt);
  }

#ifdef PARALLEL
  #pragma omp parallel for default(shared) private(i,j,k,l,found,d1,d2) schedule(dynamic,1)
#endif
  for(i=0; i<nv; ++i) {
    if (events[i].ubiquity == 1) continue;
    for(j=1+i; j<nv; ++j) {
      if (events[j].ubiquity == 1) continue;
      l = geometry->get_distance(i,j,false);
      if (l < 3.8025 || l > 4.2025) continue;
      // See if there is a third vertex that lies between these two...
      found = false;
      for(k=0; k<nv; ++k) {
        if (events[k].ubiquity == 1) continue;
        if (k == i || k == j) continue;
        d1 = geometry->get_distance(i,k,false);
        d2 = geometry->get_distance(j,k,false);
        if ((d1 > 0.81 && d1 < 1.21) && (d2 > 0.81 && d2 < 1.21)) {
          found = true;
          break;
        }
      }
      if (found) continue;
      gvalue[i] += 0.1;
#ifdef PARALLEL
      #pragma omp critical
#endif
      gvalue[j] += 0.1;
    }
  }

  for(i=0; i<nv; ++i) {
    if (events[i].ubiquity == 1) continue;
    if (events[i].neighbours.empty()) {
      R[i] = gvalue[i];
      rho[i] = events[i].energy;
      continue;
    }
    sum1 = 0.0;
    sum2 = 0.0;
    k = 0;
    Lambda[i] = 0.0;
    for(it=events[i].neighbours.begin(); it!=events[i].neighbours.end(); it++) {
      j = *it;
      l = geometry->get_distance(i,j,true);
      l_inv = 1.0/(1.0 + l);
      sum1 += gvalue[j]*l_inv;
      sum2 += events[j].energy*l_inv;
      if (l < 1.0) {
        // To handle the case of null edges...
        Lambda[i] += 10.0*(l - 1.0)*(l - 1.0);
      }
      else {
        Lambda[i] += std::log(l)*std::log(l);
      }
      k++;
      /*
      if (!(simplices[1][n].sq_volume < 0.0)) continue;
      u = simplices[1][n].ubiquity.count();
      // How to determine if j lies in the chronological future of i?
      sigma = (simplices[1][n].orientation == FUTURE) ? 1 : -1;
      if (j < i) sigma = -sigma;
      sum2 += double(sigma)*events[j].energy*l_inv;
      */
    }
    Lambda[i] = Lambda[i]/double(k);
    R[i] = gvalue[i]; // - sum1/double(k);
    rho[i] = events[i].energy; // + sum2/double(k);
  }

  // The local part of the structure equations
  for(i=0; i<nv; ++i) {
    if (events[i].ubiquity == 1) continue;
    events[i].deficiency = pfactor*(R[i] + seqn_weights[3]*events[i].obliquity + seqn_weights[2]*Lambda[i] + events[i].curvature) - lambda*rho[i];
    events[i].geometric_deficiency = pfactor*(seqn_weights[3]*events[i].obliquity + seqn_weights[2]*Lambda[i] + events[i].curvature);
  }

  // Now the chromatic energy sum...
  for(i=0; i<nv; ++i) {
    c = 0;
    for(j=0; j<nt; ++j) {
      if (NTL::divide(events[i].ubiquity,codex[j].colour) == 1) c++;
    }
    E_total += double(c)*events[i].energy;
  }
  E_total = E_total/double(nt);

  E_G = representational_energy(false);

  global_deficiency =  E_G + 2.0*M_PI*double(euler_characteristic(-1)) - E_total;

  error = 0.0;
  for(i=0; i<nv; ++i) {
    if (events[i].ubiquity == 1) continue;
    delta = events[i].deficiency;
    error += delta*delta;
  }
  error = std::sqrt(error)/na;

  // Sanity check...
  int nv_test = 0;
  for(i=0; i<nv; ++i) {
    if (events[i].ubiquity == 1) continue;
    if (std::abs(events[i].deficiency) < Spacetime::epsilon) continue;
    nv_test++;
  }
  if (nv_test == 0) assert(error < Spacetime::epsilon);

  //error += std::abs(global_deficiency)/double(na);
  delete G;
}

int Spacetime::sheet_fission(int parent)
{
  int i,j,n,l,offspring = RND.irandom(1,4);
  const int nt = (signed) codex.size();
  const int nv = (signed) events.size();

#ifdef VERBOSE
  std::cout << "Sheet " << parent << " spawning " << offspring << " daughter sheet(s)..." << std::endl;
#endif

  for(l=0; l<offspring; ++l) {
    codex.push_back(Sheet(nt+l,parent,psequence.next()));
    for(i=0; i<nv; ++i) {
      if (NTL::divide(events[i].ubiquity,codex[parent].colour) == 1) events[i].ubiquity *= codex[nt+l].colour;
    }
    for(i=1; i<=ND; ++i) {
      n = (signed) simplices[i].size();
      for(j=0; j<n; ++j) {
        if (NTL::divide(simplices[i][j].ubiquity,codex[parent].colour) == 1) simplices[i][j].ubiquity *= codex[nt+l].colour;
      }
    }
  }
  return offspring;
}

void Spacetime::sheet_dynamics()
{
  // How best to handle the issue of spinning of new sheets from existing ones according to a
  // Poisson process?
  int i,parent,nspawn = 0;
  const int nt = (signed) codex.size();
  std::set<int> candidates;
  std::set<int>::const_iterator it;

  nactive = 0;
  for(i=0; i<nt; ++i) {
    if (RND.drandom() < 0.15) codex[i].active = !codex[i].active;
    if (codex[i].active) {
      nactive++;
      candidates.insert(i);
    }
  }
  for(it=candidates.begin(); it!=candidates.end(); it++) {
    parent = *it;
    if (RND.poisson_variate()) nspawn += sheet_fission(parent);
  }
  nactive += nspawn;
}

void Spacetime::compute_global_topology(int sheet)
{
  // To calculate the global deficiency, we need to compute the Betti numbers and
  // the fundamental group, for the total spacetime, operations that are serial...
  Nexus* NX = new Nexus;
  compute_global_nexus(NX,sheet);
  if (sheet == -1) {
    // The global case...
    if (homology_method == GAP) {
      NX->compute_homology(HZ,homology_base);
    }
    else {
      // Native method...
      NX->compute_homology_native(HZ,homology_base);
    }
    // And the fundamental group...
    NX->compute_homotopy(&pi1);
    // Finally, the pseudomanifold and orientability properties
    pseudomanifold = NX->pseudomanifold(&boundary);
    if (pseudomanifold) orientable = NX->orientable();
  }
  else {
    bool bdry;
    if (homology_method == GAP) {
      NX->compute_homology(codex[sheet].HZ,homology_base);
    }
    else {
      NX->compute_homology_native(codex[sheet].HZ,homology_base);
    }
    // And the fundamental group...
    NX->compute_homotopy(&(codex[sheet].pi1));
    codex[sheet].pseudomanifold = NX->pseudomanifold(&bdry);
    codex[sheet].boundary = bdry;
    if (codex[sheet].pseudomanifold) codex[sheet].orientable = NX->orientable();
  }
  delete NX;
}

void Spacetime::write(Spacetime& state) const
{
  state.system_size = system_size;
  state.error = error;
  state.global_deficiency = global_deficiency;
  state.events = events;
  for(int i=0; i<=ND; ++i) {
    state.simplices[i] = simplices[i];
    state.index_table[i] = index_table[i];
  }
}

void Spacetime::read(const Spacetime& source)
{
  system_size = source.system_size;
  error = source.error;
  global_deficiency = source.global_deficiency;
  events = source.events;
  for(int i=0; i<=ND; ++i) {
    simplices[i] = source.simplices[i];
    index_table[i] = source.index_table[i];
  }
}

bool Spacetime::global_operations()
{
  int i,j,n,k = 0;
  bool output = false;
  NTL::ZZ chi;
  std::string filename;
  std::set<int> vmodified;
  hash_map::const_iterator qt;
  std::vector<int> fusion;
  double mu,sigma = 0.0,delta = 0.0;
  const int nd = dimension(-1);
  const int nv = (signed) events.size();
  const int ne = (signed) simplices[1].size();
#ifdef VERBOSE
  std::cout << "Main thread in global operations method..." << std::endl;
#endif
  iterations++;

  compute_delta();

  for(i=0; i<nv; ++i) {
    if (events[i].incept == -1) events[i].incept = iterations;
  }
  for(i=0; i<=ND; ++i) {
    n = (signed) simplices[i].size();
    for(j=0; j<n; ++j) {
      if (simplices[i][j].incept == -1) simplices[i][j].incept = iterations;
    }
  }

  for(i=0; i<nv; ++i) {
    if (events[i].ubiquity == 1) continue;
    if (events[i].topology_modified) events[i].global_dimension = vertex_dimension(i,-1);
  }

  // First perform the geometry and energy diffusion...
  if (adjust_dimension()) {
    compute_volume();
    compute_curvature();
    compute_obliquity();
    compute_global_topology(-1);
    for(i=0; i<(signed) codex.size(); ++i) {
      compute_global_topology(i);
    }
    structural_deficiency();
  }

  energy_diffusion();

#ifdef VERBOSE
  // Analyze the uniformity of the energy distribution...
  double E_avg = 0.0;
  double nc = double(cardinality(0,-1));
  for(i=0; i<nv; ++i) {
    if (events[i].ubiquity == 1) continue;
    E_avg += events[i].energy;
  }
  E_avg = E_avg/nc;
  sigma = 0.0;
  for(i=0; i<nv; ++i) {
    if (events[i].ubiquity == 1) continue;
    sigma += (events[i].energy - E_avg)*(events[i].energy - E_avg);
  }
  sigma = std::sqrt(sigma/nc);
  std::cout << "Standard deviation of vertex energy is " << sigma << std::endl;

  // Analyze the distribution of vertex dimensionalities...
  int histogram[ND+1],histo2[1+ND];
  for(i=0; i<=ND; ++i) {
    histogram[i] = 0;
    histo2[i] = 0;
  }
  for(i=0; i<nv; ++i) {
    if (events[i].ubiquity == 1) continue;
    histo2[vertex_dimension(i,-1)] += 1;
    if (events[i].energy < Spacetime::epsilon) continue;
    histogram[vertex_dimension(i,-1)] += 1;
  }
  for(i=0; i<=ND; ++i) {
    std::cout << "There are " << histo2[i] << " (" << histogram[i] << ") active " << i << "-dimensional vertices." << std::endl;
  }
#endif

  optimize();

#ifdef VERBOSE
  int ninitial = 0,ntouch = 0;
  for(i=0; i<nv; ++i) {
    if (events[i].incept == 0) {
      ninitial++;
      if (std::abs(events[i].deficiency) > 0.0) ntouch++;
    }
  }
  std::cout << "Percentage of perturbed initial vertices " << 100.0*double(ntouch)/double(ninitial) << std::endl;
#endif

  // Eliminate any overlapping vertices
  if (superposable) {
    superposition_fusion(vmodified);
    superposition_fission(vmodified);
  }

  if (compressible) {
    // Eliminate excessively long edges...
    // First calculate the average edge length and its variance...
    for(i=0; i<ne; ++i) {
      if (simplices[1][i].ubiquity == 1) continue;
      delta += std::abs(simplices[1][i].volume);
      k++;
    }
    mu = delta/double(k);
    for(i=0; i<ne; ++i) {
      if (simplices[1][i].ubiquity == 1) continue;
      delta = (std::abs(simplices[1][i].volume) - mu);
      sigma += delta*delta;
    }
    sigma = std::sqrt(sigma/double(k));
    // Get rid of edges that are more than one standard deviation from the mean...
    i = compression(mu+sigma,vmodified);
  }

  if (permutable) {
    // Take care of any inter-cosmic jumping...
    for(i=1; i<nd; ++i) {
      n = (signed) simplices[i].size();
      for(j=0; j<n; ++j) {
        if (simplices[i][j].ubiquity == 1) continue;
        compute_simplex_energy(i,j);
      }
    }
    delta = Spacetime::T_zero/std::sqrt(Spacetime::kappa*double(1+iterations));
    n = ubiquity_permutation(delta,vmodified);
#ifdef VERBOSE
    std::cout << "The inter-cosmic jump for this iteration is " << n << std::endl;
#endif
  }

  if (superposable || compressible || permutable) {
    compute_entourages(-1);
    compute_topological_dependency(vmodified);
    compute_geometric_dependency(vmodified);
  }

  // Calculate the local and global errors
  geometry->compute_distances(vmodified);
  compute_volume();
  compute_curvature();
  compute_obliquity();
  assert(consistent(-1));
  compute_global_topology(-1);
  for(i=0; i<(signed) codex.size(); ++i) {
    compute_global_topology(i);
  }

  structural_deficiency();

  analyze_convergence();

  if (error < Spacetime::epsilon || iterations >= max_iter) converged = true;
#ifdef VERBOSE
  std::cout << iterations << "  " << error << "  " << Spacetime::epsilon << "  " << converged << std::endl;
#endif

  if (converged) {
    output = true;
    if (!diskless) {
      filename = "data/incastrature_" + date_string + "_" + pid_string + ".dot";
      write_incastrature(filename,-1);
    }
    nactive = (signed) codex.size();
    for(i=0; i<nactive; ++i) {
      codex[i].active = true;
    }
  }
  else {
    if (foliodynamics) sheet_dynamics();
  }

  RND.increment_seed();

  // Write the current state to disk if necessary...
  if (checkpoint) {
    if (iterations >= 1 && (iterations % checkpoint_frequency == 0)) write_state();
  }

#ifdef VERBOSE
  std::cout << "At relaxation step " << iterations << " the global error is " << error << std::endl;
  std::cout << "Sheet Activity " << sheet_activity() << std::endl;
  std::cout << "Finished with global operations..." << std::endl;
#endif
  if (!diskless) write_log();
  return output;
}

void Spacetime::analyze_convergence()
{
  int i,j,m;
  hash_map::const_iterator qt;
  int tdelta = 0;
  double edelta = 0.0;
  const int nv = (signed) events.size();
  const int nva = (signed) anterior.events.size();

  for(i=0; i<nva; ++i) {
    if (events[i].ubiquity != anterior.events[i].ubiquity) {
      tdelta++;
    }
  }
  tdelta += (nv - nva);
  for(i=1; i<=ND; ++i) {
    m = (signed) simplices[i].size();
    for(j=0; j<m; ++j) {
      qt = anterior.index_table[i].find(simplices[i][j].key);
      if (qt == anterior.index_table[i].end()) {
        tdelta++;
        continue;
      }
      if (simplices[i][j].ubiquity != anterior.simplices[i][qt->second].ubiquity) {
        tdelta++;
      }
    }
  }
  topology_delta = double(tdelta);

  geometry_delta = geometry_change(geometry,&(anterior.geometry));

  for(i=0; i<nva; ++i) {
    edelta += std::abs(events[i].energy - anterior.events[i].energy);
  }
  for(i=nva; i<nv; ++i) {
    edelta += events[i].energy;
  }
  energy_delta = edelta;

  anterior.events = events;
  anterior.geometry.load(geometry);
  for(i=1; i<=ND; ++i) {
    anterior.simplices[i].clear();
    for(j=0; j<(signed) simplices[i].size(); ++j) {
      anterior.simplices[i].push_back(simplices[i][j]);
    }
  }

  for(i=1; i<=ND; ++i) {
    anterior.index_table[i].clear();
    for(j=0; j<(signed) anterior.simplices[i].size(); ++j) {
      anterior.index_table[i][anterior.simplices[i][j].key] = j;
    }
  }
}

int Spacetime::ubiquity_permutation(double temperature,std::set<int>& vmodified)
{
  int i,j,k,n,m,nd,vx[2],delta,hdistance,jz = 0;
  NTL::ZZ chi,tau;
  std::set<int>::const_iterator it;
  hash_map::const_iterator qt;
  double E,alpha;
  // This parameter must be greater than zero and determines the
  // thermo-energy scale at which intercosmic jumps become common...
  const int nv = (signed) events.size();
  const int nt = (signed) codex.size();
  const int d = dimension(-1);
  const int jlimit = nv/5;

  for(i=0; i<2*nv; ++i) {
    n = RND.irandom(1,d+1);
    m = RND.irandom(0,simplices[d].size());
    chi = simplices[n][m].ubiquity;
    if (chi == 1) continue;
    E = simplices[n][m].energy;
    if (E < Spacetime::epsilon) continue;
    tau = chi;
    alpha = std::exp(-temperature*E);
    for(j=0; j<nt; ++j) {
      if (RND.drandom() > alpha) {
        if (NTL::divide(tau,codex[j].colour) == 1) {
          tau /= codex[j].colour;
        }
        else {
          tau *= codex[j].colour;
        }
      }
    }
    hdistance = 0;
    for(j=0; j<nt; ++j) {
      if (NTL::divide(chi,codex[j].colour) != NTL::divide(tau,codex[j].colour)) hdistance++;
    }
    delta = n*hdistance;
    if (delta > 0) {
      for(it=simplices[n][m].vertices.begin(); it!=simplices[n][m].vertices.end(); it++) {
        vmodified.insert(*it);
      }
    }
    jz += delta;
    simplices[n][m].ubiquity = tau;
    if (jz >= jlimit) break;
  }
  if (jz == 0) return 0;
  // Now regularize...
  for(i=ND; i>1; i--) {
    if (simplices[i].empty()) continue;
    nd = (signed) simplices[i].size();
    for(j=0; j<nd; ++j) {
      if (simplices[i][j].ubiquity == 1) continue;
      chi = simplices[i][j].ubiquity;
      for(k=0; k<1+i; ++k) {
        qt = index_table[i-1].find(simplices[i][j].faces[k]);
        tau = simplices[i-1][qt->second].ubiquity;
        simplices[i-1][qt->second].ubiquity = (tau*chi)/NTL::GCD(tau,chi);
      }
    }
  }
  n = (signed) simplices[1].size();
  for(i=0; i<n; ++i) {
    if (simplices[1][i].ubiquity == 1) continue;
    chi = simplices[1][i].ubiquity;
    simplices[1][i].get_vertices(vx);
    tau = events[vx[0]].ubiquity;
    events[vx[0]].ubiquity = (tau*chi)/NTL::GCD(tau,chi);
    tau = events[vx[1]].ubiquity;
    events[vx[1]].ubiquity = (tau*chi)/NTL::GCD(tau,chi);
  }
  return jz;
}

void Spacetime::build_initial_state(const NTL::ZZ locale)
{
  int i,j,in1;
  Vertex vt;
  Simplex S;
#ifndef LEIBNIZ
  std::vector<double> svalue;

  for(i=0; i<Geometry::background_dimension; ++i) {
    svalue.push_back(0.0);
  }
#endif
  vt.incept = 0;
  vt.ubiquity = locale;
  S.ubiquity = locale;
  S.incept = 0;

  if (initial_state == CARTESIAN) {
    std::vector<std::pair<long,int> > factors;
    factorize(initial_size,factors);
    j = 1;
    for(i=0; i<(signed) factors.size(); ++i) {
      j *= ipow(factors[i].first,factors[i].second/Geometry::background_dimension);
    }
    const int n = j;
    const int nd = ipow(n,Geometry::background_dimension-1);
    const int nm1 = n - 1;
    //const int nperturbed = 10 + int(0.01*RND.irandom(initial_size));
    const int nperturbed = 2;
#ifndef LEIBNIZ
    const double dx = 1.0;
#endif
    int m,k,l,d,rvalue;
    hash_map::const_iterator qt;
    std::vector<int> entourage,v;
    std::vector<int>* arrangement = new std::vector<int>[initial_size];
    std::set<int> N;
    std::set<int>::const_iterator it,jt;

    for(i=0; i<Geometry::background_dimension; ++i) {
      entourage.push_back(0);
    }
    // First the vertices, making sure (if this is a non-relational
    // model) that the vertex coordinates are balanced...
    for(l=0; l<initial_size; ++l) {
      k = nd;
      in1 = l/k;
      entourage[0] = in1;
      rvalue = l - in1*nd;
#ifndef LEIBNIZ
      svalue[0] = dx*(double(in1) - 0.5*double(nm1));
#endif
      for(m=1; m<Geometry::background_dimension; ++m) {
        k /= n;
        in1 = rvalue/k;
        entourage[m] = in1;
        rvalue -= k*in1;
#ifndef LEIBNIZ
        svalue[m] = dx*(double(in1) - 0.5*double(nm1));
#endif
      }
#ifndef LEIBNIZ
      geometry->add_vertex(svalue);
#endif
      for(m=0; m<Geometry::background_dimension; ++m) {
        if (entourage[m] == 0 || entourage[m] == nm1) vt.boundary = true;
      }
      events.push_back(vt);
      vt.boundary = false;
      arrangement[l] = entourage;
    }
#ifdef LEIBNIZ
    geometry->initialize(n,"CARTESIAN");
#endif
    // Now the edges...
    for(i=0; i<initial_size; ++i) {
      v = arrangement[i];
      if (v[0] > 0) {
        in1 = (v[0]-1)*nd;
        k = nd;
        for(j=0; j<Geometry::background_dimension-1; ++j) {
          k /= n;
          in1 += v[j+1]*k;
        }
        if (in1 > i) {
          S.vertices.insert(i);
          S.vertices.insert(in1);
          S.string_assembly();
          index_table[1][S.key] = (signed) simplices[1].size();
          simplices[1].push_back(S);
          S.vertices.clear();
        }
      }
      if (v[0] < nm1) {
        in1 = (v[0]+1)*nd;
        k = nd;
        for(j=0; j<Geometry::background_dimension-1; ++j) {
          k /= n;
          in1 += v[j+1]*k;
        }
        if (in1 > i) {
          S.vertices.insert(i);
          S.vertices.insert(in1);
          S.string_assembly();
          index_table[1][S.key] = (signed) simplices[1].size();
          simplices[1].push_back(S);
          S.vertices.clear();
        }
      }
      for(j=1; j<Geometry::background_dimension; ++j) {
        if (v[j] > 0) {
          v[j] -= 1;
          in1 = v[0]*nd;
          k = nd;
          for(l=0; l<Geometry::background_dimension-1; ++l) {
            k /= n;
            in1 += v[l+1]*k;
          }
          if (in1 > i) {
            S.vertices.insert(i);
            S.vertices.insert(in1);
            S.string_assembly();
            index_table[1][S.key] = (signed) simplices[1].size();
            simplices[1].push_back(S);
            S.vertices.clear();
          }
          v[j] += 1;
        }
        if (v[j] < nm1) {
          v[j] += 1;
          in1 = v[0]*nd;
          k = nd;
          for(l=0; l<Geometry::background_dimension-1; ++l) {
            k /= n;
            in1 += v[l+1]*k;
          }
          if (in1 > i) {
            S.vertices.insert(i);
            S.vertices.insert(in1);
            S.string_assembly();
            index_table[1][S.key] = (signed) simplices[1].size();
            simplices[1].push_back(S);
            S.vertices.clear();
          }
          v[j] -= 1;
        }
      }
    }
    delete[] arrangement;
    v.clear();
    // We need to see about introducing an initial perturbation...
    if (perturb_topology) {
      k = int(double(n)/2.0);
      for(i=0; i<nperturbed; ++i) {
        j = 0;
        for(l=0; l<Geometry::background_dimension; ++l) {
          j += ipow(n,Geometry::background_dimension-1-l)*RND.irandom(k-2,k+2);
        }
        v.push_back(j);
        // Now add some edges among the neighbours of this
        // vertex
        N.insert(v[i]);
        j = RND.irandom(2,5);
        for(l=0; l<j; ++l) {
          N.insert(events.size());
          v.push_back(events.size());
          geometry->add_vertex(v[i]);
          events.push_back(vt);
        }
        d = N.size() - 1;
        S.initialize(N,locale);
        simplices[d].push_back(S);
        index_table[d][S.key] = (signed) simplices[d].size() - 1;
        N.clear();
      }
    }
    if (perturb_geometry) {
      // No changes to the spacetime topology, merely the geometry of some
      // existing vertices are altered...
      if (v.empty()) {  
        for(i=0; i<nperturbed; ++i) {
          j = RND.irandom(initial_size,v);
          v.push_back(j);
        }
      }
      for(i=0; i<(signed) v.size(); ++i) {
        geometry->perturb_vertex(v[i]);
      }
    }
    if (perturb_energy) {
      if (v.empty()) {
        // In this particular case I will simply alter the energy value of a single vertex
        // near the centre of the Cartesian network...
        j = 0;
        k = int(double(n)/2.0);
        for(i=0; i<Geometry::background_dimension; ++i) {
          j += ipow(n,Geometry::background_dimension-1-i)*RND.irandom(k-2,k+2);
        }
        events[j].energy = 1000.0*(0.5 + RND.drandom()/2.0);
      }
      else {
        for(i=0; i<(signed) v.size(); ++i) {
          events[v[i]].energy = 5.0*RND.drandom();
        }
      }
    }
  }
  else if (initial_state == SINGLETON) {
    // An initial spacetime consisting of a single isolated vertex, though with very
    // high energy
#ifdef LEIBNIZ
    geometry->initialize(0,"SINGLETON");
#else
    geometry->add_vertex(svalue);
#endif
    vt.energy = 5000.0*(0.5 + RND.drandom()/2.0);
    events.push_back(vt);
  }
  else if (initial_state == MONOPLEX) {
    // The initial spacetime is a single simplex of dimension initial_dim
    std::set<int> vx;
#ifndef LEIBNIZ
    int ulimit = Geometry::background_dimension;
    if (initial_dim > Geometry::background_dimension) ulimit = initial_dim;
#endif
    if (perturb_energy) vt.energy = 500.0 + (2000.0/double(initial_dim))*RND.drandom();

#ifndef LEIBNIZ
    svalue.clear();
    for(i=0; i<ulimit; ++i) {
      svalue.push_back(0.0);
    }
    geometry->add_vertex(svalue);
#endif
    events.push_back(vt);
    vx.insert(0);
    for(i=1; i<=initial_dim; ++i) {
#ifndef LEIBNIZ
      svalue[i-1] = 1.0;
      geometry->add_vertex(svalue);
      svalue[i-1] = 0.0;
#endif
      events.push_back(vt);
      vx.insert(i);
    }
#ifdef LEIBNIZ
    geometry->initialize(initial_dim,"MONOPLEX");
#endif
    S.vertices = vx;
    S.string_assembly();
    simplices[initial_dim].push_back(S);
    index_table[initial_dim][S.key] = 0;
  }
  else if (initial_state == RANDOM) {
    // We will use the Erdős–Rényi random graph model (the G(n,p) variant) to
    // assemble a random graph, with n = initial_size...
    int level = 2,k = 0,ulimit;
    std::set<int>* N;
    std::vector<Simplex> svector;
    std::vector<Simplex>::const_iterator vit;
    hash_map::const_iterator qt;
    bool found;
    std::set<int> v,vx,current;
    std::set<int>::const_iterator it,chk;

    RND.initialize_bernoulli(edge_probability);

    for(i=0; i<initial_size; ++i) {
      if (RND.drandom() < 0.1) {
        vt.energy = 10.0*RND.drandom();
        k++;
      }
#ifndef LEIBNIZ
      for(j=0; j<Geometry::background_dimension; ++j) {
        svalue[j] = -10.0 + 20.0*RND.drandom();
      }
      geometry->add_vertex(svalue);
#endif
      events.push_back(vt);
    }
#ifdef LEIBNIZ
    geometry->initialize(initial_size,"RANDOM");
#endif
    // A final step to ensure that at least one node has non-zero energy...
    if (k == 0) events[RND.irandom(initial_size)].energy = 20.0*RND.drandom();

    for(i=0; i<initial_size; ++i) {
      for(j=1+i; j<initial_size; ++j) {
        if (RND.bernoulli_variate() == false) continue;
        S.vertices.insert(i);
        S.vertices.insert(j);
        S.string_assembly();
        index_table[1][S.key] = (signed) simplices[1].size();
        simplices[1].push_back(S);
        events[i].neighbours.insert(j);
        events[j].neighbours.insert(i);
        S.vertices.clear();
      }
    }
    // Now we have to compute all the n-simplices (n > 1) that are implied
    // by this random graph...
    do {
      N = new std::set<int>[level];
      ulimit = (signed) simplices[level-1].size();
      for(j=0; j<ulimit; ++j) {
        v = simplices[level-1][j].vertices;
        i = 0;
        for(it=v.begin(); it!=v.end(); it++) {
          N[i] = events[*it].neighbours;
          i++;
        }
        for(it=N[0].begin(); it!=N[0].end(); it++) {
          found = true;
          for(i=1; i<level; ++i) {
            chk = std::find(N[i].begin(),N[i].end(),*it);
            if (chk == N[i].end()) {
              found = false;
              break;
            }
          }
          if (found) vx.insert(*it);
        }
        for(it=vx.begin(); it!=vx.end(); it++) {
          in1 = *it;
          current = v;
          current.insert(in1);
          S.vertices = current;
          S.string_assembly();
          svector.push_back(S);
        }
        vx.clear();
        v.clear();
      }
      delete[] N;
      if (svector.empty()) break;
      for(vit=svector.begin(); vit!=svector.end(); vit++) {
        qt = index_table[level].find(vit->key);
        if (qt == index_table[level].end()) {
          simplices[level].push_back(*vit);
          index_table[level][vit->key] = (signed) simplices[level].size() - 1;
        }
      }
      svector.clear();
#ifdef VERBOSE
      if (!diskless) std::cout << "Finished the implied " << level << "-simplices..." << std::endl;
#endif
      level += 1;
    } while(true);
  }
 
  for(i=0; i<(signed) codex.size(); ++i) {
    if (NTL::divide(locale,codex[i].colour) == 1) regularization(false,i);
  }
  regularization(false,-1);
}

void Spacetime::initialize()
{
  int i,pid;
  std::stringstream s1,s2,s3,s4;
#ifdef VERBOSE
  timeval Z;
  double t1,t2;

  gettimeofday(&Z,NULL);
  t1 = Z.tv_sec + (Z.tv_usec/1000000.0);
#endif
  if (!diskless) {
    if (std::system("mkdir -p data") < 0) {
      std::cerr << "Unable to create data directory." << std::endl;
      std::cerr << "Exiting..." << std::endl;
      std::exit(1);
    }
  }

  // Get the date and the process ID so that we can construct the
  // state_file and log_file names...
  start_time = boost::posix_time::second_clock::local_time();
  boost::gregorian::date d = start_time.date();
  boost::gregorian::date::ymd_type ymd = d.year_month_day();

  s1.width(2);
  s1 << std::setfill('0') << ymd.day.as_number();
  s2.width(2);
  s2 << std::setfill('0') << ymd.month.as_number();
  s3.width(4);
  s3 << ymd.year;
  date_string = s1.str() + s2.str() + s3.str();

  pid = getpid();
  s4.width(5);
  s4 << std::setfill('0') << pid;
  pid_string = s4.str();
  if (diskless) {
    state_file = "";
    log_file = "";
  }
  else {
    state_file += "_" + date_string + "_" + pid_string;
    log_file += "_" + date_string + "_" + pid_string + ".log";
  }

  if (initial_state == DISKFILE) {
    read_state(input_file);
    if (converged) {
      std::cout << "This spacetime complex has already converged!" << std::endl;
    }
    else if (iterations >= max_iter) {
      std::cout << "This spacetime complex has already exceeded its maximum iterations!" << std::endl;
    }
  }
  else {
    NTL::ZZ locale;

    locale = 1;
    nactive = nt_initial;
    for(i=0; i<nt_initial; ++i) {
      codex.push_back(Sheet(i,psequence.next()));
      locale *= codex[i].colour;
    }
    build_initial_state(locale);
    geometry->compute_distances();
    compute_simplicial_dimension();
    adjust_dimension();
    compute_volume();
    compute_curvature();
    compute_obliquity();
    compute_global_topology(-1);
    for(i=0; i<nt_initial; ++i) {
      compute_global_topology(i);
    }
    structural_deficiency();
  }
  RND.initialize_poisson(Spacetime::ramosity);

  assert(energy_check());

  anterior.events = events;
  for(i=1; i<=ND; ++i) {
    anterior.simplices[i] = simplices[i];
    anterior.index_table[i] = index_table[i];
  }
  if (diskless) return;

#ifdef VERBOSE
  if (pseudomanifold) {
    if (orientable) {
      if (boundary) {
        std::cout << "The spacetime complex is an orientable pseudomanifold with boundary." << std::endl;
      }
      else {
        std::cout << "The spacetime complex is an orientable pseudomanifold." << std::endl;
      }
    }
    else {
      if (boundary) {
        std::cout << "The spacetime complex is a non-orientable pseudomanifold with boundary." << std::endl;
      }
      else {
        std::cout << "The spacetime complex is a non-orientable pseudomanifold." << std::endl;
      }
    }
  }
  else {
    std::cout << "The spacetime complex is not a pseudomanifold" << std::endl;
  }
  write_topology(-1);
#endif

  write_state();
  write_log();
#ifdef VERBOSE
  std::cout << "At relaxation step " << iterations << " the global error is " << error << std::endl;
  std::cout << "Sheet activity " << sheet_activity() << std::endl;
#endif

#ifdef VERBOSE
  gettimeofday(&Z,NULL);
  t2 = Z.tv_sec + (Z.tv_usec/1000000.0);
  std::cout << "Spacetime initialization required " << t2 - t1 << " seconds to complete." << std::endl;
#endif
}
