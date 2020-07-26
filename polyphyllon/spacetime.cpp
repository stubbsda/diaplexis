#include "spacetime.h"

using namespace DIAPLEXIS;

const double Spacetime::ramosity = 0.15;
const double Spacetime::convergence_threshold = 0.00001;
const double Spacetime::T_zero = 500.0;
const double Spacetime::kappa = 1.35;
const double Spacetime::Lambda = 0.2;

Spacetime::Spacetime(bool no_disk)
{
  allocate();
  diskless = no_disk;
  if (diskless) checkpoint_frequency = 0;
  initialize();
}

Spacetime::Spacetime(const std::string& filename,bool no_disk)
{
  allocate();
  read_parameters(filename);
  diskless = no_disk;
  if (diskless) checkpoint_frequency = 0;
  initialize();
}

Spacetime::~Spacetime()
{
  codex.clear();

  delete geometry;
  delete skeleton;
}

void Spacetime::allocate()
{
  // Default geometry (Euclidean, absolute, dimensionally 
  // uniform, background dimension = 3)
  geometry = new SYNARMOSMA::Geometry;
  skeleton = new Complex;
}

void Spacetime::restart(const std::string& filename,bool save_seed)
{
  if (save_seed) {
    unsigned long n = skeleton->RND->get_seed();
    clear();
    read_parameters(filename);
    skeleton->RND->set_seed(n);
  }
  else {
    clear();
    read_parameters(filename);
  }
  initialize();
}

void Spacetime::clear()
{
  skeleton->clear();
  geometry->clear();
  iterations = 0;
  codex.clear();
}

double Spacetime::condense()
{
  // First check how many ghost events and edges there are in this spacetime....
  int n = skeleton->cardinality(0,-1),m = skeleton->cardinality(1,-1);
  const int nv = (signed) skeleton->events.size();
  const int ne = (signed) skeleton->simplices[1].size();
  double rho_v = double(n)/double(nv);
  double rho_e = double(m)/double(ne);
  double output = std::max(rho_v,rho_e);
#ifdef VERBOSE
  std::cout << "Topological density is " << rho_v << " for events and " << rho_e << " for edges." << std::endl;
#endif
  if (output > 0.5) return output;
  // So we need to condense this spacetime to reduce memory pressure...
  int i,j,offset[nv];
  std::set<int> N,S,vx,locus;
  std::set<int>::const_iterator it;
  std::vector<Event> nevents;
  std::vector<Simplex> nsimplices;

  n = 0;
  for(i=0; i<nv; ++i) {
    offset[i] = -1;
    if (!skeleton->events[i].active()) continue;
    nevents.push_back(skeleton->events[i]);
    offset[i] = n;
    n++;
  }
  skeleton->events = nevents;
  for(i=0; i<n; ++i) {
    skeleton->events[i].get_neighbours(N);
    for(it=N.begin(); it!=N.end(); ++it) {
      vx.insert(offset[*it]);
    }
    skeleton->events[i].set_neighbours(vx);
    skeleton->events[i].clear_entourage();
    vx.clear();
  }
  for(i=1; i<=Complex::ND; ++i) {
    m = (signed) skeleton->simplices[i].size();
    skeleton->index_table[i].clear();
    for(j=0; j<m; ++j) {
      if (!skeleton->simplices[i][j].active()) continue;
      skeleton->simplices[i][j].get_ubiquity(locus);
      skeleton->simplices[i][j].get_vertices(S);
      for(it=S.begin(); it!=S.end(); ++it) {
        vx.insert(offset[*it]);
      }
      skeleton->index_table[i][vx] = (signed) nsimplices.size();
      nsimplices.push_back(Simplex(vx,locus));
      vx.clear();
    }
    skeleton->simplices[i] = nsimplices;
    nsimplices.clear();
  }
  skeleton->compute_entourages(-1);  
  return 1.0;
}

bool Spacetime::clean() const
{
  unsigned int i;
  for(i=0; i<skeleton->events.size(); ++i) {
    if (skeleton->events[i].get_topology_modified()) return false;
  }
  for(int d=1; d<=Complex::ND; ++d) {
    for(i=0; i<skeleton->simplices[d].size(); ++i) {
      if (skeleton->simplices[d][i].get_modified()) return false;
    }
  }
  return true;
}

void Spacetime::fallback()
{
  if (!reversible) return;
  read_state("data/.previous_step.dat");
  reversible = false;
}

bool Spacetime::advance()
{
  int i,n = (signed) codex.size();
  std::vector<int> order;
  static double htime = 0.0;
  static double gtime = 0.0;
  boost::timer::cpu_times Z;
  boost::timer::cpu_timer t1;

  write_state("data/.previous_step.dat");
  reversible = true;

  // Begin the hyphantic phase...
  t1.start();
  skeleton->RND->shuffle(order,n);

  std::ofstream s1(hyphansis_file,std::ios::app);
  s1 << "<Iteration>" << std::endl;
  s1 << "  <Index>" << iterations << "</Index>" << std::endl;
  s1.close();
  for(i=0; i<n; ++i) {
    if (!codex[order[i]].active) continue;
    hyphansis(order[i]);
  }
  std::ofstream s2(hyphansis_file,std::ios::app);
  s2 << "</Iteration>" << std::endl;
  s2.close();

  regularization(false,-1);
  condense();
  t1.stop();
  Z = t1.elapsed();
  htime += std::stod(boost::timer::format(Z,3,"%w"));

  // Start the global operations phrase...
  t1.start();
#ifdef VERBOSE
  std::cout << "Calling global operations..." << std::endl;
#endif
  bool done = global_operations();
  t1.stop();
  Z = t1.elapsed();
  gtime += std::stod(boost::timer::format(Z,3,"%w"));

  // If the simulation is finished print out some timing data...
  if (done) {
    std::cout << "Time required for topological hyphansis was " << htime << " seconds." << std::endl;
    std::cout << "Time required for global operations was " << gtime << " seconds." << std::endl;
  }
  return done;
}

double Spacetime::evolve()
{
  for(int i=0; i<max_iter; ++i) {
    if (advance()) break;
  }

  if (diskless) return error;

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

  return error;
}

bool Spacetime::correctness()
{
  unsigned int i,j;
  double delta,original_error;
  const unsigned int nv = skeleton->events.size();

  compute_volume();
  compute_obliquity();
  structural_deficiency();
  original_error = error;

  for(i=0; i<nv; ++i) {
    skeleton->events[i].set_topology_modified(true);
    skeleton->events[i].set_geometry_modified(true);
  }
  for(i=1; i<=Complex::ND; ++i) {
    for(j=0; j<skeleton->simplices[i].size(); ++j) {
      skeleton->simplices[i][j].set_modified(true);
    }
  }

  compute_volume();
  compute_obliquity();
  structural_deficiency();
  delta = std::abs(original_error - error);
  if (delta < std::numeric_limits<double>::epsilon()) return true;
#ifdef VERBOSE
  std::cout << "Error difference in Spacetime::correctness is " << original_error << "  " << error << "  " << delta << std::endl;
#endif
  return false;
}

void Spacetime::structural_deficiency()
{
  int i,j,k,v1,v2,v3,c = 0;
  bool found;
  double alpha,sum,sum1,sum2,l,l_inv,d1,d2,delta,E_total = 0.0;
  std::set<int> S;
  std::set<int>::const_iterator it;
  SYNARMOSMA::Graph G;
  SYNARMOSMA::hash_map::const_iterator qt;
  const int na = double(skeleton->cardinality(0,-1));
  const int nv = (signed) skeleton->events.size();
  const int nt = (signed) codex.size();
  double R[nv],gvalue[nv],length_deviation[nv],rho[nv],tangle[nt];
  int avertices[na];

  for(i=0; i<nv; ++i) {
    if (skeleton->events[i].active()) {
      avertices[c] = i;
      c++;
    }
    skeleton->events[i].set_deficiency(0.0);
    skeleton->events[i].set_geometric_deficiency(0.0);
    length_deviation[i] = 0.0;
    gvalue[i] = 0.0;
    R[i] = 0.0;
    rho[i] = 0.0;
  }

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i,j,k,sum,tangle,G)
#endif
  for(i=0; i<nv; ++i) {
    for(j=0; j<nt; ++j) {
      tangle[j] = 0.0;
    }
    if (!skeleton->events[i].active()) {
      skeleton->events[i].set_entwinement(tangle,nt);
      continue;
    }
    if (!skeleton->events[i].get_topology_modified()) continue;
    for(j=0; j<nt; ++j) {
      if (!skeleton->events[i].active(j)) continue;
      tangle[j] = 0.5*double(skeleton->vertex_dimension(i,j) - 1) + compute_temporal_vorticity(i,j);
      skeleton->compute_graph(&G,i,j);
      k = G.order();
      if (k > 1) tangle[j] += G.completeness() + G.median_degree()/double(k - 1);
    }
    skeleton->events[i].set_entwinement(tangle,nt);
    skeleton->events[i].set_topological_dimension(skeleton->vertex_dimension(i,-1));
    skeleton->events[i].set_topology_modified(false);
    sum = 0.0;
    for(j=0; j<nt; ++j) {
      sum += tangle[j];
    }
    gvalue[i] = sum/double(nt);
  }

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i,j,k,v1,v2,v3,l,found,d1,d2) schedule(dynamic,1)
#endif
  for(i=0; i<na; ++i) {
    v1 = avertices[i];
    for(j=1+i; j<na; ++j) {
      v2 = avertices[j];
      l = std::abs(geometry->get_squared_distance(v1,v2,false));
      if (l < 3.8025 || l > 4.2025) continue;
      // See if there is a third vertex that lies between these two...
      found = false;
      for(k=0; k<i; ++k) {
        v3 = avertices[k];
        d1 = std::abs(geometry->get_squared_distance(v1,v3,false));
        d2 = std::abs(geometry->get_squared_distance(v2,v3,false));
        if ((d1 > 0.81 && d1 < 1.21) && (d2 > 0.81 && d2 < 1.21)) {
          found = true;
          break;
        }
      }
      if (found) continue;
      for(k=1+i; k<j; ++k) {
        v3 = avertices[k];
        d1 = std::abs(geometry->get_squared_distance(v1,v3,false));
        d2 = std::abs(geometry->get_squared_distance(v2,v3,false));
        if ((d1 > 0.81 && d1 < 1.21) && (d2 > 0.81 && d2 < 1.21)) {
          found = true;
          break;
        }
      }
      if (found) continue;
      for(k=1+j; k<na; ++k) {
        v3 = avertices[k];
        d1 = std::abs(geometry->get_squared_distance(v1,v3,false));
        d2 = std::abs(geometry->get_squared_distance(v2,v3,false));
        if ((d1 > 0.81 && d1 < 1.21) && (d2 > 0.81 && d2 < 1.21)) {
          found = true;
          break;
        }
      }
      if (found) continue;
#ifdef _OPENMP
#pragma omp critical
      {
#endif
      gvalue[v1] += 0.1;
      gvalue[v2] += 0.1;
#ifdef _OPENMP
      }
#endif
    }
  }

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i,j,k,S,v1,l,l_inv,it,sum1,sum2)
#endif
  for(i=0; i<na; ++i) {
    v1 = avertices[i];
    sum1 = 0.0;
    sum2 = 0.0;
    skeleton->events[v1].get_neighbours(S);
    k = (S.empty()) ? 1 : (signed) S.size();
    length_deviation[v1] = 0.0;
    for(it=S.begin(); it!=S.end(); ++it) {
      j = *it;
      l = std::abs(geometry->get_squared_distance(v1,j,false));
      l_inv = 1.0/(1.0 + l);
      sum1 += gvalue[j]*l_inv;
      sum2 += skeleton->events[j].get_energy()*l_inv;
      if (l < 1.0) {
        // To handle the case of null edges...
        length_deviation[v1] += 10.0*(l - 1.0)*(l - 1.0);
      }
      else {
        length_deviation[v1] += std::log(l)*std::log(l);
      }
    }
    length_deviation[v1] = length_deviation[v1]/double(k);
    R[v1] = gvalue[v1]; 
    rho[v1] = skeleton->events[v1].get_energy(); 
    // Sanity checks...
#ifdef DEBUG
    if (std::isnan(R[v1])) throw std::runtime_error("NaN detected in Spacetime::structural_deficiency for vertex " + std::to_string(v1)); 
    if (std::isnan(rho[v1])) throw std::runtime_error("NaN detected in Spacetime::structural_deficiency for vertex " + std::to_string(v1)); 
    if (std::isnan(length_deviation[v1])) throw std::runtime_error("NaN detected in Spacetime::structural_deficiency for vertex " + std::to_string(v1));
#endif
  }

  // The local part of the structure equations
  for(i=0; i<na; ++i) {
    v1 = avertices[i];
    alpha = skeleton->events[v1].get_obliquity() + length_deviation[v1];
    skeleton->events[v1].set_geometric_deficiency(alpha);
    skeleton->events[v1].set_deficiency(R[v1] + alpha - Spacetime::Lambda*rho[v1]);
  }

  // Now the chromatic energy sum...
  for(i=0; i<nv; ++i) {
    c = 0;
    for(j=0; j<nt; ++j) {
      if (skeleton->events[i].active(j)) c++;
    }
    E_total += double(c)*skeleton->events[i].get_energy();
  }
  E_total = E_total/double(nt);

  global_deficiency = representational_energy(false) + compute_temporal_nonlinearity() + 2.0*M_PI*double(skeleton->euler_characteristic(-1)) - E_total;

  error = 0.0;
  for(i=0; i<na; ++i) {
    delta = skeleton->events[avertices[i]].get_deficiency();
    error += delta*delta;
  }
  error = (std::sqrt(error) + std::abs(global_deficiency))/double(na);
#ifdef VERBOSE
  double total_error = 0.0;
  for(i=0; i<na; ++i) {
    total_error += std::abs(skeleton->events[avertices[i]].get_deficiency());
  }
  std::cout << "The total error is " << total_error << std::endl;
#endif

  // Sanity check...
#ifdef DEBUG
  int nv_test = 0;
  for(i=0; i<na; ++i) {
    if (std::abs(skeleton->events[avertices[i]].get_deficiency()) < std::numeric_limits<double>::epsilon()) continue;
    nv_test++;
  }
  if (nv_test == 0) assert(error < std::numeric_limits<double>::epsilon());
#endif
}

bool Spacetime::global_operations()
{
  int i,j,n,k = 0;
  bool output = false;
  double mu,sigma = 0.0,delta = 0.0;
  std::string filename;
  SYNARMOSMA::hash_map::const_iterator qt;
  const int nd = skeleton->dimension(-1);
  const int nv = (signed) skeleton->events.size();
  const int ne = (signed) skeleton->simplices[1].size();
#ifdef VERBOSE
  std::cout << "Main thread in global operations method with " << nv << " events." << std::endl;
#endif
  iterations++;

  skeleton->compute_modified_events();

  for(i=0; i<nv; ++i) {
    if (skeleton->events[i].get_incept() == -1) skeleton->events[i].set_incept(iterations);
  }
  for(i=0; i<=Complex::ND; ++i) {
    n = (signed) skeleton->simplices[i].size();
    for(j=0; j<n; ++j) {
      if (skeleton->simplices[i][j].get_incept() == -1) skeleton->simplices[i][j].set_incept(iterations);
    }
  }

  for(i=0; i<nv; ++i) {
    if (!skeleton->events[i].active()) continue;
    if (skeleton->events[i].get_topology_modified()) skeleton->events[i].set_topological_dimension(skeleton->vertex_dimension(i,-1));
  }

  // First perform the geometry and energy diffusion...
  skeleton->compute_parity();
  compute_lightcones();
  if (adjust_dimension()) {
    compute_volume();
    compute_obliquity();
    compute_global_topology(-1);
    for(i=0; i<(signed) codex.size(); ++i) {
      compute_global_topology(i);
    }
    structural_deficiency();
  }

  skeleton->energy_diffusion(Spacetime::Lambda);
  for(i=1; i<skeleton->dimension(-1); ++i) {
    for(j=0; j<(signed) skeleton->simplices[i].size(); ++j) {
      skeleton->compute_simplex_energy(i,j);
    }
  }

#ifdef VERBOSE
  // Analyze the uniformity of the energy distribution...
  double E_avg = 0.0;
  double nc = double(skeleton->cardinality(0,-1));
  for(i=0; i<nv; ++i) {
    if (!skeleton->events[i].active()) continue;
    E_avg += skeleton->events[i].get_energy();
  }
  E_avg = E_avg/nc;
  sigma = 0.0;
  for(i=0; i<nv; ++i) {
    if (!skeleton->events[i].active()) continue;
    sigma += (skeleton->events[i].get_energy() - E_avg)*(skeleton->events[i].get_energy() - E_avg);
  }
  sigma = std::sqrt(sigma/nc);
  std::cout << "Standard deviation of vertex energy is " << sigma << std::endl;

  // Analyze the distribution of vertex dimensionalities...
  int histogram[1 + Complex::ND],histo2[1 + Complex::ND];
  for(i=0; i<=Complex::ND; ++i) {
    histogram[i] = 0;
    histo2[i] = 0;
  }
  for(i=0; i<nv; ++i) {
    if (!skeleton->events[i].active()) continue;
    histo2[skeleton->vertex_dimension(i,-1)] += 1;
    if (skeleton->events[i].zero_energy()) continue;
    histogram[skeleton->vertex_dimension(i,-1)] += 1;
  }
  for(i=0; i<=Complex::ND; ++i) {
    std::cout << "There are " << histo2[i] << " (" << histogram[i] << ") active " << i << "-dimensional vertices." << std::endl;
  }
#endif

  optimize();

#ifdef VERBOSE
  int ninitial = 0,ntouch = 0;
  for(i=0; i<nv; ++i) {
    if (skeleton->events[i].get_incept() == 0) {
      ninitial++;
      if (std::abs(skeleton->events[i].get_deficiency()) > std::numeric_limits<double>::epsilon()) ntouch++;
    }
  }
  std::cout << "Percentage of perturbed initial events " << 100.0*double(ntouch)/double(ninitial) << std::endl;
#endif

  // Eliminate any overlapping vertices
  if (superposable) {
    superposition_fusion(0.1);
    superposition_fission(int(0.02*nv));
  }

  if (compressible) {
    // Eliminate excessively long edges...
    // First calculate the average edge length and its variance...
    for(i=0; i<ne; ++i) {
      if (!skeleton->simplices[1][i].active()) continue;
      delta += std::abs(skeleton->simplices[1][i].get_volume());
      k++;
    }
    mu = delta/double(k);
    for(i=0; i<ne; ++i) {
      if (!skeleton->simplices[1][i].active()) continue;
      delta = (std::abs(skeleton->simplices[1][i].get_volume()) - mu);
      sigma += delta*delta;
    }
    sigma = std::sqrt(sigma/double(k));
    // Get rid of edges that are more than one standard deviation from the mean...
    i = compression(mu + sigma);
  }

  if (permutable && nactive > 1) {
    // Take care of any inter-cosmic jumping...
    for(i=1; i<nd; ++i) {
      n = (signed) skeleton->simplices[i].size();
      for(j=0; j<n; ++j) {
        if (!skeleton->simplices[i][j].active()) continue;
        skeleton->compute_simplex_energy(i,j);
      }
    }
    delta = Spacetime::T_zero/std::sqrt(Spacetime::kappa*double(1+iterations));
    n = ubiquity_permutation(delta);
#ifdef VERBOSE
    std::cout << "The inter-cosmic jump for this iteration is " << n << std::endl;
#endif
  }

  if (superposable || compressible || permutable) {
    skeleton->compute_entourages(-1);
    skeleton->compute_neighbours();
    skeleton->compute_modified_events();
    n = (signed) skeleton->events.size();
    std::set<int> S;
    for(i=0; i<n; ++i) {
      if (!skeleton->events[i].active()) continue;
      if (skeleton->events[i].get_topology_modified()) S.insert(i);
    }
    geometry->compute_squared_distances(S); 
  }
  compute_volume();
  compute_obliquity();
#ifdef DEBUG
  assert(skeleton->consistent(-1));
#endif

  compute_lightcones();
  compute_global_topology(-1);
  for(i=0; i<(signed) codex.size(); ++i) {
    compute_global_topology(i);
  }

  structural_deficiency();

  if (error < Spacetime::convergence_threshold || iterations >= max_iter) converged = true;
#ifdef VERBOSE
  std::cout << iterations << "  " << error << "  " << Spacetime::convergence_threshold << "  " << converged << std::endl;
#endif

  if (converged) {
    output = true;
    if (!diskless) {
      filename = "data/incastrature_" + date_string + "_" + pid_string + ".dot";
      skeleton->write_incastrature(filename,-1);
    }
    nactive = (signed) codex.size();
    for(i=0; i<nactive; ++i) {
      codex[i].active = true;
    }
  }
  else {
    if (foliodynamics) {
      if (sheet_dynamics() > 0) {
        // If new sheets were created then we need to ensure that every event gets an entanglement 
        // property of the right length (=codex.size()), which is what the following code achieves.
        for(i=0; i<nv; ++i) {
          if (skeleton->events[i].active()) skeleton->events[i].set_topology_modified(true);
        }
        structural_deficiency();
      }
    }
  }

  skeleton->RND->increment_seed();

  // Write the current state to disk if necessary...
  if (checkpoint_frequency > 0) {
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

int Spacetime::sheet_fission(int parent)
{
  int i,j,n,l,offspring = skeleton->RND->irandom(1,4);
  const int nt = (signed) codex.size();
  const int nv = (signed) skeleton->events.size();

#ifdef VERBOSE
  std::cout << "Sheet " << parent << " spawning " << offspring << " daughter sheet(s)..." << std::endl;
#endif

  for(l=0; l<offspring; ++l) {
    codex.push_back(Sheet(nt+l,parent,skeleton->H->get_field(),skeleton->H->get_method()));
    for(i=0; i<nv; ++i) {
      if (skeleton->events[i].active(parent)) {
        skeleton->events[i].activate(nt+l); 
        skeleton->events[i].set_topology_modified(true);
      }
    }
    for(i=1; i<=Complex::ND; ++i) {
      n = (signed) skeleton->simplices[i].size();
      for(j=0; j<n; ++j) {
        if (skeleton->simplices[i][j].active(parent)) skeleton->simplices[i][j].activate(nt+l);
      }
    }
  }
  return offspring;
}

int Spacetime::sheet_dynamics()
{
  // How best to handle the issue of spinning of new sheets from existing ones according to a
  // Poisson process?
  int i,parent,nspawn = 0;
  const int nt = (signed) codex.size();
  std::set<int> candidates;
  std::set<int>::const_iterator it;

  nactive = 0;
  for(i=0; i<nt; ++i) {
    if (skeleton->RND->drandom() < 0.15) codex[i].active = !codex[i].active;
    if (codex[i].active) {
      nactive++;
      candidates.insert(i);
    }
  }
  for(it=candidates.begin(); it!=candidates.end(); ++it) {
    parent = *it;
    if (skeleton->RND->poisson_variate()) nspawn += sheet_fission(parent);
  }
  nactive += nspawn;
  return nspawn;
}

void Spacetime::compute_global_topology(int sheet)
{
  // To calculate the global deficiency, we need to compute the Betti numbers and
  // the fundamental group, for the total spacetime, operations that are serial...
  if (skeleton->cardinality(0,sheet) < 2) return;

  if (sheet == -1) {
    // The global case...
    skeleton->compute_global_topology(high_memory);
  }
  else {
    SYNARMOSMA::Nexus* NX = new SYNARMOSMA::Nexus;
    bool bdry;

    skeleton->compute_global_nexus(NX,sheet);
    codex[sheet].H->compute(NX);
    if (high_memory) codex[sheet].pi1->compute(NX);
    codex[sheet].pseudomanifold = NX->pseudomanifold(&bdry);
    codex[sheet].boundary = bdry;
    if (codex[sheet].pseudomanifold) codex[sheet].orientable = NX->orientable();
    delete NX;
  }
}

void Spacetime::get_ubiquity(int d,int n,std::string& output) const
{
  const int nt = (signed) codex.size();

  if (d == 0) {
    output = "{";
    for(int i=0; i<nt-1; ++i) {
      if (skeleton->events[n].active(i)) {
        output += "1,";
      }
      else {
        output += "0,";
      }
    }
    if (skeleton->events[n].active(nt-1)) {
      output += "1}";
    }
    else {
      output += "0}";
    }
  }
  else {
    for(int i=0; i<nt-1; ++i) {
      if (skeleton->simplices[d][n].active(i)) {
        output += "1,";
      }
      else {
        output += "0,";
      }
    }
    if (skeleton->simplices[d][n].active(nt-1)) {
      output += "1}";
    }
    else {
      output += "0}";
    }
  }
}

void Spacetime::get_ubiquity_vector(std::vector<int>& output) const
{
  int i,j,n;
  const int nv = (signed) skeleton->events.size();
  const int ne = (signed) skeleton->simplices[1].size();
  const int nt = (signed) codex.size();

  output.clear();

  for(i=0; i<nv; ++i) {
    if (!skeleton->events[i].active()) continue;
    for(j=0; j<nt; ++j) {
      n = skeleton->events[i].active(j) ? 1 : 0;
      output.push_back(n);
    }
  }
  for(i=0; i<ne; ++i) {
    if (!skeleton->simplices[1][i].active()) continue;
    for(j=0; j<nt; ++j) {
      n = skeleton->simplices[1][i].active(j) ? 1 : 0;
      output.push_back(n);
    }
  }
}

int Spacetime::ubiquity_permutation(double temperature)
{
  int i,j,k,n,m,nd,vx[2],delta,hdistance,jz = 0;
  std::set<int> S,ubiquity,mutated;
  std::set<int>::const_iterator it;
  std::vector<std::set<int> > F;
  SYNARMOSMA::hash_map::const_iterator qt;
  double E,alpha;
  // This parameter must be greater than zero and determines the
  // thermo-energy scale at which intercosmic jumps become common...
  const int nv = (signed) skeleton->events.size();
  const int nt = (signed) codex.size();
  const int d = skeleton->dimension(-1);
  const int jlimit = nv/5;

  for(i=0; i<2*nv; ++i) {
    n = skeleton->RND->irandom(1,d+1);
    m = skeleton->RND->irandom(0,skeleton->simplices[d].size());
    if (!skeleton->simplices[n][m].active()) continue;
    E = skeleton->simplices[n][m].get_energy();
    if (E < std::numeric_limits<double>::epsilon()) continue;
    skeleton->simplices[n][m].get_ubiquity(ubiquity);
    mutated = ubiquity;
    alpha = std::exp(-temperature*E);
    hdistance = 0;
    for(j=0; j<nt; ++j) {
      if (!codex[j].active) continue;
      if (skeleton->RND->drandom() > alpha) {
        if (ubiquity.count(j) > 0) {
          mutated.erase(j);
        }
        else {
          mutated.insert(j);
        }
        hdistance++;
      }
    }
    delta = n*hdistance;
    if (delta > 0) {
      skeleton->simplices[n][m].get_vertices(S);
      for(it=S.begin(); it!=S.end(); ++it) {
        skeleton->events[*it].set_topology_modified(true);
      }
      skeleton->simplices[n][m].set_ubiquity(mutated);
    }
    jz += delta;
    if (jz >= jlimit) break;
  }
  if (jz == 0) return 0;
  // Now regularize...
  for(i=Complex::ND; i>1; i--) {
    if (skeleton->simplices[i].empty()) continue;
    nd = (signed) skeleton->simplices[i].size();
    for(j=0; j<nd; ++j) {
      if (!skeleton->simplices[i][j].active()) continue;
      skeleton->simplices[i][j].get_ubiquity(ubiquity);
      skeleton->simplices[i][j].get_faces(F);
      for(k=0; k<1+i; ++k) {
        qt = skeleton->index_table[i-1].find(F[k]);
        for(it=ubiquity.begin(); it!=ubiquity.end(); ++it) {
          skeleton->simplices[i-1][qt->second].activate(*it);
        }
      }
    }
  }
  n = (signed) skeleton->simplices[1].size();
  for(i=0; i<n; ++i) {
    if (!skeleton->simplices[1][i].active()) continue;
    skeleton->simplices[1][i].get_ubiquity(ubiquity);
    skeleton->simplices[1][i].get_vertices(vx);
    for(it=ubiquity.begin(); it!=ubiquity.end(); ++it) {
      skeleton->events[vx[0]].activate(*it);
      skeleton->events[vx[1]].activate(*it);
    }
  }
  return jz;
}
