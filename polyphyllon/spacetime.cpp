#include "spacetime.h"

using namespace DIAPLEXIS;

const double Spacetime::ramosity = 0.15;
const double Spacetime::convergence_threshold = 0.00001;
const double Spacetime::T_zero = 500.0;
const double Spacetime::kappa = 1.35;
const double Spacetime::Lambda = 0.2;

Spacetime::Spacetime()
{
  allocate();
  initialize();
}

Spacetime::Spacetime(bool no_disk)
{
  allocate();
  diskless = no_disk;
  if (diskless) checkpoint_frequency = 0;
  initialize();
}

Spacetime::Spacetime(const std::string& filename)
{
  allocate();
  read_parameters(filename);
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
  codex.clear();
  skeleton->clear();
  iterations = 0;
}

void Spacetime::condense()
{
  // First check how many ghost vertices and edges there are in this spacetime....
  int i,n = 0,m = 0;
  const int nv = (signed) skeleton->events.size();
  const int ne = (signed) skeleton->simplices[1].size();
  for(i=0; i<nv; ++i) {
    if (skeleton->events[i].active()) n++;
  }
  for(i=0; i<ne; ++i) {
    if (skeleton->simplices[1][i].active()) m++;
  }
  double rho_v = double(n)/double(nv);
  double rho_e = double(m)/double(ne);
#ifdef VERBOSE
  std::cout << "Topological density is for vertices " << rho_v << " and for edges " << rho_e << std::endl;
#endif
  if (rho_v > 0.5 || rho_e > 0.5) return;
  // So we need to condense this spacetime to reduce memory pressure...
  int j,offset[nv];
  Simplex S;
  std::set<int> N,vx;
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
    skeleton->events[i].entourage.clear();
    vx.clear();
  }
  for(i=1; i<=Complex::ND; ++i) {
    m = (signed) skeleton->simplices[i].size();
    skeleton->index_table[i].clear();
    for(j=0; j<m; ++j) {
      if (!skeleton->simplices[i][j].active()) continue;
      S = skeleton->simplices[i][j];
      for(it=S.vertices.begin(); it!=S.vertices.end(); ++it) {
        vx.insert(offset[*it]);
      }
      S.vertices = vx;
      S.entourage.clear();
      S.calculate_faces();
      skeleton->index_table[i][vx] = (signed) nsimplices.size();
      nsimplices.push_back(S);
      vx.clear();
    }
    skeleton->simplices[i] = nsimplices;
    nsimplices.clear();
  }
  skeleton->compute_entourages(-1);  
}

void Spacetime::clean() const
{
  int i,j;
  for(i=0; i<(signed) skeleton->events.size(); ++i) {
    assert(!skeleton->events[i].topology_modified);
  }
  for(i=1; i<=Complex::ND; ++i) {
    for(j=0; j<(signed) skeleton->simplices[i].size(); ++j) {
      assert(!skeleton->simplices[i][j].modified);
    }
  }
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
  RND->shuffle(order,n);

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
  htime += boost::lexical_cast<double>(boost::timer::format(Z,3,"%w"));

  // Start the global operations phrase...
  t1.start();
#ifdef VERBOSE
  std::cout << "Calling global operations..." << std::endl;
#endif
  bool done = global_operations();
  t1.stop();
  Z = t1.elapsed();
  gtime += boost::lexical_cast<double>(boost::timer::format(Z,3,"%w"));

  // If the simulation is finished print out some timing data...
  if (done) {
    std::cout << "Time required for topological hyphansis was " << htime << " seconds." << std::endl;
    std::cout << "Time required for global operations was " << gtime << " seconds." << std::endl;
  }
  return done;
}

void Spacetime::evolve()
{
  for(int i=0; i<max_iter; ++i) {
    if (step_forwards()) break;
  }

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
  const int nv = (signed) skeleton->events.size();

  compute_volume();
  compute_curvature();
  compute_obliquity();
  structural_deficiency();
  double anterior_error = error;

  for(i=0; i<nv; ++i) {
    skeleton->events[i].topology_modified = true;
    skeleton->events[i].geometry_modified = true;
  }
  for(i=1; i<=Complex::ND; ++i) {
    for(j=0; j<(signed) skeleton->simplices[i].size(); ++j) {
      skeleton->simplices[i][j].modified = true;
    }
  }

  compute_volume();
  compute_curvature();
  compute_obliquity();
  structural_deficiency();
  if (std::abs(anterior_error - error) < std::numeric_limits<double>::epsilon()) return true;
  std::cout << "Error difference is " << anterior_error << "  " << error << "  " << std::abs(anterior_error - error) << std::endl;
  return false;
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
  int i,j,k,v1,v2,v3,c = 0;
  double sum,sum1,sum2,l,l_inv,d1,d2,delta,E_G,E_total = 0.0;
  bool found;
  std::set<int>::const_iterator it;
  SYNARMOSMA::hash_map::const_iterator qt;
  SYNARMOSMA::Graph G;
  const int na = double(cardinality(0,-1));
  const int nv = (signed) skeleton->events.size();
  const int nt = (signed) codex.size();
  double R[nv],gvalue[nv],length_deviation[nv],rho[nv];
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
#pragma omp parallel default(shared) private(i,j,G)
  {
#endif
  std::vector<double> tangle;
  for(i=0; i<nt; ++i) {
    tangle.push_back(0.0);
  }
#ifdef _OPENMP
#pragma omp for
#endif
  for(i=0; i<nv; ++i) {
    for(j=0; j<nt; ++j) {
      tangle[j] = 0.0;
    }
    if (!skeleton->events[i].active()) {
      skeleton->events[i].set_entwinement(tangle);
      continue;
    }
    if (!skeleton->events[i].get_topology_modified()) continue;
    for(j=0; j<nt; ++j) {
      if (!skeleton->events[i].active(j)) continue;
      tangle[j] = 0.5*double(vertex_dimension(i,j) - 1) + compute_temporal_vorticity(i,j);
      compute_graph(&G,i,j);
      if (G.order() < 2) continue;
      tangle[j] += G.completeness();
      tangle[j] += G.entwinement()/double(G.order() - 1);
    }
    skeleton->events[i].set_entwinement(tangle);
    skeleton->events[i].set_topological_dimension(vertex_dimension(i,-1));
    skeleton->events[i].set_topology_modified(false);
  }
#ifdef _OPENMP
  }
#endif

  for(i=0; i<nv; ++i) {
    sum = 0.0;
    for(j=0; j<nt; ++j) {
      sum += skeleton->events[i].entwinement[j];
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
    assert(!std::isnan(R[v1])); 
    assert(!std::isnan(rho[v1])); 
    assert(!std::isnan(length_deviation[v1]));
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

  E_G = representational_energy(false);

  global_deficiency =  E_G + 2.0*M_PI*double(euler_characteristic(-1)) - E_total;

  error = 0.0;
  for(i=0; i<na; ++i) {
    delta = skeleton->events[avertices[i]].get_deficiency();
    error += delta*delta;
  }
  error = std::sqrt(error)/double(na);
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
  error += std::abs(global_deficiency)/double(na);
}

bool Spacetime::global_operations()
{
  int i,j,n,k = 0;
  bool output = false;
  double mu,sigma = 0.0,delta = 0.0;
  std::string filename;
  SYNARMOSMA::hash_map::const_iterator qt;
  const int nd = dimension(-1);
  const int nv = (signed) skeleton->events.size();
  const int ne = (signed) skeleton->simplices[1].size();
#ifdef VERBOSE
  std::cout << "Main thread in global operations method with " << nv << " vertices." << std::endl;
#endif
  iterations++;

  skeleton->compute_modified_vertices();

  for(i=0; i<nv; ++i) {
    if (skeleton->events[i].incept == -1) skeleton->events[i].set_incept(iterations);
  }
  for(i=0; i<=Complex::ND; ++i) {
    n = (signed) skeleton->simplices[i].size();
    for(j=0; j<n; ++j) {
      if (skeleton->simplices[i][j].incept == -1) skeleton->simplices[i][j].set_incept(iterations);
    }
  }

  for(i=0; i<nv; ++i) {
    if (!skeleton->events[i].active()) continue;
    if (skeleton->events[i].topology_modified) skeleton->events[i].set_topological_dimension(vertex_dimension(i,-1));
  }

  // First perform the geometry and energy diffusion...
  compute_parity();
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
  for(i=1; i<dimension(-1); ++i) {
    for(j=0; j<(signed) skeleton->simplices[i].size(); ++j) {
      compute_simplex_energy(i,j);
    }
  }

#ifdef VERBOSE
  // Analyze the uniformity of the energy distribution...
  double E_avg = 0.0;
  double nc = double(cardinality(0,-1));
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
    histo2[vertex_dimension(i,-1)] += 1;
    if (skeleton->events[i].zero_energy()) continue;
    histogram[vertex_dimension(i,-1)] += 1;
  }
  for(i=0; i<=Complex::ND; ++i) {
    std::cout << "There are " << histo2[i] << " (" << histogram[i] << ") active " << i << "-dimensional vertices." << std::endl;
  }
#endif

  optimize();

#ifdef VERBOSE
  int ninitial = 0,ntouch = 0;
  for(i=0; i<nv; ++i) {
    if (skeleton->events[i].incept == 0) {
      ninitial++;
      if (std::abs(skeleton->events[i].deficiency) > 0.0) ntouch++;
    }
  }
  std::cout << "Percentage of perturbed initial vertices " << 100.0*double(ntouch)/double(ninitial) << std::endl;
#endif

  // Eliminate any overlapping vertices
  if (superposable) {
    superposition_fusion(0.1); assert(skeleton->consistent());
    superposition_fission(int(0.02*nv)); assert(skeleton->consistent());
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
    i = compression(mu + sigma); assert(skeleton->consistent());
  }

  if (permutable) {
    // Take care of any inter-cosmic jumping...
    for(i=1; i<nd; ++i) {
      n = (signed) skeleton->simplices[i].size();
      for(j=0; j<n; ++j) {
        if (!skeleton->simplices[i][j].active()) continue;
        compute_simplex_energy(i,j);
      }
    }
    delta = Spacetime::T_zero/std::sqrt(Spacetime::kappa*double(1+iterations));
    n = ubiquity_permutation(delta); assert(skeleton->consistent());
#ifdef VERBOSE
    std::cout << "The inter-cosmic jump for this iteration is " << n << std::endl;
#endif
  }

  if (superposable || compressible || permutable) {
    skeleton->compute_entourages(-1);
    skeleton->compute_neighbours();
    skeleton->compute_modified_vertices();
    nv = (signed) skeleton->events.size();
    std::set<int> S;
    for(i=0; i<nv; ++i) {
      if (!skeleton->active_event(i)) continue;
      if (skeleton->events[i].get_topology_modified()) S.insert(i);
    }
    geometry->compute_squared_distances(S); 
  }
  compute_volume();
  compute_obliquity();
#ifdef DEBUG
  assert(consistent(-1));
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
      write_incastrature(filename,-1);
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
          if (skeleton->events[i].active()) skeleton->events[i].topology_modified = true;
        }
        structural_deficiency();
      }
    }
  }

  RND->increment_seed();

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
  int i,j,n,l,offspring = RND->irandom(1,4);
  const int nt = (signed) codex.size();
  const int nv = (signed) skeleton->events.size();

#ifdef VERBOSE
  std::cout << "Sheet " << parent << " spawning " << offspring << " daughter sheet(s)..." << std::endl;
#endif

  for(l=0; l<offspring; ++l) {
    codex.push_back(Sheet(nt+l,parent,H->get_field(),H->get_method()));
    for(i=0; i<nv; ++i) {
      if (skeleton->events[i].active(parent)) {
        skeleton->events[i].set_active(nt+l); 
        skeleton->events[i].set_topology_modified(true);
      }
    }
    for(i=1; i<=Complex::ND; ++i) {
      n = (signed) skeleton->simplices[i].size();
      for(j=0; j<n; ++j) {
        if (skeleton->simplices[i][j].active(parent)) skeleton->simplices[i][j].set_active(nt+l);
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
    if (RND->drandom() < 0.15) codex[i].active = !codex[i].active;
    if (codex[i].active) {
      nactive++;
      candidates.insert(i);
    }
  }
  for(it=candidates.begin(); it!=candidates.end(); ++it) {
    parent = *it;
    if (RND->poisson_variate()) nspawn += sheet_fission(parent);
  }
  nactive += nspawn;
  return nspawn;
}

void Spacetime::compute_global_topology(int sheet)
{
  // To calculate the global deficiency, we need to compute the Betti numbers and
  // the fundamental group, for the total spacetime, operations that are serial...
  SYNARMOSMA::Nexus* NX = new SYNARMOSMA::Nexus;

  skeleton->compute_global_nexus(NX,sheet);

  if (sheet == -1) {
    // The global case...
    skeleton->compute_global_topology(NX,high_memory);
  }
  else {
    bool bdry;
    codex[sheet].H->compute(NX);
    if (high_memory) codex[sheet].pi1->compute(NX);
    codex[sheet].pseudomanifold = NX->pseudomanifold(&bdry);
    codex[sheet].boundary = bdry;
    if (codex[sheet].pseudomanifold) codex[sheet].orientable = NX->orientable();
  }

  delete NX;
}

void Spacetime::get_ubiquity(int d,int n,std::string& output) const
{
  const int nt = (signed) codex.size();

  if (d == 0) {
    output = "{";
    for(int i=0; i<nt-1; ++i) {
      if (skeleton->skeleton->events[n].active(i)) {
        output += "1,";
      }
      else {
        output += "0,";
      }
    }
    if (skeleton->skeleton->events[n].active(nt-1)) {
      output += "1}";
    }
    else {
      output += "0}";
    }
  }
  else {
    for(int i=0; i<nt-1; ++i) {
      if (skeleton->skeleton->simplices[d][n].active(i)) {
        output += "1,";
      }
      else {
        output += "0,";
      }
    }
    if (skeleton->skeleton->simplices[d][n].active(nt-1)) {
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

int Spacetime::ubiquity_permutation(double temperature,std::set<int>& vmodified)
{
  int i,j,k,n,m,nd,vx[2],delta,hdistance,jz = 0;
  std::set<int> ubiquity,mutated;
  std::set<int>::const_iterator it;
  SYNARMOSMA::hash_map::const_iterator qt;
  double E,alpha;
  // This parameter must be greater than zero and determines the
  // thermo-energy scale at which intercosmic jumps become common...
  const int nv = (signed) skeleton->events.size();
  const int nt = (signed) codex.size();
  const int d = dimension(-1);
  const int jlimit = nv/5;

  for(i=0; i<2*nv; ++i) {
    n = RND->irandom(1,d+1);
    m = RND->irandom(0,skeleton->simplices[d].size());
    if (!skeleton->simplices[n][m].active()) continue;
    E = skeleton->simplices[n][m].energy;
    if (E < std::numeric_limits<double>::epsilon()) continue;
    skeleton->simplices[n][m].get_ubiquity(ubiquity);
    mutated = ubiquity;
    alpha = std::exp(-temperature*E);
    hdistance = 0;
    for(j=0; j<nt; ++j) {
      if (RND->drandom() > alpha) {
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
      for(it=skeleton->simplices[n][m].vertices.begin(); it!=skeleton->simplices[n][m].vertices.end(); ++it) {
        vmodified.insert(*it);
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
      for(k=0; k<1+i; ++k) {
        qt = index_table[i-1].find(skeleton->simplices[i][j].faces[k]);
        for(it=ubiquity.begin(); it!=ubiquity.end(); ++it) {
          skeleton->simplices[i-1][qt->second].set_active(*it);
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
      skeleton->events[vx[0]].set_active(*it);
      skeleton->events[vx[1]].set_active(*it);
    }
  }
  return jz;
}
