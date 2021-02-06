#include "spacetime.h"

using namespace DIAPLEXIS;

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
  delete skeleton;
  delete geometry;
  if (score_allocated) delete[] hyphantic_notes;
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
  hyphantic_ops = "";
}

double Spacetime::condense()
{
  // First check how many ghost events and edges there are in this spacetime....
  int n = skeleton->cardinality(0),m = skeleton->cardinality(1);
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
  std::set<int> N,S,vx;
  std::vector<double> xc;
  std::vector<std::vector<double> > svalues;
  std::set<int>::const_iterator it;
  std::vector<Event> nevents;
  std::vector<Simplex> nsimplices;

  n = 0;
  for(i=0; i<nv; ++i) {
    offset[i] = -1;
    if (!skeleton->active_event(i)) continue;
    geometry->get_coordinates(i,xc);
    svalues.push_back(xc);
    nevents.push_back(skeleton->events[i]);
    offset[i] = n;
    n++;
  }
  skeleton->events = nevents;
  geometry->multiple_vertex_addition(svalues);
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
      if (!skeleton->active_simplex(i,j)) continue;
      skeleton->simplices[i][j].get_vertices(S);
      for(it=S.begin(); it!=S.end(); ++it) {
        vx.insert(offset[*it]);
      }
      skeleton->index_table[i][vx] = (signed) nsimplices.size();
      nsimplices.push_back(Simplex(vx));
      vx.clear();
    }
    skeleton->simplices[i] = nsimplices;
    nsimplices.clear();
  }
  skeleton->compute_entourages();
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
  static double htime = 0.0;
  static double gtime = 0.0;
  std::chrono::duration<double> elapsed_seconds;

  write_state("data/.previous_step.dat");
  reversible = true;

  // Begin the hyphantic phase...
  auto t1 = std::chrono::steady_clock::now();
  std::ofstream s1(hyphansis_file,std::ios::app);
  s1 << "<Iteration>" << std::endl;
  s1 << "  <Index>" << iterations << "</Index>" << std::endl;
  s1.close();
 
  hyphansis();

  std::ofstream s2(hyphansis_file,std::ios::app);
  s2 << "</Iteration>" << std::endl;
  s2.close();

  regularization(false);
  condense();
  auto t2 = std::chrono::steady_clock::now();
  elapsed_seconds = t2 - t1;
  htime += elapsed_seconds.count();

  // Start the global operations phrase...
  t1 = std::chrono::steady_clock::now();
#ifdef VERBOSE
  std::cout << "Calling global operations..." << std::endl;
#endif
  bool done = global_operations();
  t2 = std::chrono::steady_clock::now();
  elapsed_seconds = t2 - t1;
  gtime += elapsed_seconds.count();
#ifdef DEBUG
  assert(consistent());
#endif
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
  time_t now = std::time(NULL);
  std:: string cdate = std::string(std::ctime(&now));
  cdate.erase(std::remove(cdate.begin(),cdate.end(),'\n'),cdate.end());  

  pugi::xml_document logfile;
  pugi::xml_node rstep;

  logfile.load_file(log_file.c_str());
  rstep = logfile.child("LogFile").insert_child_after("FinishTime",logfile.child("LogFile").child("StartTime"));
  rstep.append_child(pugi::node_pcdata).set_value(cdate.c_str());
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
  std::cout << "Error difference in Spacetime::correctness method is " << original_error << "  " << error << "  " << delta << std::endl;
#endif
  return false;
}

void Spacetime::structural_deficiency()
{
  int i,j,v1,v2,v3,k = 0;
  double sum1,sum2,l,l_inv,alpha,d1,d2,delta,E_total = 0.0;
  bool found;
  std::set<int> S;
  std::set<int>::const_iterator it;
  SYNARMOSMA::hash_map::const_iterator qt;
  SYNARMOSMA::Graph G;
  const int na = skeleton->cardinality(0);
  const int nv = (signed) skeleton->events.size();
  const double LL = (2.0 - 0.5*abnormality_threshold)*(2.0 - 0.5*abnormality_threshold);
  const double UL = (2.0 + 0.5*abnormality_threshold)*(2.0 + 0.5*abnormality_threshold);
  const double ulimit = (1.0 + abnormality_threshold)*(1.0 + abnormality_threshold);
  const double llimit = (1.0 - abnormality_threshold)*(1.0 - abnormality_threshold);
  double R[nv],gvalue[nv],length_deviation[nv],rho[nv];
  int avertices[na];

  for(i=0; i<nv; ++i) {
    if (skeleton->active_event(i)) {
      avertices[k] = i;
      k++;
    }
    else {
      skeleton->events[i].set_entwinement(0.0);
    }
    skeleton->events[i].set_deficiency(0.0);
    skeleton->events[i].set_geometric_deficiency(0.0);
    length_deviation[i] = 0.0;
    gvalue[i] = 0.0;
    R[i] = 0.0;
    rho[i] = 0.0;
  }

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i,j,k,alpha,v1,G)
#endif
  for(i=0; i<na; ++i) {
    v1 = avertices[i];
    if (!skeleton->events[v1].get_topology_modified()) continue;
    j = skeleton->vertex_dimension(v1);
    alpha = 0.5*double(j - 1) + compute_temporal_vorticity(v1);
    skeleton->compute_graph(&G,v1);
    k = G.order();
    if (k > 1) alpha += G.completeness() + G.median_degree()/double(k - 1);
    skeleton->events[v1].set_entwinement(alpha);
    skeleton->events[v1].set_topological_dimension(j);
    skeleton->events[v1].set_topology_modified(false);
  }

  for(i=0; i<nv; ++i) {
    gvalue[i] = skeleton->events[i].get_entwinement();
  }

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i,j,k,v1,v2,v3,l,found,d1,d2) schedule(dynamic,1)
#endif
  for(i=0; i<na; ++i) {
    v1 = avertices[i];
    for(j=1+i; j<na; ++j) {
      v2 = avertices[j];
      l = std::abs(geometry->get_squared_distance(v1,v2,false));
      if (l < LL || l > UL) continue;
      // See if there is a third event that lies between these two...
      found = false;
      for(k=0; k<i; ++k) {
        v3 = avertices[k];
        d1 = std::abs(geometry->get_squared_distance(v1,v3,false));
        d2 = std::abs(geometry->get_squared_distance(v2,v3,false));
        if ((d1 > llimit && d1 < ulimit) && (d2 > llimit && d2 < ulimit)) {
          found = true;
          break;
        }
      }
      if (found) continue;
      for(k=1+i; k<j; ++k) {
        v3 = avertices[k];
        d1 = std::abs(geometry->get_squared_distance(v1,v3,false));
        d2 = std::abs(geometry->get_squared_distance(v2,v3,false));
        if ((d1 > llimit && d1 < ulimit) && (d2 > llimit && d2 < ulimit)) {
          found = true;
          break;
        }
      }
      if (found) continue;
      for(k=1+j; k<na; ++k) {
        v3 = avertices[k];
        d1 = std::abs(geometry->get_squared_distance(v1,v3,false));
        d2 = std::abs(geometry->get_squared_distance(v2,v3,false));
        if ((d1 > llimit && d1 < ulimit) && (d2 > llimit && d2 < ulimit)) {
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
    skeleton->events[v1].set_deficiency(R[v1] + alpha - coupling_constant*rho[v1]);
  }

  // Now the chromatic energy sum...
  E_total = 0.0;
  for(i=0; i<nv; ++i) {
    E_total += skeleton->events[i].get_energy();
  }
#ifdef VERBOSE
  if (std::isnan(E_total)) throw std::runtime_error("NaN in total energy!");
#endif
  global_deficiency =  representational_energy(false) + compute_temporal_nonlinearity() + 2.0*M_PI*double(skeleton->euler_characteristic()) - E_total;

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
  if (std::isnan(total_error)) throw std::runtime_error("NaN in total error!");
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
  std::string filename;
  SYNARMOSMA::hash_map::const_iterator qt;
  double mu,sigma = 0.0,delta = 0.0;
  const int nv = (signed) skeleton->events.size();
  const int ne = (signed) skeleton->simplices[1].size();
#ifdef VERBOSE
  std::cout << "Main thread in global operations method..." << std::endl;
#endif
  iterations++;

  skeleton->compute_modified_events();

  for(i=0; i<nv; ++i) {
    if (skeleton->events[i].get_incept() == -1) skeleton->events[i].set_incept(iterations);
  }
  for(i=1; i<=Complex::ND; ++i) {
    n = (signed) skeleton->simplices[i].size();
    for(j=0; j<n; ++j) {
      if (skeleton->simplices[i][j].get_incept() == -1) skeleton->simplices[i][j].set_incept(iterations);
    }
  }

  for(i=0; i<nv; ++i) {
    if (!skeleton->active_event(i)) continue;
    if (skeleton->events[i].get_topology_modified()) skeleton->events[i].set_topological_dimension(skeleton->vertex_dimension(i));
  }

  // First perform the geometry and energy diffusion...
  skeleton->compute_parity();
  compute_lightcones();
  adjust_dimension();
  compute_volume();
  compute_obliquity();
  skeleton->compute_global_topology(geometry->get_memory_type());
  structural_deficiency();

  skeleton->energy_diffusion(coupling_constant,convergence_threshold);
  for(i=1; i<=skeleton->dimension(); ++i) {
    for(j=0; j<(signed) skeleton->simplices[i].size(); ++j) {
      skeleton->compute_simplex_energy(i,j);
    }
  }
  
#ifdef VERBOSE
  // Analyze the uniformity of the energy distribution...
  double E_avg = 0.0;
  double nc = double(skeleton->cardinality(0));
  for(i=0; i<nv; ++i) {
    if (!skeleton->active_event(i)) continue;
    E_avg += skeleton->events[i].get_energy();
  }
  E_avg = E_avg/nc;
  sigma = 0.0;
  for(i=0; i<nv; ++i) {
    if (!skeleton->active_event(i)) continue;
    sigma += (skeleton->events[i].get_energy() - E_avg)*(skeleton->events[i].get_energy() - E_avg);
  }
  sigma = std::sqrt(sigma/nc);
  std::cout << "Standard deviation of event energy is " << sigma << std::endl;

  // Analyze the distribution of event dimensionalities...
  int histogram[1 + Complex::ND],histo2[1 + Complex::ND];
  for(i=0; i<=Complex::ND; ++i) {
    histogram[i] = 0;
    histo2[i] = 0;
  }
  for(i=0; i<nv; ++i) {
    if (!skeleton->active_event(i)) continue;
    histo2[skeleton->vertex_dimension(i)] += 1;
    if (skeleton->events[i].zero_energy()) continue;
    histogram[skeleton->vertex_dimension(i)] += 1;
  }
  for(i=0; i<=Complex::ND; ++i) {
    std::cout << "There are " << histo2[i] << " (" << histogram[i] << ") active " << i << "-dimensional events." << std::endl;
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

  // Eliminate any overlapping events
  if (superposable) {
    superposition_fusion(); 
    superposition_fission(int(0.02*nv)); 
  }

  if (compressible) {
    // Eliminate excessively long edges...
    // First calculate the average edge length and its variance...
    for(i=0; i<ne; ++i) {
      if (!skeleton->active_simplex(1,i)) continue;
      delta += std::abs(skeleton->simplices[1][i].get_volume());
      k++;
    }
    mu = delta/double(k);
    for(i=0; i<ne; ++i) {
      if (!skeleton->active_simplex(1,i)) continue;
      delta = std::abs(skeleton->simplices[1][i].get_volume()) - mu;
      sigma += delta*delta;
    }
    sigma = std::sqrt(sigma/double(k));
    // Get rid of edges that are more than one standard deviation from the mean...
    compression(mu + sigma); 
  }

  if (superposable || compressible) {
    skeleton->compute_entourages();
    skeleton->compute_neighbours();
    skeleton->compute_modified_events();
    n = (signed) skeleton->events.size();
    std::set<int> S;
    for(i=0; i<n; ++i) {
      if (!skeleton->active_event(i)) continue;
      if (skeleton->events[i].get_topology_modified()) S.insert(i);
    }
    geometry->compute_squared_distances(S);
  }
  compute_volume();
  compute_obliquity();
#ifdef DEBUG
  assert(skeleton->consistent());
#endif

  compute_lightcones();
  skeleton->compute_global_topology(geometry->get_memory_type());

  structural_deficiency();

  if (error < Spacetime::convergence_threshold || iterations >= max_iter) converged = true;
#ifdef VERBOSE
  std::cout << iterations << "  " << error << "  " << Spacetime::convergence_threshold << "  " << converged << std::endl;
#endif

  if (converged) {
    output = true;
    if (!diskless) {
      filename = "data/incastrature_" + date_string + "_" + pid_string + ".dot";
      skeleton->write_incastrature(filename);
    }
  }

  skeleton->RND->increment_seed();

  // Write the current state to disk if necessary...
  if (checkpoint_frequency > 0) {
    if (iterations >= 1 && (iterations % checkpoint_frequency == 0)) write_state();
  }

#ifdef VERBOSE
  std::cout << "At relaxation step " << iterations << " the global error is " << error << std::endl;
  std::cout << "Finished with global operations..." << std::endl;
#endif
  if (!diskless) write_log();
  return output;
}

bool Spacetime::consistent() const
{
  if (!skeleton->consistent()) return false;

  const int nv = (signed) skeleton->events.size();

  for(int i=0; i<nv; ++i) {
    if (std::isnan(skeleton->events[i].get_deficiency())) {
      std::cout << "NaN in deficiency for event " << i << std::endl;
      return false;
    }
    if (std::isnan(skeleton->events[i].get_energy())) {
      std::cout << "NaN in energy for event " << i << std::endl;
      return false;
    }
    if (std::isnan(skeleton->events[i].get_geometric_deficiency())) {
      std::cout << "NaN in geometric_deficiency for event " << i << std::endl;
      return false;
    }
    if (std::isnan(skeleton->events[i].get_obliquity())) {
      std::cout << "NaN in obliquity for event " << i << std::endl;
      return false;
    }
    if (std::isnan(skeleton->events[i].get_entwinement())) {
      std::cout << "NaN in entwinement for event " << i << std::endl;
      return false;
    }

    // Energy must be non-negative...
    if (skeleton->events[i].get_energy() < -std::numeric_limits<double>::epsilon()) {
      std::cout << "Negative energy for event " << i << std::endl;
      return false;
    }
    // Inactive events should have zero energy...
    if (!skeleton->active_event(i)) {
      if (!skeleton->events[i].zero_energy()) {
        std::cout << "Positive energy for inactive event " << i << std::endl;
        return false;
      }
    }
  }
  return true;
}
