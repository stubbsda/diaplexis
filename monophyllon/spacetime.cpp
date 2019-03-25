#include "spacetime.h"

using namespace DIAPLEXIS;

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
  delete skeleton;
  delete geometry;
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
  hyphantic_ops = "";
  iterations = 0;
}

void Spacetime::condense()
{
  // First check how many ghost vertices and edges there are in this spacetime....
  int i,n = 0,m = 0;
  const int nv = (signed) skeleton->events.size();
  const int ne = (signed) skeleton->simplices[1].size();
  for(i=0; i<nv; ++i) {
    if (skeleton->active_event(i)) n++;
  }
  for(i=0; i<ne; ++i) {
    if (skeleton->active_simplex(1,i)) m++;
  }
  double rho_v = double(n)/double(nv);
  double rho_e = double(m)/double(ne);
#ifdef VERBOSE
  std::cout << "Topological density is " << rho_v << " for vertices and " << rho_e << " for edges." << std::endl;
#endif
  if (rho_v > 0.5 || rho_e > 0.5) return;
  // So we need to condense this spacetime to reduce memory pressure...
  int j,offset[nv];
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
}

void Spacetime::clean() const
{
  int i,j;
  for(i=0; i<(signed) skeleton->events.size(); ++i) {
    assert(!skeleton->events[i].get_topology_modified());
  }
  for(i=1; i<=Complex::ND; ++i) {
    for(j=0; j<(signed) skeleton->simplices[i].size(); ++j) {
      assert(!skeleton->simplices[i][j].get_modified());
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
  static double htime = 0.0;
  static double gtime = 0.0;
  boost::timer::cpu_times Z;
  boost::timer::cpu_timer t1;

  write_state("data/.previous_step.dat");
  reversible = true;

  // Begin the hyphantic phase...
  std::ofstream s1(hyphansis_file,std::ios::app);
  s1 << "<Iteration>" << std::endl;
  s1.close();
 
  hyphansis();

  std::ofstream s2(hyphansis_file,std::ios::app);
  s2 << "</Iteration>" << std::endl;
  s2.close();

  regularization(false);
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
    if (advance()) break;
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
  compute_obliquity();
  structural_deficiency();
  double anterior_error = error;

  for(i=0; i<nv; ++i) {
    skeleton->events[i].set_topology_modified(true);
    skeleton->events[i].set_geometry_modified(true);
  }
  for(i=1; i<=Complex::ND; ++i) {
    for(j=0; j<(signed) skeleton->simplices[i].size(); ++j) {
      skeleton->simplices[i][j].set_modified(true);
    }
  }

  compute_volume();
  compute_obliquity();
  structural_deficiency();
  if (std::abs(anterior_error - error) < std::numeric_limits<double>::epsilon()) return true;
  std::cout << "Error difference is " << anterior_error << "  " << error << "  " << std::abs(anterior_error - error) << std::endl;
  return false;
}

void Spacetime::structural_deficiency()
{
  int i,j,v1,v2,v3,k = 0;
  double sum1,sum2,l,l_inv,alpha,d1,d2,delta,E_G,E_total = 0.0;
  bool found;
  std::set<int> S;
  std::set<int>::const_iterator it;
  SYNARMOSMA::hash_map::const_iterator qt;
  SYNARMOSMA::Graph G;
  const int na = skeleton->cardinality(0);
  const int nv = (signed) skeleton->events.size();
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
#pragma omp parallel for default(shared) private(i,j,alpha,v1,G)
#endif
  for(i=0; i<na; ++i) {
    v1 = avertices[i];
    if (!skeleton->events[v1].get_topology_modified()) continue;
    j = skeleton->vertex_dimension(v1);
    alpha = 0.5*double(j - 1) + compute_temporal_vorticity(v1);
    skeleton->compute_graph(&G,v1);
    alpha += G.completeness();
    if (G.order() > 1) alpha += G.entwinement()/double(G.order() - 1); 
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
  E_total = 0.0;
  for(i=0; i<nv; ++i) {
    E_total += skeleton->events[i].get_energy();
  }

  E_G = representational_energy(false);

  global_deficiency =  E_G + 2.0*M_PI*double(skeleton->euler_characteristic()) - E_total;

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
  std::string filename;
  SYNARMOSMA::hash_map::const_iterator qt;
  double mu,sigma = 0.0,delta = 0.0;
  const int nv = (signed) skeleton->events.size();
  const int ne = (signed) skeleton->simplices[1].size();
#ifdef VERBOSE
  std::cout << "Main thread in global operations method..." << std::endl;
#endif
  iterations++;

  skeleton->compute_modified_vertices();

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

  skeleton->energy_diffusion(Spacetime::Lambda);
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
  std::cout << "Standard deviation of vertex energy is " << sigma << std::endl;

  // Analyze the distribution of vertex dimensionalities...
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
  std::cout << "Percentage of perturbed initial vertices " << 100.0*double(ntouch)/double(ninitial) << std::endl;
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
    skeleton->compute_modified_vertices();
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
