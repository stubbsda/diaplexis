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

void Spacetime::set_checkpoint_frequency(int a)
{
  checkpoint_frequency = a;
}

void Spacetime::restart(const std::string& filename,bool save_seed)
{
  if (save_seed) {
    unsigned int n = skeleton->RND->get_seed();
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
    if (skeleton->events[i].active) n++;
  }
  for(i=0; i<ne; ++i) {
    if (skeleton->simplices[1][i].active) m++;
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
  std::set<int> vx;
  std::vector<double> xc;
  std::vector<std::vector<double> > svalues;
  std::set<int>::const_iterator it;
  std::vector<Event> nevents;
  std::vector<Simplex> nsimplices;

  n = 0;
  for(i=0; i<nv; ++i) {
    offset[i] = -1;
    if (!skeleton->events[i].active) continue;
    geometry->get_coordinates(i,xc);
    svalues.push_back(xc);
    nevents.push_back(skeleton->events[i]);
    offset[i] = n;
    n++;
  }
  skeleton->events = nevents;
  geometry->multiple_vertex_addition(svalues);
  for(i=0; i<n; ++i) {
    for(it=skeleton->events[i].neighbours.begin(); it!=skeleton->events[i].neighbours.end(); ++it) {
      vx.insert(offset[*it]);
    }
    skeleton->events[i].neighbours = vx;
    skeleton->events[i].entourage.clear();
    vx.clear();
  }
  for(i=1; i<=Complex::ND; ++i) {
    m = (signed) skeleton->simplices[i].size();
    skeleton->index_table[i].clear();
    for(j=0; j<m; ++j) {
      if (!skeleton->simplices[i][j].active) continue;
      S = skeleton->simplices[i][j];
      for(it=S.vertices.begin(); it!=S.vertices.end(); ++it) {
        vx.insert(offset[*it]);
      }
      S.vertices = vx;
      S.entourage.clear();
      S.calculate_faces();
      skeleton->index_table[i][S.vertices] = (signed) nsimplices.size();
      nsimplices.push_back(S);
      vx.clear();
    }
    skeleton->simplices[i] = nsimplices;
    nsimplices.clear();
  }
  skeleton->compute_entourages();
}

void Spacetime::distribute(int nprocs) const
{
  assert(nprocs > 0);
  int i,j,k,n,p,p_old,ecount,bcount,its = 0,volume[nprocs],max_dim = 0,current = -1,cproc = 0,nreal = 0;
  int cneighbour;
  bool done,bdry;
  double cost,current_cost;
  std::vector<int> affinity;
  std::vector<std::pair<int,int> > candidates;
  std::set<int> vx,neg,next;
  std::set<int>::const_iterator it;
  const int nv = (signed) skeleton->events.size();
  const int D = skeleton->dimension();

  // The algorithm should begin by making current equal to the 
  // vertex with the highest simplicial dimension...
  for(i=0; i<nv; ++i) {
    affinity.push_back(-1);
    if (!skeleton->events[i].active) continue;
    nreal++;
    if (skeleton->events[i].topological_dimension > max_dim) {
      current = i;
      max_dim = skeleton->events[i].topological_dimension;
    }
  }
  for(i=0; i<nprocs; ++i) {
    volume[i] = 0;
  }
#ifdef VERBOSE
  std::cout << "Distributing " << nreal << " vertices across " << nprocs << " processor elements, with maximum simplicial dimension of " << max_dim << "." << std::endl;
#endif
  // First do the higher-dimensional vertices...
  do {
    // So take all of the n-skeleton->simplices (n > 1) that are active and which 
    // contain the current vertex and assign all of their vertices to the 
    // current processor...
    if (skeleton->events[current].topological_dimension > 1) {
      affinity[current] = cproc;
      volume[cproc] += 1;
      // Need to loop over all d-skeleton->simplices, d > 1
      for(i=2; i<=skeleton->events[current].topological_dimension; ++i) {
        for(j=0; j<(signed) skeleton->simplices[i].size(); ++j) {
          if (!skeleton->simplices[i][j].active) continue;
          if (skeleton->simplices[i][j].contains(current)) {
            skeleton->simplices[i][j].get_vertices(vx);
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
      // Now handle the neighbouring d-skeleton->simplices (d > 1)...
      do {
        for(i=2; i<=D; ++i) {
          for(j=0; j<(signed) skeleton->simplices[i].size(); ++j) {
            if (!skeleton->simplices[i][j].active) continue;
            skeleton->simplices[i][j].get_vertices(vx);
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
      if (!skeleton->events[i].active) continue;
      if (affinity[i] == -1 && skeleton->events[i].topological_dimension > 1) {
        done = false;
        current = i;
        break;
      }
    }
    cproc = (cproc + 1)%nprocs;
  } while(!done);
  // Now everything else, adding vertices to equilibrate the population counts 
  // and minimize the number of boundary edges...
#ifdef VERBOSE
  std::cout << "Done handling higher-dimensional vertices, now doing edges and vertices..." << std::endl;
#endif  
  for(i=0; i<nprocs; ++i) {
    if (volume[i] == 0) {
      // Find an initial vertex, ideally far from any existing vertices that have 
      // been assigned to a processor...
      candidates.clear();
      for(j=0; j<nv; ++j) {
        if (!skeleton->events[j].active) continue;
        if (affinity[j] > -1) continue;
        cneighbour = 0;
        for(it=skeleton->events[j].neighbours.begin(); it!=skeleton->events[j].neighbours.end(); ++it) {
          if (affinity[*it] > -1) cneighbour++; 
        }
        candidates.push_back(std::pair<int,int>(j,cneighbour));         
      }
      if (candidates.empty()) break;
      std::sort(candidates.begin(),candidates.end(),SYNARMOSMA::pair_predicate_int);
      n = candidates[0].first;
      affinity[n] = i;
      volume[i] += 1;
    }
    done = false;
    do {
      next.clear();
      for(j=0; j<nv; ++j) {
        if (!skeleton->events[j].active) continue;
        if (affinity[j] > -1) continue;
        for(it=skeleton->events[j].neighbours.begin(); it!=skeleton->events[j].neighbours.end(); ++it) {
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
    } while(!done);
  }
  for(i=0; i<nv; ++i) {
    if (!skeleton->events[i].active) continue;
    if (affinity[i] == -1) {
      for(it=skeleton->events[i].neighbours.begin(); it!=skeleton->events[i].neighbours.end(); ++it) {
        p = affinity[*it];
        if (p > -1) {
          affinity[i] = p;
          volume[p] += 1;
          break;
        }
      }
    }
  }
  current_cost = skeleton->distribution_fitness(volume,affinity,nprocs);
#ifdef VERBOSE
  std::cout << "At iteration 0 the cost is " << current_cost << std::endl;
#endif
  do {
    do {
      n = skeleton->RND->irandom(nv);
      if (!skeleton->events[n].active) continue;
      if (skeleton->events[n].topological_dimension > 1) continue;
      break;
    } while(true);
    p_old = affinity[n];
    do {
      p = skeleton->RND->irandom(nprocs);
      if (p != p_old) break;
    } while(true);
    affinity[n] = p;
    volume[p] += 1;
    volume[p_old] -= 1;
    cost = skeleton->distribution_fitness(volume,affinity,nprocs);
    if (cost < current_cost) {
      current_cost = cost;
    }
    else {
      affinity[n] = p_old;
      volume[p] -= 1;
      volume[p_old] += 1;
    }
    its++;
#ifdef VERBOSE
    if (its % 250 == 0) std::cout << "At iteration " << its << " the cost is " << current_cost << std::endl; 
#endif
  } while(its < 10000);
  // Some basic sanity checks: every vertex has a processor...
  n = 0;
  for(i=0; i<nv; ++i) {
    if (!skeleton->events[i].active) continue;
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
      if (!skeleton->events[j].active) continue;
      if (affinity[j] != i) continue;
      for(it=skeleton->events[j].neighbours.begin(); it!=skeleton->events[j].neighbours.end(); ++it) {
        k = *it;
        if (affinity[k] == i) {
          ecount++;
        }
        else {
          bcount++;
        }
      }
    }
#ifdef VERBOSE 
    std::cout << "Processor " << i << " owns " << volume[i] << " vertices and " << ecount << " edges, with " << bcount << " boundary edges." << std::endl;
#endif
  }
  write_distribution(affinity);
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
    skeleton->events[i].topology_modified = true;
    skeleton->events[i].geometry_modified = true;
  }
  for(i=1; i<=Complex::ND; ++i) {
    for(j=0; j<(signed) skeleton->simplices[i].size(); ++j) {
      skeleton->simplices[i][j].modified = true;
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
  int i,j,v1,v2,v3,k = 0;
  double sum1,sum2,l,l_inv,d1,d2,delta,E_G,E_total = 0.0;
  bool found;
  std::set<int>::const_iterator it;
  SYNARMOSMA::hash_map::const_iterator qt;
  SYNARMOSMA::Graph G;
  const int na = skeleton->cardinality(0);
  const int nv = (signed) skeleton->events.size();
  double R[nv],gvalue[nv],length_deviation[nv],rho[nv];
  int avertices[na];

  for(i=0; i<nv; ++i) {
    if (skeleton->events[i].active) {
      avertices[k] = i;
      k++;
    }
    else {
      skeleton->events[i].entwinement = 0.0;
    }
    skeleton->events[i].deficiency = 0.0;
    skeleton->events[i].geometric_deficiency = 0.0;
    length_deviation[i] = 0.0;
    gvalue[i] = 0.0;
    R[i] = 0.0;
    rho[i] = 0.0;
  }

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i,j,v1,G)
#endif
  for(i=0; i<na; ++i) {
    v1 = avertices[i];
    if (!skeleton->events[v1].topology_modified) continue;
    j = skeleton->vertex_dimension(v1);
    skeleton->events[v1].entwinement = 0.5*double(j - 1) + compute_temporal_vorticity(i);
    skeleton->compute_graph(&G,v1);
    skeleton->events[v1].entwinement += G.completeness();
    if (G.order() > 1) skeleton->events[v1].entwinement += G.entwinement()/double(G.order() - 1); 
    skeleton->events[v1].topological_dimension = j;
    skeleton->events[v1].topology_modified = false;
  }

  for(i=0; i<nv; ++i) {
    gvalue[i] = skeleton->events[i].entwinement;
  }

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i,j,k,v1,v2,v3,l,found,d1,d2) schedule(dynamic,1)
#endif
  for(i=0; i<na; ++i) {
    v1 = avertices[i];
    for(j=1+i; j<na; ++j) {
      v2 = avertices[j];
      l = geometry->get_squared_distance(v1,v2,false);
      if (l < 3.8025 || l > 4.2025) continue;
      // See if there is a third vertex that lies between these two...
      found = false;
      for(k=0; k<i; ++k) {
        v3 = avertices[k];
        d1 = geometry->get_squared_distance(v1,v3,false);
        d2 = geometry->get_squared_distance(v2,v3,false);
        if ((d1 > 0.81 && d1 < 1.21) && (d2 > 0.81 && d2 < 1.21)) {
          found = true;
          break;
        }
      }
      if (found) continue;
      for(k=1+i; k<j; ++k) {
        v3 = avertices[k];
        d1 = geometry->get_squared_distance(v1,v3,false);
        d2 = geometry->get_squared_distance(v2,v3,false);
        if ((d1 > 0.81 && d1 < 1.21) && (d2 > 0.81 && d2 < 1.21)) {
          found = true;
          break;
        }
      }
      if (found) continue;
      for(k=1+j; k<na; ++k) {
        v3 = avertices[k];
        d1 = geometry->get_squared_distance(v1,v3,false);
        d2 = geometry->get_squared_distance(v2,v3,false);
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
#pragma omp parallel for default(shared) private(i,j,k,v1,l,l_inv,it,sum1,sum2) 
#endif
  for(i=0; i<na; ++i) {
    v1 = avertices[i];
    sum1 = 0.0;
    sum2 = 0.0;
    k = (skeleton->events[v1].neighbours.empty()) ? 1 : (signed) skeleton->events[v1].neighbours.size();
    length_deviation[v1] = 0.0;
    for(it=skeleton->events[v1].neighbours.begin(); it!=skeleton->events[v1].neighbours.end(); ++it) {
      j = *it;
      l = geometry->get_squared_distance(v1,j,false);
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
      //if (skeleton->simplices[1][n].orientation == SYNARMOSMA::DISPARATE) continue;
      // How to determine if j lies in the chronological future of i?
      //sigma = (skeleton->simplices[1][n].orientation == SYNARMOSMA::AFTER) ? 1 : -1;
      //if (j < i) sigma = -sigma;
      //sum2 += double(sigma)*skeleton->events[j].get_energy()*l_inv;
    }
    length_deviation[v1] = length_deviation[v1]/double(k);
    R[v1] = gvalue[v1]; // - sum1/double(k);
    rho[v1] = skeleton->events[v1].get_energy(); // + sum2/double(k);
  }
  
  // The local part of the structure equations
  for(i=0; i<na; ++i) {
    v1 = avertices[i];
    skeleton->events[v1].deficiency = R[v1] + skeleton->events[v1].obliquity + length_deviation[v1] - Spacetime::Lambda*rho[v1];
    skeleton->events[v1].geometric_deficiency = skeleton->events[v1].obliquity + length_deviation[v1];
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
    delta = skeleton->events[avertices[i]].deficiency;
    error += delta*delta;
  }
  error = std::sqrt(error)/double(na);
#ifdef VERBOSE
  double total_error = 0.0;
  for(i=0; i<na; ++i) {
    total_error += std::abs(skeleton->events[avertices[i]].deficiency);
  }
  std::cout << "The total error is " << total_error << std::endl;
#endif

  // Sanity check...
#ifdef DEBUG
  int nv_test = 0;
  for(i=0; i<na; ++i) {
    if (std::abs(skeleton->events[avertices[i]].deficiency) < std::numeric_limits<double>::epsilon()) continue;
    nv_test++;
  }
  if (nv_test == 0) assert(error < std::numeric_limits<double>::epsilon());
#endif
  //error += std::abs(global_deficiency)/double(na);
}

bool Spacetime::global_operations()
{
  int i,j,n,k = 0;
  bool output = false;
  std::string filename;
  std::set<int> vmodified;
  SYNARMOSMA::hash_map::const_iterator qt;
  std::vector<int> fusion;
  double mu,sigma = 0.0,delta = 0.0;
  const int nv = (signed) skeleton->events.size();
  const int ne = (signed) skeleton->simplices[1].size();
#ifdef VERBOSE
  std::cout << "Main thread in global operations method..." << std::endl;
#endif
  iterations++;

  skeleton->compute_delta(modified_vertices);

  for(i=0; i<nv; ++i) {
    if (skeleton->events[i].incept == -1) skeleton->events[i].incept = iterations;
  }
  for(i=0; i<=Complex::ND; ++i) {
    n = (signed) skeleton->simplices[i].size();
    for(j=0; j<n; ++j) {
      if (skeleton->simplices[i][j].incept == -1) skeleton->simplices[i][j].incept = iterations;
    }
  }

  for(i=0; i<nv; ++i) {
    if (!skeleton->events[i].active) continue;
    if (skeleton->events[i].topology_modified) skeleton->events[i].topological_dimension = skeleton->vertex_dimension(i);
  }

  // First perform the geometry and energy diffusion...
  skeleton->compute_parity();
  compute_lightcones();
  if (adjust_dimension()) {
    compute_volume();
    compute_obliquity();
    skeleton->compute_global_topology(geometry->get_memory_type());
    structural_deficiency();
  }

  skeleton->energy_diffusion(Spacetime::Lambda);
  for(i=1; i<skeleton->dimension(); ++i) {
    for(j=0; j<(signed) skeleton->simplices[i].size(); ++j) {
      skeleton->compute_simplex_energy(i,j);
    }
  }
  
#ifdef VERBOSE
  // Analyze the uniformity of the energy distribution...
  double E_avg = 0.0;
  double nc = double(skeleton->cardinality(0));
  for(i=0; i<nv; ++i) {
    if (!skeleton->events[i].active) continue;
    E_avg += skeleton->events[i].get_energy();
  }
  E_avg = E_avg/nc;
  sigma = 0.0;
  for(i=0; i<nv; ++i) {
    if (!skeleton->events[i].active) continue;
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
    if (!skeleton->events[i].active) continue;
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
    if (skeleton->events[i].incept == 0) {
      ninitial++;
      if (std::abs(skeleton->events[i].deficiency) > 0.0) ntouch++;
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
      if (!skeleton->simplices[1][i].active) continue;
      delta += std::abs(skeleton->simplices[1][i].volume);
      k++;
    }
    mu = delta/double(k);
    for(i=0; i<ne; ++i) {
      if (!skeleton->simplices[1][i].active) continue;
      delta = (std::abs(skeleton->simplices[1][i].volume) - mu);
      sigma += delta*delta;
    }
    sigma = std::sqrt(sigma/double(k));
    // Get rid of edges that are more than one standard deviation from the mean...
    i = compression(mu+sigma,vmodified);
  }

  if (superposable || compressible) {
    skeleton->compute_entourages();
    skeleton->compute_topological_dependency(vmodified);
    skeleton->compute_geometric_dependency(vmodified);
  }

  // Calculate the local and global errors
  geometry->compute_squared_distances(vmodified);
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
