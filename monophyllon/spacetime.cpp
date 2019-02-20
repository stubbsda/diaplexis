#include "spacetime.h"

using namespace DIAPLEXIS;

const double Spacetime::convergence_threshold = 0.00001;
const double Spacetime::T_zero = 500.0;
const double Spacetime::kappa = 1.35;
const double Spacetime::Lambda = 0.2;
const int Spacetime::N_EXP;
const int Spacetime::N_IMP;
const std::string Spacetime::EXP_OP[] = {"D","Ux","Ox","R","C","N","A","G","Sg","Sm","Y"};
const std::string Spacetime::IMP_OP[] = {"I","Um","Om","E","F","P","V","Δ"};

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
  delete RND;
  delete skeleton;
  delete geometry;
}

void Spacetime::allocate()
{
  // Default geometry (Euclidean, absolute, dimensionally 
  // uniform, background dimension = 3)
  geometry = new SYNARMOSMA::Geometry;
  skeleton = new Complex;
  RND = new SYNARMOSMA::Random;
}

void Spacetime::set_checkpoint_frequency(int a)
{
  checkpoint_frequency = a;
}

void Spacetime::restart(const std::string& filename,bool save_seed)
{
  if (save_seed) {
    unsigned int n = RND->get_seed();
    clear();
    read_parameters(filename);
    RND->set_seed(n);
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
  compute_entourages();
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
      n = RND->irandom(nv);
      if (!skeleton->events[n].active) continue;
      if (skeleton->events[n].topological_dimension > 1) continue;
      break;
    } while(true);
    p_old = affinity[n];
    do {
      p = RND->irandom(nprocs);
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

void Spacetime::get_energy_values(std::vector<double>& output) const
{
  int i;
  const int nv = (signed) skeleton->events.size();

  output.clear();

  for(i=0; i<nv; ++i) {
    if (!skeleton->events[i].active) continue;
    output.push_back(skeleton->events[i].get_energy());
  }
}

void Spacetime::get_deficiency_values(std::vector<double>& output) const
{
  int i;
  const int nv = (signed) skeleton->events.size();

  output.clear();

  for(i=0; i<nv; ++i) {
    if (!skeleton->events[i].active) continue;
    output.push_back(skeleton->events[i].deficiency);
  }
}

std::string Spacetime::implicative_scale(int key,std::vector<double>& parameters) const
{
  // This consists of twelve "notes", eight of which belong to the scale itself 
  // (diatonic notes) and four chromatic notes
  // The implicative scale is the treble clef in the F major scale, so the following 
  // twelve piano keys in ascending pitch, 
  // F4, G4, G4 sharp, A4, A4 sharp, C5, C5 sharp, D5, E5, F5, F5 sharp, G5
  // The four chromatic notes are G4 sharp, C5 sharp, F5 sharp and G5
  // The twelve implicative operations are Um, Om, V, P, I1, I2, E1, E2, E3, F1, F2 and F3
  // V <=> G5*
  // P <=> F5 sharp*
  // E3 <=> F5
  // F3 <=> E5
  // Om <=> D5
  // F2 <=> C5 sharp*
  // F1 <=> C5 
  // E2 <=> A4 sharp
  // Um <=> A4
  // I2 <=> G4 sharp*
  // E1 <=> G4
  // I1 <=> F4 
  std::string output = "NULL";
  parameters.clear();
  switch (key) {
    case 45:
      output = "I";
      parameters.push_back(0.25);
      break;
    case 47:
      output = "E";
      parameters.push_back(0.25);
      break;
    case 48:
      output = "I";
      parameters.push_back(0.75);
      break;
    case 49:
      output = "Um";
      break;
    case 50:
      output = "E";
      parameters.push_back(0.5);
      break;
    case 52:
      output = "F";
      parameters.push_back(0.25);
      break;
    case 53:
      output = "F";
      parameters.push_back(0.5);
      break;
    case 54:
      output = "Om";
      break;
    case 56:
      output = "F";
      parameters.push_back(0.75);
      break;
    case 57:
      output = "E";
      parameters.push_back(0.75);
      break;
    case 58:
      output = "P";
      parameters.push_back(double(3));
      break;
    case 59:
      output = "V";
      break;
    default:
      throw std::invalid_argument("Illegal key value in implicative scale!"); 
      break;
  }
  return output;
}

std::string Spacetime::explicative_scale(int key,std::vector<double>& parameters) const
{
  // This consists of twelve "notes", eight of which belong to the scale itself 
  // (diatonic notes) and four chromatic notes
  // The explicative scale is the bass clef in the F major scale, so the following 
  // twelve piano keys in descending pitch, 
  // A3, G3, F3 sharp, F3, E3, D3, C3 sharp, C3, A2 sharp, A2, G2, F2 
  // The four chromatic notes are A3, G3, F3 sharp and C3 sharp
  // The twelve explicative operations are D, Ox, R, C, G, Sg, Sm, A, N1, N2, Ux1 and Ux2
  // G <=> F2
  // A <=> G2
  // D <=> A2
  // C <=> A2 sharp
  // Sg <=> C3
  // N2 <=> C3 sharp*
  // Ux2 <=> D3 
  // Sm <=> E3
  // R <=> F3
  // Ox <=> F3 sharp*
  // N1 <=> G3*
  // Ux1 <=> A3*
  std::string output = "NULL";
  parameters.clear();
  switch (key) {
    case 37:
      output = "Ux";
      parameters.push_back(0.2);
      break;
    case 35:
      output = "N";
      parameters.push_back(2.5);
      break;
    case 34:
      output = "Ox";
      break;
    case 33:
      output = "R";
      break;
    case 32:
      output = "Sm";
      break;
    case 30:
      output = "Ux";
      parameters.push_back(0.5);
      break;
    case 29:
      output = "N";
      parameters.push_back(1.2);
      break;
    case 28:
      output = "Sg";
      break;
    case 26:
      output = "C";
      break;
    case 25:
      output = "D";
      break;
    case 23:
      output = "A";
      break;
    case 21:
      output = "G";
      break;
    default:
      throw std::invalid_argument("Illegal key value in explicative scale!");
      break;
  }
  return output;
}

void Spacetime::implication(std::string& output) const
{
  // Should return one of the implicative operators: {F,Um,Om,E,I,P,V,Δ}
  double alpha;
  if (iterations < 50) {
    alpha = RND->drandom();
    if (alpha < 0.3) {
      output = "F";
    }
    else if (alpha < 0.6) {
      if (RND->drandom() < 0.5) {
        output = "E";
      }
      else {
        output = "I";
      }
    }
    else if (alpha < 0.75) {
      output = "Om";
    }
    else if (alpha < 0.9) {
      output = "Um";
    }
    else {
      if (RND->drandom() < 0.33) {
        output = "P";
      }
      else {
        output = "V";
      }
    }
  }
  else {
    if (RND->drandom() < 0.5) {
      if (RND->drandom() < 0.5) {
        output = "P";
      }
      else {
        output = "V";
      }
    }
    else {
      alpha = RND->drandom();
      if (alpha < 0.4) {
        output = "Om";
      }
      else if (alpha < 0.6) {
        output = "Um";
      }
      else if (alpha < 0.8) {
        output = "F";
      }
      else {
        if (RND->drandom() < 0.67) {
          output = "E";
        }
        else {
          output = "I";
        }
      }
    }
  }
}

void Spacetime::explication(std::string& output) const
{
  // Should return one of the explicative operators: {G,C,A,Sg,Sm,D,N,Y,R,Ox,Ux}
  //if (RND->drandom() < 0.1) return 'G';
  double alpha;
  if (iterations < 50) {
    if (RND->drandom() < (0.35 + 1.0/(1+iterations/2))) {
      output = "Sg";
    }
    else {
      if (iterations <= 10) {
        output = "C";
      }
      else {
        if (RND->drandom() < 0.5) {
          output = "C";
        }
        else {
          output = "A";
        }
      }
    }
  }
  else {
    alpha = RND->drandom();
    if (alpha < 0.25) {
      output = "C";
    }
    else if (alpha < 0.5) {
      output = "Sg";
    }
    else if (alpha < 0.75) {
      output = "N";
    }
    else {
      output = "Ux";
    }
  }
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

}

bool Spacetime::advance()
{
  static double htime = 0.0;
  static double gtime = 0.0;
  boost::timer::cpu_times Z;
  boost::timer::cpu_timer t1;

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

void Spacetime::write_vertex_data(int v) const
{
  if (v < 0 || v >= (signed) skeleton->events.size()) return;
  std::cout << "For vertex " << v << " we have:" << std::endl;
  std::cout << "    Incept = " << skeleton->events[v].incept << std::endl;
  std::cout << "    Structural deficiency = " << skeleton->events[v].deficiency << std::endl;
  std::cout << "    Energy = " << skeleton->events[v].get_energy() << std::endl;
  std::cout << "    Orthogonality = " << skeleton->events[v].obliquity << std::endl;
  std::cout << std::endl;
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
    skeleton->events[v1].deficiency = R[v1] + skeleton->events[v1].obliquity + length_deviation[v1] + skeleton->events[v1].curvature - Spacetime::Lambda*rho[v1];
    skeleton->events[v1].geometric_deficiency = skeleton->events[v1].obliquity + length_deviation[v1] + skeleton->events[v1].curvature;
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

  compute_delta();

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
    compute_curvature();
    compute_obliquity();
    skeleton->compute_global_topology();
    structural_deficiency();
  }

  energy_diffusion();
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
    compute_entourages();
    compute_topological_dependency(vmodified);
    compute_geometric_dependency(vmodified);
  }

  // Calculate the local and global errors
  geometry->compute_squared_distances(vmodified);
  compute_volume();
  compute_curvature();
  compute_obliquity();
#ifdef DEBUG
  assert(skeleton->consistent());
#endif

  compute_lightcones();
  skeleton->compute_global_topology();

  structural_deficiency();

  if (instrument_convergence) analyze_convergence();

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

  RND->increment_seed();

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

void Spacetime::analyze_convergence()
{
  /*
  int i,j,m;
  SYNARMOSMA::hash_map::const_iterator qt;
  int tdelta = 0;
  double edelta = 0.0;
  const int nv = (signed) skeleton->events.size();
  const int nva = (signed) anterior.skeleton->events.size();

  for(i=0; i<nva; ++i) {
    if (skeleton->events[i].active != anterior.skeleton->events[i].active) {
      tdelta++;
    }
  }
  tdelta += (nv - nva);
  for(i=1; i<=Complex::ND; ++i) {
    m = (signed) skeleton->simplices[i].size();
    for(j=0; j<m; ++j) {
      qt = anterior.skeleton->index_table[i].find(skeleton->simplices[i][j].vertices);
      if (qt == anterior.skeleton->index_table[i].end()) {
        tdelta++;
        continue;
      }
      if (skeleton->simplices[i][j].active != anterior.skeleton->simplices[i][qt->second].active) {
        tdelta++;
      }
    }
  }
  topology_delta = double(tdelta);

  geometry_delta = geometry_change(geometry,&(anterior.geometry));

  for(i=0; i<nva; ++i) {
    edelta += std::abs(skeleton->events[i].get_energy() - anterior.skeleton->events[i].get_energy());
  }
  for(i=nva; i<nv; ++i) {
    edelta += skeleton->events[i].get_energy();
  }
  energy_delta = edelta;

  anterior.skeleton->events = skeleton->events;
  anterior.geometry.load(geometry);
  for(i=1; i<=Complex::ND; ++i) {
    anterior.skeleton->simplices[i].clear();
    for(j=0; j<(signed) skeleton->simplices[i].size(); ++j) {
      anterior.skeleton->simplices[i].push_back(skeleton->simplices[i][j]);
    }
  }

  for(i=1; i<=Complex::ND; ++i) {
    anterior.skeleton->index_table[i].clear();
    for(j=0; j<(signed) anterior.skeleton->simplices[i].size(); ++j) {
      anterior.skeleton->index_table[i][anterior.skeleton->simplices[i][j].vertices] = j;
    }
  }
  */
}

void Spacetime::build_initial_state()
{
  int i,j,in1;
  unsigned int m;
  std::string geometry_type;
  Event vt;
  Simplex S;
  const bool relational = geometry->get_relational();
  std::vector<double> svalue;

  for(m=0; m<geometry->dimension(); ++m) {
    svalue.push_back(0.0);
  }
  vt.incept = 0;
  S.incept = 0;

  if (initial_state == Initial_Topology::cartesian) {
    std::vector<std::pair<long,int> > factors;
    geometry_type = "CARTESIAN";
    SYNARMOSMA::factorize(initial_size,factors);
    j = 1;
    for(m=0; m<factors.size(); ++m) {
      j *= SYNARMOSMA::ipow(factors[m].first,factors[m].second/geometry->dimension());
    }
    const int n = j;
    const int nd = SYNARMOSMA::ipow(n,geometry->dimension()-1);
    const int nm1 = n - 1;
    const int nperturbed = 10 + int(0.01*RND->irandom(initial_size));
    const double dx = 1.0;
    int k,l,d,rvalue;
    unsigned int r;
    SYNARMOSMA::hash_map::const_iterator qt;
    std::vector<int> entourage,v;
    std::vector<int>* arrangement = new std::vector<int>[initial_size];
    std::set<int> N;
    std::set<int>::const_iterator it,jt;
    std::vector<double> vpoints;

    for(m=0; m<geometry->dimension(); ++m) {
      entourage.push_back(0);
    }
    
    // First the vertices, making sure (if this is a non-relational
    // model) that the vertex coordinates are balanced...
    for(l=0; l<initial_size; ++l) {
      k = nd;
      in1 = l/k;
      entourage[0] = in1;
      rvalue = l - in1*nd;
      vpoints.push_back(dx*(double(in1) - 0.5*double(nm1)));
      for(m=1; m<geometry->dimension(); ++m) {
        k /= n;
        in1 = rvalue/k;
        entourage[m] = in1;
        rvalue -= k*in1;
        vpoints.push_back(dx*(double(in1) - 0.5*double(nm1)));
      }
      for(m=0; m<geometry->dimension(); ++m) {
        if (entourage[m] == 0 || entourage[m] == nm1) vt.boundary = true;
      }
      skeleton->events.push_back(vt);
      vt.boundary = false;
      arrangement[l] = entourage;
    }

    if (relational) {
      geometry->create(n,geometry_type);
    }
    else {
      geometry->multiple_vertex_addition(initial_size,false,vpoints);
    }

    // Now the edges...
    for(i=0; i<initial_size; ++i) {
      v = arrangement[i];
      if (v[0] > 0) {
        in1 = (v[0]-1)*nd;
        k = nd;
        for(m=0; m<geometry->dimension()-1; ++m) {
          k /= n;
          in1 += v[m+1]*k;
        }
        if (in1 > i) {
          S.initialize(i,in1);
          skeleton->index_table[1][S.vertices] = (signed) skeleton->simplices[1].size();
          skeleton->simplices[1].push_back(S);
          S.vertices.clear();
        }
      }
      if (v[0] < nm1) {
        in1 = (v[0]+1)*nd;
        k = nd;
        for(m=0; m<geometry->dimension()-1; ++m) {
          k /= n;
          in1 += v[m+1]*k;
        }
        if (in1 > i) {
          S.initialize(i,in1);
          skeleton->index_table[1][S.vertices] = (signed) skeleton->simplices[1].size();
          skeleton->simplices[1].push_back(S);
          S.vertices.clear();
        }
      }
      for(m=1; m<geometry->dimension(); ++m) {
        if (v[m] > 0) {
          v[m] -= 1;
          in1 = v[0]*nd;
          k = nd;
          for(r=0; r<geometry->dimension()-1; ++r) {
            k /= n;
            in1 += v[r+1]*k;
          }
          if (in1 > i) {
            S.initialize(i,in1);
            skeleton->index_table[1][S.vertices] = (signed) skeleton->simplices[1].size();
            skeleton->simplices[1].push_back(S);
            S.vertices.clear();
          }
          v[m] += 1;
        }
        if (v[m] < nm1) {
          v[m] += 1;
          in1 = v[0]*nd;
          k = nd;
          for(r=0; r<geometry->dimension()-1; ++r) {
            k /= n;
            in1 += v[r+1]*k;
          }
          if (in1 > i) {
            S.initialize(i,in1);
            skeleton->index_table[1][S.vertices] = (signed) skeleton->simplices[1].size();
            skeleton->simplices[1].push_back(S);
            S.vertices.clear();
          }
          v[m] -= 1;
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
        for(m=0; m<geometry->dimension(); ++m) {
          j += SYNARMOSMA::ipow(n,geometry->dimension()-1-m)*RND->irandom(k-2,k+2);
        }
        v.push_back(j);
        // Now add some edges among the neighbours of this
        // vertex
        N.insert(v[i]);
        j = RND->irandom(2,5);
        for(l=0; l<j; ++l) {
          N.insert(skeleton->events.size());
          v.push_back(skeleton->events.size());
          geometry->vertex_addition(v[i]);
          skeleton->events.push_back(vt);
        }
        d = N.size() - 1;
        S.initialize(N);
        skeleton->simplices[d].push_back(S);
        skeleton->index_table[d][S.vertices] = (signed) skeleton->simplices[d].size() - 1;
        N.clear();
      }
    }
    if (perturb_geometry) {
      // No changes to the spacetime topology, merely the geometry of some
      // existing vertices are altered...
      if (v.empty()) {  
        for(i=0; i<nperturbed; ++i) {
          j = RND->irandom(initial_size,v);
          v.push_back(j);
        }
      }
      for(m=0; m<v.size(); ++m) {
        geometry->mutation(v[m],true,false,0.5);
      }
    }
    if (perturb_energy) {
      if (v.empty()) {
        // In this particular case I will simply alter the energy value of a single vertex
        // near the centre of the Cartesian network...
        j = 0;
        k = int(double(n)/2.0);
        for(m=0; m<geometry->dimension(); ++m) {
          j += SYNARMOSMA::ipow(n,geometry->dimension()-1-m)*RND->irandom(k-2,k+2);
        }
        skeleton->events[j].set_energy(1000.0*(0.5 + RND->drandom()/2.0));
      }
      else {
        for(m=0; m<v.size(); ++m) {
          skeleton->events[v[m]].set_energy(5.0*RND->drandom());
        }
      }
    }
  }
  else if (initial_state == Initial_Topology::singleton) {
    // An initial spacetime consisting of a single isolated vertex, though with very
    // high energy
    geometry_type = "SINGLETON";
    if (relational) {
      geometry->create(0,geometry_type);
    }
    else {
      geometry->vertex_addition(svalue);
    }
    vt.set_energy(5000.0*(0.5 + RND->drandom()/2.0));
    skeleton->events.push_back(vt);
  }
  else if (initial_state == Initial_Topology::monoplex) {
    // The initial spacetime is a single simplex of dimension initial_dim
    std::set<int> vx;
    geometry_type = "MONOPLEX";
    int ulimit = geometry->dimension();
    if (initial_dim > geometry->dimension()) ulimit = initial_dim;
    if (perturb_energy) vt.set_energy(500.0 + (2000.0/double(initial_dim))*RND->drandom());

    if (!relational) {
      svalue.clear();
      for(i=0; i<ulimit; ++i) {
        svalue.push_back(0.0);
      }
      geometry->vertex_addition(svalue);
    }
    skeleton->events.push_back(vt);
    vx.insert(0);

    for(m=1; m<=initial_dim; ++m) {
      if (!relational) {
        svalue[m-1] = 1.0;
        geometry->vertex_addition(svalue);
        svalue[m-1] = 0.0;
      }
      skeleton->events.push_back(vt);
      vx.insert(m);
    }

    if (relational) geometry->create(initial_dim,geometry_type);

    S.initialize(vx);
    skeleton->simplices[initial_dim].push_back(S);
    skeleton->index_table[initial_dim][S.vertices] = 0;
  }
  else if (initial_state == Initial_Topology::random) {
    // We will use the Erdős–Rényi random graph model (the G(n,p) variant) to
    // assemble a random graph, with n = initial_size
    int level = 2,k = 0,ulimit;
    double percent = 0.0;
    bool found;
    std::set<int>* N;
    std::vector<Simplex> svector;
    std::vector<Simplex>::const_iterator vit;
    SYNARMOSMA::hash_map::const_iterator qt;
    std::set<int> v,vx,current;
    std::set<int>::const_iterator it,chk;
    const double nv = double(initial_size);

    geometry_type = "RANDOM";
    // Add the vertices
    for(i=0; i<initial_size; ++i) {
      skeleton->events.push_back(vt);
    }
    
    // Next distribute energy among the vertices, ensuring at least one vertex 
    // has non-zero energy
    do {
      i = RND->irandom(initial_size);
      if (!skeleton->events[i].zero_energy()) continue;
      skeleton->events[i].set_energy(10.0*RND->drandom());
      k += 1;
      percent = double(k)/nv;
    } while(percent < 0.1);

    // Next initialize the geometry
    if (relational) {
      geometry->create(initial_size,geometry_type);
    }
    else {
      std::vector<double> climits;
      for(m=0; m<geometry->dimension(); ++m) {
        climits.push_back(-10.0);
        climits.push_back(10.0);
      }
      geometry->multiple_vertex_addition(initial_size,true,climits);
    }

    // Now initialize the probability distribution for creating edges...
    RND->initialize_bernoulli(edge_probability);
    // and create the edges
    for(i=0; i<initial_size; ++i) {
      for(j=1+i; j<initial_size; ++j) {
        if (RND->bernoulli_variate() == false) continue;
        S.initialize(i,j);
        skeleton->index_table[1][S.vertices] = (signed) skeleton->simplices[1].size();
        skeleton->simplices[1].push_back(S);
        skeleton->events[i].neighbours.insert(j);
        skeleton->events[j].neighbours.insert(i);
        S.vertices.clear();
      }
    }

    // Finally we have to compute all the n-skeleton->simplices (n > 1) that are implied
    // by this random graph
    do {
      N = new std::set<int>[level];
      ulimit = (signed) skeleton->simplices[level-1].size();
      for(j=0; j<ulimit; ++j) {
        v = skeleton->simplices[level-1][j].vertices;
        i = 0;
        for(it=v.begin(); it!=v.end(); ++it) {
          N[i] = skeleton->events[*it].neighbours;
          i++;
        }
        for(it=N[0].begin(); it!=N[0].end(); ++it) {
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
        for(it=vx.begin(); it!=vx.end(); ++it) {
          in1 = *it;
          current = v;
          current.insert(in1);
          svector.push_back(Simplex(current));
        }
        vx.clear();
        v.clear();
      }
      delete[] N;
      if (svector.empty()) break;
      for(vit=svector.begin(); vit!=svector.end(); ++vit) {
        qt = skeleton->index_table[level].find(vit->vertices);
        if (qt == skeleton->index_table[level].end()) {
          skeleton->simplices[level].push_back(*vit);
          skeleton->index_table[level][vit->vertices] = (signed) skeleton->simplices[level].size() - 1;
        }
      }
      svector.clear();
#ifdef VERBOSE
      if (!diskless) std::cout << "Finished the implied " << level << "-skeleton->simplices..." << std::endl;
#endif
      level += 1;
    } while(true);
  }
 
  regularization(false);
}

void Spacetime::initialize()
{
  std::stringstream day,month,year,pid;
#ifdef VERBOSE
  boost::timer::cpu_timer t1;
#endif

  if (!diskless) {
    if (std::system("mkdir -p data") < 0) throw std::runtime_error("Unable to create data directory!");
  }

  // Get the date and the process ID so that we can construct the
  // state_file and log_file names...
  start_time = boost::posix_time::second_clock::local_time();
  boost::gregorian::date d = start_time.date();
  boost::gregorian::date::ymd_type ymd = d.year_month_day();

  day.width(2);
  day << std::setfill('0') << ymd.day.as_number();
  month.width(2);
  month << std::setfill('0') << ymd.month.as_number();
  year.width(4);
  year << ymd.year;
  // Using the ISO 8601 norm
  date_string = year.str() + "-" + month.str() + "-" + day.str();

  pid.width(5);
  pid << std::setfill('0') << getpid();
  pid_string = pid.str();
  if (diskless) {
    state_file = "";
    log_file = "";
    hyphansis_file = "";
  }
  else {
    state_file += "_" + date_string + "_" + pid_string;
    log_file += "_" + date_string + "_" + pid_string + ".log";
    hyphansis_file += "_" + date_string + "_" + pid_string + ".xml";
  }

  if (initial_state == Initial_Topology::diskfile) {
    read_state(input_file);
    if (converged) {
      std::cout << "This spacetime complex has already converged!" << std::endl;
    }
    else if (iterations >= max_iter) {
      std::cout << "This spacetime complex has already exceeded its maximum iterations!" << std::endl;
    }
  }
  else {
    build_initial_state();
    geometry->compute_squared_distances();
    skeleton->compute_simplicial_dimension();
    adjust_dimension();
    skeleton->compute_parity();
    compute_lightcones();
    compute_volume();
    compute_curvature();
    compute_obliquity();
    compute_lightcones();
    skeleton->compute_global_topology();
    structural_deficiency();
  }
#ifdef DEBUG
  assert(skeleton->energy_check());
#endif

  if (instrument_convergence) {
    /*
    anterior.skeleton->events = skeleton->events;
    for(int i=1; i<=Complex::ND; ++i) {
      anterior.skeleton->simplices[i] = skeleton->simplices[i];
      anterior.skeleton->index_table[i] = skeleton->index_table[i];
    }
    */
  }
  if (diskless) return;

#ifdef VERBOSE
  if (skeleton->pseudomanifold) {
    if (skeleton->orientable) {
      if (skeleton->boundary) {
        std::cout << "The spacetime complex is an orientable pseudomanifold with boundary." << std::endl;
      }
      else {
        std::cout << "The spacetime complex is an orientable pseudomanifold." << std::endl;
      }
    }
    else {
      if (skeleton->boundary) {
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
  skeleton->write_topology();
#endif

  write_state();
  write_log();

  if (iterations == 0) {
    std::ofstream s(hyphansis_file,std::ios::trunc);
    s << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << std::endl;
    s.close();
  }
#ifdef VERBOSE
  std::cout << "At relaxation step " << iterations << " the global error is " << error << std::endl;
#endif

#ifdef VERBOSE
  t1.stop();
  std::cout << "Spacetime initialization required " << boost::timer::format(t1.elapsed(),3,"%w") << " seconds to complete." << std::endl;
#endif
}
