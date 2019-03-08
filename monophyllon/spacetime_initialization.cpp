#include "spacetime.h"

using namespace DIAPLEXIS;

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
    const int nperturbed = 10 + int(0.01*skeleton->RND->irandom(initial_size));
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
          j += SYNARMOSMA::ipow(n,geometry->dimension()-1-m)*skeleton->RND->irandom(k-2,k+2);
        }
        v.push_back(j);
        // Now add some edges among the neighbours of this
        // vertex
        N.insert(v[i]);
        j = skeleton->RND->irandom(2,5);
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
          j = skeleton->RND->irandom(initial_size,v);
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
          j += SYNARMOSMA::ipow(n,geometry->dimension()-1-m)*skeleton->RND->irandom(k-2,k+2);
        }
        skeleton->events[j].set_energy(1000.0*(0.5 + skeleton->RND->drandom()/2.0));
      }
      else {
        for(m=0; m<v.size(); ++m) {
          skeleton->events[v[m]].set_energy(5.0*skeleton->RND->drandom());
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
    vt.set_energy(5000.0*(0.5 + skeleton->RND->drandom()/2.0));
    skeleton->events.push_back(vt);
  }
  else if (initial_state == Initial_Topology::monoplex) {
    // The initial spacetime is a single simplex of dimension initial_dim
    std::set<int> vx;
    geometry_type = "MONOPLEX";
    int ulimit = geometry->dimension();
    if (initial_dim > geometry->dimension()) ulimit = initial_dim;
    if (perturb_energy) vt.set_energy(500.0 + (2000.0/double(initial_dim))*skeleton->RND->drandom());

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
      i = skeleton->RND->irandom(initial_size);
      if (!skeleton->events[i].zero_energy()) continue;
      skeleton->events[i].set_energy(10.0*skeleton->RND->drandom());
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
    skeleton->RND->initialize_bernoulli(edge_probability);
    // and create the edges
    for(i=0; i<initial_size; ++i) {
      for(j=1+i; j<initial_size; ++j) {
        if (skeleton->RND->bernoulli_variate() == false) continue;
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
    compute_obliquity();
    compute_lightcones();
    skeleton->compute_global_topology(geometry->get_memory_type());
    structural_deficiency();
  }
#ifdef DEBUG
  assert(skeleton->energy_check());
#endif

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
