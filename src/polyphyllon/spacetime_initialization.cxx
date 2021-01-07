#include "spacetime.h"

using namespace DIAPLEXIS;

void Spacetime::build_initial_state(const std::set<int>& locus)
{
  int i,j,m,in1;
  std::string geometry_type;
  Event vt;
  const int D = geometry->dimension();
  const bool relational = geometry->get_relational();
  std::vector<double> svalue;

  for(m=0; m<D; ++m) {
    svalue.push_back(0.0);
  }
  vt.set_incept(0);
  vt.set_ubiquity(locus);

  if (initial_state == Initial_Topology::cartesian) {
    std::vector<std::pair<long,int> > factors;
    geometry_type = "CARTESIAN";
    SYNARMOSMA::factorize(initial_size,factors);
    j = 1;
    for(m=0; m<(signed) factors.size(); ++m) {
      j *= SYNARMOSMA::ipow(factors[m].first,factors[m].second/D);
    }
    const int n = j;
    const int nd = SYNARMOSMA::ipow(n,D-1);
    const int nm1 = n - 1;
    const int nperturbed = 10 + int(0.01*skeleton->RND->irandom(initial_size));
    const double dx = 1.0;
    int k,l,r,rvalue;
    SYNARMOSMA::hash_map::const_iterator qt;
    std::vector<int> entourage,v;
    std::vector<int>* arrangement = new std::vector<int>[initial_size];
    std::set<int> N;
    std::set<int>::const_iterator it,jt;
    std::vector<double> vpoints;

    for(m=0; m<D; ++m) {
      entourage.push_back(0);
    }
    
    // First the events, making sure (if this is a non-relational
    // model) that the event coordinates are balanced...
    for(l=0; l<initial_size; ++l) {
      k = nd;
      in1 = l/k;
      entourage[0] = in1;
      rvalue = l - in1*nd;
      vpoints.push_back(dx*(double(in1) - 0.5*double(nm1)));
      for(m=1; m<D; ++m) {
        k /= n;
        in1 = rvalue/k;
        entourage[m] = in1;
        rvalue -= k*in1;
        vpoints.push_back(dx*(double(in1) - 0.5*double(nm1)));
      }
      for(m=0; m<D; ++m) {
        if (entourage[m] == 0 || entourage[m] == nm1) vt.set_boundary(true);
      }
      skeleton->events.push_back(vt);
      vt.set_boundary(false);
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
        for(m=0; m<D-1; ++m) {
          k /= n;
          in1 += v[m+1]*k;
        }
        if (in1 > i) skeleton->simplex_addition(i,in1,locus,0);
      }
      if (v[0] < nm1) {
        in1 = (v[0]+1)*nd;
        k = nd;
        for(m=0; m<D-1; ++m) {
          k /= n;
          in1 += v[m+1]*k;
        }
        if (in1 > i) skeleton->simplex_addition(i,in1,locus,0);
      }
      for(m=1; m<D; ++m) {
        if (v[m] > 0) {
          v[m] -= 1;
          in1 = v[0]*nd;
          k = nd;
          for(r=0; r<D-1; ++r) {
            k /= n;
            in1 += v[r+1]*k;
          }
          if (in1 > i) skeleton->simplex_addition(i,in1,locus,0);
          v[m] += 1;
        }
        if (v[m] < nm1) {
          v[m] += 1;
          in1 = v[0]*nd;
          k = nd;
          for(r=0; r<D-1; ++r) {
            k /= n;
            in1 += v[r+1]*k;
          }
          if (in1 > i) skeleton->simplex_addition(i,in1,locus,0);
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
        for(m=0; m<D; ++m) {
          j += SYNARMOSMA::ipow(n,D-1-m)*(skeleton->RND->irandom(k-2,k+2));
        }
        if (j < 0) j = skeleton->RND->irandom(initial_size,v);
        v.push_back(j);
        // Now add some edges among the neighbours of this
        // event
        N.insert(v[i]);
        j = skeleton->RND->irandom(2,5);
        for(l=0; l<j; ++l) {
          N.insert(skeleton->events.size());
          v.push_back(skeleton->events.size());
          geometry->vertex_addition(v[i]);
          skeleton->events.push_back(vt);
        }
        skeleton->simplex_addition(N,locus,0);
        N.clear();
      }
    }
    if (perturb_geometry) {
      // No changes to the spacetime topology, merely the geometry of some
      // existing events are altered...
      if (v.empty()) {  
        for(i=0; i<nperturbed; ++i) {
          j = skeleton->RND->irandom(initial_size,v);
          v.push_back(j);
        }
      }
      for(m=0; m<(signed) v.size(); ++m) {
        geometry->mutation(v[m],true,false,0.5);
      }
    }
    if (perturb_energy) {
      if (v.empty()) {
        // In this particular case I will simply alter the energy value of a single event
        // near the centre of the Cartesian network...
        j = 0;
        k = int(double(n)/2.0);
        for(m=0; m<D; ++m) {
          j += SYNARMOSMA::ipow(n,D-1-m)*skeleton->RND->irandom(k-2,k+2);
        }
        if (j < 0) j = skeleton->RND->irandom(initial_size);
        skeleton->events[j].set_energy(1000.0*(0.5 + skeleton->RND->drandom()/2.0));
      }
      else {
        for(m=0; m<(signed) v.size(); ++m) {
          skeleton->events[v[m]].set_energy(5.0*skeleton->RND->drandom());
        }
      }
    }
  }
  else if (initial_state == Initial_Topology::singleton) {
    // An initial spacetime consisting of a single isolated event, though with very
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
    int ulimit = D;
    if (initial_dim > D) ulimit = initial_dim;
    if (perturb_energy) vt.set_energy(500.0 + (2000.0/double(initial_dim))*skeleton->RND->drandom());

    if (!relational) {
      svalue.clear();
      for(m=0; m<ulimit; ++m) {
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

    skeleton->simplex_addition(vx,locus,0);
  }
  else if (initial_state == Initial_Topology::random) {
    // We will use the Erdős–Rényi random graph model (the G(n,p) variant) to
    // assemble a random graph, with n = initial_size
    int n = 0,level = 2,k = 0,ulimit;
    double percent = 0.0;
    bool found;
    std::set<int>* N;
    std::vector<std::set<int> > svector;
    SYNARMOSMA::hash_map::const_iterator qt;
    std::set<int> v,vx,current;
    std::set<int>::const_iterator it,chk;
    const double nv = double(initial_size);

    geometry_type = "RANDOM";
    // Add the events
    for(i=0; i<initial_size; ++i) {
      skeleton->events.push_back(vt);
    }

    // Next distribute energy among the events, ensuring at least one event
    // has non-zero energy
    do {
      i = skeleton->RND->irandom(initial_size);
      if (!skeleton->events[i].zero_energy()) continue;
      skeleton->events[i].set_energy(0.01 + 10.0*skeleton->RND->drandom());
      k += 1;
      percent = double(k)/nv;
    } while(percent < 0.1);

    // Next initialize the geometry
    if (relational) {
      geometry->create(initial_size,geometry_type);
    }
    else {
      std::vector<double> climits;
      for(m=0; m<D; ++m) {
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
        skeleton->simplex_addition(i,j,locus,0); n++;
      }
    }
#ifdef VERBOSE
    std::cout << "Added " << n << " edges among the " << initial_size << " events..." << std::endl;
#endif
    // Finally we have to compute all the n-simplices (n > 1) that are implied
    // by this random graph...
    do {
      N = new std::set<int>[level];
      ulimit = (signed) skeleton->simplices[level-1].size();
      for(j=0; j<ulimit; ++j) {
        skeleton->simplices[level-1][j].get_vertices(v);
        i = 0;
        for(it=v.begin(); it!=v.end(); ++it) {
          skeleton->events[*it].get_neighbours(N[i]);
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
          svector.push_back(current);
        }
        vx.clear();
        v.clear();
      }
      delete[] N;
      if (svector.empty()) break;
      for(j=0; j<(signed) svector.size(); ++j) {
        skeleton->simplex_addition(svector[j],locus,0);
      }
      svector.clear();
#ifdef VERBOSE
      if (!diskless) std::cout << "Finished the implied " << level << "-simplices..." << std::endl;
#endif
      level += 1;
    } while(true);
  }
 
  for(i=0; i<(signed) codex.size(); ++i) {
    if (locus.count(i) > 0) regularization(false,i);
  }
  regularization(false,-1);
#ifdef DEBUG
  if (initial_state == Initial_Topology::singleton) {
    assert(skeleton->dimension(-1) == 0);
  }
  else if (initial_state == Initial_Topology::monoplex) {
    assert(skeleton->dimension(-1) == initial_dim);
  }
#endif
}

void Spacetime::initialize()
{
  int i;
  std::stringstream day,month,year,pid;
#ifdef VERBOSE
  auto t1 = std::chrono::steady_clock::now();
#endif

  if (!diskless) {
    if (std::system("mkdir -p data") < 0) throw std::runtime_error("Unable to create data directory!");
  }

  // Get the date and the process ID so that we can construct the
  // state_file and log_file names...
  start_time = std::time(NULL); 
  tm* ltime = localtime(&start_time);

  day.width(2);
  day << std::setfill('0') << ltime->tm_mday;
  month.width(2);
  month << std::setfill('0') << 1 + ltime->tm_mon;
  year.width(4);
  year << 1900 + ltime->tm_year;
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
    std::set<int> locus;

    nactive = nt_initial;
    for(i=0; i<nt_initial; ++i) {
      codex.push_back(Sheet(i,skeleton->get_homology_field(),skeleton->get_homology_method()));
      locus.insert(i);
    }
    if (weaving == Hyphansis::musical) {
      std::set<int>::const_iterator it;

      i = 0;
      for(it=voices.begin(); it!=voices.end(); ++it) {
        codex[i].set_voice(*it);
        codex[i].parse_music_score(max_iter,hyphansis_score);
        i++;
      }
    }
    build_initial_state(locus);
    geometry->compute_squared_distances();
    skeleton->compute_simplicial_dimension();
    adjust_dimension();
#ifdef DEBUG
    assert(skeleton->consistent(-1));
#endif
    skeleton->compute_parity();
    compute_lightcones();
    compute_volume();
    compute_obliquity();
    compute_lightcones();
    compute_global_topology(-1);
    for(i=0; i<nt_initial; ++i) {
      codex[i].set_topology(skeleton->H,skeleton->pi1,skeleton->pseudomanifold,skeleton->boundary,skeleton->orientable);  
    }
    structural_deficiency();
  }
  skeleton->RND->initialize_poisson(Spacetime::ramosity);
#ifdef DEBUG
  assert(consistent());
  assert(skeleton->energy_check());
#endif

  if (weaving == Hyphansis::musical) {
    for(int i=0; i<60; ++i) {
      key_mapping[i] = -1;
    }
    key_mapping[21] = 0;
    key_mapping[23] = 1;
    key_mapping[25] = 2;
    key_mapping[26] = 3;
    key_mapping[28] = 4;
    key_mapping[29] = 5;
    key_mapping[30] = 6;
    key_mapping[32] = 7;
    key_mapping[33] = 8;
    key_mapping[34] = 9;
    key_mapping[35] = 10;
    key_mapping[37] = 11;
    key_mapping[40] = 12;
    key_mapping[45] = 13;
    key_mapping[47] = 14;
    key_mapping[48] = 15;
    key_mapping[49] = 16;
    key_mapping[50] = 17;
    key_mapping[52] = 18;
    key_mapping[53] = 19;
    key_mapping[54] = 20;
    key_mapping[56] = 21;
    key_mapping[57] = 22;
    key_mapping[58] = 23;
    key_mapping[59] = 24;
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
  skeleton->write_topology(-1);
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
  std::cout << "Sheet activity " << sheet_activity() << std::endl;
#endif

#ifdef VERBOSE
  auto t2 = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed_seconds = t2 - t1;
  std::cout << "Spacetime initialization required " << elapsed_seconds.count() << " seconds to complete." << std::endl;
#endif
}