#include "spacetime.h"

extern Random RND;

void Spacetime::mechanical_force(const std::vector<int>& offset,const std::vector<double>& y,double* output) const
{
  int i,j,l,m,k;
  double delta,r_true;
  std::set<int>::const_iterator it;
  const int nv = (signed) events.size();
  const int nreal = cardinality(0,-1);
  const int D = Geometry::background_dimension;
  const double pfactor = (2.0/M_PI)*5.0;
  double force[D*nreal];

#ifdef PARALLEL
#pragma omp parallel for default(shared) private(i,j,k,m,l,r_true,it,delta)
#endif 
  for(i=0; i<nv; ++i) {
    if (events[i].ubiquity == 1) continue;
    k = offset[i];
    // This vertex's force begins at zero...
    for(j=0; j<D; ++j) {
      force[D*k+j] = 0.0;
    }
    // Vertex-Vertex Repulsion...
    for(j=0; j<nv; ++j) {
      if (events[j].ubiquity == 1) continue;
      m = offset[j];
      if (k == m) continue;
      r_true = repulsion_constant/(1.0 + pfactor*std::atan(0.5*(events[i].energy + events[j].energy)));
      delta = 0.0;
      for(l=0; l<D; ++l) {
        delta += (y[D*k+l] - y[D*m+l])*(y[D*k+l] - y[D*m+l]);
      }
      for(l=0; l<D; ++l) {
        force[D*k+l] += r_true*(y[D*k+l] - y[D*m+l])/delta;
      }
    }
    // Edge-Edge Repulsion?
    // This means two edges, so four vertices would have their force 
    // modified 
    // Spring-like attraction between connected vertices...
    for(it=events[i].neighbours.begin(); it!=events[i].neighbours.end(); ++it) {
      j = offset[*it];
      for(l=0; l<D; ++l) {
        force[D*k+l] += spring_constant*(y[D*k+l] - y[D*j+l]); 
      } 
    }
  }
  // Position loop...
  for(i=0; i<D*nreal; ++i) {
    output[i] = y[i+D*nreal];
  }
  // Velocity loop...
  for(i=D*nreal; i<2*D*nreal; ++i) {
    output[i] = force[i-D*nreal] - damping_constant*y[i];
  }
}

void Spacetime::compute_delta()
{
  int i,j,m,n,l,nhop;
  std::set<int> vx,current,next;
  std::set<int>::const_iterator it,jt,kt;
  const int nv = (signed) events.size();
  const int nt = (signed) codex.size();
  int done[nv];

  for(i=0; i<nv; ++i) {
    done[i] = 0;
    events[i].topology_modified = false;
  }

  // All of the new vertices...
  for(i=0; i<nv; ++i) {
    if (events[i].incept == -1) vx.insert(i);
  }
  for(i=1; i<=ND; ++i) {
    n = (signed) simplices[i].size();
    // And all the new d-simplices (d >= 1)...
    for(j=0; j<n; ++j) {
      if (simplices[i][j].incept >= 0) continue;
      for(it=simplices[i][j].vertices.begin(); it!=simplices[i][j].vertices.end(); it++) {
        vx.insert(*it);
      }
    }
  }

  for(i=0; i<nt; ++i) {
    if (!codex[i].active) continue;
    for(it=codex[i].vx_delta.begin(); it!=codex[i].vx_delta.end(); it++) {
      vx.insert(*it);
    }
  }
  for(i=0; i<nt; ++i) {
    codex[i].vx_delta.clear();
  }
  // Now with this know we can calculate which events need to have their
  // entwinement and/or dimensional stress recalculated
#ifdef VERBOSE
  std::cout << "There are " << vx.size() << " vertices directly implicated" << std::endl;
#endif
  for(it=vx.begin(); it!=vx.end(); it++) {
    n = *it;
    nhop = 0;
    // Every vertex within topological_radius hops of n is labelled as
    // modified
    current.insert(n);
    done[n] = 1;
    do {
      for(jt=current.begin(); jt!=current.end(); jt++) {
        m = *jt;
        for(kt=events[m].neighbours.begin(); kt!=events[m].neighbours.end(); kt++) {
          l = *kt;
          if (done[l] == 0) next.insert(l);
        }
      }
      if (next.empty()) break;
      for(jt=next.begin(); jt!=next.end(); jt++) {
        done[*jt] = 1;
      }
      current = next;
      nhop++;
      next.clear();
    } while(nhop < Spacetime::topological_radius);
    current.clear();
    next.clear();
    for(i=0; i<nv; ++i) {
      if (done[i] == 1) events[i].topology_modified = true;
      done[i] = 0;
    }
  }
  int nmod = 0;
  for(i=0; i<nv; ++i) {
    if (events[i].topology_modified) nmod++;
  }
#ifdef VERBOSE
  std::cout << "There are " << nmod << " modified vertices out of " << nv << std::endl;
#endif
}

void Spacetime::energy_diffusion()
{
  const int nv = (signed) events.size();
  double Enew[nv],E,Z;
  int i,v,n,m,d1,d2;
  std::set<int>::const_iterator it;
  std::vector<int> vx,dimensions;
  std::vector<double> vdimension;
  std::vector<std::pair<int,double> > candidates;

  for(i=0; i<nv; ++i) {
    Enew[i] = -1.0;
    vdimension.push_back(double(events[i].global_dimension));
    if (events[i].energy < Spacetime::epsilon) {
      candidates.push_back(std::pair<int,double>(i,-1.0));
      continue;
    }
    if (events[i].ubiquity == 1) {
      candidates.push_back(std::pair<int,double>(i,-1.0));
    }
    else {
      candidates.push_back(std::pair<int,double>(i,std::abs(events[i].deficiency)));
    }
  }
  std::sort(candidates.begin(),candidates.end(),pair_predicate_dbl);
  for(i=nv-1; i>0; --i) {
    vx.clear();
    dimensions.clear();
    v = candidates[i].first;
    if (Enew[v] > 0.0) continue;
    E = events[v].energy;
    d1 = events[v].global_dimension;
    dimensions.push_back(d1);
    vx.push_back(v);
    // We need to reject neighbours if their dimensionality is different and/or
    // they are far away
    for(it=events[v].neighbours.begin(); it!=events[v].neighbours.end(); it++) {
      m = *it;
      if (Enew[m] > 0.0) continue;
      // How to handle the dimensionality issue?
      d2 = events[m].global_dimension;
      // 1) Forbid all inter-dimensional energy transfer...
      //if (d1 != d2) continue;
      // 2) Forbid it across certain "special" frontiers like 2/1...
      //if ((d1 > 1 && d2 == 1) || (d1 == 1 && d2 > 1)) continue;
      // 3) Forbid it according to a Boltzmann criterion...
      Z = 0.25*double(std::abs(d1 - d2));
      if (RND.drandom() > std::exp(-Z)) continue;
      // Now the geometric criterion: the farther away the vertex, the less
      // likely that it participates in energy transfer...
      Z = 0.5*(geometry->get_distance(v,m,false) - 1.0);
      if (RND.drandom() > std::exp(-Z)) continue;
      dimensions.push_back(d2);
      E += events[m].energy;
      vx.push_back(m);
    }
    if (vx.size() < 2) continue;
    m = (signed) vx.size();
    Z = 0.0;
    for(n=0; n<m; ++n) {
      Z += double(dimensions[n]);
    }
    for(n=0; n<m; ++n) {
      Enew[vx[n]] = double(dimensions[n])*E/Z;
    }
  }
  for(i=0; i<nv; ++i) {
    if (Enew[i] > 0.0) events[i].energy = Enew[i];
  }
}

bool Spacetime::adjust_dimension()
{
  // This method handles the geometry and energy changes to
  // minimize the structural deficiency...
  int i,n;
  bool modified;
  std::vector<int> vdimension;
  const int nvertex = (signed) events.size();

  system_size = 0;

  for(i=0; i<nvertex; ++i) {
    if (events[i].ubiquity == 1) {
      vdimension.push_back(-1);
      continue;
    }
    n = events[i].global_dimension;
#ifdef FLAT
    system_size += Geometry::background_dimension;
#else
    system_size += (n <= Geometry::background_dimension) ? Geometry::background_dimension : n;
#endif
    vdimension.push_back(n);
  }

  modified = geometry->adjust_dimension(vdimension);
  return modified;
}

void Spacetime::optimize()
{
  int i,n;
  bool deficient = false;
  std::vector<double> v1;
  const int nvertex = (signed) events.size();

  for(i=0; i<nvertex; ++i) {
    if (events[i].ubiquity == 1) continue;
    if (std::abs(events[i].deficiency) >= Spacetime::epsilon) {
      deficient = true;
      break;
    }
  }
  if (!deficient) return;
#ifdef VERBOSE
  std::cout << "Calling the geometric minimization routine with system dimension " << system_size << std::endl;
#endif
  if (solver == MINIMAL) {
    int its = 0;
    bool found;
    std::set<int> vmodified,candidates,bbarrel;
    std::set<int>::const_iterator it;
    double old_error,sigma;

    structural_deficiency();

    for(i=0; i<nvertex; ++i) {
      if (events[i].ubiquity == 1) continue;
      if (std::abs(events[i].geometric_deficiency) < Spacetime::epsilon) continue;
      bbarrel.insert(i);
      // Now check to see if all of this vertex's neighbours have a geometric deficiency > 0
      found = false;
      for(it=events[i].neighbours.begin(); it!=events[i].neighbours.end(); it++) {
        if (std::abs(events[*it].geometric_deficiency) < Spacetime::epsilon) {
          found = true;
          break;
        }
      }
      if (found) continue;
      candidates.insert(i);
    }
    if (candidates.empty() && bbarrel.empty()) return;
    if (candidates.empty()) candidates = bbarrel;
#ifdef VERBOSE
    std::cout << "There are " << candidates.size() << " candidate vertices for geometric optimization." << std::endl;
    std::cout << "In the geometry solver the initial error is " << error << std::endl;
#endif
    old_error = error;
    do {
      // Do something clever here to minimize the structural deficiency
      its++;
      // This is decidedly not very clever but it is quick!
      // First the geometry...
      n = RND.irandom(candidates);
      // In fact the variance we use should grow smaller as the geometric deficiency of the
      // vertex diminishes with 0.05 the maximum
      sigma = 0.1*events[n].geometric_deficiency;
      if (sigma > 0.05) sigma = 0.05;
      geometry->multiplicative_modification(n,true,1.0,sigma);
      vmodified.insert(n);
      compute_geometric_dependency(vmodified);
      geometry->compute_distances(vmodified);
      compute_volume();
      compute_curvature();
      compute_obliquity();

      // Finally, see if it's an improvement...
      structural_deficiency();
#ifdef VERBOSE
      std::cout << its << "  " << n << "  " << old_error << "  " << error << std::endl;
#endif
      vmodified.clear();
      if (error > old_error) {
        geometry->rollback();
        vmodified.insert(n);
        continue;
      }
      old_error = error;
    } while(its < solver_its);
  }
  else if (solver == EVOLUTIONARY) {
    int j,n,in1,vc,joust,generation = 1;
    bool viable,mtype = geometry->euclidean;
    double p,f1,f2,severity,fbest,fitness[2*pool_size],ftemp[pool_size];
    std::vector<std::pair<int,int> > vcount;
    std::set<int> vmodified,vx;
    std::set<int>::const_iterator it;
    Geometry* pool = new Geometry[2*pool_size];
    Geometry* ptemp = new Geometry[pool_size];
    Geometry* optimal = new Geometry(mtype);
    Geometry* initial_state = new Geometry(mtype);
    const int pmagnitude = int(0.15*system_size);
    const double initial_error = error;

#ifdef VERBOSE
    std::cout << "In evolutionary solver with initial error " << error << std::endl;
#endif

    for(i=0; i<2*pool_size; ++i) {
      pool[i].set_metric(mtype);
      vcount.push_back(std::pair<int,int>(0,0));
    }
    for(i=0; i<pool_size; ++i) {
      ptemp[i].set_metric(mtype);
    }

    geometry->store(initial_state);
    geometry->store(optimal);
    fbest = initial_error;

    // Initialize the population based on some Gaussian flux around the
    // current spacetime geometry
    for(i=0; i<pool_size; ++i) {
      geometry->store(&pool[i]);
      for(j=0; j<pmagnitude; ++j) {
        n = RND.irandom(system_size);
        pool[i].get_implied_vertices(n,vx);
        viable = false;
        for(it=vx.begin(); it!=vx.end(); ++it) {
          if (std::abs(events[*it].geometric_deficiency) > Spacetime::epsilon) viable = true;
        }
        if (!viable) continue;
        pool[i].geometry_modification(n,0.0,0.05);
        for(it=vx.begin(); it!=vx.end(); it++) {
          vmodified.insert(*it);
        }
        vx.clear();
      }
      geometry->load(&pool[i]);
      compute_geometric_dependency(vmodified);
      geometry->compute_distances(vmodified);
      compute_volume();
      compute_curvature();
      compute_obliquity();
      structural_deficiency();
      fitness[i] = error;
      vmodified.clear();
      geometry->load(initial_state);
    }

    // During each generation, a single individual creates one child
    // which is then mutated
    do {
      // Create the offspring
      p = 0.0;
      f1 = fitness[0];
      for(i=0; i<pool_size; ++i) {
        pool[pool_size+i] = pool[i];
        p += fitness[i];
        if (fitness[i] < f1) f1 = fitness[i];
        if (fitness[i] < fbest) {
          optimal = &pool[i];
          fbest = f1;
        }
      }
#ifdef VERBOSE
      std::cout << "At generation " << generation << " the average population fitness is " << p/double(pool_size) << " and minimum fitness is " << f1 << " with global optimum at " << fbest << std::endl;
#endif
      // Now the mutations
      RND.initialize_beta(1.0,double(generation)/10.0);
      for(i=pool_size; i<2*pool_size; ++i) {
        severity = RND.beta_variate();
        for(j=0; j<system_size; ++j) {
          if (RND.drandom() >= severity) continue;
          pool[i].get_implied_vertices(j,vx);
          viable = false;
          for(it=vx.begin(); it!=vx.end(); ++it) {
            if (std::abs(events[*it].geometric_deficiency) > Spacetime::epsilon) viable = true;
          }
          if (!viable) continue;
          pool[i].geometry_modification(j,0.0,severity);
          for(it=vx.begin(); it!=vx.end(); it++) {
            vmodified.insert(*it);
          }
          vx.clear();
        }
        geometry->load(&pool[i]);
        compute_geometric_dependency(vmodified);
        geometry->compute_distances(vmodified);
        compute_volume();
        compute_curvature();
        compute_obliquity();
        structural_deficiency();
        fitness[i] = error;
        vmodified.clear();
      }
      // Now we need to see among the total population (2*npop) which are the
      // top npop individuals via a stochastic tournament
      for(i=0; i<2*pool_size; ++i) {
        joust = 0;
        f1 = fitness[i];
        vc = 0;
        do {
          in1 = RND.irandom(2*pool_size);
          if (in1 == i) continue;
          f2 = fitness[in1];
          // Calculate the probability that f1 emerges triumphant...
          p = 0.5*(1.0 + std::tanh(10.0*(f2 - f1)));
          if (p > RND.drandom()) vc++;
          joust += 1;
        } while(joust < njousts);
        vcount[i] = std::pair<int,int>(i,vc);
      }
      std::sort(vcount.begin(),vcount.end(),pair_predicate_int);
      for(i=0; i<pool_size; ++i) {
        ptemp[i] = pool[vcount[i].first];
        ftemp[i] = fitness[vcount[i].first];
      }
      for(i=0; i<pool_size; ++i) {
        pool[i] = ptemp[i];
        fitness[i] = ftemp[i];
      }
      generation += 1;
    } while(generation <= ngenerations);
    // We need to pick the best individual in the pool in order to
    // write its "genome" into the spacetime geometry...
    // With the above algorithm, the fittest individual should always
    // be the first element in the array "pool", in the sense that this
    // individual won the most combats in the tournament
    if (fbest < initial_error) {
#ifdef VERBOSE
      std::cout << "Loading optimized geometry with " << 100.0*(initial_error - fbest)/initial_error << " percent improvement" << std::endl;
#endif
      geometry->load(optimal);
    }
    else {
#ifdef VERBOSE
      std::cout << "Geometry optimization failed, reverting to original geometry." << std::endl;
#endif
      geometry->load(initial_state);
    }
    geometry->compute_distances();
    compute_volume();
    compute_curvature();
    compute_obliquity();
    structural_deficiency();
    delete optimal;
    delete initial_state;
    delete[] ptemp;
    delete[] pool;
  }
  else if (solver == ANNEALING) {
    int j,m,step,naccept,nim,nwo_a,nwo_r;
    double S,q,t,E_old,paccept,sigma,E_avg,E_best,delta,temperature = 100.0;
    bool done,viable;
    std::set<int> vmodified;
    std::set<int>::const_iterator it;
    std::vector<double> E,lengths,base,output,output1,output2;
    Geometry* optimal = new Geometry(geometry->euclidean);
    const double initial_error = error;

    geometry->store(optimal);

    // The goal here should be to find the temperature at which some 80% of
    // random changes are accepted - this will be the melting point of the
    // system
    do {
      naccept = 0;
      E_old = initial_error;
      for(i=0; i<1000; ++i) {
        do {
          m = RND.irandom(system_size);
          geometry->get_implied_vertices(m,vmodified);
          viable = false;
          for(it=vmodified.begin(); it!=vmodified.end(); ++it) {
            if (std::abs(events[*it].geometric_deficiency) > Spacetime::epsilon) viable = true;
          }
          if (viable) break;
          vmodified.clear();
        } while(true);
        geometry->geometry_modification(m,0.0,thermal_variance);
        compute_geometric_dependency(vmodified);
        geometry->compute_distances(vmodified);
        compute_volume();
        compute_curvature();
        compute_obliquity();
        structural_deficiency();
        if (error < E_old) {
          naccept++;
          E_old = error;
        }
        else {
          q = RND.drandom();
          sigma = error - E_old;
          if (q < std::exp(-sigma/temperature)) {
            naccept++;
            E_old = error;
          }
          else {
            geometry->geometry_restoration();
            compute_geometric_dependency(vmodified);
            geometry->compute_distances(vmodified);
            compute_volume();
            compute_curvature();
            compute_obliquity();
            structural_deficiency();
          }
        }
        vmodified.clear();
      }
      paccept = double(naccept)/1000.0;
#ifdef VERBOSE
      std::cout << "At temperature " << temperature << ", the acceptance percentage is " << 100.0*paccept << std::endl;
#endif
      geometry->load(optimal);
      if (paccept > 0.85) {
        temperature *= 0.9;
      }
      else if (paccept < 0.8) {
        temperature *= 1.1;
      }
      else {
        break;
      }
    } while(true);
#ifdef VERBOSE
    std::cout << "Beginning annealing process with initial temperature of " << temperature << std::endl;
#endif
    E_best = initial_error;
    // The main annealing loop...
    for(i=0; i<annealing_steps; ++i) {
      done = false;
      step = 1;
      do {
        if (step == 1) {
          geometry->compute_distances();
          compute_volume();
          compute_curvature();
          compute_obliquity();
          structural_deficiency();
          E_old = error;
          output = base;
        }
        nim = 0;
        nwo_a = 0;
        nwo_r = 0;
        for(j=0; j<thermal_sweep; ++j) {
          do {
            n = RND.irandom(system_size);
            geometry->get_implied_vertices(n,vmodified);
            viable = false;
            for(it=vmodified.begin(); it!=vmodified.end(); ++it) {
              if (std::abs(events[*it].geometric_deficiency) > Spacetime::epsilon) viable = true;
            }
            if (viable) break;
            vmodified.clear();
          } while(true);
          geometry->geometry_modification(n,0.0,thermal_variance);
          // Now we need to recaculate some edge lengths and
          // then various defect angles associated with the
          // triangles in their entourage...
          compute_geometric_dependency(vmodified);
          geometry->compute_distances(vmodified);
          compute_volume();
          compute_curvature();
          compute_obliquity();
          structural_deficiency();
          q = error;
          if (q < E_old) {
            E_old = q;
            E.push_back(q);
            if (q < E_best) {
              E_best = q;
              geometry->store(optimal);
            }
            nim++;
          }
          else {
            t = std::exp((E_old - q)/temperature);
            sigma = RND.drandom();
            if (sigma < t) {
              E_old = q;
              E.push_back(q);
              nwo_a++;
            }
            else {
              E.push_back(E_old);
              geometry->geometry_restoration();
              compute_geometric_dependency(vmodified);
              geometry->compute_distances(vmodified);
              compute_volume();
              compute_curvature();
              compute_obliquity();
              structural_deficiency();
              nwo_r++;
            }
          }
          vmodified.clear();
        }
        E_avg = 0.0;
        for(j=0; j<step*thermal_sweep; ++j) {
          E_avg += E[j];
        }
        E_avg /= double(step*thermal_sweep);
        sigma = 0.0;
        for(j=0; j<step*thermal_sweep; ++j) {
          delta = E[j] - E_avg;
          sigma += delta*delta;
        }
        S = std::sqrt(sigma/double(step*thermal_sweep));
        if (S < thermalization) done = true;
#ifdef VERBOSE
        std::cout << step << "  " << E_best << "  " << nim << "  " << nwo_a << "  " << nwo_r << "  " << double(nim+nwo_a)/double(thermal_sweep) << "  " << E_avg << "  " << S << std::endl;
#endif
        step += 1;
      } while(!done);
      E.clear();
      temperature *= 0.9;
    }
#ifdef VERBOSE
    if (E_best < initial_error) {
      std::cout << "Loading optimized geometry with " << 100.0*(initial_error - E_best)/initial_error << " percent improvement" << std::endl;
    }
    else {
      std::cout << "Geometry optimization failed, reverting to original geometry." << std::endl;
    }
#endif
    geometry->load(optimal);
    geometry->compute_distances();
    compute_volume();
    compute_curvature();
    compute_obliquity();
    structural_deficiency();
    delete optimal;
  }
  else if (solver == MECHANICAL) {
    // Method to use an effective force algorithm to calculate the best 
    // geometric configuration for the spacetime complex, with the possibility 
    // of a final conjugate gradient optimization stage to attempt to reduce 
    // the geometry's abnormality. 
    int j,m = 0,k = 0,its = 1;
    double vnorm;
    std::vector<int> offset;
    std::vector<double> y,ynew;
    const int nreal = cardinality(0,-1);
    const int D = Geometry::background_dimension;

#ifdef VERBOSE
    std::cout << "Using effective force geometry solver with " << nreal << " active vertices and background dimension = " << D << "." << std::endl;
#endif

    for(i=0; i<nvertex; ++i) {
      if (events[i].ubiquity == 1) continue;
      offset.push_back(m);
      m++;
    }
    // Initial position from current geometry...
    for(i=0; i<D*nreal; ++i) {
      y.push_back(geometry->get_element(i));
    }

    // Initial velocity is zero...
    for(i=D*nreal; i<2*D*nreal; ++i) {
      y.push_back(0.0); 
    }
    for(i=0; i<2*D*nreal; ++i) {
      ynew.push_back(0.0);
    }
    if (int_engine == "EULER") {
      double F[2*D*nreal];
      do {
        mechanical_force(offset,y,F);
        // This just uses the Euler method, not very sophisticated...
        for(i=0; i<2*D*nreal; ++i) {
          ynew[i] = y[i] + step_size*F[i];
        }
        vnorm = 0.0;
        for(i=D*nreal; i<2*D*nreal; ++i) {
          vnorm += ynew[i]*ynew[i];
        }
        vnorm = std::sqrt(vnorm);
        if (vnorm < geometry_cutoff || its > max_int_steps) break;
        y = ynew; 
        its++;
      } while(true);
    }
    else if (int_engine == "RK4") {
      double k1[2*D*nreal],k2[2*D*nreal],k3[2*D*nreal],k4[2*D*nreal];
      std::vector<double> temp;
      for(i=0; i<2*D*nreal; ++i) {
        temp.push_back(0.0);
      }
      do {
        mechanical_force(offset,y,k1);
        for(i=0; i<2*D*nreal; ++i) {
          temp[i] = y[i] + 0.5*step_size*k1[i];
        }
        mechanical_force(offset,temp,k2);
        for(i=0; i<2*D*nreal; ++i) {
          temp[i] = y[i] + 0.5*step_size*k2[i];
        }
        mechanical_force(offset,temp,k3);
        for(i=0; i<2*D*nreal; ++i) {
          temp[i] = y[i] + step_size*k3[i];
        }
        mechanical_force(offset,temp,k4);    
        for(i=0; i<2*D*nreal; ++i) {
          ynew[i] = y[i] + step_size/6.0*(k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]); 
        }
        vnorm = 0.0;
        for(i=D*nreal; i<2*D*nreal; ++i) {
          vnorm += ynew[i]*ynew[i];
        }
        vnorm = std::sqrt(vnorm);
        if (vnorm < geometry_cutoff || its > max_int_steps) break;
        y = ynew; 
        its++;
      } while(true);
    }
    for(i=0; i<D*nreal; ++i) {
      geometry->set_element(i,ynew[i]);
    }
    geometry->compute_distances();
    compute_volume();
    if (!cgradient_refinement) return;
    double d,q,E,prior,alpha,beta,E_initial,n = 0.0,sigma = 0.1;
    Geometry* initial_state = new Geometry(geometry->euclidean);
    std::vector<double> s,snew,dx,dy,dx_old,x,c,fx;

    determine_flexible_edges();
    E_initial = compute_abnormality();
#ifdef VERBOSE
    std::cout << "Initial error is " << E_initial << std::endl;
#endif
    geometry->store(initial_state);
    for(i=0; i<system_size; ++i) {
      x.push_back(geometry->get_element(i));
      snew.push_back(0.0);
    }

    compute_geometric_gradient(dx,true);
    for(i=0; i<system_size; ++i) {
      fx.push_back(-dx[i]);
      c.push_back(dx[i]*fx[i]);
      n += c[i];
    }
    for(i=0; i<max_LS_steps; ++i) {
      for(j=0; j<system_size; ++j) {
        geometry->set_element(j,x[j] + sigma*dx[j]);
      }
      geometry->compute_distances();
      compute_geometric_gradient(dy,false);
      d = 0.0;
      for(j=0; j<system_size; ++j) {
        d += dy[j]*dx[j] - c[j];
      }
      q = n/d;
      alpha = -sigma*q;
      sigma = -alpha;
    }
    for(i=0; i<system_size; ++i) {
      geometry->set_element(i,x[i] + alpha*dx[i]);
    }
    geometry->compute_distances();
    E = compute_abnormality();

    prior = E;
    s = dx;
    dx_old = dx;
    for(i=0; i<system_size; ++i) {
      x[i] = geometry->get_element(i);
    }
#ifdef VERBOSE
    std::cout << "At conjugate gradient iteration 1, the error is " << E << std::endl;
#endif
    for(i=1; i<max_CG_steps; ++i) {
      compute_geometric_gradient(dx,true);
      n = 0.0;
      d = 0.0;
      for(j=0; j<system_size; ++j) {
        n += dx[j]*(dx[j] - dx_old[j]);
        d += dx_old[j]*dx_old[j];
      }
      beta = std::max(0.0,n/d);
      for(j=0; j<system_size; ++j) {
        snew[j] = dx[j] + beta*s[j];
      }
      // Now begin the line search...
      n = 0.0;
      for(j=0; j<system_size; ++j) {
        fx[j] = -dx[j];
        c[j] = snew[j]*fx[j];
        n += c[j];
      }
      sigma = 0.1;
      for(j=0; j<max_LS_steps; ++j) {
        for(k=0; k<system_size; ++k) {
          geometry->set_element(k,x[k] + sigma*snew[k]);
        }
        geometry->compute_distances();
        compute_geometric_gradient(dy,false);
        d = 0.0;
        for(k=0; k<system_size; ++k) {
          d += dy[k]*snew[k] - c[k];
        }
        q = n/d;
        alpha = -sigma*q;
        sigma = -alpha;
      }
      for(j=0; j<system_size; ++j) {
        geometry->set_element(j,x[j] + alpha*snew[j]);
      }
      geometry->compute_distances();
      E = compute_abnormality();
#ifdef VERBOSE
      std::cout << "At conjugate gradient iteration " << i+1 << ", the error is " << E << std::endl;
#endif
      if (E > prior) {
        for(j=0; j<system_size; ++j) {
          geometry->set_element(j,x[j]);
        }
        geometry->compute_distances();
        break;
      }
      // Update the base values...
      prior = E;
      s = snew;
      dx_old = dx;
      for(j=0; j<system_size; ++j) {
        x[j] = geometry->get_element(j);
      }
      if (E < geometry_cutoff) break;
    }
    if (E > E_initial) {
#ifdef VERBOSE
      std::cout << "Geometry optimization failed, reverting to original geometry." << std::endl;
#endif
      geometry->load(initial_state);
      geometry->compute_distances();
      for(i=0; i<nvertex; ++i) {
        if (events[i].ubiquity == 1) continue;
        events[i].geometry_modified = true;
      }
      compute_volume();
    }
    delete initial_state;
  }
  else if (solver == SIMPLEX) {
    int in1,j,k,bindex,windex,ntrans = 0;
    double f,q,centroid[system_size];
    Geometry SR,SE;
    Geometry* initial_state = new Geometry(geometry->euclidean);
    std::vector<std::pair<int,double> > fitness;
    std::set<int> vmodified;
    std::vector<Geometry> S;
    const double initial_error = error;

#ifdef VERBOSE
    std::cout << "Using Nelder-Mead method with dimension " << system_size << std::endl;
#endif
    geometry->store(initial_state);
    SR.set_metric(geometry->euclidean);
    SE.set_metric(geometry->euclidean);

    // We begin by creating a simplex of size system_size....
    geometry->store(&SR);
    S.push_back(SR);
    fitness.push_back(std::pair<int,double>(0,error));
    for(i=0; i<system_size; ++i) {
      geometry->get_implied_vertices(i,vmodified);
      compute_geometric_dependency(vmodified);
      geometry->geometry_modification(i,0.0,0.1);
      geometry->compute_distances(vmodified);
      geometry->store(&SR);
      S.push_back(SR);
      geometry->geometry_restoration();
      vmodified.clear();
      for(j=0; j<nvertex; ++j) {
        events[j].geometry_modified = false;
      }
      for(j=1; j<=ND; ++j) {
        for(k=0; k<(signed) simplices[j].size(); ++k) {
          simplices[j][k].modified = false;
        }
      }
    }
    // Compute the error for the new vertices...
    for(i=1; i<1+system_size; ++i) {
      geometry->load(&S[i]);
      geometry->compute_distances();
      compute_volume();
      compute_curvature();
      compute_obliquity();
      structural_deficiency();
      fitness.push_back(std::pair<int,double>(i,error));
    }
    // Sort the simplex vertices from lowest to highest error...
    std::sort(fitness.begin(),fitness.end(),pair_predicate_dbl);
    do {
      ntrans++;
      bindex = fitness[0].first;
      windex = fitness[system_size].first;
      // Compute the centroid of the simplex...
      for(i=0; i<system_size; ++i) {
        centroid[i] = 0.0;
      }
      for(i=0; i<system_size; ++i) {
        for(j=0; j<system_size; ++j) {
          centroid[j] += S[i].get_element(j);
        }
      }
      for(i=0; i<system_size; ++i) {
        centroid[i] /= double(system_size);
      }
      // Now calculate the reflection...
      for(i=0; i<system_size; ++i) {
        f = (1.0 + simplex_alpha)*centroid[k] - simplex_alpha*S[windex].get_element(i);
        SR.set_element(i,f);
      }
      geometry->load(&SR);
      geometry->compute_distances();
      for(i=0; i<nvertex; ++i) {
        if (events[i].ubiquity > 1) events[i].geometry_modified = true;
      }
      compute_volume();
      compute_curvature();
      compute_obliquity();
      structural_deficiency();
      f = error;
      if (f < fitness[0].second) {
        // Try expansion in this case...
        SE = SR;
        for(i=0; i<system_size; ++i) {
          q = (1.0 + simplex_gamma*simplex_alpha)*centroid[i] - simplex_gamma*simplex_alpha*S[windex].get_element(i);
          SE.set_element(i,q);
        }
        geometry->load(&SE);
        geometry->compute_distances();
        for(i=0; i<nvertex; ++i) {
          if (events[i].ubiquity > 1) events[i].geometry_modified = true;
        }
        compute_volume();
        compute_curvature();
        compute_obliquity();
        structural_deficiency();
        if (error < f) {
          S[windex] = SE;
          fitness[system_size].second = error;
        }
        else {
          S[windex] = SR;
          fitness[system_size].second = f;
        }
      }
      else if (f >= fitness[0].second && f < fitness[system_size-1].second) {
        S[windex] = SR;
        fitness[system_size].second = f;
      }
      else if (f >= fitness[system_size-1].second && f < fitness[system_size].second) {
        // Try external contraction...
        SE = SR;
        for(i=0; i<system_size; ++i) {
          q = (1.0 + simplex_alpha*simplex_rho)*centroid[i] - simplex_alpha*simplex_rho*S[windex].get_element(i);
          SE.set_element(i,q);
        }
        geometry->load(&SE);
        geometry->compute_distances();
        for(i=0; i<nvertex; ++i) {
          if (events[i].ubiquity > 1) events[i].geometry_modified = true;
        }
        compute_volume();
        compute_curvature();
        compute_obliquity();
        structural_deficiency();
        if (error < f) {
          S[windex] = SE;
          fitness[system_size].second = error;
        }
        else {
          for(i=1; i<=system_size; ++i) {
            in1 = fitness[i].first;
            for(j=0; j<system_size; ++j) {
              q = S[bindex].get_element(j) + simplex_sigma*(S[in1].get_element(j) - S[bindex].get_element(j));
              S[in1].set_element(j,q);
            }
            geometry->load(&S[in1]);
            geometry->compute_distances();
            for(j=0; j<nvertex; ++j) {
              if (events[j].ubiquity > 1) events[j].geometry_modified = true;
            }
            compute_volume();
            compute_curvature();
            compute_obliquity();
            structural_deficiency();
            fitness[i].second = error;
          }
        }
      }
      else {
        // Try internal contraction...
        SE = SR;
        for(i=0; i<system_size; ++i) {
          q = (1.0 - simplex_rho)*centroid[i] + simplex_rho*S[windex].get_element(i);
          SE.set_element(i,q);
        }
        geometry->load(&SE);
        geometry->compute_distances();
        for(i=0; i<nvertex; ++i) {
          if (events[i].ubiquity > 1) events[i].geometry_modified = true;
        }
        compute_volume();
        compute_curvature();
        compute_obliquity();
        structural_deficiency();
        if (error < fitness[system_size].second) {
          S[windex] = SE;
          fitness[system_size].second = error;
        }
        else {
          for(i=1; i<=system_size; ++i) {
            in1 = fitness[i].first;
            for(j=0; j<system_size; ++j) {
              q = S[bindex].get_element(j) + simplex_sigma*(S[in1].get_element(j) - S[bindex].get_element(j));
              S[in1].set_element(j,q);
            }
            geometry->load(&S[in1]);
            geometry->compute_distances();
            for(j=0; j<nvertex; ++j) {
              if (events[j].ubiquity > 1) events[j].geometry_modified = true;
            }
            compute_volume();
            compute_curvature();
            compute_obliquity();
            structural_deficiency();
            fitness[i].second = error;
          }
        }
      }
      // Now resort the simplex vertices...
      std::sort(fitness.begin(),fitness.end(),pair_predicate_dbl);
#ifdef VERBOSE
      std::cout << "At simplex transformation step " << ntrans << " the error is " << fitness[0].second << std::endl;
#endif
    } while(ntrans <= 10*system_size && fitness[0].second > geometry_cutoff);
    if (fitness[0].first < initial_error) {
#ifdef VERBOSE
      std::cout << "Loading optimized geometry with " << 100.0*(initial_error - fitness[0].first)/initial_error << " percent improvement" << std::endl;
#endif
      geometry->load(&S[fitness[0].first]);
    }
    else {
#ifdef VERBOSE
      std::cout << "Geometry optimization failed, reverting to original geometry." << std::endl;
#endif
      geometry->load(initial_state);
    }
    geometry->compute_distances();
    compute_volume();
    compute_curvature();
    compute_obliquity();
    structural_deficiency();
    delete initial_state;
  }
}





