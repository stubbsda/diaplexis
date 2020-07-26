#include "spacetime.h"

using namespace DIAPLEXIS;

void Spacetime::mechanical_force(const std::vector<int>& offset,const std::vector<int>& vx_map,const std::vector<double>& energy,const std::vector<double>& y,double* output) const
{
  int i,j,k;
  double delta,r_true;
  std::set<int> S;
  std::set<int>::const_iterator it;
  const int nreal = skeleton->cardinality(0,-1);
  const int D = geometry->dimension();
  const double pfactor = (2.0/M_PI)*5.0;
  double force[D*nreal];

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i,j,k,delta,r_true,S,it)
#endif 
  for(i=0; i<nreal; ++i) {
    // This event's force begins at zero...
    for(j=0; j<D; ++j) {
      force[D*i+j] = 0.0;
    }
    // Event-Event Repulsion...
    for(j=0; j<nreal; ++j) {
      if (i == j) continue;
      r_true = repulsion_constant/(1.0 + pfactor*std::atan(0.5*(energy[i] + energy[j])));
      delta = 0.0;
      for(k=0; k<D; ++k) {
        delta += (y[D*i+k] - y[D*j+k])*(y[D*i+k] - y[D*j+k]);
      }
      for(k=0; k<D; ++k) {
        force[D*i+k] += r_true*(y[D*i+k] - y[D*j+k])/delta;
      }
    }
    // Edge-Edge Repulsion?
    // This means two edges, so four events would have their force 
    // modified 
    // Spring-like attraction between connected events...
    skeleton->events[offset[i]].get_neighbours(S);
    for(it=S.begin(); it!=S.end(); ++it) {
      j = vx_map[*it];
      for(k=0; k<D; ++k) {
        force[D*i+k] += spring_constant*(y[D*i+k] - y[D*j+k]); 
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

void Spacetime::mechanical_solver()
{
  // Method to use an effective force algorithm to calculate the best 
  // geometric configuration for the spacetime complex, with the possibility 
  // of a final conjugate gradient optimization stage to attempt to reduce 
  // the geometry's abnormality. 
  int i,j,m = 0,k = 0,its = 1;
  double vnorm;
  std::vector<int> offset,vx_map;
  std::vector<double> xc,y,ynew,energy;
  const int nv = (signed) skeleton->events.size();
  const int nreal = skeleton->cardinality(0,-1);
  const int D = geometry->dimension();

#ifdef VERBOSE
  std::cout << "Using effective force geometry solver with " << nreal << " active events and background dimension = " << D << "." << std::endl;
#endif

  for(i=0; i<nv; ++i) {
    if (!skeleton->events[i].active()) {
      vx_map.push_back(-1);
      continue;
    }
    vx_map.push_back(m);
    offset.push_back(i);
    energy.push_back(skeleton->events[i].get_energy());
    m++;
  }
  // Initial position from current geometry...
  for(i=0; i<nreal; ++i) {
    k = offset[i];
    geometry->get_coordinates(k,xc);
    for(j=0; j<D; ++j) {
      y.push_back(xc[j]);
    }
  }

  // Initial velocity is zero...
  for(i=D*nreal; i<2*D*nreal; ++i) {
    y.push_back(0.0); 
  }
  for(i=0; i<2*D*nreal; ++i) {
    ynew.push_back(0.0);
  }
  if (engine == Integrator::euler) {
    double F[2*D*nreal];
    do {
      mechanical_force(offset,vx_map,energy,y,F);
      // This just uses the Euler method, not very sophisticated...
      for(i=0; i<2*D*nreal; ++i) {
        ynew[i] = y[i] + step_size*F[i];
      }
      vnorm = 0.0;
      for(i=D*nreal; i<2*D*nreal; ++i) {
        vnorm += ynew[i]*ynew[i];
      }
      vnorm = std::sqrt(vnorm);
      if (vnorm < geometry_tolerance) {
#ifdef VERBOSE
        std::cout << "Mechanical solver converged, exiting loop..." << std::endl;
#endif
        break;
      }
      else if (vnorm > 100.0*double(nreal)) {
#ifdef VERBOSE
        std::cout << "Mechanical solver has excessive velocity norm " << vnorm << " after " << its << " iterations, exiting loop..." << std::endl;
#endif
        break;
      } 
      else if (its > max_int_steps) {
#ifdef VERBOSE
        std::cout << "Mechanical solver reached maximum iterations with velocity norm " << vnorm << ", exiting loop..." << std::endl;
#endif
        break;
      }
#ifdef DEBUG
      for(i=0; i<2*D*nreal; ++i) {
        if (std::isnan(ynew[i])) throw std::runtime_error("NaN detected in mechanical solver at iteration " + std::to_string(its));
      }
#endif
      y = ynew; 
      its++;
    } while(true);
  }
  else if (engine == Integrator::rk4) {
    double k1[2*D*nreal],k2[2*D*nreal],k3[2*D*nreal],k4[2*D*nreal];
    std::vector<double> temp;
    for(i=0; i<2*D*nreal; ++i) {
      temp.push_back(0.0);
    }
    do {
      mechanical_force(offset,vx_map,energy,y,k1);
      for(i=0; i<2*D*nreal; ++i) {
        temp[i] = y[i] + 0.5*step_size*k1[i];
      }
      mechanical_force(offset,vx_map,energy,temp,k2);
      for(i=0; i<2*D*nreal; ++i) {
        temp[i] = y[i] + 0.5*step_size*k2[i];
      }
      mechanical_force(offset,vx_map,energy,temp,k3);
      for(i=0; i<2*D*nreal; ++i) {
        temp[i] = y[i] + step_size*k3[i];
      }
      mechanical_force(offset,vx_map,energy,temp,k4);    
      for(i=0; i<2*D*nreal; ++i) {
        ynew[i] = y[i] + step_size/6.0*(k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]); 
      }
      vnorm = 0.0;
      for(i=D*nreal; i<2*D*nreal; ++i) {
        vnorm += ynew[i]*ynew[i];
      }
      vnorm = std::sqrt(vnorm);
      if (vnorm < geometry_tolerance) {
#ifdef VERBOSE
        std::cout << "Mechanical solver converged, exiting loop..." << std::endl;
#endif
        break;
      }
      else if (vnorm > 100.0*double(nreal)) {
#ifdef VERBOSE
        std::cout << "Mechanical solver has excessive velocity norm " << vnorm << " after " << its << " iterations, exiting loop..." << std::endl;
#endif
        break;
      } 
      else if (its > max_int_steps) {
#ifdef VERBOSE
        std::cout << "Mechanical solver reached maximum iterations with velocity norm " << vnorm << ", exiting loop..." << std::endl;
#endif
        break;
      }
#ifdef DEBUG
      for(i=0; i<2*D*nreal; ++i) {
        if (std::isnan(ynew[i])) throw std::runtime_error("NaN detected in mechanical solver at iteration " + std::to_string(its));
      }
#endif
      y = ynew; 
      its++;
    } while(true);
  }
  for(i=0; i<nreal; ++i) {
    for(j=0; j<D; ++j) {
      xc[j] = ynew[D*i+j];
    }
    k = offset[i];
    geometry->set_coordinates(k,xc);
  }
  geometry->compute_squared_distances();
  compute_volume();
  if (!cgradient_refinement || geometry->get_euclidean() == false) return;
  double d,q,E,prior,alpha = 0.0,beta,E_initial,nx = 0.0,sigma = 0.1;
  SYNARMOSMA::Geometry initial_state; 
  std::vector<int> flexible_edge;
  std::vector<double> s,snew,dx,dy,dx_old,x,c,fx;

  skeleton->determine_flexible_edges(flexible_edge);
  E_initial = compute_abnormality(flexible_edge);
#ifdef VERBOSE
  std::cout << "Initial error is " << E_initial << std::endl;
#endif
  geometry->store(&initial_state);
  for(i=0; i<system_size; ++i) {
    x.push_back(geometry->get_element(i));
    snew.push_back(0.0);
  }

  compute_geometric_gradient(dx,true,flexible_edge);
  for(i=0; i<system_size; ++i) {
    fx.push_back(-dx[i]);
    c.push_back(dx[i]*fx[i]);
    nx += c[i];
  }
  for(i=0; i<max_LS_steps; ++i) {
    for(j=0; j<system_size; ++j) {
      geometry->set_element(j,x[j] + sigma*dx[j]);
    }
    geometry->compute_squared_distances();
    compute_geometric_gradient(dy,false,flexible_edge);
    d = 0.0;
    for(j=0; j<system_size; ++j) {
      d += dy[j]*dx[j] - c[j];
    }
    q = nx/d;
    alpha = -sigma*q;
    sigma = -alpha;
  }
  for(i=0; i<system_size; ++i) {
    geometry->set_element(i,x[i] + alpha*dx[i]);
  }
  geometry->compute_squared_distances();
  E = compute_abnormality(flexible_edge);

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
    compute_geometric_gradient(dx,true,flexible_edge);
    nx = 0.0;
    d = 0.0;
    for(j=0; j<system_size; ++j) {
      nx += dx[j]*(dx[j] - dx_old[j]);
      d += dx_old[j]*dx_old[j];
    }
    beta = std::max(0.0,nx/d);
    for(j=0; j<system_size; ++j) {
      snew[j] = dx[j] + beta*s[j];
    }
    // Now begin the line search...
    nx = 0.0;
    for(j=0; j<system_size; ++j) {
      fx[j] = -dx[j];
      c[j] = snew[j]*fx[j];
      nx += c[j];
    }
    sigma = 0.1;
    for(j=0; j<max_LS_steps; ++j) {
      for(k=0; k<system_size; ++k) {
        geometry->set_element(k,x[k] + sigma*snew[k]);
      }
      geometry->compute_squared_distances();
      compute_geometric_gradient(dy,false,flexible_edge);
      d = 0.0;
      for(k=0; k<system_size; ++k) {
        d += dy[k]*snew[k] - c[k];
      }
      q = nx/d;
      alpha = -sigma*q;
      sigma = -alpha;
    }
    for(j=0; j<system_size; ++j) {
      geometry->set_element(j,x[j] + alpha*snew[j]);
    }
    geometry->compute_squared_distances();
    E = compute_abnormality(flexible_edge);
#ifdef VERBOSE
    std::cout << "At conjugate gradient iteration " << i+1 << ", the error is " << E << std::endl;
#endif
    if (E > prior) {
      for(j=0; j<system_size; ++j) {
        geometry->set_element(j,x[j]);
      }
      geometry->compute_squared_distances();
      break;
    }
    // Update the base values...
    prior = E;
    s = snew;
    dx_old = dx;
    for(j=0; j<system_size; ++j) {
      x[j] = geometry->get_element(j);
    }
    if (E < geometry_tolerance) break;
  }
  if (E > E_initial) {
#ifdef VERBOSE
    std::cout << "Geometry optimization failed, reverting to original geometry." << std::endl;
#endif
    geometry->load(&initial_state);
    geometry->compute_squared_distances();
    for(i=0; i<nv; ++i) {
      if (!skeleton->events[i].active()) continue;
      skeleton->events[i].set_geometry_modified(true);
    }
    compute_volume();
  }
}

void Spacetime::evolutionary_solver()
{
  int i,j,n,in1,vc,joust,generation = 1;
  bool viable;
  double p,f1,f2,severity,fbest,fitness[2*pool_size],ftemp[pool_size];
  std::vector<std::pair<int,int> > vcount;
  std::set<int> modified_vertices,vx;
  std::set<int>::const_iterator it;
  std::vector<SYNARMOSMA::Geometry> pool,ptemp;
  SYNARMOSMA::Geometry optimal; 
  SYNARMOSMA::Geometry initial_state; 
  const int pmagnitude = int(0.15*system_size);
  const double initial_error = error;

#ifdef VERBOSE
  std::cout << "In evolutionary solver with initial error " << error << std::endl;
#endif

  for(i=0; i<2*pool_size; ++i) {
    pool.push_back(SYNARMOSMA::Geometry(*geometry));
    vcount.push_back(std::pair<int,int>(0,0));
  }
  for(i=0; i<pool_size; ++i) {
    ptemp.push_back(SYNARMOSMA::Geometry(*geometry));
  }

  geometry->store(&initial_state);
  geometry->store(&optimal);
  fbest = initial_error;

  // Initialize the population based on some Gaussian flux around the
  // current spacetime geometry
  for(i=0; i<pool_size; ++i) {
    geometry->store(&pool[i]);
    for(j=0; j<pmagnitude; ++j) {
      n = skeleton->RND->irandom(system_size);
      pool[i].get_implied_vertices(n,vx);
      viable = false;
      for(it=vx.begin(); it!=vx.end(); ++it) {
        if (std::abs(skeleton->events[*it].get_geometric_deficiency()) > std::numeric_limits<double>::epsilon()) viable = true;
      }
      if (!viable) continue;
      pool[i].mutation(n,false,false,0.05);
      for(it=vx.begin(); it!=vx.end(); ++it) {
        modified_vertices.insert(*it);
      }
      vx.clear();
    }
    geometry->load(&pool[i]);
    skeleton->compute_dependent_simplices(modified_vertices);
    geometry->compute_squared_distances(modified_vertices);
    compute_volume();
    compute_obliquity();
    structural_deficiency();
    fitness[i] = error;
    modified_vertices.clear();
    geometry->load(&initial_state);
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
        optimal = pool[i];
        fbest = f1;
      }
    }
#ifdef VERBOSE
    std::cout << "At generation " << generation << " the average population fitness is " << p/double(pool_size) << " and minimum fitness is " << f1 << " with global optimum at " << fbest << std::endl;
#endif
    // Now the mutations
    skeleton->RND->initialize_beta(1.0,double(generation)/10.0);
    for(i=pool_size; i<2*pool_size; ++i) {
      severity = skeleton->RND->beta_variate();
      for(j=0; j<system_size; ++j) {
        if (skeleton->RND->drandom() >= severity) continue;
        pool[i].get_implied_vertices(j,vx);
        viable = false;
        for(it=vx.begin(); it!=vx.end(); ++it) {
          if (std::abs(skeleton->events[*it].get_geometric_deficiency()) > std::numeric_limits<double>::epsilon()) viable = true;
        }
        if (!viable) continue;
        pool[i].mutation(j,false,false,severity);
        for(it=vx.begin(); it!=vx.end(); ++it) {
          modified_vertices.insert(*it);
        }
        vx.clear();
      }
      geometry->load(&pool[i]);
      skeleton->compute_dependent_simplices(modified_vertices);
      geometry->compute_squared_distances(modified_vertices);
      compute_volume();
      compute_obliquity();
      structural_deficiency();
      fitness[i] = error;
      modified_vertices.clear();
    }
    // Now we need to see among the total population (2*npop) which are the
    // top npop individuals via a stochastic tournament
    for(i=0; i<2*pool_size; ++i) {
      joust = 0;
      f1 = fitness[i];
      vc = 0;
      do {
        in1 = skeleton->RND->irandom(2*pool_size);
        if (in1 == i) continue;
        f2 = fitness[in1];
        // Calculate the probability that f1 emerges triumphant...
        p = 0.5*(1.0 + std::tanh(10.0*(f2 - f1)));
        if (p > skeleton->RND->drandom()) vc++;
        joust += 1;
      } while(joust < njousts);
      vcount[i] = std::pair<int,int>(i,vc);
    }
    std::sort(vcount.begin(),vcount.end(),SYNARMOSMA::pair_predicate_int);
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
    geometry->load(&optimal);
  }
  else {
#ifdef VERBOSE
    std::cout << "Geometry optimization failed, reverting to original geometry." << std::endl;
#endif
    geometry->load(&initial_state);
  }
  geometry->compute_squared_distances();
  compute_volume();
  compute_obliquity();
  structural_deficiency();
}

void Spacetime::annealing_solver()
{
  int i,j,n,m,step,naccept,nim,nwo_a,nwo_r;
  double S,q,t,E_old,paccept,sigma,E_avg,E_best,delta,temperature = 100.0;
  bool done,viable;
  std::set<int> modified_vertices;
  std::set<int>::const_iterator it;
  std::vector<double> E,lengths,base,output,output1,output2;
  SYNARMOSMA::Geometry optimal; 
  const double initial_error = error;

  geometry->store(&optimal);

  // The goal here should be to find the temperature at which some 80% of
  // random changes are accepted - this will be the melting point of the
  // system
  do {
    naccept = 0;
    E_old = initial_error;
    for(i=0; i<1000; ++i) {
      do {
        m = skeleton->RND->irandom(system_size);
        geometry->get_implied_vertices(m,modified_vertices);
        viable = false;
        for(it=modified_vertices.begin(); it!=modified_vertices.end(); ++it) {
          if (std::abs(skeleton->events[*it].get_geometric_deficiency()) > std::numeric_limits<double>::epsilon()) viable = true;
        }
        if (viable) break;
        modified_vertices.clear();
      } while(true);
      geometry->mutation(m,false,false,thermal_variance);
      skeleton->compute_dependent_simplices(modified_vertices);
      geometry->compute_squared_distances(modified_vertices);
      compute_volume();
      compute_obliquity();
      structural_deficiency();
      if (error < E_old) {
        naccept++;
        E_old = error;
      }
      else {
        q = skeleton->RND->drandom();
        sigma = error - E_old;
        if (q < std::exp(-sigma/temperature)) {
          naccept++;
          E_old = error;
        }
        else {
          geometry->rollback(true);
          skeleton->compute_dependent_simplices(modified_vertices);
          geometry->compute_squared_distances(modified_vertices);
          compute_volume();
          compute_obliquity();
          structural_deficiency();
        }
      }
      modified_vertices.clear();
    }
    paccept = double(naccept)/1000.0;
#ifdef VERBOSE
    std::cout << "At temperature " << temperature << ", the acceptance percentage is " << 100.0*paccept << std::endl;
#endif
    geometry->load(&optimal);
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
        geometry->compute_squared_distances();
        compute_volume();
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
          n = skeleton->RND->irandom(system_size);
          geometry->get_implied_vertices(n,modified_vertices);
          viable = false;
          for(it=modified_vertices.begin(); it!=modified_vertices.end(); ++it) {
            if (std::abs(skeleton->events[*it].get_geometric_deficiency()) > std::numeric_limits<double>::epsilon()) viable = true;
          }
          if (viable) break;
          modified_vertices.clear();
        } while(true);
        geometry->mutation(n,false,false,thermal_variance);
        // Now we need to recaculate some edge lengths and
        // then various defect angles associated with the
        // triangles in their entourage...
        skeleton->compute_dependent_simplices(modified_vertices);
        geometry->compute_squared_distances(modified_vertices);
        compute_volume();
        compute_obliquity();
        structural_deficiency();
        q = error;
        if (q < E_old) {
          E_old = q;
          E.push_back(q);
          if (q < E_best) {
            E_best = q;
            geometry->store(&optimal);
          }
          nim++;
        }
        else {
          t = std::exp((E_old - q)/temperature);
          sigma = skeleton->RND->drandom();
          if (sigma < t) {
            E_old = q;
            E.push_back(q);
            nwo_a++;
          }
          else {
            E.push_back(E_old);
            geometry->rollback(true);
            skeleton->compute_dependent_simplices(modified_vertices);
            geometry->compute_squared_distances(modified_vertices);
            compute_volume();
            compute_obliquity();
            structural_deficiency();
            nwo_r++;
          }
        }
        modified_vertices.clear();
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
  geometry->load(&optimal);
  geometry->compute_squared_distances();
  compute_volume();
  compute_obliquity();
  structural_deficiency();
}

void Spacetime::simplex_solver()
{
  int i,j,k = 0,in1,bindex,windex,ntrans = 0;
  double f,q,centroid[system_size];
  SYNARMOSMA::Geometry SR(*geometry),SE(*geometry); 
  SYNARMOSMA::Geometry initial_state; 
  std::vector<std::pair<int,double> > fitness;
  std::set<int> modified_vertices;
  std::vector<SYNARMOSMA::Geometry> S;
  const int nv = (signed) skeleton->events.size();
  const double initial_error = error;

#ifdef VERBOSE
  std::cout << "Using Nelder-Mead method with dimension " << system_size << std::endl;
#endif
  geometry->store(&initial_state);

  // We begin by creating a simplex of size system_size....
  geometry->store(&SR);
  S.push_back(SR);
  fitness.push_back(std::pair<int,double>(0,error));
  for(i=0; i<system_size; ++i) {
    geometry->get_implied_vertices(i,modified_vertices);
    skeleton->compute_dependent_simplices(modified_vertices);
    geometry->mutation(i,false,false,0.1);
    geometry->compute_squared_distances(modified_vertices);
    geometry->store(&SR);
    S.push_back(SR);
    geometry->rollback(true);
    modified_vertices.clear();
    for(j=0; j<nv; ++j) {
      skeleton->events[j].set_geometry_modified(false);
    }
    for(j=1; j<=Complex::ND; ++j) {
      for(k=0; k<(signed) skeleton->simplices[j].size(); ++k) {
        skeleton->simplices[j][k].set_modified(false);
      }
    }
  }
  // Compute the error for the new events...
  for(i=1; i<1+system_size; ++i) {
    geometry->load(&S[i]);
    geometry->compute_squared_distances();
    compute_volume();
    compute_obliquity();
    structural_deficiency();
    fitness.push_back(std::pair<int,double>(i,error));
  }
  // Sort the simplex vertices from lowest to highest error...
  std::sort(fitness.begin(),fitness.end(),SYNARMOSMA::pair_predicate_dbl);
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
    geometry->compute_squared_distances();
    for(i=0; i<nv; ++i) {
      if (skeleton->events[i].active()) skeleton->events[i].set_geometry_modified(true);
    }
    compute_volume();
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
      geometry->compute_squared_distances();
      for(i=0; i<nv; ++i) {
        if (skeleton->events[i].active()) skeleton->events[i].set_geometry_modified(true);
      }
      compute_volume();
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
      geometry->compute_squared_distances();
      for(i=0; i<nv; ++i) {
        if (skeleton->events[i].active()) skeleton->events[i].set_geometry_modified(true);
      }
      compute_volume();
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
          geometry->compute_squared_distances();
          for(j=0; j<nv; ++j) {
            if (skeleton->events[j].active()) skeleton->events[j].set_geometry_modified(true);
          }
          compute_volume();
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
      geometry->compute_squared_distances();
      for(i=0; i<nv; ++i) {
        if (skeleton->events[i].active()) skeleton->events[i].set_geometry_modified(true);
      }
      compute_volume();
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
          geometry->compute_squared_distances();
          for(j=0; j<nv; ++j) {
            if (skeleton->events[j].active()) skeleton->events[j].set_geometry_modified(true);
          }
          compute_volume();
          compute_obliquity();
          structural_deficiency();
          fitness[i].second = error;
        }
      }
    }
    // Now resort the simplex vertices...
    std::sort(fitness.begin(),fitness.end(),SYNARMOSMA::pair_predicate_dbl);
#ifdef VERBOSE
    std::cout << "At simplex transformation step " << ntrans << " the error is " << fitness[0].second << std::endl;
#endif
  } while(ntrans <= 10*system_size && fitness[0].second > geometry_tolerance);
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
    geometry->load(&initial_state);
  }
  geometry->compute_squared_distances();
  compute_volume();
  compute_obliquity();
  structural_deficiency();
}

void Spacetime::optimize()
{
  int i,n;
  bool deficient = false;
  std::vector<double> v1;
  const int nv = (signed) skeleton->events.size();

  for(i=0; i<nv; ++i) {
    if (!skeleton->events[i].active()) continue;
    if (std::abs(skeleton->events[i].get_deficiency()) > std::numeric_limits<double>::epsilon()) {
      deficient = true;
      break;
    }
  }
  if (!deficient) return;
#ifdef VERBOSE
  std::cout << "Calling the geometric minimization routine with system dimension " << system_size << std::endl;
#endif
  if (solver == Geometry_Solver::minimal) {
    int its = 0;
    bool found;
    std::set<int> N,modified_vertices,candidates,bbarrel;
    std::set<int>::const_iterator it;
    double old_error,sigma;

    structural_deficiency();

    for(i=0; i<nv; ++i) {
      if (!skeleton->events[i].active()) continue;
      if (std::abs(skeleton->events[i].get_geometric_deficiency()) < std::numeric_limits<double>::epsilon()) continue;
      bbarrel.insert(i);
      // Now check to see if all of this vertex's neighbours have a geometric deficiency > 0
      found = false;
      skeleton->events[i].get_neighbours(N);
      for(it=N.begin(); it!=N.end(); ++it) {
        if (std::abs(skeleton->events[*it].get_geometric_deficiency()) < std::numeric_limits<double>::epsilon()) {
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
    std::cout << "There are " << candidates.size() << " candidate events for geometric optimization." << std::endl;
    std::cout << "In the geometry solver the initial error is " << error << std::endl;
#endif
    old_error = error;
    do {
      // Do something clever here to minimize the structural deficiency
      its++;
      // This is decidedly not very clever but it is quick!
      // First the geometry...
      n = skeleton->RND->irandom(candidates);
      // In fact the variance we use should grow smaller as the geometric deficiency of the
      // vertex diminishes with 0.05 the maximum
      sigma = 0.1*skeleton->events[n].get_geometric_deficiency();
      if (sigma > 0.05) sigma = 0.05;
      geometry->mutation(n,true,true,sigma);
      modified_vertices.insert(n);
      skeleton->compute_dependent_simplices(modified_vertices);
      geometry->compute_squared_distances(modified_vertices);
      compute_volume();
      compute_obliquity();

      // Finally, see if it's an improvement...
      structural_deficiency();
#ifdef VERBOSE
      std::cout << its << "  " << n << "  " << old_error << "  " << error << std::endl;
#endif
      modified_vertices.clear();
      if (error > old_error) {
        geometry->rollback(false);
        modified_vertices.insert(n);
        continue;
      }
      old_error = error;
    } while(its < solver_its);
  }
  else if (solver == Geometry_Solver::evolutionary) {
    evolutionary_solver();
  }
  else if (solver == Geometry_Solver::annealing) {
    annealing_solver();
  }
  else if (solver == Geometry_Solver::mechanical) {
    mechanical_solver();
  }
  else if (solver == Geometry_Solver::simplex) {
    simplex_solver();
  }
}
