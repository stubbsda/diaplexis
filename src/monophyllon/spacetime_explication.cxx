#include "spacetime.h"

using namespace DIAPLEXIS;

template<class kind1,class kind2>
bool Spacetime<kind1,kind2>::stellar_addition(int base)
{
  if (skeleton->vertex_dimension(base) < 2) return false;
  // Take one of these 2-simplices and eliminate it in favour of a
  // new event
  int i,m = (signed) skeleton->simplices[2].size();
  unsigned int l;
  std::vector<int> sx;
  std::vector<double> xc,xtemp;
  std::set<int> candidates,S;
  SYNARMOSMA::hash_map::const_iterator qt;

  for(i=0; i<m; ++i) {
    if (!skeleton->active_simplex(2,i)) continue;
    if (skeleton->simplices[2][i].contains(base)) candidates.insert(i);
  }
  // Choose a 2-simplex at random, get its three events and then
  // eliminate each of the three edges in succession...
  m = skeleton->RND->irandom(candidates);
  skeleton->simplices[2][m].get_vertices(sx);

  S.clear();
  S.insert(sx[0]); S.insert(sx[1]);
  qt = skeleton->index_table[1].find(S);
  skeleton->simplex_deletion(1,qt->second);

  S.clear();
  S.insert(sx[0]); S.insert(sx[2]);
  qt = skeleton->index_table[1].find(S);
  skeleton->simplex_deletion(1,qt->second);

  S.clear();
  S.insert(sx[1]); S.insert(sx[2]);
  qt = skeleton->index_table[1].find(S);
  skeleton->simplex_deletion(1,qt->second);    
  // Now add the new event and the three edges...
  for(l=0; l<geometry->dimension(); ++l) {
    xc.push_back(0.0);
  }
  for(i=0; i<3; ++i) {
    geometry->get_coordinates(sx[i],xtemp);
    for(l=0; l<geometry->dimension(); ++l) {
      xc[l] += xtemp[l]; 
    }
  }
  for(l=0; l<geometry->dimension(); ++l) {
    xc[l] = xc[l]/3.0;
  }
  m = event_addition(xc);
  for(i=0; i<3; ++i) {
    S.clear();
    S.insert(m);
    S.insert(sx[i]);
    skeleton->events[sx[i]].set_topology_modified(true);
    skeleton->simplex_addition(S,-1);
  }
  regularization(true);
  return true;    
}

template<class kind1,class kind2>
bool Spacetime<kind1,kind2>::correction(int base)
{
  // This method needs to loop over all active events and then find those which
  // are capable of adding another event at a distance of (roughly) one and which is
  // orthogonal to the event's current set of edges.
  int i,j,n,m,in1;
  unsigned int r;
  bool active,modified = false;
  double l,d1,d2,dbest;
  SYNARMOSMA::hash_map::const_iterator qt;
  std::set<int> candidates;
  Simplex s1;
  const int nv = (signed) skeleton->events.size();
  const double LL = (2.0 - 0.5*abnormality_threshold)*(2.0 - 0.5*abnormality_threshold);
  const double UL = (2.0 + 0.5*abnormality_threshold)*(2.0 + 0.5*abnormality_threshold);
  const double ulimit = (1.0 + abnormality_threshold)*(1.0 + abnormality_threshold);
  const double llimit = (1.0 - abnormality_threshold)*(1.0 - abnormality_threshold);

  for(i=0; i<nv; ++i) {
    if (i == base) continue;
    if (!skeleton->active_event(i)) continue;
    if (std::abs(skeleton->events[i].get_deficiency()) < std::numeric_limits<double>::epsilon()) continue;
    l = std::abs(geometry->get_squared_distance(base,i,true));
    if (l < LL || l > UL) continue;
    // See if there is a third event that lies between these two...
    in1 = -1;
    active = false;
    dbest = 100.0;
    for(j=0; j<nv; ++j) {
      if (j == base || j == i) continue;
      d1 = std::abs(geometry->get_squared_distance(base,j,true));
      d2 = std::abs(geometry->get_squared_distance(i,j,true));
      if ((d1 > llimit && d1 < ulimit) && (d2 > llimit && d2 < ulimit)) {
        l = (d1 - 1.0)*(d1 - 1.0) + (d2 - 1.0)*(d2 - 1.0);
        if (l < dbest) {
          dbest = l;
          in1 = j;
        }
        if (skeleton->active_event(j)) {
          active = true;
          break;
        }
      }
    }
    if (active) continue;
    if (in1 >= 0) {
      if (!skeleton->active_event(in1)) {
#ifdef VERBOSE
        std::cout << "Restoring event " << in1 << " in orthonormal fission between " << base << "  " << i << std::endl;
#endif
        skeleton->events[in1].activate();
        modified = true;
      }
      // Connect this event to v and i if necessary
      modified = skeleton->simplex_addition(base,in1,-1);
      modified = skeleton->simplex_addition(i,in1,-1);
      continue;
    }
    candidates.insert(i);
  }
  if (modified) regularization(true);
  if (candidates.empty()) return false;
  n = skeleton->RND->irandom(candidates);
  std::vector<double> xc,x1,x2;
  // Perform the event fission on the n'th event pair...
  geometry->get_coordinates(base,x1);
  geometry->get_coordinates(n,x2);
  for(r=0; r<geometry->dimension(); ++r) {
    l = 0.5*(x1[r] + x2[r]);
    xc.push_back(l);
  }
  // Add something to check here if this new event is within 0.5 of an existing event...
  for(i=0; i<nv; ++i) {
    if (!skeleton->active_event(i)) continue;
    if (std::abs(geometry->get_squared_distance(i,xc)) <= 0.5) return false;
  }
#ifdef VERBOSE
  std::cout << "Adding event between " << base << " and " << n << std::endl;
#endif
  m = event_addition(xc);
  skeleton->simplex_addition(base,m,-1);
  skeleton->simplex_addition(n,m,-1);

  return true;
}

template<class kind1,class kind2>
bool Spacetime<kind1,kind2>::contraction(int base,double l)
{
  int i,u,vx[2];
  std::set<int> pool;
  const double L2 = l*l;
  const int ne = (signed) skeleton->simplices[1].size();

  for(i=0; i<ne; ++i) {
    if (!skeleton->active_simplex(1,i)) continue;
    skeleton->simplices[1][i].get_vertices(vx);
    if (vx[0] != base && vx[1] != base) continue;
    u = (vx[0] == base) ? vx[1] : vx[0];
    if (std::abs(geometry->get_squared_distance(base,u,true)) > L2) pool.insert(i);
  }
  if (pool.empty()) return false;
  i = skeleton->RND->irandom(pool);
  skeleton->simplex_deletion(1,i);
  return true;
}

template<class kind1,class kind2>
bool Spacetime<kind1,kind2>::compensation_m(int base)
{
  if (skeleton->vertex_dimension(base) < 2) return false;

  int i,j,vx[2];
  double l;
  std::set<int> candidates,S;
  SYNARMOSMA::hash_map::const_iterator qt;
  const int ne = (signed) skeleton->simplices[1].size();
  const double ulimit = (1.0 + abnormality_threshold)*(1.0 + abnormality_threshold);
  const double llimit = (1.0 - abnormality_threshold)*(1.0 - abnormality_threshold);

#ifdef VERBOSE
  std::cout << "Compensation with " << skeleton->dimension() << std::endl;
#endif

  // Remove an edge from this event: ideally one that connects with a event whose degree
  // is also excessive *and* which is relatively far away (edge length >> 1)
  for(i=0; i<ne; ++i) {
    if (!skeleton->active_simplex(1,i)) continue;
    skeleton->simplices[1][i].get_vertices(vx);
    if (vx[0] != base && vx[1] != base) continue;
    j = (vx[0] == base) ? vx[1] : vx[0];
    if (skeleton->vertex_dimension(j) < 2) continue;
    if (skeleton->events[j].get_deficiency() < -std::numeric_limits<double>::epsilon()) continue;
    l = std::abs(geometry->get_squared_distance(base,j,true));
    if (l > llimit && l < ulimit) continue;
    S.clear();
    S.insert(base);
    S.insert(j);
    qt = skeleton->index_table[1].find(S);
    candidates.insert(qt->second);
  }

  if (candidates.empty()) {
#ifdef VERBOSE
    std::cout << "Unable to perform dimensional compensation, no viable candidates..." << std::endl;
#endif
    return false;
  }
  j = skeleton->RND->irandom(candidates);
#ifdef VERBOSE
  std::cout << "Deleting edge with key " << skeleton->get_simplex_key(1,j) << " to reduce the complex's dimensionality" << std::endl;
#endif
  skeleton->simplex_deletion(1,j);
  return true;
}

template<class kind1,class kind2>
bool Spacetime<kind1,kind2>::compensation_g(int base)
{
  int i,j,vx[2],sdegree = 0;
  double l;
  std::set<int> candidates,S;
  SYNARMOSMA::hash_map::const_iterator qt;
  const int D = (signed) geometry->dimension();
  const int ne = (signed) skeleton->simplices[1].size();
  const double ulimit = (1.0 + abnormality_threshold)*(1.0 + abnormality_threshold);
  const double llimit = (1.0 - abnormality_threshold)*(1.0 - abnormality_threshold);

  for(i=0; i<ne; ++i) {
    if (!skeleton->active_simplex(1,i)) continue;
    skeleton->simplices[1][i].get_vertices(vx);
    if (base != vx[0] && base != vx[1]) continue;
    sdegree++;
  }
  if (sdegree == 2*D) return false;
  const int idegree = 2*D;
  const int nv = (signed) skeleton->events.size();

  bool up = (sdegree < idegree) ? true : false;
#ifdef VERBOSE
  std::cout << "Compensation with " << base << "  " << up << "  " << sdegree - idegree << std::endl;
#endif
  if (up) {
    // Add an edge to the event base: ideally one that connects with a event whose degree is
    // also too low and which is relatively close at hand. If the only candidates are distant
    // then do nothing and return false!
    for(i=0; i<nv; ++i) {
      if (i == base) continue;
      if (!skeleton->active_event(i)) continue;
      if (skeleton->edge_exists(base,i)) continue;
      l = std::abs(geometry->get_squared_distance(base,i,true));
      if (l > llimit && l < ulimit) candidates.insert(i);
    }

    if (candidates.empty()) {
#ifdef VERBOSE
      std::cout << "Failure in positive compensation..." << std::endl;
#endif
      return false;
    }
    j = skeleton->RND->irandom(candidates);
    skeleton->events[base].set_topology_modified(true);
    skeleton->events[j].set_topology_modified(true);
    skeleton->simplex_addition(base,j,-1);
  }
  else {
    // Remove an edge from this event: ideally one that connects with a event whose degree
    // is also excessive *and* which is relatively far away (edge length >> 1)
    for(i=0; i<ne; ++i) {
      if (!skeleton->active_simplex(1,i)) continue;
      skeleton->simplices[1][i].get_vertices(vx);
      if (vx[0] != base && vx[1] != base) continue;
      j = (vx[0] == base) ? vx[1] : vx[0];
      l = std::abs(geometry->get_squared_distance(base,j,true));
      if (l > llimit && l < ulimit) continue;
      S.clear();
      S.insert(base);
      S.insert(j);
      qt = skeleton->index_table[1].find(S);
      candidates.insert(qt->second);
    }
    if (candidates.empty()) {
#ifdef VERBOSE
      std::cout << "Failure in negative compensation..." << std::endl;
#endif
      return false;
    }
    i = skeleton->RND->irandom(candidates);
#ifdef VERBOSE
    std::cout << "Deleting edge with key " << skeleton->get_simplex_key(1,i) << " in negative compensation" << std::endl;
#endif
    skeleton->simplex_deletion(1,i);
  }
  return true;
}

template<class kind1,class kind2>
bool Spacetime<kind1,kind2>::foliation_x(int base)
{
  int i,p,n1,n2,vx[2];
  std::set<int> candidates,S;
  SYNARMOSMA::hash_map::iterator qt;
  const int ne = (signed) skeleton->simplices[1].size();

  for(i=0; i<ne; ++i) {
    if (!skeleton->active_simplex(1,i)) continue;
    skeleton->simplices[1][i].get_vertices(vx);
    if (vx[0] != base && vx[1] != base) continue;
    p = (vx[0] == base) ? vx[1] : vx[0];
    if (std::abs(skeleton->events[p].get_deficiency()) < std::numeric_limits<double>::epsilon()) continue;
    candidates.insert(p);
  }
  if (candidates.size() < 2) return false;
  n1 = skeleton->RND->irandom(candidates);
  do {
    n2 = skeleton->RND->irandom(candidates);
    if (n2 != n1) break;
  } while(true);
  S.insert(n1); S.insert(n2);
  qt = skeleton->index_table[1].find(S);
  if (qt == skeleton->index_table[1].end()) return false;
  if (!skeleton->active_simplex(1,qt->second)) return false;
  skeleton->simplex_deletion(1,qt->second);
  return true;
}

template<class kind1,class kind2>
bool Spacetime<kind1,kind2>::deflation(int base)
{
  int d = skeleton->vertex_dimension(base);
  if (d < 2) return false;
  int i,n,dw = 1;
  std::set<int> candidates;
  if (d > 2) dw = skeleton->RND->irandom(1,d);
  const int m = (signed) skeleton->simplices[dw].size();
  for(i=0; i<m; ++i) {
    if (!skeleton->active_simplex(dw,i)) continue;
    if (skeleton->simplices[dw][i].contains(base)) candidates.insert(i);
  }
  n = skeleton->RND->irandom(candidates);
#ifdef VERBOSE
  std::cout << "Deflating a " << d << "-simplex with base " << base << " to a " << dw << "-simplex." << std::endl;
#endif
  skeleton->simplex_deletion(dw,n);
  return true;
}

template<class kind1,class kind2>
bool Spacetime<kind1,kind2>::reduction(int base)
{
  const int d = skeleton->vertex_dimension(base);
  if (d < 2) return false;
  int i,n,m,vx[1+d];
  std::set<int> candidates,s1;
  SYNARMOSMA::hash_map::const_iterator qt;
  const int N = (signed) skeleton->simplices[d].size();

  // Find which edges this event possesses are used in n-simplices (n > 1) and
  // eliminate one of them...
  for(i=0; i<N; ++i) {
    if (!skeleton->active_simplex(d,i)) continue;
    if (skeleton->simplices[d][i].contains(base)) candidates.insert(i);
  }
  m = skeleton->RND->irandom(candidates);
  skeleton->simplices[d][m].get_vertices(vx);
  do {
    n = skeleton->RND->irandom(1+d);
    if (vx[n] != base) break;
  } while(true);
  // Delete the edge v:vx[n]...
#ifdef VERBOSE
  std::cout << "Deleting edge with key " << base << ":" << vx[n] << std::endl;
#endif
  s1.insert(base); s1.insert(vx[n]);
  qt = skeleton->index_table[1].find(s1);
  skeleton->simplex_deletion(1,qt->second);
  return true;
}

template<class kind1,class kind2>
bool Spacetime<kind1,kind2>::amputation(int base,double tolerance)
{
  int i,j,n,p,vx[2];
  std::set<int> candidates;
  const int ne = (signed) skeleton->simplices[1].size();
  const int ulimit = skeleton->dimension();

  if (skeleton->events[base].get_deficiency() > tolerance) candidates.insert(base);
  for(i=0; i<ne; ++i) {
    if (!skeleton->active_simplex(1,i)) continue;
    skeleton->simplices[1][i].get_vertices(vx);
    if (base != vx[0] && base != vx[1]) continue;
    n = (vx[0] == base) ? vx[1] : vx[0];
    if (skeleton->events[n].get_deficiency() > tolerance) candidates.insert(n);
  }
  
  if (candidates.empty()) return false;
  n = skeleton->RND->irandom(candidates);
 
#ifdef VERBOSE
  std::cout << "Amputating event " << n << " and all its dependent simplices" << std::endl;
#endif
  // Delete this event and all its edges...
  skeleton->events[n].deactivate();
  skeleton->events[n].set_topology_modified(true);
  for(i=ulimit; i>=1; --i) {
    p = (signed) skeleton->simplices[i].size();
    for(j=0; j<p; ++j) {
      if (!skeleton->active_simplex(i,j)) continue;
      if (skeleton->simplices[i][j].contains(n)) skeleton->simplex_deletion(i,j);
    }
  }
  if (!skeleton->active_event(n)) {
    if (!skeleton->events[n].zero_energy()) {
      skeleton->events[base].increment_energy(skeleton->events[n].get_energy());
      skeleton->events[n].nullify_energy();
    }
  }
  
  return true;
}

template<class kind1,class kind2>
bool Spacetime<kind1,kind2>::fusion_x(int base,double tolerance)
{
  int i,u;
  std::set<int> candidates;
  const int nv = (signed) skeleton->events.size();

  for(i=0; i<nv; ++i) {
    if (i == base) continue;
    if (!skeleton->active_event(i)) continue;
    if (std::abs(skeleton->events[i].get_deficiency()) < std::numeric_limits<double>::epsilon()) continue;
    if (std::abs(geometry->get_squared_distance(base,i,true)) > tolerance) continue;
    candidates.insert(i);
  }
  if (candidates.empty()) return false;
  u = skeleton->RND->irandom(candidates);
#ifdef VERBOSE
  std::cout << "Fusing deficient events: " << u << " => " << base << std::endl;
#endif
  event_fusion(base,u);
  return true;
}

template<class kind1,class kind2>
bool Spacetime<kind1,kind2>::germination(int base)
{
  // This method constructs new neighbour events w_i for the event base which are
  // unit distance from the base and orthogonal to the base's existing edges, if 
  // possible.
  if (skeleton->events[base].get_boundary()) return false;

  bool good,modified = false;
  double a,b,d_min,delta;
  std::set<int> N,Dm2,Dm1,current,tset,free_dims;
  std::set<int>::const_iterator it,jt;
  std::vector<double> x,y,z,xc,bvector;
  SYNARMOSMA::hash_map::const_iterator qt;
  int i,j,vx[2],D1 = -1,D2 = -1,m,in1,n = -1,mi = 0;
  const int nv = (signed) skeleton->events.size();
  const int ne = (signed) skeleton->simplices[1].size();
  const int D = (signed) geometry->dimension();

  for(i=0; i<ne; ++i) {
    if (!skeleton->active_simplex(1,i)) continue;
    skeleton->simplices[1][i].get_vertices(vx);
    if (base != vx[0] && base != vx[1]) continue;
    j = (vx[0] == base) ? vx[1] : vx[0];
    N.insert(j);
  }
  for(it=N.begin(); it!=N.end(); ++it) {
    n = *it;
    if (skeleton->active_event(n)) break;
  }
#ifdef VERBOSE
  std::cout << "Germination with " << base << "  " << n << std::endl;
#endif
  Dm2.insert(n);

  // The problem here is that this assumes that a^2 + b^2 != 0, which we are
  // by no means guaranteed if the background dimension > 2
  geometry->vertex_difference(n,base,z);
  geometry->get_coordinates(base,bvector);
  xc = bvector;
  for(i=0; i<D; ++i) {
    x.push_back(0.0);
  }
  for(i=0; i<D; ++i) {
    for(j=1+i; j<D; ++j) {
      a = z[i];
      b = z[j];
      delta = a*a + b*b;
      if (delta < convergence_threshold) continue;
      delta = std::sqrt(delta);
      D1 = i;
      D2 = j;
      x[D1] = a/delta;
      x[D2] = b/delta;
    }
    if (D1 > -1) break;
  }
  if (D1 == -1) return false;

  // Check the three possible locations for a event orthogonal to this edge:
  // First, a 90 degree rotation...
  d_min = 0.5;
  in1 = -1;
  xc[D1] = -x[D2] + bvector[D1];
  xc[D2] =  x[D1] + bvector[D2];
  for(i=0; i<nv; ++i) {
    if (i == base) continue;
    delta = std::abs(geometry->get_squared_distance(i,xc));
    delta = std::sqrt(delta);
    if (delta < d_min) {
      in1 = i;
      d_min = delta;
    }
  }
  if (in1 >= 0) {
    if (!skeleton->active_event(in1)) {
      modified = true;
      skeleton->events[in1].activate();
    }
    modified = modified || skeleton->simplex_addition(base,in1,-1);
  }
  else {
    modified = true;
    m = event_addition(base);
    Dm1.insert(m);
    geometry->set_coordinates(m,xc);
    skeleton->simplex_addition(base,m,-1);
  }

  // Second, a 180 degree rotation...
  d_min = 0.5;
  in1 = -1;
  xc[D1] = -x[D1] + bvector[D1];
  xc[D2] = -x[D2] + bvector[D2];
  for(i=0; i<nv; ++i) {
    if (i == base) continue;
    delta = std::abs(geometry->get_squared_distance(i,xc));
    delta = std::sqrt(delta);
    if (delta < d_min) {
      in1 = i;
      d_min = delta;
    }
  }
  if (in1 >= 0) {
    if (!skeleton->active_event(in1)) {
      modified = true;
      skeleton->events[in1].activate();
    }
    modified = modified || skeleton->simplex_addition(base,in1,-1);
  }
  else {
    modified = true;
    m = event_addition(base);
    geometry->set_coordinates(m,xc);
    skeleton->simplex_addition(base,m,-1);
  }

  // And finally the 270 degree rotation...
  d_min = 0.5;
  in1 = -1;
  xc[D1] =  x[D2] + bvector[D1];
  xc[D2] = -x[D1] + bvector[D2];
  for(i=0; i<nv; ++i) {
    if (i == base) continue;
    delta = std::abs(geometry->get_squared_distance(i,xc));
    delta = std::sqrt(delta);
    if (delta < d_min) {
      in1 = i;
      d_min = delta;
    }
  }
  if (in1 >= 0) {
    if (!skeleton->active_event(in1)) {
      modified = true;
      skeleton->events[in1].activate();
    }
    modified = modified || skeleton->simplex_addition(base,in1,-1);
  }
  else {
    modified = true;
    m = event_addition(base);
    geometry->set_coordinates(m,xc);
    Dm1.insert(m);
    skeleton->simplex_addition(base,m,-1);
  }

  if (modified) regularization(false);

  // Now we use rolling cross products to branch out into higher
  // dimensions...
  for(i=0; i<D; ++i) {
    if (i == D1 || i == D2) continue;
    free_dims.insert(i);
  }

  do {
    if (Dm1.empty() || Dm2.empty()) break;

    good = false;
    for(it=Dm2.begin(); it!=Dm2.end(); ++it) {
      n = *it;
      geometry->vertex_difference(n,base,y);
      x.clear();
      x.push_back(y[D1]);
      x.push_back(y[D2]);
      x.push_back(0.0);
      for(jt=free_dims.begin(); jt!=free_dims.end(); ++jt) {
        mi = *jt;
        x[2] = y[mi];
        delta = SYNARMOSMA::norm(x);
        if (delta > convergence_threshold) {
          good = true;
          break;
        }
      }
      if (good) break;
    }
    if (!good) break;
    for(i=0; i<3; ++i) {
      x[i] = x[i]/delta;
    }

    // What is m equal to?
    good = false;
    for(it=Dm1.begin(); it!=Dm1.end(); ++it) {
      m = *it;
      geometry->vertex_difference(m,base,z);
      y.clear();
      y.push_back(z[D1]);
      y.push_back(z[D2]);
      y.push_back(z[mi]);
      delta = SYNARMOSMA::norm(y);
      if (delta > convergence_threshold) {
        good = true;
        break;
      }
    }
    if (!good) break;
    for(i=0; i<3; ++i) {
      y[i] = y[i]/delta;
    }

    SYNARMOSMA::cross_product(x,y,z);

    // Now look to see if there are any existing events near xc
    // and -xc...
    xc.clear();
    for(i=0; i<D; ++i) {
      xc.push_back(bvector[i]);
    }
    xc[D1] = z[0] + bvector[D1];
    xc[D2] = z[1] + bvector[D2];
    xc[mi] = z[2] + bvector[mi];

    d_min = 0.5;
    in1 = -1;
    for(i=0; i<nv; ++i) {
      if (i == base) continue;
      delta = std::abs(geometry->get_squared_distance(i,xc));
      delta = std::sqrt(delta);
      if (delta < d_min) {
        in1 = i;
        d_min = delta;
      }
    }
    if (in1 >= 0) {
      if (!skeleton->active_event(in1)) {
        modified = true;
        skeleton->events[in1].activate();
      }
      modified = modified || skeleton->simplex_addition(base,in1,-1);
    }
    else {
      modified = true;
      m = event_addition(base);
      current.insert(m);
      geometry->set_coordinates(m,xc);
      skeleton->simplex_addition(base,m,-1);
    }

    // Now the mirror image in a three-dimensional subspace...
    xc[D1] = -z[0] + bvector[D1];
    xc[D2] = -z[1] + bvector[D2];
    xc[mi] = -z[2] + bvector[mi];
    d_min = 0.5;
    in1 = -1;
    for(i=0; i<nv; ++i) {
      if (i == base) continue;
      delta = std::abs(geometry->get_squared_distance(i,xc));
      delta = std::sqrt(delta);
      if (delta < d_min) {
        in1 = i;
        d_min = delta;
      }
    }
    if (in1 >= 0) {
      if (!skeleton->active_event(in1)) {
        modified = true;
        skeleton->events[in1].activate();
      }
      modified = modified || skeleton->simplex_addition(base,in1,-1);
    }
    else {
      modified = true;
      m = event_addition(base);
      current.insert(m);
      geometry->set_coordinates(m,xc);
      skeleton->simplex_addition(base,m,-1);
    }

    for(it=free_dims.begin(); it!=free_dims.end(); ++it) {
      if (*it == mi) continue;
      tset.insert(*it);
    }
    free_dims = tset;
    tset.clear();

    Dm2 = Dm1;
    Dm1 = current;
    current.clear();

    D1 = D2;
    D2 = mi;
    mi = -1;
  } while(!free_dims.empty());

  return modified;
}
