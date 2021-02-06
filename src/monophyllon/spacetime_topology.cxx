#include "spacetime.h"

using namespace DIAPLEXIS;

int Spacetime::superposition_fusion()
{
  int i,j,nf,v1,v2,nfail = 0,nfused = 0;
  double delta,pfusion = 0.0;
  std::vector<int> modified;
  std::set<int> vx;
  std::vector<std::pair<int,int> > candidates;
  const int nv = (signed) skeleton->events.size();
  const double na = double(skeleton->cardinality(0));
  const double ulimit = superposition_threshold*superposition_threshold;

  // Event fusion - if two events are close enough together, they
  // should coalesce.
  for(i=0; i<nv; ++i) {
    modified.push_back(0);
    if (!skeleton->active_event(i)) continue;
    for(j=1+i; j<nv; ++j) {
      if (!skeleton->active_event(j)) continue;
      delta = std::abs(geometry->get_squared_distance(i,j,false));
      if (delta < ulimit) candidates.push_back(std::pair<int,int>(i,j));
    }
  }
  if (candidates.empty()) return 0;

  // Now we need to carry out the fusions...
  nf = (signed) candidates.size();
  do {
    i = skeleton->RND->irandom(nf);
    v1 = candidates[i].first;
    v2 = candidates[i].second;
    if (modified[v1] == 1 || modified[v2] == 1) {
      nfail++;
      continue;
    }
    modified[v1] = 1;
    modified[v2] = 1;
    if (event_fusion(v1,v2)) {
#ifdef VERBOSE
      std::cout << "Fusing events: " << v2 << " => " << v1 << " via superposition" << std::endl;
#endif
      nfused++;
    }
#ifdef DEBUG
    assert(consistent());
#endif
    pfusion = double(nfused)/na;
  } while(pfusion < 0.05 && nfused < nf && nfail < 20);

  return nfused;
}

void Spacetime::superposition_fission(int ulimit)
{
  int n,nc = 0;

  // Finally, the opposite possibility - that a given event might undergo
  // spontaneous fission, creating a new event in its immediate vicinity...
  do {
    n = skeleton->RND->irandom(skeleton->events.size());
    if (!skeleton->active_event(n)) continue;
    if (skeleton->RND->drandom() > 0.01) continue;
    if (fission(n,1.0)) {
#ifdef VERBOSE
      std::cout << "Event " << n << " undergoing spontaneous fission during superposition." << std::endl;
#endif
      nc++;
    }
#ifdef DEBUG
    assert(consistent());
#endif
  } while(nc < ulimit);
}

int Spacetime::compression(double threshold)
{
  int i,n,nc,vx[2],output = 0;
  std::set<int>::const_iterator it;
  double s,alpha;
  std::vector<int> candidates;
  const int nv = (signed) skeleton->events.size();
  const int ne = (signed) skeleton->simplices[1].size();

  for(i=0; i<ne; ++i) {
    if (!skeleton->active_simplex(1,i)) continue;
    alpha = std::abs(skeleton->simplices[1][i].get_volume());
    if (alpha > threshold) {
      s = std::exp(-(alpha - threshold));
      if (skeleton->RND->drandom() > s) candidates.push_back(i);
    }
  }
  nc = (signed) candidates.size();
  for(i=0; i<nc; ++i) {
    n = candidates[i];
    if (skeleton->simplex_deletion(1,n)) output++;
  }
  if (!skeleton->connected()) {
    // We need to add enough (short!) edges to keep the spacetime connected...
    int j,k,nsedge,ncomp;
    double l,tlength;
    std::set<int>::const_iterator jt,kt;
    std::set<int> linked,S;
    std::vector<int> component,sedge;
    std::vector<std::tuple<int,int,double> > connect;

    ncomp = skeleton->component_analysis(component);
#ifdef VERBOSE
    std::cout << "There are " << ncomp << " components in this spacetime." << std::endl;
#endif
    std::vector<int>* cvertex = new std::vector<int>[ncomp];
    for(i=0; i<nv; ++i) {
      if (!skeleton->active_event(i)) continue;
      cvertex[component[i]].push_back(i);
    }

    // Now we want to find the shortest possible edge that connects all of these
    // components together, first the brute force solution
    for(i=0; i<ncomp; ++i) {
      for(j=1+i; j<ncomp; ++j) {
        l = minimize_lengths(cvertex[i],cvertex[j],vx);
        connect.push_back(std::tuple<int,int,double>(vx[0],vx[1],l));
      }
    }
    n = (signed) connect.size();
    bool used[n];
    for(i=0; i<n; ++i) {
      used[i] = false;
    }
    std::sort(connect.begin(),connect.end(),SYNARMOSMA::tuple_predicate);
    // Assuming the array "connect" has been sorted in ascending order for the last
    // element...
    j = std::get<0>(connect[0]);
    k = std::get<1>(connect[0]);
    sedge.push_back(j);
    sedge.push_back(k);
    linked.insert(component[j]);
    linked.insert(component[k]);
    tlength = std::get<2>(connect[0]);
    used[0] = true;
    do {
      for(i=1; i<n; ++i) {
        if (used[i]) continue;
        j = component[std::get<0>(connect[i])];
        k = component[std::get<1>(connect[i])];
        jt = std::find(linked.begin(),linked.end(),j);
        kt = std::find(linked.begin(),linked.end(),k);
        if (jt != linked.end() && kt != linked.end()) continue;
        if (jt == linked.end() && kt == linked.end()) continue;
        sedge.push_back(std::get<0>(connect[i]));
        sedge.push_back(std::get<1>(connect[i]));
        tlength += std::get<2>(connect[i]);
        linked.insert(j);
        linked.insert(k);
        break;
      }
    } while((signed) linked.size() < ncomp);
    nsedge = (signed) sedge.size()/2;
#ifdef DEBUG
    assert((nsedge+1) == ncomp);
#endif
#ifdef VERBOSE
    std::cout << "Adding " << nsedge << " edge(s) to reconnect the spacetime complex..." << std::endl;
#endif
    // Now add the edges in sedge...
    for(i=0; i<nsedge; ++i) {
      j = sedge[2*i];
      k = sedge[2*i+1];
      S.insert(j); S.insert(k);
      skeleton->simplex_addition(S,iterations);
      S.clear();
    }
    nc -= nsedge;
#ifdef DEBUG
    assert(skeleton->connected());
#endif
    delete[] cvertex;
  }
  return output;
}

bool Spacetime::interplication(int centre,double size,int D)
{
  // Method to construct a highly-entwined knot after carving out a hole for 
  // in the spacetime complex - it really only makes sense if D is large 
  // enough...  
  if (D <= 3) return false;

  const double dm = double(D);

  int i,j,k,d,p,q,m,nc,nbridge,its;
  double l,alpha,width = size/dm;
  std::vector<double> xc,cvertex;
  std::vector<std::pair<int,double> > ambient;
  std::set<int> S,base,nvertex,bvertex,kvertex;
  std::set<int>::const_iterator it;
  const int g_dim = (signed) geometry->dimension();

  // We begin by eliminating the central event along with any events on this 
  // sheet that are within the a sphere of radius "size"
  geometry->get_coordinates(centre,cvertex); 
  if (!event_deletion(centre)) throw std::runtime_error("Failed event deletion in Spacetime::interplication method!");
  for(i=0; i<(signed) skeleton->events.size(); ++i) {
    if (!skeleton->active_event(i)) continue;
    l = std::abs(geometry->get_squared_distance(centre,i,false));
    if (l < size) {
      if (!event_deletion(i)) throw std::runtime_error("Failed event deletion in Spacetime::interplication method!");
      continue;
    }
    ambient.push_back(std::pair<int,double>(i,l - size));
  }
  std::sort(ambient.begin(),ambient.end(),SYNARMOSMA::pair_predicate_dbl);

  // Ideally the number of boundary events should be the surface area of a d-sphere
  d = g_dim;
  alpha = double(d)/2.0;
  l = std::pow(M_PI,alpha)*std::pow(size,double(d - 1))/std::tgamma(1.0 + alpha);  
  d *= int(l);
  q = (signed) ambient.size();
  if (d > q) d = q;
  for(i=0; i<d; ++i) {
    // If we can't get enough boundary events that are close enough, exit
    if (ambient[i].second > 1.5*size) break;
    bvertex.insert(ambient[i].first);
  }
#ifdef DEBUG
  assert(bvertex.size() > 1);
#endif

  // First add the central element of the knot, a D-simplex
  for(i=0; i<1+D; ++i) {
    for(j=0; j<g_dim; ++j) {
      xc.push_back(cvertex[j] + width*(skeleton->RND->drandom() - 0.5));
    }
    if (!geometry->get_uniform()) {
      for(j=g_dim; j<D; ++j) {
        xc.push_back(-dm + dm/2.0*skeleton->RND->drandom()); 
      }
    }
    S.insert(event_addition(xc));
    xc.clear();
  }
  skeleton->simplex_addition(S,-1);

  // Now add a set of lower-dimensional simplices to this knot and also tie it 
  // in to the existing spacetime complex...
  base = S;
  kvertex = S;
  S.clear();
  for(i=D-1; i>=3; --i) {
    // As the dimension shrinks, the number of simplices should grow
    l = double(i);
    m = int(skeleton->RND->nrandom(dm - l,2.0));
#ifdef VERBOSE
    std::cout << "Adding " << m << " " << i << "-simplices to this complex" << std::endl;
#endif
    if (m <= 0) continue;
    nc = 0;
    width = 2.0*size/l;
    if (width > size) width = size;
    do {    
      for(j=0; j<1+i; ++j) {
        q = -1;
        if (skeleton->RND->drandom() < (0.5*l/dm)) {
          for(k=0; k<g_dim; ++k) {
            xc.push_back(cvertex[k] + width*(skeleton->RND->drandom() - 0.5));
          }
          if (!geometry->get_uniform()) {
            for(k=g_dim; k<i; ++k) {
              alpha = skeleton->RND->nrandom(-l,dm);
              if (alpha < -dm) alpha = -dm + 0.1 + 1.5*skeleton->RND->drandom();
              if (alpha > -0.1) alpha = -0.1 - skeleton->RND->drandom();
              xc.push_back(alpha); 
            }
          }
          q = event_addition(xc);
          xc.clear();
          nvertex.insert(q);
          kvertex.insert(q);
        }
        else {
          its = 0;
          do {
            its++;
            q = skeleton->RND->irandom(base);
            // Make sure it isn't already in S,
            if (S.count(q) == 0) break;
          } while(its < 2*((signed) base.size()));
        }
        if (q == -1) {
          // Create a new event from scratch
          for(k=0; k<g_dim; ++k) {
            xc.push_back(cvertex[k] + width*(skeleton->RND->drandom() - 0.5));
          }
          if (!geometry->get_uniform()) {
            for(k=g_dim; k<i; ++k) {
              alpha = skeleton->RND->nrandom(-l,dm);
              if (alpha < -dm) alpha = -dm + 0.1 + 1.5*skeleton->RND->drandom();
              if (alpha > -0.1) alpha = -0.1 - skeleton->RND->drandom();
              xc.push_back(alpha); 
            }
          }
          q = event_addition(xc);
          xc.clear();
          nvertex.insert(q);
          kvertex.insert(q);          
        }
        S.insert(q);
      }
#ifdef VERBOSE
      std::cout << "Created " << i << "-simplex with " << nvertex.size() << " new events" << std::endl;
#endif
      skeleton->simplex_addition(S,-1);
      S.clear();
      nc++;
    } while(nc < m);
    base = nvertex;
    nvertex.clear();
  }
  // Finally the 2-simplexes to bridge the knot with the ambient spacetime complex
  d = bvertex.size()*(bvertex.size() - 1)/2;
  nbridge = int((0.1 + 0.15*skeleton->RND->drandom())*d);
#ifdef VERBOSE
  std::cout << "Adding " << nbridge << " 2-simplices" << std::endl;
#endif
  d = 0;
  do {
    q = skeleton->RND->irandom(bvertex);
    if (skeleton->RND->drandom() < 0.33) {
      do { 
        p = skeleton->RND->irandom(bvertex);
        if (p != q) break;
      } while(true);
    }
    else {
      alpha = 10.0*size;
      for(it=kvertex.begin(); it!=kvertex.end(); ++it) {
        l = std::abs(geometry->get_squared_distance(q,*it,false));
        if (l < alpha) {
          p = *it;
          alpha = l;
        }
      }
    }
    S.insert(q);
    S.insert(p);
    // Find the closest new event...
    ambient.clear();
    for(it=kvertex.begin(); it!=kvertex.end(); ++it) {
      if (p == *it) continue;
      l = std::abs(geometry->get_squared_distance(q,*it,false));
      alpha = std::abs(geometry->get_squared_distance(p,*it,false));
      ambient.push_back(std::pair<int,double>(*it,l + alpha));
    }
    std::sort(ambient.begin(),ambient.end(),SYNARMOSMA::pair_predicate_dbl);
    if (skeleton->RND->drandom() < 0.5) {
      S.insert(ambient[0].first);
    }
    else {
      S.insert(ambient[1].first);
    }
    if (skeleton->simplex_addition(S,-1)) d++;
    S.clear();
  } while(d < nbridge);

  regularization(true);
  return true;
}

void Spacetime::regularization(bool minimal)
{
  if (!minimal) skeleton->simplicial_implication();

  skeleton->compute_neighbours();
  if (skeleton->connected()) {
    skeleton->compute_entourages();
    return;
  }

  int i,j,k,v1,v2,n1,n2,nc,loc = -1;
  unsigned int max = 0;
  double l,mdelta;
  std::vector<int> component;
  const int nv = (signed) skeleton->events.size();

  nc = skeleton->component_analysis(component);
#ifdef VERBOSE
  std::cout << "There are " << nc << " components in this spacetime." << std::endl;
#endif
  // So there are nc (>1) components in this graph, we will unite them in the
  // most minimal fashion, using (nc-1) edges
  // November 2, 2012: We need to consider ensuring that the minimalist linking of these
  // distinct components is also minimalist in a geometric sense (i.e. the shortest possible
  // edge) and if a "component" is nothing more than an isolated event with zero energy
  // then drop the event (set ubiquity = 1).
  std::vector<int>* cvertex = new std::vector<int>[nc];

  for(i=0; i<nv; ++i) {
    if (!skeleton->active_event(i)) continue;
    cvertex[component[i]].push_back(i);
  }

  for(i=0; i<nc; ++i) {
    if (cvertex[i].size() > max) {
      max = cvertex[i].size();
      loc = i;
    }
  }
  for(i=0; i<nc; ++i) {
    if (cvertex[i].size() == 1) {
      // Isolated event, let's see what its energy is...
      if (skeleton->events[cvertex[i][0]].zero_energy()) {
        // Eliminate this event altogether...
#ifdef VERBOSE
        std::cout << "Eliminating isolated event " << cvertex[i][0] << std::endl;
#endif
        skeleton->events[cvertex[i][0]].deactivate();
        cvertex[i].erase(cvertex[i].begin());
      }
    }
  }
  n1 = (signed) cvertex[loc].size();

  for(i=0; i<nc; ++i) {
    if (i == loc) continue;
    if (cvertex[i].empty()) continue;
    // We need to choose the closest pair of events, not just a random pair
    n2 = (signed) cvertex[i].size();
    mdelta = std::numeric_limits<double>::infinity();
    v1 = -1;
    v2 = -1;
    for(j=0; j<n1; ++j) {
      for(k=0; k<n2; ++k) {
        l = std::abs(geometry->get_squared_distance(cvertex[loc][j],cvertex[i][k],true));
        if (l < mdelta) {
          mdelta = l;
          v1 = cvertex[loc][j];
          v2 = cvertex[i][k];
        }
      }
    }
#ifdef DEBUG
    assert(v1 > -1);
    assert(v2 > -1);
#endif
#ifdef VERBOSE
    std::cout << "Linking together " << v1 << " and " << v2 << " at an absolute distance of " << mdelta << std::endl;
#endif
    skeleton->simplex_addition(v1,v2,-1);
    skeleton->events[v1].set_topology_modified(true);
    skeleton->events[v2].set_topology_modified(true);
  }
  skeleton->compute_neighbours();
#ifdef DEBUG
  assert(skeleton->connected());
#endif
  skeleton->compute_entourages();
  delete[] cvertex;
}
