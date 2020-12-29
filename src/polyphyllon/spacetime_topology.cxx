#include "spacetime.h"

using namespace DIAPLEXIS;

int Spacetime::superposition_fusion(double threshold)
{
  int i,j,m,nf,v1,v2,nfail = 0,nfused = 0;
  double delta,pfusion = 0.0;
  std::vector<int> modified;
  std::set<int> vx;
  std::vector<std::pair<int,int> > candidates;
  const double na = double(skeleton->cardinality(0,-1));
  const int nv = (signed) skeleton->events.size();
  const int ulimit = skeleton->dimension(-1);

  // Event fusion - if two events are close enough together, they
  // should coalesce.
  for(i=0; i<nv; ++i) {
    modified.push_back(0);
    if (!skeleton->events[i].active()) continue;
    for(j=1+i; j<nv; ++j) {
      if (!skeleton->events[j].active()) continue;
      delta = std::abs(geometry-> get_squared_distance(i,j,false));
      if (delta < threshold) candidates.push_back(std::pair<int,int>(i,j));
    }
  }
  if (candidates.empty()) return 0;

  // Now we need to carry out the fusions...
  nf = (signed) candidates.size();
  do {
    i = skeleton->RND->irandom(nf);
    v1 = candidates[i].first;
    v2 = candidates[i].second;
    if (modified[v1] || modified[v2]) {
      nfail++;
      continue;
    }
    modified[v1] = 1;
    modified[v2] = 1;
#ifdef VERBOSE
    std::cout << "Fusing events: " << v2 << " => " << v1 << " via superposition" << std::endl;
#endif
    if (event_fusion(v1,v2,-1)) nfused++;
#ifdef DEBUG
    assert(consistent());
#endif
    pfusion = double(nfused)/na;
  } while(pfusion < 0.05 && nfused < nf && nfail < 20);

  return nfused;
}

void Spacetime::superposition_fission(int ulimit)
{
  int n,sheet,nc = 0;
  std::set<int> locus;

  // Finally, the opposite possibility - that a given event might undergo
  // spontaneous fission, creating a new event in its immediate vicinity...
  do {
    n = skeleton->RND->irandom(skeleton->events.size());
    if (!skeleton->events[n].active()) continue;
    if (skeleton->RND->drandom() > 0.01) continue;
    skeleton->events[n].get_ubiquity(locus);
    sheet = skeleton->RND->irandom(locus);
    if (fission(n,1.0,sheet)) nc++;
  } while(nc < ulimit);
}

int Spacetime::compression(double threshold)
{
  int i,n,nc,v[2],output = 0;
  std::set<int> locus;
  std::set<int>::const_iterator it;
  double s,alpha;
  std::vector<int> candidates;
  const int nv = (signed) skeleton->events.size();
  const int ne = (signed) skeleton->simplices[1].size();
  const int nt = (signed) codex.size();

  for(i=0; i<ne; ++i) {
    if (!skeleton->simplices[1][i].active()) continue;
    alpha = std::abs(skeleton->simplices[1][i].get_volume());
    if (alpha > threshold) {
      s = std::exp(-(alpha - threshold));
      if (skeleton->RND->drandom() > s) candidates.push_back(i);
    }
  }
  nc = (signed) candidates.size();
  for(i=0; i<nc; ++i) {
    n = candidates[i];
    if (skeleton->simplex_deletion(1,n,-1)) output++;
  }
  if (!skeleton->connected(-1)) {
    // We need to add enough (short!) edges to keep the spacetime connected...
    int j,k,nsedge,ncomp;
    double l,tlength;
    std::set<int>::const_iterator jt,kt;
    std::set<int> linked;
    std::vector<int> component,sedge;
    std::vector<std::tuple<int,int,double> > connect;
    Simplex S;

    ncomp = skeleton->component_analysis(component,-1);
#ifdef VERBOSE
    std::cout << "There are " << ncomp << " components in this spacetime complex." << std::endl;
#endif
    std::vector<int>* cvertex = new std::vector<int>[ncomp];
    for(i=0; i<nv; ++i) {
      if (!skeleton->events[i].active()) continue;
      cvertex[component[i]].push_back(i);
    }

    // Now we want to find the shortest possible edge that connects all of these
    // components together, first the brute force solution
    for(i=0; i<ncomp; ++i) {
      for(j=1+i; j<ncomp; ++j) {
        l = minimize_lengths(cvertex[i],cvertex[j],v);
        connect.push_back(std::tuple<int,int,double>(v[0],v[1],l));
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
      l = skeleton->RND->irandom(nt);
      locus.insert(l);
      skeleton->simplex_addition(j,k,locus);
      locus.erase(l);
    }
    nc -= nsedge;
#ifdef DEBUG
    assert(skeleton->connected(-1));
#endif
    delete[] cvertex;
  }
  return output;
}

bool Spacetime::interplication(int centre,double size,int D,int sheet)
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
  std::set<int> S,locus,base,nvertex,bvertex,kvertex;
  std::set<int>::const_iterator it;
  const int g_dim = (signed) geometry->dimension();

  locus.insert(sheet);

  // We begin by eliminating the central event along with any events on this 
  // sheet that are within the a sphere of radius "size"
  geometry->get_coordinates(centre,cvertex); 
  if (!event_deletion(centre,sheet)) throw std::runtime_error("Failed event deletion in Spacetime::interplication method!");
  for(i=0; i<(signed) skeleton->events.size(); ++i) {
    if (!skeleton->events[i].active(sheet)) continue;
    l = std::abs(geometry->get_squared_distance(centre,i,false));
    if (l < size) {
      if (!event_deletion(i,sheet)) throw std::runtime_error("Failed event deletion in Spacetime::interplication method!");
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
    S.insert(event_addition(xc,sheet));
    xc.clear();
  }
  skeleton->simplex_addition(S,locus);

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
          q = event_addition(xc,sheet);
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
          // Create a new vertex from scratch
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
          q = event_addition(xc,sheet);
          xc.clear();
          nvertex.insert(q);
          kvertex.insert(q);          
        }
        S.insert(q);
      }
#ifdef VERBOSE
      std::cout << "Created " << i << "-simplex with " << nvertex.size() << " new events" << std::endl;
#endif
      skeleton->simplex_addition(S,locus);
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
    if (skeleton->simplex_addition(S,locus)) d++;
    S.clear();
  } while(d < nbridge);

  regularization(true,sheet);
  return true;
}

void Spacetime::regularization(bool minimal,int sheet)
{
  if (!minimal) skeleton->simplicial_implication(sheet);

  skeleton->compute_neighbours();
  if (skeleton->connected(sheet)) {
    skeleton->compute_entourages(sheet);
    return;
  }

  int i,j,k,nc,v1,v2,n1,n2,loc = -1;
  double l,mdelta;
  unsigned int max = 0;
  Simplex S;
  std::set<int> locus,colours;
  std::vector<int> component;
  SYNARMOSMA::hash_map::const_iterator qt;
  const int nv = (signed) skeleton->events.size();
  const int nt = (signed) codex.size();

  for(i=0; i<nt; ++i) {
    if (codex[i].active) colours.insert(i);
  }
  nc = skeleton->component_analysis(component,sheet);
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
  if (sheet == -1) {
    for(i=0; i<nv; ++i) {
      if (!skeleton->events[i].active()) continue;
      cvertex[component[i]].push_back(i);
    }
  }
  else {
    for(i=0; i<nv; ++i) {
      if (!skeleton->events[i].active(sheet)) continue;
      cvertex[component[i]].push_back(i);
    }
  }

  for(i=0; i<nc; ++i) {
    if (cvertex[i].size() > max) {
      max = cvertex[i].size();
      loc = i;
    }
  }
  for(i=0; i<nc; ++i) {
    if (cvertex[i].size() == 1) {
      // Isolated vertex, let's see what its energy is...
      if (skeleton->events[cvertex[i][0]].zero_energy()) {
        // Eliminate this vertex altogether...
#ifdef VERBOSE
        std::cout << "Eliminating isolated vertex " << cvertex[i][0] << std::endl;
#endif
        skeleton->events[cvertex[i][0]].deactivate();
        cvertex[i].erase(cvertex[i].begin());
      }
    }
  }
  n1 = (signed) cvertex[loc].size();
  if (sheet == -1) {
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
      j = skeleton->RND->irandom(colours);
      locus.insert(j);
      skeleton->simplex_addition(v1,v2,locus);
      locus.erase(j);
      // Need to ensure that v1 and v2 are also coloured by sheet...
      if (!skeleton->events[v1].active(j)) skeleton->events[v1].activate(j); 
      if (!skeleton->events[v2].active(j)) skeleton->events[v2].activate(j); 
    }
  }
  else {
    locus.insert(sheet);
    for(i=0; i<nc; ++i) {
      if (i == loc) continue;
      if (cvertex[i].empty()) continue;
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
      skeleton->simplex_addition(v1,v2,locus);
    }
  }
  skeleton->compute_neighbours();
#ifdef DEBUG
  assert(skeleton->connected(sheet));
#endif
  skeleton->compute_entourages(sheet);
  delete[] cvertex;
}
