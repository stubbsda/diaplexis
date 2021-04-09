#include "spacetime.h"

using namespace DIAPLEXIS;

template<class kind1,class kind2>
bool Spacetime<kind1,kind2>::stellar_deletion(int base)
{
  // If this event is already part of a d-simplex, d >= 2, then
  // this operation is pointless...
  if (skeleton->vertex_dimension(base) >= 2) return false;
  int m,vx[2];
  std::set<int> S,nset;
  std::set<int>::const_iterator it;

  skeleton->events[base].get_entourage(S);
  for(it=S.begin(); it!=S.end(); ++it) {
    if (skeleton->active_simplex(1,*it)) {
      skeleton->simplices[1][*it].get_vertices(vx);
      m = (vx[0] == base) ? vx[1] : vx[0];
      nset.insert(m);
    }
  }
  if (nset.size() != 3) return false;
  event_deletion(base);
  skeleton->simplex_addition(nset,-1);
  regularization(true);
  return true;
}

template<class kind1,class kind2>
bool Spacetime<kind1,kind2>::foliation_m(int base)
{
  int i,p,n1,n2,vx[2];
  std::set<int> candidates,s1;
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

  // Clearly adding this edge creates at least a 2-simplex (v,n1,n2) which should be
  // added to the spacetime complex!
  std::set<int> leaf;
  leaf.insert(base);
  leaf.insert(n1);
  leaf.insert(n2);
  if (skeleton->simplex_addition(leaf,-1)) {
    skeleton->events[n1].set_topology_modified(true);
    skeleton->events[n2].set_topology_modified(true);
    return true;
  } 
  return false;
}

template<class kind1,class kind2>
bool Spacetime<kind1,kind2>::fusion_m(int base)
{
  int i,u,vx[2];
  std::set<int> candidates;
  const int ne = (signed) skeleton->simplices[1].size();

  for(i=0; i<ne; ++i) {
    if (!skeleton->active_simplex(1,i)) continue;
    skeleton->simplices[1][i].get_vertices(vx);
    if (vx[0] != base && vx[1] != base) continue;
    u = (vx[0] == base) ? vx[1] : vx[0];
    candidates.insert(u);
  }
  if (candidates.empty()) return false;
  u = skeleton->RND->irandom(candidates);
#ifdef VERBOSE
  std::cout << "Fusing events: " << u << " => " << base << std::endl;
#endif
  event_fusion(base,u);
  return true;
}

template<class kind1,class kind2>
bool Spacetime<kind1,kind2>::fission(int base,double density)
{
  int i,p,q,n,vx[2];
  const int ne = (signed) skeleton->simplices[1].size();

  if (base >= 0) {
    n = 0;

    p = event_addition(base);

    for(i=0; i<ne; ++i) {
      if (!skeleton->active_simplex(1,i)) continue;
      skeleton->simplices[1][i].get_vertices(vx);
      if (vx[0] != base && vx[1] != base) continue;
      q = (vx[0] == base) ? vx[1] : vx[0];
      if (skeleton->RND->drandom() < density) {
        skeleton->simplex_addition(q,p,-1);
        n++;
      }
    }
    skeleton->simplex_addition(base,p,-1);
#ifdef VERBOSE
    std::cout << "Added " << 1+n << " edges to the complex after fission of " << base << " to " << p << std::endl;
#endif
    return true;
  }
  else {
    // This method takes a low-dimensional simplex (d = 0, 1 or 2) and causes
    // a new d-simplex to be created that shares the same neighbour tables and
    // includes between 1 and (1+d)^2 links the new d-simplex and its ancestor.
    std::set<int> nsimplex;

    if (skeleton->dimension() < 1) return false;

    if (skeleton->RND->drandom() < 0.5) {
      // The simplest case: we just need to add a new event, clone the antecedent
      // event's one-dimensional entourage and lastly create a 1-simplex joining
      // the two events.
      int d;

      do {
        p = skeleton->RND->irandom(skeleton->events.size());
        if (skeleton->active_event(p)) break;
      } while(true);
      q = event_addition(p);
      nsimplex.insert(p);
      nsimplex.insert(q);

      for(i=0; i<ne; ++i) {
        if (!skeleton->active_simplex(1,i)) continue;
        skeleton->simplices[1][i].get_vertices(vx);
        if (vx[0] != p && vx[1] != p) continue;
        n = (vx[0] == p) ? vx[1] : vx[0];
        nsimplex.insert(n);
      }
      d = (signed) nsimplex.size() - 1;
      if (d >= Complex<kind1>::ND) {
        std::set<int> s1;
        d = skeleton->RND->irandom(2,Complex<kind1>::ND);
        do {
          n = skeleton->RND->irandom(nsimplex);
          s1.insert(n);
        } while((signed) s1.size() < (1+d));
        nsimplex = s1;
      }
#ifdef VERBOSE
      std::cout << "Event fission on " << p << " with new event " << q << " and simplex dimension " << d << std::endl;
#endif
      skeleton->simplex_addition(nsimplex,-1);
    }
    else {
      // The plot thickens: we must create a new 1-simplex out of an existing one, perhaps by forming
      // a pair of 2-simplices: (w1,u_new,v_new) and (w2,u_new,v_new)
      std::set<int> antecedent;

      do {
        n = skeleton->RND->irandom(skeleton->simplices[1].size());
        if (skeleton->active_simplex(1,n)) break;
      } while(true);
      skeleton->simplices[1][n].get_vertices(vx);
      antecedent.insert(vx[0]);
      antecedent.insert(vx[1]);
      p = event_addition(antecedent);
      antecedent.insert(p);
      q = event_addition(antecedent);

      nsimplex.insert(vx[0]);
      nsimplex.insert(p);
      nsimplex.insert(q);
      skeleton->simplex_addition(nsimplex,-1);

      nsimplex.clear();
      nsimplex.insert(vx[1]);
      nsimplex.insert(p);
      nsimplex.insert(q);
      skeleton->simplex_addition(nsimplex,-1);
    }
    return true;
  }
}

template<class kind1,class kind2>
bool Spacetime<kind1,kind2>::circumvolution()
{
  // This method attempts a circumvolution using boundary edges, it is
  // normally only called when the global energy reaches a critical threshold.
  int i,j,k,l,u[2],w[2];
  std::set<int> edge_set,S,N;
  std::set<int>::const_iterator it,jt;
  std::vector<std::set<int> > F;
  SYNARMOSMA::hash_map::const_iterator qt;
  const int D = skeleton->dimension();
  int vx[D];

  for(i=0; i<(signed) skeleton->simplices[D].size(); ++i) {
    if (!skeleton->active_simplex(D,i)) continue;
    skeleton->simplices[D][i].get_faces(F);
    for(j=0; j<=D; ++j) {
      qt = skeleton->index_table[D-1].find(F[j]);
      skeleton->simplices[D-1][qt->second].get_entourage(N);
      if (N.size() == 1) {
        skeleton->simplices[D-1][qt->second].get_vertices(vx);
        for(k=0; k<D; ++k) {
          for(l=k+1; l<D; ++l) {
            S.clear();
            S.insert(vx[k]);
            S.insert(vx[l]);
            qt = skeleton->index_table[1].find(S);
            edge_set.insert(qt->second);
          }
        }
      }
    }
  }
  if (edge_set.empty()) return false;

  int d,n1,n2,np;
  std::vector<int> candidates;
#ifdef VERBOSE
  std::cout << "There are " << edge_set.size() << " boundary edges." << std::endl;
#endif
  for(it=edge_set.begin(); it!=edge_set.end(); ++it) {
    n1 = *it;
    for(jt=edge_set.begin(); jt!=edge_set.end(); ++jt) {
      n2 = *jt;
      if (n1 <= n2) continue;

      skeleton->simplices[1][n1].get_vertices(w);
      skeleton->simplices[1][n2].get_vertices(u);

      d = skeleton->combinatorial_distance(w[0],u[0]);
      if (d <= 2) continue;

      d = skeleton->combinatorial_distance(w[1],u[0]);
      if (d <= 2) continue;

      d = skeleton->combinatorial_distance(w[0],u[1]);
      if (d <= 2) continue;

      d = skeleton->combinatorial_distance(w[1],u[1]);
      if (d <= 2) continue;

      candidates.push_back(n1);
      candidates.push_back(n2);
    }
  }
  if (candidates.empty()) return false;
  np = candidates.size()/2;
  i = skeleton->RND->irandom(np);
  n1 = candidates[2*i];
  n2 = candidates[2*i+1];
  skeleton->simplices[1][n1].get_vertices(w);
  skeleton->simplices[1][n2].get_vertices(u);
  skeleton->events[u[0]].set_topology_modified(true);
  skeleton->events[u[1]].set_topology_modified(true);

  if (skeleton->RND->drandom() < 0.5) {
    event_fusion(w[0],u[0]);
    event_fusion(w[1],u[1]);
  }
  else {
    event_fusion(w[0],u[1]);
    event_fusion(w[1],u[0]);
  }
  return true;
}

template<class kind1,class kind2>
bool Spacetime<kind1,kind2>::circumvolution(int base)
{
  // This method seeks to fuse together two d-simplices, one of
  // which contains the event base
  int i,d,nd,s1,s2;
  Simplex S;
  std::set<int> candidates;
  std::vector<int> order;

  if (base >= 0) {
    nd = skeleton->vertex_dimension(base);
    if (nd < 1) return false;
    d = skeleton->RND->irandom(1,nd);
    nd = (signed) skeleton->simplices[d].size();
    for(i=0; i<nd; ++i) {
      if (!skeleton->active_simplex(d,i)) continue;
      if (skeleton->simplices[d][i].contains(base)) candidates.insert(i);
    }
    s1 = skeleton->RND->irandom(candidates);
    candidates.clear();
    for(i=0; i<nd; ++i) {
      if (!skeleton->active_simplex(d,i)) continue;
      if (i == s1) continue;
      S = skeleton->simplices[d][s1] ^ skeleton->simplices[d][i];
      if (S.empty()) candidates.insert(i);
    }
    if (candidates.empty()) return false;
    // Now find which of these candidates is the closest to the d-simplex s1, where the
    // distance between two simplices is defined to be the maximum inter-event distance
    int j,vx1[1+d],vx2[1+d];
    double D1,delta;
    std::set<int> slist;
    std::set<int>::const_iterator it;

    skeleton->simplices[d][s1].get_vertices(vx1);
    for(it=candidates.begin(); it!=candidates.end(); ++it) {
      skeleton->simplices[d][*it].get_vertices(vx2);
      D1 = 0.0;
      for(i=0; i<=d; ++i) {
        for(j=1+i; j<=d; ++j) {
          delta = std::abs(geometry->get_squared_distance(vx1[i],vx2[j],true));
          if (delta > D1) D1 = delta;
        }
      }
      if (D1 < 2.5) slist.insert(*it);
    }
    if (slist.empty()) return false;
    s2 = skeleton->RND->irandom(slist);
  }
  else {
    do {
      d = skeleton->RND->irandom(1,skeleton->dimension());
      for(i=0; i<(signed) skeleton->simplices[d].size(); ++i) {
        if (skeleton->active_simplex(d,i)) candidates.insert(i);
      }
      if (candidates.size() > 3) break;
    } while(true);
    do {
      s1 = skeleton->RND->irandom(candidates);
      s2 = skeleton->RND->irandom(candidates);
      if (s1 == s2) continue;
      // Now check to make sure that the intersection of
      // these two d-simplices is null...
      S = skeleton->simplices[d][s1] ^ skeleton->simplices[d][s2];
      if (S.empty()) break;
    } while(true);
  }
#ifdef VERBOSE
  std::cout << "Circumvolving the " << d << "-simplices: " << skeleton->get_simplex_key(d,s2) << " => " << skeleton->get_simplex_key(d,s1) << std::endl;
#endif
  int v1[d+1],v2[d+1];
  for(i=0; i<=d; ++i) {
    order.push_back(i);
  }
  std::random_shuffle(order.begin(),order.end());
  skeleton->simplices[d][s1].get_vertices(v1);
  skeleton->simplices[d][s2].get_vertices(v2);

  for(i=0; i<=d; ++i) {
    event_fusion(v1[i],v2[order[i]]);
  }
  return true;
}

template<class kind1,class kind2>
bool Spacetime<kind1,kind2>::expansion(int base,double creativity)
{
  // Create an entirely new d-simplex
  int i,d,k,novum = 0;
  std::set<int> vx;

  if (base >= 0) {
    int n1 = skeleton->vertex_dimension(base);
    if (n1 == Complex<kind1>::ND) return false;
    int u,vtx[2],its = 0;
    std::set<int> M;
    double tau = std::abs(skeleton->events[base].get_deficiency()) + 25.0;
    const int ne = (signed) skeleton->simplices[1].size();

    d = int(double(Complex<kind1>::ND-1-n1)*1.0/(1.0 + std::exp(tau)) + double(1 + n1));

    for(i=0; i<ne; ++i) {
      if (!skeleton->active_simplex(1,i)) continue;
      skeleton->simplices[1][i].get_vertices(vtx);
      if (vtx[0] != base && vtx[1] != base) continue;
      u = (vtx[0] == base) ? vtx[1] : vtx[0];
      if (skeleton->events[u].get_deficiency() < -std::numeric_limits<double>::epsilon()) M.insert(u);
    }
    if (M.empty()) creativity = 1.0;
    vx.insert(base);
    do {
      its++;
      if (its == Complex<kind1>::ND) creativity = 1.0;
      if (skeleton->RND->drandom() < creativity) {
        // Create a new event...
        k = event_addition(base);
        novum++;
      }
      else {
        k = skeleton->RND->irandom(M);
      }
      vx.insert(k);
    } while((signed) vx.size() < (1+d));
  }
  else {
    int j,nv;
    bool success;
    std::set<int>::const_iterator it;

    i = skeleton->dimension();
    d = (i <= 2) ? 3 : skeleton->RND->irandom(2,i);
    if (double(d+1)/double(skeleton->events.size()) > 0.8) creativity = 1.0;
    do {
      if (skeleton->RND->drandom() < creativity) {
        // Create a new event...
        k = event_addition(vx);
        novum++;
      }
      else {
        // Grab an existing event...
        j = 0;
        nv = (signed) skeleton->events.size();
        success = false;
        do {
          k = skeleton->RND->irandom(nv);
          it = std::find(vx.begin(),vx.end(),k);
          if (it == vx.end()) {
            skeleton->events[k].activate();
            success = true;
            break;
          }
          ++j;
        } while(j < nv);
        if (success == false) {
          k = event_addition(vx);
          novum++;
        }
      }
      vx.insert(k);
    } while((signed) vx.size() < (1+d));
  }
#ifdef VERBOSE
  std::cout << "Created a " << d << "-simplex with " << novum << " new events." << std::endl;
#endif
  skeleton->simplex_addition(vx,-1);
  return true;
}

template<class kind1,class kind2>
bool Spacetime<kind1,kind2>::expansion(int base)
{
  int n = skeleton->vertex_dimension(base);
  if (n == Complex<kind1>::ND) return false;
  int i,u,d,m,vtx[2];
  std::set<int> vx,s,N;
  std::set<int>::const_iterator it;
  SYNARMOSMA::hash_map::const_iterator qt;
  Simplex S;
  const int ne = (signed) skeleton->simplices[1].size();
  double tau = std::abs(skeleton->events[base].get_deficiency()) + 25.0;

  d = int(double(Complex<kind1>::ND-1-n)*1.0/(1.0 + std::exp(tau)) + double(1 + n));

  vx.insert(base);
  for(i=0; i<d; ++i) {
    // Create a new event...
    m = event_addition(base);
    vx.insert(m);
  }
#ifdef VERBOSE
  std::cout << "Created a " << d << "-simplex with base " << base << std::endl;
#endif
  skeleton->simplex_addition(vx,-1);
  for(i=0; i<ne; ++i) {
    if (!skeleton->active_simplex(1,i)) continue;
    skeleton->simplices[1][i].get_vertices(vtx);
    if (vtx[0] != base && vtx[1] != base) continue;
    u = (vtx[0] == base) ? vtx[1] : vtx[0];
    if (skeleton->events[u].get_deficiency() < -std::numeric_limits<double>::epsilon()) N.insert(u);
  }

  for(it=N.begin(); it!=N.end(); ++it) {
    n = *it;
    m = skeleton->RND->irandom(vx);
    s.clear();
    s.insert(base);
    s.insert(n);
    qt = skeleton->index_table[1].find(s);
    skeleton->simplices[1][qt->second].deactivate();
    skeleton->simplex_addition(m,n,-1);
  }
  return true;
}

template<class kind1,class kind2>
bool Spacetime<kind1,kind2>::perforation(int base,int d)
{
  int i,j,k,n,nd;
  bool good,found;
  std::set<int> candidates;
  std::vector<std::set<int> > S;

  if (base >= 0) {
    // This call of the perforation operator is localized
    int n1 = skeleton->vertex_dimension(base);
    if (n1 < 2) return false;
    d = skeleton->RND->irandom(2,n1);
    nd = (signed) skeleton->simplices[d].size();
    for(i=0; i<nd; ++i) {
      if (!skeleton->active_simplex(d,i)) continue;
      if (!skeleton->simplices[d][i].contains(base)) continue;
      // We need to check now if each of this simplex's 1+d faces is
      // also the face of another d-simplex
      good = true;
      skeleton->simplices[d][i].get_faces(S);
      for(j=0; j<1+d; ++j) {
        found = false;
        for(k=0; k<nd; ++k) {
          if (k == i || !skeleton->active_simplex(d,k)) continue;
          if (skeleton->simplices[d][k].face(S[j])) {
            found = true;
            break;
          }
        }
        if (!found) {
          good = false;
          break;
        }
      }
      if (good) candidates.insert(i);
    }
  }
  else {
    if (d < 2 || d >= Complex<kind1>::ND) return false;
    if (skeleton->simplices[d].empty()) return false;
    nd = (signed) skeleton->simplices[d].size();
    for(i=0; i<nd; ++i) {
      if (!skeleton->active_simplex(d,i)) continue;
      // We need to check now if each of this simplex's 1+d faces is
      // also the face of another d-simplex
      good = true;
      skeleton->simplices[d][i].get_faces(S);
      for(j=0; j<1+d; ++j) {
        found = false;
        for(k=0; k<nd; ++k) {
          if (k == i || !skeleton->active_simplex(d,k)) continue;
          if (skeleton->simplices[d][k].face(S[j])) {
            found = true;
            break;
          }
        }
        if (!found) {
          good = false;
          break;
        }
      }
      if (good) candidates.insert(i);
    }
  }
  if (candidates.empty()) return false;
#ifdef VERBOSE
  std::cout << "Perforating a " << d << "-simplex." << std::endl;
#endif
  n = skeleton->RND->irandom(candidates);
  skeleton->simplex_deletion(d,n);
  return true;
}

template<class kind1,class kind2>
bool Spacetime<kind1,kind2>::inflation(int base,double creativity)
{
  // Performs an inflation on the event base 
  int i,k,n1,na,delta,its = 0;
  Event<kind1> vt;
  bool success;
  std::set<int>::const_iterator it;
  SYNARMOSMA::hash_map::const_iterator qt;
  std::set<int> vx,candidates;

  if (base >= 0) {
    n1 = skeleton->vertex_dimension(base);
    if (n1 == Complex<kind1>::ND) return false;
    int j,vtx[2];
    std::set<int>::const_iterator jt;
    std::set<int> M,N;
    double tau;
    const int ne = (signed) skeleton->simplices[1].size();

    // We want to inflate this to a d-simplex, where d is between
    // 1+n1 and ND - the greater the magnitude of v's deficiency,
    // the higher the dimension d should be...
    tau = std::abs(skeleton->events[base].get_deficiency()) + 25.0;
    delta = int(double(Complex<kind1>::ND-1-n1)*1.0/(1.0 + std::exp(tau)) + double(1 + n1));

    for(i=0; i<(signed) skeleton->simplices[n1].size(); ++i) {
      if (!skeleton->active_simplex(n1,i)) continue;
      if (skeleton->simplices[n1][i].contains(base)) candidates.insert(i);
    }
    if (candidates.empty()) return false;

    for(i=0; i<ne; ++i) {
      if (!skeleton->active_simplex(1,i)) continue;
      skeleton->simplices[1][i].get_vertices(vtx);
      if (base != vtx[0] && base != vtx[1]) continue;
      j = (vtx[0] == base) ? vtx[1] : vtx[0];
      N.insert(j);
    }
    skeleton->simplices[n1][skeleton->RND->irandom(candidates)].get_vertices(vx);
    for(it=N.begin(); it!=N.end(); ++it) {
      jt = std::find(vx.begin(),vx.end(),*it);
      if (jt == vx.end()) {
        if (skeleton->events[*it].get_deficiency() < -std::numeric_limits<double>::epsilon()) M.insert(*it);
      }
    }
    if (M.empty()) creativity = 1.0;
    do {
      its++;
      if (its == Complex<kind1>::ND) creativity = 1.0;
      if (skeleton->RND->drandom() < creativity) {
        // Create a new event...
        k = event_addition(base);
      }
      else {
        k = skeleton->RND->irandom(M);
      }
      vx.insert(k);
    } while((signed) vx.size() < (1 + delta));
#ifdef VERBOSE
    std::cout << "Inflated a " << n1 << "-simplex into a " << delta << "-simplex based on event " << base << std::endl;
#endif
  }
  else {
    if (skeleton->dimension() == 0) {
      n1 = 0;
    }
    else if (skeleton->dimension() < (signed) geometry->dimension()) {
      n1 = 1;
    }
    else {
      n1 = 1 + skeleton->RND->irandom(skeleton->dimension());
    }
    for(i=0; i<(signed) skeleton->events.size(); ++i) {
      if (skeleton->active_event(i)) candidates.insert(i);
    }
    na = (signed) candidates.size();
    if (n1 == 0) {
      delta = skeleton->RND->irandom(2,2*geometry->dimension());
      vx.insert(skeleton->RND->irandom(candidates));
    }
    else {
      candidates.clear();
      for(i=0; i<(signed) skeleton->simplices[n1].size(); ++i) {
        if (skeleton->active_simplex(n1,i)) candidates.insert(i);
      }
      if (candidates.empty()) return false;
      skeleton->simplices[n1][skeleton->RND->irandom(candidates)].get_vertices(vx);
      delta = n1 + skeleton->RND->irandom(1,geometry->dimension());
    }
    if (double(delta)/double(na) > 0.25) creativity = 1.0;
    if (double(vx.size())/double(na) > 0.9) creativity = 1.0;
    do {
      if (skeleton->RND->drandom() < creativity) {
        // Create a new event...
        k = event_addition(vx);
      }
      else {
        // Grab an existing event...
        its = 0;
        success = false;
        do {
          // We should perhaps alter this to favour events that are
          // few hops away
          k = skeleton->RND->irandom(skeleton->events.size());
          if (!skeleton->active_event(k)) continue;
          it = std::find(vx.begin(),vx.end(),k);
          if (it == vx.end()) {
            success = true;
            break;
          }
          its++;
        } while(its < (signed) skeleton->events.size());
        if (success == false) {
          k = (signed) skeleton->events.size();
          skeleton->events.push_back(vt);
        }
      }
      vx.insert(k);
    } while((signed) vx.size() < (1+delta));
#ifdef VERBOSE
    std::cout << "Inflated a " << n1 << "-simplex into a " << delta << "-simplex" << std::endl;
#endif
  }
  skeleton->simplex_addition(vx,-1);
  return true;
}
