#include "spacetime.h"

using namespace DIAPLEXIS;

bool Spacetime::stellar_deletion(int base,int sheet)
{
  // If this vertex is already part of a d-simplex, d >= 2, then 
  // this operation is pointless... 
  if (vertex_dimension(base,sheet) >= 2) return false;
  int m,vx[2];
  std::set<int> nset;
  std::set<int>::const_iterator it;

  for(it=events[base].entourage.begin(); it!=events[base].entourage.end(); ++it) {
    if (simplices[1][*it].active(sheet)) {
      simplices[1][*it].get_vertices(vx);
      m = (vx[0] == base) ? vx[1] : vx[0];
      nset.insert(m);
    }
  }
  if (nset.size() != 3) return false;    
  vertex_deletion(base,sheet);
  simplex_addition(nset,sheet);
  regularization(true,sheet);
  return true;
}

bool Spacetime::foliation_m(int base,int sheet)
{
  int i,p,n1,n2,vx[2];
  std::set<int> locus,candidates,s1;
  SYNARMOSMA::hash_map::iterator qt;
  const int ne = (signed) simplices[1].size();

  locus.insert(sheet);

  for(i=0; i<ne; ++i) {
    if (!simplices[1][i].active(sheet)) continue;
    simplices[1][i].get_vertices(vx);
    if (vx[0] != base && vx[1] != base) continue;
    p = (vx[0] == base) ? vx[1] : vx[0];
    if (std::abs(events[p].deficiency) < std::numeric_limits<double>::epsilon()) continue;
    candidates.insert(p);
  }
  if (candidates.size() < 2) return false;
  n1 = RND->irandom(candidates);
  do {
    n2 = RND->irandom(candidates);
    if (n2 != n1) break;
  } while(true);


  std::set<int> leaf;
  leaf.insert(base);
  leaf.insert(n1);
  leaf.insert(n2);
  Simplex S(leaf,locus);
  qt = index_table[2].find(S.vertices);
  if (qt == index_table[2].end()) {
    simplices[2].push_back(S);
    index_table[2][S.vertices] = (signed) simplices[2].size() - 1;
    return true;
  }
  if (simplices[2][qt->second].active(sheet)) return false;
  simplices[2][qt->second].set_active(sheet);
  s1.insert(n1); s1.insert(n2);
  qt = index_table[1].find(s1);
  if (!simplices[1][qt->second].active(sheet)) {
    simplices[1][qt->second].set_active(sheet);
  }
  codex[sheet].vx_delta.insert(n1);
  codex[sheet].vx_delta.insert(n2);
  return true;
}

bool Spacetime::fusion_m(int base,int sheet)
{
  int i,u,vx[2];
  std::set<int> candidates;
  const int ne = (signed) simplices[1].size();

  for(i=0; i<ne; ++i) {
    if (!simplices[1][i].active(sheet)) continue;
    simplices[1][i].get_vertices(vx);
    if (vx[0] != base && vx[1] != base) continue;
    u = (vx[0] == base) ? vx[1] : vx[0];
    candidates.insert(u);
  }
  if (candidates.empty()) return false;
  u = RND->irandom(candidates);
#ifdef VERBOSE
  std::cout << "Fusing vertices: " << u << " => " << base << std::endl;
#endif
  vertex_fusion(base,u,sheet);
  return true;
}

bool Spacetime::fission(int base,double density,int sheet)
{
  int i,p,q,n,vx[2];
  std::set<int> locus;
  Simplex S;
  const int ne = (signed) simplices[1].size();

  locus.insert(sheet);

  if (base >= 0) {
    n = 0;

    p = vertex_addition(base,sheet);

    for(i=0; i<ne; ++i) {
      if (!simplices[1][i].active(sheet)) continue;
      simplices[1][i].get_vertices(vx);
      if (vx[0] != base && vx[1] != base) continue;
      q = (vx[0] == base) ? vx[1] : vx[0];
      if (RND->drandom() < density) {
        S.initialize(q,p,locus);
        simplices[1].push_back(S);
        index_table[1][S.vertices] = (signed) simplices[1].size() - 1;
        n++;
      }
    }
    S.initialize(base,p,locus);
    simplices[1].push_back(S);
    index_table[1][S.vertices] = (signed) simplices[1].size() - 1;
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

    if (dimension(sheet) < 1) return false;

    if (RND->drandom() < 0.5) {
      // The simplest case: we just need to add a new vertex, clone the antecedent
      // vertex's one-dimensional entourage and lastly create a 1-simplex joining
      // the two vertices.
      int d;

      do {
        p = RND->irandom(events.size());
        if (events[p].active(sheet)) break;
      } while(true);
      q = vertex_addition(p,sheet);
      nsimplex.insert(p);
      nsimplex.insert(q);

      for(i=0; i<ne; ++i) {
        if (!simplices[1][i].active(sheet)) continue;
        simplices[1][i].get_vertices(vx);
        if (vx[0] != p && vx[1] != p) continue;
        n = (vx[0] == p) ? vx[1] : vx[0];
        nsimplex.insert(n);
      }
      d = (signed) nsimplex.size() - 1;
      if (d >= Spacetime::ND) {
        std::set<int> s1;
        d = RND->irandom(2,Spacetime::ND);
        do {
          n = RND->irandom(nsimplex);
          s1.insert(n);
        } while((signed) s1.size() < (1+d));
        nsimplex = s1;
      }
      S.initialize(nsimplex,locus);
#ifdef VERBOSE
      std::cout << "Vertex fission on " << p << " with new vertex " << q << " and simplex dimension " << d << std::endl;
#endif
      simplices[d].push_back(S);
      index_table[d][S.vertices] = simplices[d].size() - 1;
    }
    else {
      // The plot thickens: we must create a new 1-simplex out of an existing one, perhaps by forming
      // a pair of 2-simplices: (w1,u_new,v_new) and (w2,u_new,v_new)
      std::set<int> antecedent;

      do {
        n = RND->irandom(simplices[1].size());
        if (simplices[1][n].active(sheet)) break;
      } while(true);
      simplices[1][n].get_vertices(vx);
      antecedent.insert(vx[0]);
      antecedent.insert(vx[1]);
      p = vertex_addition(antecedent,sheet);
      antecedent.insert(p);
      q = vertex_addition(antecedent,sheet);

      nsimplex.insert(vx[0]);
      nsimplex.insert(p);
      nsimplex.insert(q);
      S.initialize(nsimplex,locus);
      simplices[2].push_back(S);
      index_table[2][S.vertices] = simplices[2].size() - 1;
      nsimplex.clear();

      nsimplex.insert(vx[1]);
      nsimplex.insert(p);
      nsimplex.insert(q);
      S.initialize(nsimplex,locus);
      simplices[2].push_back(S);
      index_table[2][S.vertices] = simplices[2].size() - 1;
    }
    return true;
  }
}

bool Spacetime::circumvolution(int sheet)
{
  // This method attempts a circumvolution using boundary edges, it is
  // normally only called when the global energy reaches a critical threshold.
  int i,j,k,l,u[2],w[2];
  std::set<int> edge_set,S;
  std::set<int>::const_iterator it,jt;
  SYNARMOSMA::hash_map::const_iterator qt;
  const int D = dimension(sheet);
  int vx[D];

  for(i=0; i<(signed) simplices[D].size(); ++i) {
    if (!simplices[D][i].active(sheet)) continue;
    for(j=0; j<=D; ++j) {
      qt = index_table[D-1].find(simplices[D][i].faces[j]);
      if (simplices[D-1][qt->second].entourage.size() == 1) {
        simplices[D-1][qt->second].get_vertices(vx);
        for(k=0; k<D; ++k) {
          for(l=k+1; l<D; ++l) {
            S.clear();
            S.insert(vx[k]);
            S.insert(vx[l]);
            qt = index_table[1].find(S);
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

      simplices[1][n1].get_vertices(w);
      simplices[1][n2].get_vertices(u);

      d = combinatorial_distance(w[0],u[0],-1);
      if (d <= 2) continue;

      d = combinatorial_distance(w[1],u[0],-1);
      if (d <= 2) continue;

      d = combinatorial_distance(w[0],u[1],-1);
      if (d <= 2) continue;

      d = combinatorial_distance(w[1],u[1],-1);
      if (d <= 2) continue;

      candidates.push_back(n1);
      candidates.push_back(n2);
    }
  }
  if (candidates.empty()) return false;
  np = candidates.size()/2;
  i = RND->irandom(np);
  n1 = candidates[2*i];
  n2 = candidates[2*i+1];
  simplices[1][n1].get_vertices(w);
  simplices[1][n2].get_vertices(u);
  codex[sheet].vx_delta.insert(u[0]);
  codex[sheet].vx_delta.insert(u[1]);
  if (RND->drandom() < 0.5) {
    vertex_fusion(w[0],u[0],sheet);
    vertex_fusion(w[1],u[1],sheet);
  }
  else {
    vertex_fusion(w[0],u[1],sheet);
    vertex_fusion(w[1],u[0],sheet);
  }
  regularization(false,sheet);
  return true;
}

bool Spacetime::circumvolution(int base,int sheet)
{
  // This method seeks to fuse together two d-simplices, one of
  // which contains the vertex base
  int i,d,nd,s1,s2;
  Simplex S;
  std::set<int> candidates;
  std::vector<int> order;

  if (base >= 0) {
    nd = vertex_dimension(base,sheet);
    if (nd < 1) return false;
    d = RND->irandom(1,nd);
    nd = (signed) simplices[d].size();
    for(i=0; i<nd; ++i) {
      if (!simplices[d][i].active(sheet)) continue;
      if (simplices[d][i].contains(base)) candidates.insert(i);
    }
    s1 = RND->irandom(candidates);
    candidates.clear();
    for(i=0; i<nd; ++i) {
      if (!simplices[d][i].active(sheet)) continue;
      if (i == s1) continue;
      S = simplices[d][s1] ^ simplices[d][i];
      if (S.empty()) candidates.insert(i);
    }
    if (candidates.empty()) return false;
    int j,vx1[1+d],vx2[1+d];
    double D1,delta;
    std::set<int> slist;
    std::set<int>::const_iterator it;

    simplices[d][s1].get_vertices(vx1);
    for(it=candidates.begin(); it!=candidates.end(); ++it) {
      simplices[d][*it].get_vertices(vx2);
      D1 = 0.0;
      for(i=0; i<=d; ++i) {
        for(j=1+i; j<=d; ++j) {
          delta = geometry->get_squared_distance(vx1[i],vx2[j],true);
          if (delta > D1) D1 = delta;
        }
      }
      if (D1 < 2.5) slist.insert(*it);
    }
    if (slist.empty()) return false;
    s2 = RND->irandom(slist);
  }
  else {
    do {
      d = RND->irandom(1,dimension(sheet));
      for(i=0; i<(signed) simplices[d].size(); ++i) {
        if (simplices[d][i].active(sheet)) candidates.insert(i);
      }
      if (candidates.size() > 3) break;
    } while(true);
    do {
      s1 = RND->irandom(candidates);
      s2 = RND->irandom(candidates);
      if (s1 == s2) continue;
      // Now check to make sure that the intersection of
      // these two d-simplices is null...
      S = simplices[d][s1] ^ simplices[d][s2];
      if (S.empty()) break;
    } while(true);
  }
#ifdef VERBOSE
  std::cout << "Circumvolving the " << d << "-simplices: " << SYNARMOSMA::make_key(simplices[d][s2].vertices) << " => " << SYNARMOSMA::make_key(simplices[d][s1].vertices) << std::endl;
#endif
  int v1[d+1],v2[d+1];
  for(i=0; i<=d; ++i) {
    order.push_back(i);
  }
  std::random_shuffle(order.begin(),order.end());
  simplices[d][s1].get_vertices(v1);
  simplices[d][s2].get_vertices(v2);

  for(i=0; i<=d; ++i) {
    codex[sheet].vx_delta.insert(v2[order[i]]);
    vertex_fusion(v1[i],v2[order[i]],sheet);
  }
  return true;
}

bool Spacetime::expansion(int base,double creativity,int sheet)
{
  // Create an entirely new d-simplex
  int i,j,d,k,nv,novum = 0;
  bool success;
  std::set<int> vx;
  std::set<int>::const_iterator it;
  SYNARMOSMA::hash_map::const_iterator qt;

  if (base >= 0) {
    int n1 = vertex_dimension(base,sheet);
    if (n1 == Spacetime::ND) return false;
    int u,vtx[2],its = 0;
    std::set<int> M;
    double tau = std::abs(events[base].deficiency) + 25.0;
    const int ne = (signed) simplices[1].size();

    d = int(double(Spacetime::ND-1-n1)*1.0/(1.0 + std::exp(tau)) + double(1 + n1));

    for(i=0; i<ne; ++i) {
      if (!simplices[1][i].active(sheet)) continue;
      simplices[1][i].get_vertices(vtx);
      if (vtx[0] != base && vtx[1] != base) continue;
      u = (vtx[0] == base) ? vtx[1] : vtx[0];
      if (events[u].deficiency < -std::numeric_limits<double>::epsilon()) M.insert(u);
    }
    if (M.empty()) creativity = 1.0;
    vx.insert(base);

    do {
      its++;
      if (its == Spacetime::ND) creativity = 1.0;
      if (RND->drandom() < creativity) {
        // Create a new vertex...
        k = vertex_addition(base,sheet);
        novum++;
      }
      else {
        k = RND->irandom(M);
      }
      vx.insert(k);
    } while((signed) vx.size() < (1+d));
  }
  else {
    i = dimension(sheet);
    d = (i <= 2) ? 3 : RND->irandom(2,i);
    if (double(d+1)/double(events.size()) > 0.8) creativity = 1.0;
    do {
      if (RND->drandom() < creativity) {
        // Create a new vertex...
        k = vertex_addition(vx,sheet);
        novum++;
      }
      else {
        // Grab an existing vertex...
        j = 0;
        nv = (signed) events.size();
        success = false;
        do {
          k = RND->irandom(nv);
          it = std::find(vx.begin(),vx.end(),k);
          if (it == vx.end()) {
            events[k].set_active(sheet);
            codex[sheet].vx_delta.insert(k);
            success = true;
            break;
          }
          ++j;
        } while(j < nv);
        if (success == false) {
          k = vertex_addition(vx,sheet);
          novum++;
        }
      }
      vx.insert(k);
    } while((signed) vx.size() < (1+d));
  }
#ifdef VERBOSE
  std::cout << "Created a " << d << "-simplex with " << novum << " new vertices." << std::endl;
#endif
  simplex_addition(vx,sheet);
  return true;
}

bool Spacetime::expansion(int base,int sheet)
{
  int i,u,d,m,vtx[2],n = vertex_dimension(base,sheet);
  if (n == Spacetime::ND) return false;
  std::set<int> vx,s,N,locus;
  std::set<int>::const_iterator it;
  SYNARMOSMA::hash_map::const_iterator qt;
  Simplex S;
  const int ne = (signed) simplices[1].size();
  double tau = std::abs(events[base].deficiency) + 25.0;

  locus.insert(sheet);

  d = int(double(Spacetime::ND-1-n)*1.0/(1.0 + std::exp(tau)) + double(1 + n));

  vx.insert(base);
  for(i=0; i<d; ++i) {
    // Create a new vertex...
    m = vertex_addition(base,sheet);
    vx.insert(m);
  }
#ifdef VERBOSE
  std::cout << "Created a " << d << "-simplex with base " << base << std::endl;
#endif
  simplex_addition(vx,sheet);
  for(i=0; i<ne; ++i) {
    if (!simplices[1][i].active(sheet)) continue;
    simplices[1][i].get_vertices(vtx);
    if (vtx[0] != base && vtx[1] != base) continue;
    u = (vtx[0] == base) ? vtx[1] : vtx[0];
    if (events[u].deficiency < -std::numeric_limits<double>::epsilon()) N.insert(u);
  }
  for(it=N.begin(); it!=N.end(); ++it) {
    n = *it;
    m = RND->irandom(vx);
    s.clear();
    s.insert(base);
    s.insert(n);
    qt = index_table[1].find(s);
    simplices[1][qt->second].set_inactive(sheet);
    S.initialize(m,n,locus);
    qt = index_table[1].find(S.vertices);
    if (qt == index_table[1].end()) {
      simplices[1].push_back(S);
      index_table[1][S.vertices] = simplices[1].size() - 1;
    }
    else {
      simplices[1][qt->second].set_active(sheet); 
    }
  }
  return true;
}

bool Spacetime::perforation(int base,int d,int sheet)
{
  int i,j,k,n,nd;
  bool good,found;
  std::set<int> candidates,S;

  if (base >= 0) {
    // This call of the perforation operator is localized
    int n1 = vertex_dimension(base,sheet);
    if (n1 < 2) return false;
    d = RND->irandom(2,n1);
    nd = (signed) simplices[d].size();
    for(i=0; i<nd; ++i) {
      if (!simplices[d][i].active(sheet)) continue;
      if (!simplices[d][i].contains(base)) continue;
      // We need to check now if each of this simplex's 1+d faces is
      // also the face of another d-simplex
      good = true;
      for(j=0; j<1+d; ++j) {
        found = false;
        S = simplices[d][i].faces[j];
        for(k=0; k<nd; ++k) {
          if (k == i || !simplices[d][k].active(sheet)) continue;
          if (simplices[d][k].face(S)) {
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
    if (d < 2 || d >= Spacetime::ND) return false;
    if (simplices[d].empty()) return false;
    nd = (signed) simplices[d].size();
    for(i=0; i<nd; ++i) {
      if (!simplices[d][i].active(sheet)) continue;
      // We need to check now if each of this simplex's 1+d faces is
      // also the face of another d-simplex
      good = true;
      for(j=0; j<1+d; ++j) {
        found = false;
        S = simplices[d][i].faces[j];
        for(k=0; k<nd; ++k) {
          if (k == i || !simplices[d][k].active(sheet)) continue;
          if (simplices[d][k].face(S)) {
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
  n = RND->irandom(candidates);
  simplex_deletion(d,n,sheet);
  return true;
}

bool Spacetime::inflation(int base,double creativity,int sheet)
{
  // Performs an inflation on the vertex base with colour sheet
  int i,k,n1,na,delta,its = 0;
  bool success;
  std::set<int> vx,candidates;
  std::set<int>::const_iterator it;
  Event vt;
  SYNARMOSMA::hash_map::const_iterator qt;

  vt.set_active(sheet);

  if (base >= 0) {
    n1 = vertex_dimension(base,sheet);
    if (n1 == Spacetime::ND) return false;
    int j,vtx[2];
    std::set<int>::const_iterator jt;
    std::set<int> M,N;
    double tau;
    const int ne = (signed) simplices[1].size();

    // We want to inflate this to a d-simplex, where d is between
    // 1+n1 and ND - the greater the magnitude of v's deficiency,
    // the higher the dimension d should be...
    tau = std::abs(events[base].deficiency) + 25.0;
    delta = int(double(Spacetime::ND-1-n1)*1.0/(1.0 + std::exp(tau)) + double(1 + n1));

    for(i=0; i<(signed) simplices[n1].size(); ++i) {
      if (!simplices[n1][i].active(sheet)) continue;
      if (simplices[n1][i].contains(base)) candidates.insert(i);
    }
    if (candidates.empty()) return false;

    for(i=0; i<ne; ++i) {
      if (!simplices[1][i].active(sheet)) continue;
      simplices[1][i].get_vertices(vtx);
      if (base != vtx[0] && base != vtx[1]) continue;
      j = (base == vtx[0]) ? vtx[1] : vtx[0];
      N.insert(j);
    }
    vx = simplices[n1][RND->irandom(candidates)].vertices;
    for(it=N.begin(); it!=N.end(); ++it) {
      jt = std::find(vx.begin(),vx.end(),*it);
      if (jt == vx.end()) {
        if (events[*it].deficiency < -std::numeric_limits<double>::epsilon()) M.insert(*it);
      }
    }
    if (M.empty()) creativity = 1.0;
    do {
      its++;
      if (its == ND) creativity = 1.0;
      if (RND->drandom() < creativity) {
        // Create a new vertex...
        k = vertex_addition(base,sheet);
      }
      else {
        k = RND->irandom(M);
      }
      vx.insert(k);
    } while((signed) vx.size() < (1 + delta));
#ifdef VERBOSE
    std::cout << "Inflated a " << n1 << "-simplex into a " << delta << "-simplex based on vertex " << base << std::endl;
#endif
  }
  else {
    if (dimension(sheet) == 0) {
      n1 = 0;
    }
    else if (dimension(sheet) < (signed) geometry->dimension()) {
      n1 = 1;
    }
    else {
      n1 = 1 + RND->irandom(dimension(sheet));
    }
    for(i=0; i<(signed) events.size(); ++i) {
      if (events[i].active(sheet)) candidates.insert(i);
    }
    na = (signed) candidates.size();
    if (n1 == 0) {
      delta = RND->irandom(2,2*geometry->dimension());
      vx.insert(RND->irandom(candidates));
    }
    else {
      candidates.clear();
      for(i=0; i<(signed) simplices[n1].size(); ++i) {
        if (simplices[n1][i].active(sheet)) candidates.insert(i);
      }
      if (candidates.empty()) return false;
      vx = simplices[n1][RND->irandom(candidates)].vertices;
      delta = n1 + RND->irandom(1,geometry->dimension());
    }
    if (double(delta)/double(na) > 0.25) creativity = 1.0;
    if (double(vx.size())/double(na) > 0.9) creativity = 1.0;
    do {
      if (RND->drandom() < creativity) {
        // Create a new vertex...
        k = vertex_addition(vx,sheet);
      }
      else {
        // Grab an existing vertex...
        its = 0;
        success = false;
        do {
          // We should perhaps alter this to favour vertices that are
          // few hops away
          k = RND->irandom(events.size());
          if (!events[k].active(sheet)) continue;
          it = std::find(vx.begin(),vx.end(),k);
          if (it == vx.end()) {
            success = true;
            break;
          }
          its++;
        } while(its < (signed) events.size());
        if (success == false) {
          k = (signed) events.size();
          events.push_back(vt);
        }
      }
      vx.insert(k);
    } while((signed) vx.size() < (1+delta));
#ifdef VERBOSE
    std::cout << "Inflated a " << n1 << "-simplex into a " << delta << "-simplex" << std::endl;
#endif
  }
  simplex_addition(vx,sheet);
  return true;
}

