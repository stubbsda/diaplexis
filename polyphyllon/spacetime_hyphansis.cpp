#include "spacetime.h"

using namespace DIAPLEXIS;

const int Spacetime::N_EXP;
const int Spacetime::N_IMP;
const std::string Spacetime::EXP_OP[] = {"D","Ux","Ox","R","C","N","A","G","Sg","Sm","Y"};
const std::string Spacetime::IMP_OP[] = {"I","Um","Om","E","F","P","V","Δ"};

bool Spacetime::edge_parity_mutation(int base,int sheet)
{
  int n;
  std::set<int> candidates;
  std::set<int>::const_iterator it;

  for(it=events[base].entourage.begin(); it!=events[base].entourage.end(); ++it) {
    if (simplices[1][*it].active(sheet)) {
      candidates.insert(*it);
    }
  }
  if (candidates.empty()) return false;      
  n = RND->irandom(candidates);
  if (simplices[1][n].parity == 0) {
    simplices[1][n].parity = (RND->irandom(2) == 0) ? 1 : -1;
  }
  else {
    simplices[1][n].parity *= -1;
  }
  recompute_parity(n);
  return true;
}

bool Spacetime::edge_parity_mutation(int u,int v,int sheet)
{
  // This is the method used for dynamic hyphansis where there is a single call 
  // to recompute the orientation of all the higher-dimensional simplices, so no 
  // need to include one at the method's end. 
  int n;
  std::set<int> S;
  SYNARMOSMA::hash_map::const_iterator qt;

  S.insert(u); S.insert(v);
  qt = index_table[1].find(S);
  if (qt == index_table[1].end()) return false;
  if (!simplices[1][qt->second].active(sheet)) return false;
  n = qt->second;
  if (simplices[1][n].parity == 0) {
    simplices[1][n].parity = (RND->irandom(2) == 0) ? 1 : -1;
  }
  else {
    simplices[1][n].parity *= -1;
  }
  return true;
}

void Spacetime::simplicial_implication(int base,int sheet) const
{
  // This method will calculate all of the n-simplices (n > 1) that are "implied" by
  // the vertex base and its neighbours (via their mutual edges) and list them, as well
  // as checking to see if they already exist in the spacetime complex.
  // One wrinkle with this method is that it can only be used when the "entourage" and
  // "neighbours" properties of the events are well-defined... so after the calls to
  // "regularize".
  int i,j,k,n,l,d,M,nsimp,nfound,w1,w2;
  std::vector<int> C,vx;
  std::vector<std::set<int> >* implied_simplex;
  std::vector<std::set<int> >::const_iterator it;
  bool failure;
  SYNARMOSMA::hash_map::const_iterator qt;
  std::set<int> S,sv;

  if (sheet == -1) {
    d = (signed) events[base].entourage.size();
    M = 1 + d;
    implied_simplex = new std::vector<std::set<int> >[M+1];

    S = events[base].neighbours;
    S.insert(base);

    for(i=M; i>2; --i) {
      // See how many i-dimensional simplices exist among the relations between v and its
      // neighbours, so there should be (M choose i) such possible i-simplices
      n = SYNARMOSMA::combinations(S,i,C);
      for(l=0; l<n; ++l) {
        // Grab the i elements from S and put them into the vector vx...
        for(j=0; j<i; ++j) {
          vx.push_back(C[l*i+j]);
        }
        failure = false;
        for(j=0; j<i; ++j) {
          for(k=1+j; k<i; ++k) {
            w1 = vx[j];
            w2 = vx[k];
            sv.clear();
            sv.insert(w1);
            sv.insert(w2);
            qt = index_table[1].find(sv);
            if (qt == index_table[1].end()) {
              failure = true;
              break;
            }
            else {
              if (!simplices[1][qt->second].active()) {
                failure = true;
                break;
              }
            }
          }
          if (failure) break;
        }
        if (!failure) {
          sv.clear();
          for(j=0; j<i; ++j) {
            sv.insert(vx[j]);
          }
          implied_simplex[i].push_back(sv);
        }
        vx.clear();
      }
    }
  }
  else {
#ifdef DEBUG
    assert(events[base].active(sheet));
#endif
    std::set<int>::const_iterator jt;

    d = 0;
    for(jt=events[base].entourage.begin(); jt!=events[base].entourage.end(); ++jt) {
      if (!simplices[1][*jt].active(sheet)) continue;
      d++;
    }
    M = 1 + d;
    implied_simplex = new std::vector<std::set<int> >[M+1];

    for(jt=events[base].neighbours.begin(); jt!=events[base].neighbours.end(); ++jt) {
      if (!events[*jt].active(sheet)) continue;
      S.insert(*jt);
    }
    S.insert(base);

    for(i=M; i>2; --i) {
      // See how many i-dimensional simplices exist among the relations between v and its
      // neighbours, so there should be (M choose i) such possible i-simplices
      n = SYNARMOSMA::combinations(S,i,C);
      for(l=0; l<n; ++l) {
        // Grab the i elements from S and put them into the vector vx...
        for(j=0; j<i; ++j) {
          vx.push_back(C[l*i+j]);
        }
        failure = false;
        for(j=0; j<i; ++j) {
          for(k=1+j; k<i; ++k) {
            w1 = vx[j];
            w2 = vx[k];
            sv.clear();
            sv.insert(w1);
            sv.insert(w2);
            qt = index_table[1].find(sv);
            if (qt == index_table[1].end()) {
              failure = true;
              break;
            }
            else {
              if (!simplices[1][qt->second].active(sheet)) {
                failure = true;
                break;
              }
            }
          }
          if (failure) break;
        }
        if (!failure) {
          sv.clear();
          for(j=0; j<i; ++j) {
            sv.insert(vx[j]);
          }
          implied_simplex[i].push_back(sv);
        }
        vx.clear();
      }
    }
  }
  // So how many of these implied simplices actually exist? And having found one or more which officially don't
  // exist, what to do with them? If deleting a single edge can eliminate an entire tower of dependent n-simplices
  // (n > 1) then shouldn't recreating this same edge cause the tower of n-simplices to be restored?
  for(i=M; i>2; --i) {
    nsimp = 0;
    nfound = 0;
    for(it=implied_simplex[i].begin(); it!=implied_simplex[i].end(); ++it) {
      qt = index_table[i-1].find(*it);
      if (qt != index_table[i-1].end()) {
        if (simplices[i-1][qt->second].active()) nfound++;
      }
      nsimp++;
    }
#ifdef VERBOSE
    std::cout << "There are " << nsimp << " implied " << i-1 << "-simplices of which " << nfound << " already exist in the spacetime complex." << std::endl;
#endif
  }
  delete[] implied_simplex;
}

bool Spacetime::unravel(int base,int sheet)
{
  int i,j,c,vx[2],in_max = -1;
  unsigned int n1,n2;
  std::set<int>::const_iterator it;
  double vf,vf_max = 0.0;

  if (base >= 0) {
    n1 = (signed) events[base].entourage.size();
    for(it=events[base].entourage.begin(); it!=events[base].entourage.end(); ++it) {
      j = *it;
      if (!simplices[1][j].active(sheet)) continue;
      vf = double(simplices[1][j].entourage.size());
      if (vf > vf_max) {
        vf_max = vf;
        in_max = j;
      }
    }
  }
  else {
    for(i=0; i<(signed) simplices[1].size(); ++i) {
      if (!simplices[1][i].active(sheet)) continue;
      vf = double(simplices[1][i].entourage.size());
      if (vf > vf_max) {
        vf_max = vf;
        in_max = i;
      }
    }
  }
  if (in_max > -1) {
    // Eliminate the edge in_max...
    simplex_deletion(1,in_max,sheet);
    return true;
  }
  if (base >= 0) {
    n1 = events[base].entourage.size();
    if (n1 <= 2*geometry->dimension()) return false;
    for(it=events[base].entourage.begin(); it!=events[base].entourage.end(); ++it) {
      j = *it;
      if (!simplices[1][j].active(sheet)) continue;
      simplices[1][j].get_vertices(vx);
      c = (vx[0] == base) ? vx[1] : vx[0];
      n2 = events[c].entourage.size();
      if (n2 <= 2*geometry->dimension()) continue;
      vf = 0.5*double(n1 + n2 - 4*geometry->dimension());
      if (vf > vf_max) {
        vf_max = vf;
        in_max = j;
      }
    }
  }
  else {
    for(i=0; i<(signed) simplices[1].size(); ++i) {
      if (!simplices[1][i].active(sheet)) continue;
      simplices[1][i].get_vertices(vx);
      n1 = events[vx[0]].entourage.size();
      if (n1 <= 2*geometry->dimension()) continue;
      n2 = events[vx[1]].entourage.size();
      if (n2 <= 2*geometry->dimension()) continue;
      vf = 0.5*double(n1 + n2 - 4*geometry->dimension());
      if (vf > vf_max) {
        vf_max = vf;
        in_max = i;
      }
    }
  }
  if (in_max == -1) return false;
  // Eliminate the edge in_max...
  simplex_deletion(1,in_max,sheet);
  return true;
}

void Spacetime::compute_geometric_dependency(const std::set<int>& vx)
{
  // A method that takes a set of vertices whose coordinates have
  // been modified and determines which simplices will have their
  // modified property set to true...
  if (vx.empty()) return;
  int i,m,n;
  std::set<int> current,next;
  std::set<int>::const_iterator it,jt,kt;

#ifdef VERBOSE
  std::cout << "There are " << vx.size() << " vertices directly implicated." << std::endl;
#endif
  // We assume that the cardinality of vmodified is small relative
  // to the total number of vertices in the spacetime complex
  for(it=vx.begin(); it!=vx.end(); ++it) {
    n = *it;
    current = events[n].entourage;
    for(i=1; i<=Spacetime::ND; ++i) {
      for(jt=current.begin(); jt!=current.end(); ++jt) {
        m = *jt;
        simplices[i][m].modified = true;
        for(kt=simplices[i][m].entourage.begin(); kt!=simplices[i][m].entourage.end(); ++kt) {
          next.insert(*kt);
        }
      }
      if (next.empty()) break;
      current = next;
      next.clear();
    }
  }
}

void Spacetime::compute_topological_dependency(const std::set<int>& vx)
{
  int i,n,m,l,nhop;
  std::set<int> current,next;
  std::set<int>::const_iterator it,jt,kt;
  const int nv = (signed) events.size();
  int done[nv];

  for(i=0; i<nv; ++i) {
    done[i] = 0;
  }
  for(it=vx.begin(); it!=vx.end(); ++it) {
    n = *it;
    nhop = 0;
    // Every vertex within topological_radius hops of n is labelled as
    // modified
    current.insert(n);
    done[n] = 1;
    do {
      for(jt=current.begin(); jt!=current.end(); ++jt) {
        m = *jt;
        for(kt=events[m].neighbours.begin(); kt!=events[m].neighbours.end(); ++kt) {
          l = *kt;
          if (done[l] == 0) next.insert(l);
        }
      }
      if (next.empty()) break;
      for(jt=next.begin(); jt!=next.end(); ++jt) {
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
}

void Spacetime::compute_entourages(int sheet)
{
  int i,j,k,ns,vx[2];
  std::set<int> s;
  std::set<int>::const_iterator it;
  std::string fx;
  SYNARMOSMA::hash_map::const_iterator qt;
  const int ulimit = dimension(sheet);

  // What about removing items from the entourage of a d-simplex, when this item has
  // changed its ubiquity?
  for(i=1; i<=Spacetime::ND; ++i) {
    for(j=0; j<(signed) simplices[i].size(); ++j) {
      if (!simplices[i][j].active()) continue;
      for(it=simplices[i][j].entourage.begin(); it!=simplices[i][j].entourage.end(); ++it) {
        if (!simplices[i+1][*it].active()) continue;
        s.insert(*it);
      }
      simplices[i][j].entourage = s;
      s.clear();
    }
  }
  for(i=0; i<(signed) events.size(); ++i) {
    if (!events[i].active()) continue;
    for(it=events[i].entourage.begin(); it!=events[i].entourage.end(); ++it) {
      if (!simplices[1][*it].active()) continue;
      s.insert(*it);
    }
    events[i].entourage = s;
    s.clear();
  }

  if (sheet == -1) {
    std::set<int> ubiquity;

    for(i=ulimit; i>=2; i--) {
      ns = (signed) simplices[i].size();
      for(j=0; j<ns; ++j) {
        if (!simplices[i][j].active()) continue;
        simplices[i][j].get_ubiquity(ubiquity);
        for(k=0; k<1+i; ++k) {
          s = simplices[i][j].faces[k];
          qt = index_table[i-1].find(s);
          if (qt == index_table[i-1].end()) throw std::runtime_error("Missing entourage element!");
          for(it=ubiquity.begin(); it!=ubiquity.end(); ++it) {
            simplices[i-1][qt->second].set_active(*it);
          }
          simplices[i-1][qt->second].entourage.insert(j);
        }
      }
    }
    // Now the edges...
    ns = (signed) simplices[1].size();
    for(i=0; i<ns; ++i) {
      if (!simplices[1][i].active()) continue;
      simplices[1][i].get_ubiquity(ubiquity);
      simplices[1][i].get_vertices(vx);
      for(it=ubiquity.begin(); it!=ubiquity.end(); ++it) {
        events[vx[0]].set_active(*it);
        events[vx[1]].set_active(*it);
      }
      events[vx[0]].entourage.insert(i);
      events[vx[1]].entourage.insert(i);
    }
  }
  else {
    for(i=ulimit; i>=2; i--) {
      ns = (signed) simplices[i].size();
      for(j=0; j<ns; ++j) {
        if (!simplices[i][j].active(sheet)) continue;
        for(k=0; k<1+i; ++k) {
          s = simplices[i][j].faces[k];
          qt = index_table[i-1].find(s);
          if (qt == index_table[i-1].end()) throw std::runtime_error("Missing entourage element!");
          if (!simplices[i-1][qt->second].active(sheet)) simplices[i-1][qt->second].set_active(sheet);
          simplices[i-1][qt->second].entourage.insert(j);
        }
      }
    }
    // Now the edges...
    ns = (signed) simplices[1].size();
    for(i=0; i<ns; ++i) {
      if (!simplices[1][i].active(sheet)) continue;
      simplices[1][i].get_vertices(vx);

      if (!events[vx[0]].active(sheet)) events[vx[0]].set_active(sheet);
      events[vx[0]].entourage.insert(i);

      if (!events[vx[1]].active(sheet)) events[vx[1]].set_active(sheet);
      events[vx[1]].entourage.insert(i);
    }
  }
}

void Spacetime::compute_neighbours()
{
  int i,vx[2];
  std::set<int> empty;

  for(i=0; i<(signed) events.size(); ++i) {
    events[i].neighbours.clear();
  }
  for(i=0; i<(signed) simplices[1].size(); ++i) {
    if (!simplices[1][i].active()) continue;
    simplices[1][i].get_vertices(vx);
    events[vx[0]].neighbours.insert(vx[1]);
    events[vx[1]].neighbours.insert(vx[0]);
  }
}

void Spacetime::inversion()
{
  int i,j,p;
  std::set<int> locus,hold;
  std::set<int>::const_iterator it;
  Simplex S;
  const int nv = (signed) events.size();
  const int nt = (signed) codex.size();

  for(i=0; i<nv; ++i) {
    // hold = V / events[i].neighbours
    for(j=0; j<nv; ++j) {
      if (!events[j].active()) continue;
      if (j == i) continue;
      it = std::find(events[i].neighbours.begin(),events[i].neighbours.end(),j);
      if (it == events[i].neighbours.end()) hold.insert(j);
    }
    events[i].neighbours = hold;
    hold.clear();
  }
  simplices[1].clear();
  index_table[1].clear();
  for(i=0; i<nv; ++i) {
    for(it=events[i].neighbours.begin(); it!=events[i].neighbours.end(); ++it) {
      j = *it;
      if (i < j) {
        p = RND->irandom(nt);
        locus.insert(p);
        S.initialize(i,j,locus);
        simplices[1].push_back(S);
        index_table[1][S.vertices] = (signed) simplices[1].size() - 1;
        locus.erase(p);
      }
    }
  }
}

bool Spacetime::simplex_addition(const std::set<int>& S,int sheet)
{
  int i,j;
  std::set<int> locus,fc;
  std::set<int>::const_iterator it;
  std::vector<int> vec,vx;
  SYNARMOSMA::hash_map::const_iterator qt;
  const int d = (signed) S.size() - 1;

  locus.insert(sheet);

  Simplex s(S,locus);
  if (d == 1) {
    s.parity = 0;
    if (RND->drandom() < 0.2) {
      s.parity = (RND->irandom(2) == 0) ? 1 : -1;
    }
  }

  qt = index_table[d].find(S);
  if (qt == index_table[d].end()) {
    simplices[d].push_back(s);
    index_table[d][S] = simplices[d].size() - 1;
  }
  else {
    if (simplices[d][qt->second].active(sheet)) {
      return false;
    }
    else {
      simplices[d][qt->second].set_active(sheet);
    }
  }

#ifdef VERBOSE
  std::cout << "Adding a " << d << "-simplex to the spacetime complex..." << std::endl;
#endif

  for(it=S.begin(); it!=S.end(); ++it) {
    codex[sheet].vx_delta.insert(*it);
  }

  if (d == 1) {
    int vn[2];
    s.get_vertices(vn);
    events[vn[0]].neighbours.insert(vn[1]);
    events[vn[1]].neighbours.insert(vn[0]);
    return true;
  }
  for(it=S.begin(); it!=S.end(); ++it) {
    vx.push_back(*it);
  }

  for(i=d-1; i>=1; i--) {
    for(j=0; j<=i; ++j) {
      vec.push_back(j);
      fc.insert(vx[j]);
    }
    // Add this simplex...
    qt = index_table[i].find(fc);
    if (qt == index_table[i].end()) {
      simplices[i].push_back(Simplex(fc,locus));
      index_table[i][fc] = simplices[i].size() - 1;
    }
    else {
      simplices[i][qt->second].set_active(sheet);
    }
    fc.clear();
    while(SYNARMOSMA::next_combination(vec,1+d)) {
      for(j=0; j<=i; ++j) {
        fc.insert(vx[vec[j]]);
      }
      qt = index_table[i].find(fc);
      if (qt == index_table[i].end()) {
        simplices[i].push_back(Simplex(fc,locus));
        index_table[i][fc] = simplices[i].size() - 1;
      }
      else {
        simplices[i][qt->second].set_active(sheet);
      }
      fc.clear();
    }
    vec.clear();
  }
  simplicial_implication(sheet);
  return true;
}

void Spacetime::simplex_deletion(int d,int n,int sheet)
{
  std::set<int>::const_iterator it;
  std::set<int> parents;
  int i,dp1 = d + 1;

  if (sheet == -1) {
    for(i=0; i<(signed) codex.size(); ++i) {
      simplices[d][n].set_inactive(i);
    }    
  } 
  else {
    if (!simplices[d][n].active(sheet)) return;
    simplices[d][n].set_inactive(sheet);
    for(it=simplices[d][n].vertices.begin(); it!=simplices[d][n].vertices.end(); ++it) {
      codex[sheet].vx_delta.insert(*it);
    }
  }
  parents = simplices[d][n].entourage;
  for(it=parents.begin(); it!=parents.end(); ++it) {
    i = *it;
    simplex_deletion(dp1,i,sheet);
  }
}

void Spacetime::simplicial_implication(int sheet)
{
  int i,j,k,n,m,vx[2];
  Simplex S;
  std::string sx;
  std::set<int> colours;
  std::set<int>::const_iterator it;
  SYNARMOSMA::hash_map::const_iterator qt;
  const int ulimit = dimension(sheet);

  if (sheet == -1) {
    std::set<int> ubiquity;

    for(i=0; i<(signed) codex.size(); ++i) {
      if (dimension(i) < 0) continue;
      colours.insert(i);
    }
    for(i=ulimit; i>=2; i--) {
      n = (signed) simplices[i].size();
      m = (signed) simplices[i-1].size();
      for(j=0; j<n; ++j) {
        if (!simplices[i][j].active()) continue;
        simplices[i][j].get_ubiquity(ubiquity);
        for(k=0; k<1+i; ++k) {
          qt = index_table[i-1].find(simplices[i][j].faces[k]);
          if (qt == index_table[i-1].end()) {
            S.initialize(simplices[i][j].faces[k],ubiquity);
            simplices[i-1].push_back(S);
            index_table[i-1][S.vertices] = m;
            m++;
          }
          else {
            for(it=ubiquity.begin(); it!=ubiquity.end(); ++it) {
              simplices[i-1][qt->second].set_active(*it);
            }
          }
        }
      }
    }
    n = (signed) simplices[1].size();
    for(i=0; i<n; ++i) {
      if (!simplices[1][i].active()) continue;
      simplices[1][i].get_ubiquity(ubiquity);
      simplices[1][i].get_vertices(vx);
      for(it=ubiquity.begin(); it!=ubiquity.end(); ++it) {
        events[vx[0]].set_active(*it);
        events[vx[1]].set_active(*it);
      }
    }
  }
  else {
    std::set<int> locus;

    locus.insert(sheet);
    for(i=ulimit; i>=2; i--) {
      n = (signed) simplices[i].size();
      m = (signed) simplices[i-1].size();
      for(j=0; j<n; ++j) {
        if (!simplices[i][j].active(sheet)) continue;
        for(k=0; k<1+i; ++k) {
          qt = index_table[i-1].find(simplices[i][j].faces[k]);
          if (qt == index_table[i-1].end()) {
#ifdef VERBOSE
            std::cout << "Adding simplex with key " << SYNARMOSMA::make_key(simplices[i][j].faces[k]) << " to regularize the complex" << std::endl;
#endif
            S.initialize(simplices[i][j].faces[k],locus);
            simplices[i-1].push_back(S);
            index_table[i-1][S.vertices] = m;
            m++;
          }
          else {
            if (!simplices[i-1][qt->second].active(sheet)) {
#ifdef VERBOSE
              std::cout << "Restoring simplex with key " << SYNARMOSMA::make_key(simplices[i-1][qt->second].vertices) << " to regularize the complex" << std::endl;
#endif
              simplices[i-1][qt->second].set_active(sheet);
            }
          }
        }
      }
    }
    n = (signed) simplices[1].size();
    for(i=0; i<n; ++i) {
      if (!simplices[1][i].active(sheet)) continue;
      simplices[1][i].get_vertices(vx);
      if (!events[vx[0]].active(sheet)) {
#ifdef VERBOSE
        std::cout << "Restoring vertex " << vx[0] << " to regularize the complex" << std::endl;
#endif
        events[vx[0]].set_active(sheet);
      }
      if (!events[vx[1]].active(sheet)) {
#ifdef VERBOSE
        std::cout << "Restoring vertex " << vx[1] << " to regularize the complex" << std::endl;
#endif
        events[vx[1]].set_active(sheet);
      }
    }
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

int Spacetime::select_vertex(const std::vector<int>& candidates,double intensity,int sheet) const
{
  if (candidates.empty()) return -1;
  // The closer the intensity is to unity, the more we should try to choose an element 
  // of candidates close to the beginning
  int i,output,n = (signed) candidates.size();
  double cdeficit,tdeficit;
  std::vector<int> vcandidates;

  for(i=0; i<n; ++i) {
    if (!events[candidates[i]].active(sheet)) continue;
    vcandidates.push_back(candidates[i]);
  }
  if (vcandidates.empty()) return -1;
  if (vcandidates.size() == 1) return vcandidates[0];
  n = (signed) vcandidates.size();
  tdeficit = std::abs(events[vcandidates[0]].deficiency - events[vcandidates[n-1]].deficiency);
  for(i=0; i<n; ++i) {
    output = vcandidates[i];
    cdeficit = std::abs(events[output].deficiency - events[vcandidates[0]].deficiency);
    if (intensity <= cdeficit/tdeficit) break;
  }
  return output;
}

void Spacetime::musical_hyphansis(const std::vector<std::pair<int,double> >& candidates,int sheet)
{
  int i,j,v,its,opcount;
  double m_width = 1.0,x_width = 1.0;
  bool success = false;
  std::string line,op;
  std::stringstream opstring;
  std::vector<int> key_list,m_vertices,x_vertices,neutral_vertices,m_keys,x_keys;
  std::vector<std::string> elements;
  std::vector<double> pvalues;
  boost::char_separator<char> sp("/");
  const int nc = (signed) candidates.size();

  // Open the file containing the hyphantic score 
  std::ifstream mscore;
  mscore.exceptions(std::ifstream::badbit);
  try {
    mscore.open(hyphansis_score);
    // Now read the measure that corresponds to this iteration and sheet...
    while(mscore.good()) {
      getline(mscore,line);
      // Break the line up at the forward slash
      elements.clear();
      boost::tokenizer<boost::char_separator<char> > tok(line,sp);
      for(boost::tokenizer<boost::char_separator<char> >::iterator beg=tok.begin(); beg!=tok.end(); beg++) {
        elements.push_back(*beg);
      }
      if (elements.empty()) continue;
#ifdef DEBUG
      assert(elements.size() == 3);
#endif
      its = boost::lexical_cast<int>(elements[0]) - 1;
      if (its < iterations) continue;
      if (its > iterations) break;
      // So this is a line for this relaxation step, check if it is the right sheet/voice...
      v = boost::lexical_cast<int>(elements[1]);
      if (v != sheet) continue;
      // So, grab the piano key...
      key_list.push_back(boost::lexical_cast<int>(elements[2]));
    }
  }
  catch (const std::ifstream::failure& e) {
    std::cout << "Error in opening or reading the " << hyphansis_score << " file!" << std::endl;
  }
  // Close the score file
  mscore.close();

  // Open the hyphantic log file
  std::ofstream s(hyphansis_file,std::ios::app);

  if (key_list.empty()) {
    // We're done!
    s << "  </Sheet>" << std::endl;
    s.close();
    return; 
  }

  // Start "playing" the notes for this voice - our instrument is the topology of spacetime...
  opcount = (signed) key_list.size();

  // This is a fairly complicated operation - we need to play the notes the in the right order, 
  // while at the same time highest pitched key above 44 is assigned to the vertex with the most 
  // negative deficiency, the next highest pitched key above 44 acts upon the vertex with the  
  // second most negative deficiency etc.
  // We need to create a list of the distinct explicative and implicative piano keys in this measure, 
  // paired to the appropriate vertex
  for(i=0; i<opcount; ++i) {
#ifdef DEBUG
    assert(key_list[i] >= 1 && key_list[i] <= 88);
#endif
    if (key_list[i] > 40) {
      if (std::count(m_keys.begin(),m_keys.end(),key_list[i]) == 0) m_keys.push_back(key_list[i]);
    }
    else if (key_list[i] < 40) {
      if (std::count(x_keys.begin(),x_keys.end(),key_list[i]) == 0) x_keys.push_back(key_list[i]);
    }
  }
  // Now sort these piano key values in the correct order, meaning ascending 
  // for explicative (1 to 39) and descending for implicative (88 to 41)
  std::sort(m_keys.begin(),m_keys.end(),std::greater<int>());
  std::sort(x_keys.begin(),x_keys.end());

  for(i=nc-1; i>0; --i) {
    v = candidates[i].first;
    if (!events[v].active(sheet)) continue;
    neutral_vertices.push_back(v);
    if (events[v].deficiency < -std::numeric_limits<double>::epsilon()) {
      m_vertices.push_back(v);
    }
    else if (events[v].deficiency > std::numeric_limits<double>::epsilon()) {
      x_vertices.push_back(v);
    }
  }

  if (!m_keys.empty()) m_width = double(m_keys[0] - m_keys.back());
  if (!x_keys.empty()) x_width = double(x_keys.back() - x_keys[0]);
  for(i=0; i<opcount; ++i) {
    j = key_list[i];
    if (j > 40) {
      v = select_vertex(m_vertices,double(j - m_keys.back())/m_width,sheet);
    }
    else if (j < 40) {
      v = select_vertex(x_vertices,double(j - x_keys[0])/x_width,sheet);
    }
    else {
      // Any active vertex will do for a neutral operator
      v = RND->irandom(neutral_vertices);
    }
    if (v == -1) continue;
    // Now we have the base vertex v, next we need to get the operator and 
    // parameters for this piano key
    if (j > 40) {
      op = implicative_scale(j,pvalues);
    }
    else if (j < 40) {
      op = explicative_scale(j,pvalues);
    }
    else {
      // Edge reorientation
      op = "T";
    }
    opstring << op << "," << v;
    if (op == "F") {
      success = fission(v,pvalues[0],sheet);
      opstring << "," << pvalues[0]; 
    }
    else if (op == "Um") {
      success = fusion_m(v,sheet);
    }
    else if (op == "Om") {
      success = foliation_m(v,sheet);
    }
    else if (op == "E") {
      success = expansion(v,pvalues[0],sheet);
      opstring << "," << pvalues[0];
    }
    else if (op == "I") {
      success = inflation(v,pvalues[0],sheet);
      opstring << "," << pvalues[0]; 
    }
    else if (op == "P") {
      success = perforation(v,0,sheet);
      opstring << ",0"; 
    }
    else if (op == "V") {
      success = circumvolution(v,sheet);
    }
    else if (op == "D") {
      success = deflation(v,sheet);
    }
    else if (op == "Ux") {
      success = fusion_x(v,pvalues[0],sheet);
      opstring << "," << pvalues[0]; 
    }
    else if (op == "Ox") {
      success = foliation_x(v,sheet);
    }
    else if (op == "Sg") {
      success = compensation_g(v,sheet);
    }
    else if (op == "Sm") {
      success = compensation_m(v,sheet);
    }
    else if (op == "R") {
      success = reduction(v,sheet);
    }
    else if (op == "C") {
      success = correction(v,sheet);
    }
    else if (op == "Δ") {
      success = stellar_deletion(v,sheet);
    }
    else if (op == "Y") {
      success = stellar_addition(v,sheet);
    }
    else if (op == "N") {
      success = contraction(v,pvalues[0],sheet);
      opstring << "," << pvalues[0]; 
    }
    else if (op == "A") {
      success = amputation(v,10.0,sheet);
      opstring << ",10.0"; 
    }
    else if (op == "G") {
      success = germination(v,sheet);
    }
    else if (op == "T") {
      success = edge_parity_mutation(v,sheet);
    }
    if (success) {
      s << "    <Operation>" << opstring.str() << "</Operation>" << std::endl;
      regularization(false,sheet);
      codex[sheet].ops += op;
    }
    opstring.str("");
  }

  // We're done, so close the hyphantic log file and return
  s << "  </Sheet>" << std::endl;
  s.close(); 
}

void Spacetime::dynamic_hyphansis(const std::vector<std::pair<int,double> >& candidates,int sheet)
{
  int i,v,n,vx[2],nsuccess = 0;
  double alpha;
  std::string op;
  std::stringstream opstring;
  std::set<int> r_edges;
  bool success = false;
  const int nc = (signed) candidates.size();

  std::ofstream s(hyphansis_file,std::ios::app);

  for(i=nc-1; i>0; --i) {
    v = candidates[i].first;
    // An earlier hyphantic operation may have rendered this
    // vertex inactive...
    if (!events[v].active(sheet)) continue;
    alpha = events[v].deficiency;
    if (alpha < -std::numeric_limits<double>::epsilon()) {
      implication(op);
      opstring << op << "," << v;
      if (op == "F") {
        success = fission(v,0.4,sheet);
        opstring << ",0.4";
      }
      else if (op == "Um") {
        success = fusion_m(v,sheet);
      }
      else if (op == "Om") {
        success = foliation_m(v,sheet);
      }
      else if (op == "E") {
        success = expansion(v,0.15,sheet);
        opstring << ",0.15";
      }
      else if (op == "I") {
        success = inflation(v,0.25,sheet);
        opstring << ",0.25";
      }
      else if (op == "P") {
        success = perforation(v,0,sheet);
        opstring << ",0";
      }
      else if (op == "V") {
        success = circumvolution(v,sheet);
      }
      else if (op == "Δ") {
        success = stellar_deletion(v,sheet);
      }
    }
    else if (alpha > std::numeric_limits<double>::epsilon()) {
      if (vertex_dimension(v,sheet) > 1) {
        if (RND->drandom() < alpha/10.0) {
          op = "D";
          success = deflation(v,sheet);
          opstring << "D," << v;
        }
        else {
          op = "R";
          success = reduction(v,sheet);
          opstring << "R," << v;
        }
      }
      else {
        explication(op);
        opstring << op << "," << v;
        if (op == "C") {
          success = correction(v,sheet);
        }
        else if (op == "N") {
          success = contraction(v,1.2,sheet);
          opstring << ",1.2";
        }
        else if (op == "Ux") {
          success = fusion_x(v,0.5,sheet);
          opstring << "0.5";
        }
        else if (op == "Sg") {
          success = compensation_g(v,sheet);
        }
        else if (op == "Sm") {
          success = compensation_m(v,sheet);
        }
        else if (op == "G") {
          success = germination(v,sheet);
        }
        else if (op == "Y") {
          success = stellar_addition(v,sheet);
        }
        else if (op == "A") {
          success = amputation(v,10.0,sheet);
          opstring << ",10.0";
        }
        if (!success && op == "R") {
          opstring.str("");
          if (RND->irandom(2) == 0) {
            op = "N";
            success = contraction(v,1.2,sheet);
            opstring << "N," << v << ",1.2"; 
          }
          else {
            op = "Ux";
            success = fusion_x(v,0.5,sheet);
            opstring << "Ux," << v << ",0.5";
          }
        }
        if (!success) {
          op = "Sg";
          success = compensation_g(v,sheet);
          opstring.str("");
          opstring << "Sg," << v; 
        }
      }
    }
    if (success) {
      s << "    <Operation>" << opstring.str() << "</Operation>" << std::endl;
      regularization(false,sheet);
      codex[sheet].ops += op;
      nsuccess++;
    }
    if (double(nsuccess)/nactive > 0.1) break;
  }
  n = (signed) simplices[1].size();
  for(i=0; i<n; ++i) {
    if (!simplices[1][i].active(sheet)) continue;
    if (RND->drandom() < parity_mutation) {
      simplices[1][i].get_vertices(vx);
      if (edge_parity_mutation(vx[0],vx[1],sheet)) {
        s << "    <Operation>T," << vx[0] << "</Operation>" << std::endl;
        r_edges.insert(i);
      }
    }
  }
  recompute_parity(r_edges);
  s << "  </Sheet>" << std::endl;
  s.close();
}

void Spacetime::hyphansis(int sheet)
{
  int i,v;
  double alpha;
  std::vector<std::pair<int,double> > candidates;
  const int nv = (signed) events.size();

  codex[sheet].ops = "";

  std::ofstream s(hyphansis_file,std::ios::app);
  s << "  <Sheet>" << std::endl;
  s << "    <Index>" << sheet << "</Index>" << std::endl;

#ifdef VERBOSE
  int npos = 0,nneg = 0,nze = 0;
#endif
  for(i=0; i<nv; ++i) {
    if (!events[i].active(sheet)) continue;
#ifdef VERBOSE
    if (!events[i].zero_energy()) nze++;
    if (events[i].deficiency > std::numeric_limits<double>::epsilon()) npos++;
    if (events[i].deficiency < -std::numeric_limits<double>::epsilon()) nneg++;
#endif
    // Eliminate vertices whose structural deficiency is close to zero...
    alpha = std::abs(events[i].deficiency);
    if (alpha < std::numeric_limits<double>::epsilon()) continue;
    candidates.push_back(std::pair<int,double>(i,alpha));
  }
#ifdef VERBOSE
  std::cout << "There are " << nze << " vertices with positive energy or " << 100.0*double(nze)/double(cardinality(0,sheet)) << " percent of the total." << std::endl;
  std::cout << "There are " << npos << " positive vertices and " << nneg << " negative vertices in the spacetime complex." << std::endl;
#endif
  if (candidates.empty()) {
    s << "  </Sheet>" << std::endl;
    s.close();
    return;
  }

  if (candidates.size() == 1) {
    v = candidates[0].first;
    if (events[v].deficiency < -std::numeric_limits<double>::epsilon()) {
      assert(expansion(v,sheet));
      codex[sheet].ops += 'E';
      regularization(false,sheet);
      s << "    <Operation>E," << v << "</Operation>" << std::endl;
      s << "  </Sheet>" << std::endl;
      s.close();
      return;
    }
  }
  s.close();

  std::sort(candidates.begin(),candidates.end(),SYNARMOSMA::pair_predicate_dbl);

  if (weaving == Hyphansis::dynamic) {
    dynamic_hyphansis(candidates,sheet);
  }
  else {
    musical_hyphansis(candidates,sheet);
  }
}

void Spacetime::vertex_fusion(int n1,int n2,int sheet)
{
  int i,j,l;
  std::set<int> vx,locus;
  const int ulimit = dimension(sheet);

  if (sheet == -1) {
    // A method that converts all occurences of n2 into n1...
    int k,m,n,im1,in1;
    std::set<int> duplicate;
    std::set<int>::const_iterator jt;
    std::set<int>::reverse_iterator it;
    std::vector<int> mindex;
    std::vector<int>* mutation = new std::vector<int>[ulimit+1];

    // Do the actual vertex swap...
    for(i=1; i<=ulimit; ++i) {
      m = (signed) simplices[i].size();
      for(j=0; j<m; ++j) {
        if (simplices[i][j].exchange(n1,n2)) mutation[i].push_back(j);
      }
    }

    // Place simplices in the right dimension...
    for(i=ulimit; i>1; i--) {
      im1 = i - 1;
      n = (signed) simplices[i-1].size();
      for(j=0; j<(signed) mutation[i].size(); ++j) {
        in1 = mutation[i][j];
        if (simplices[i][in1].dimension() == im1) duplicate.insert(in1);
      }
      for(it=duplicate.rbegin(); it!=duplicate.rend(); ++it) {
        in1 = *it;
        simplices[im1].push_back(simplices[i][in1]);
        mutation[im1].push_back(n);
        n++;
        simplices[i].erase(simplices[i].begin() + in1);
        for(j=0; j<(signed) mutation[i].size(); ++j) {
          k = mutation[i][j];
          if (k == in1) continue;
          if (k > in1) {
            mindex.push_back(k-1);
          }
          else {
            mindex.push_back(k);
          }
        }
        mutation[i] = mindex;
        mindex.clear();
      }
      duplicate.clear();
    }
    m = (signed) mutation[1].size();
    for(i=0; i<m; ++i) {
      in1 = mutation[1][i];
      if (simplices[1][in1].dimension() == 0) duplicate.insert(in1);
    }
    for(it=duplicate.rbegin(); it!=duplicate.rend(); ++it) {
      in1 = *it;
      simplices[1].erase(simplices[1].begin() + in1);
      for(j=0; j<m; ++j) {
        k = mutation[1][j];
        if (k == in1) continue;
        if (k > in1) {
          mindex.push_back(k-1);
        }
        else {
          mindex.push_back(k);
        }
      }
      mutation[1] = mindex;
    }
    duplicate.clear();

    // Now, check for duplicates at each dimension and delete them...
    for(i=1; i<=ulimit; ++i) {
      n = (signed) simplices[i].size();
      for(j=0; j<(signed) mutation[i].size(); ++j) {
        in1 = mutation[i][j];
#ifdef DEBUG
        assert(in1 < n);
#endif
        vx = simplices[i][in1].vertices;
        for(l=0; l<n; ++l) {
          if (l == in1) continue;
          if (simplices[i][l].vertices == vx) {
            duplicate.insert(in1);
            simplices[i][in1].get_ubiquity(locus);
            for(jt=locus.begin(); jt!=locus.end(); ++jt) {
              simplices[i][l].set_active(*jt);
            }
          }
        }
      }
      for(it=duplicate.rbegin(); it!=duplicate.rend(); ++it) {
        in1 = *it;
        simplices[i].erase(simplices[i].begin() + in1);
      }
      duplicate.clear();
    }
    events[n2].get_ubiquity(locus);
    for(jt=locus.begin(); jt!=locus.end(); ++jt) {
      events[n1].set_active(*jt);
    }
    events[n2].ubiquity.clear();
    delete[] mutation;

    // Then recalculate the index_table hash map..
    for(i=1; i<=ulimit; ++i) {
      index_table[i].clear();
      m = (signed) simplices[i].size();
      for(j=0; j<m; ++j) {
        simplices[i][j].entourage.clear();
        index_table[i][simplices[i][j].vertices] = j;
      }
    }
    for(i=0; i<(signed) events.size(); ++i) {
      events[i].entourage.clear();
    }

    // Then recalculate the entourages...
    compute_entourages(-1);
    compute_neighbours();
  }
  else {
    bool found;
    std::set<int> v;
    std::set<int>::const_iterator it;
    SYNARMOSMA::hash_map::const_iterator qt;

    locus.insert(sheet);

    codex[sheet].vx_delta.insert(n1);
    codex[sheet].vx_delta.insert(n2);
    for(i=1; i<=Spacetime::ND; ++i) {
      for(j=0; j<(signed) simplices[i].size(); ++j) {
        if (!simplices[i][j].active(sheet)) continue;
        // Check of this simplex contains the vertex "n2", if so
        // divide its colour by sheet and create a new simplex with the
        // right vertices. Check to see if it already exists, if so
        // multiply its colour by sheet otherwise add the new simplex.
        simplices[i][j].get_vertices(v);
        found = false;
        for(it=v.begin(); it!=v.end(); ++it) {
          if (*it == n2) {
            found = true;
            vx.insert(n1);
            continue;
          }
          vx.insert(*it);
        }
        if (!found) {
          vx.clear();
          continue;
        }
        simplices[i][j].set_inactive(sheet);
        for(it=simplices[i][j].vertices.begin(); it!=simplices[i][j].vertices.end(); ++it) {
          codex[sheet].vx_delta.insert(*it);
        }
        l = (signed) vx.size() - 1;
        if (l >= 1) {
          qt = index_table[l].find(vx);
          if (qt != index_table[l].end()) {
            simplices[l][qt->second].set_active(sheet); 
          }
          else {
            simplices[l].push_back(Simplex(vx,locus));
            index_table[l][vx] = (signed) simplices[l].size() - 1;
          }
        }
        else {
          l = SYNARMOSMA::element(vx);
          events[l].set_active(sheet); 
        }
        vx.clear();
      }
    }
    events[n2].set_inactive(sheet); 
  }
  if (!events[n2].active()) {
    if (!events[n2].zero_energy()) {
      events[n1].increment_energy(events[n2].get_energy());
      events[n2].nullify_energy();
    }
  }
}

bool Spacetime::vertex_twist(int sheet)
{
  // This method will fuse two 0-simplices with each other, so as to twist the
  // complex's topology, creating non-orientability and torsion groups in the
  // simplicial homology.
  int i,j,n1,n2,nc,delta,tolerance,n = (signed) events.size();
  std::vector<std::string> ekeys;
  SYNARMOSMA::hash_map::const_iterator qt;
  std::set<int> parents;
  std::set<int>::const_iterator it;
  std::vector<int> candidates,used;
  std::vector<std::pair<int,int> > slist;
  std::stringstream s;
  static int first = 1;

#ifdef DEBUG
  assert(consistent(sheet));
#endif

  tolerance = (first == 1) ? 2 : 1;

  for(i=0; i<n; ++i) {
    if (events[i].active(sheet)) candidates.push_back(i);
  }
  nc = (signed) candidates.size();
  if (nc < 3) return false;
#ifdef VERBOSE
  std::cout << "Vertex twist with " << nc << " candidates." << std::endl;
#endif
  for(i=0; i<nc; ++i) {
    used.push_back(0);
  }
  // Grab a random pair (v1,v2) of elements from the vector candidates, set v1
  // to be inactive and then go through all d-simplices, d > 0, setting every
  // occurrence of v1 to v2.
  for(i=0; i<nc; ++i) {
    for(j=1+i; j<nc; ++j) {
      delta = combinatorial_distance(candidates[i],candidates[j],-1);
      if (delta > tolerance) {
        slist.push_back(std::pair<int,int>(candidates[i],candidates[j]));
      }
    }
  } 
  if (slist.empty()) return false;
  n1 = RND->irandom(slist.size());
  n1 = slist[i].first;
  n2 = slist[i].second;
#ifdef VERBOSE
  std::cout << "Twist is " << n1 << " => " << n2 << std::endl;
#endif
  events[n2].set_inactive(sheet);
  codex[sheet].vx_delta.insert(n2);
  if (!events[n2].active()) {
    if (!events[n2].zero_energy()) {
      events[n1].increment_energy(events[n2].get_energy());
      events[n2].nullify_energy();
    }
  }
  vertex_fusion(n1,n2,sheet);
  if (first == 1) first = 0;
  return true;
}

bool Spacetime::vertex_deletion(int n,int sheet)
{
  if (!events[n].active(sheet)) return false;
  int i,vx[2],ne = (signed) simplices[1].size();
  events[n].set_inactive(sheet);
  for(i=0; i<ne; ++i) {
    if (!simplices[1][i].active(sheet)) continue;
    simplices[1][i].get_vertices(vx);
    if (vx[0] == n || vx[1] == n) simplex_deletion(1,i,sheet);
  }
  return true;  
}

int Spacetime::vertex_addition(const std::vector<double>& xc,int sheet)
{
  int n = (signed) events.size();
  Event vt;

  geometry->vertex_addition(xc);
  vt.set_active(sheet); 
  events.push_back(vt);

  return n;
}

int Spacetime::vertex_addition(const std::set<int>& antecedents,int sheet)
{
  int n = (signed) events.size();
  Event vt;

  geometry->vertex_addition(antecedents);
  vt.set_active(sheet);
  events.push_back(vt);

  return n;
}

int Spacetime::vertex_addition(int base,int sheet)
{
  int n = (signed) events.size();
  Event vt;

  geometry->vertex_addition(base);
  vt.set_active(sheet);
  events.push_back(vt);

  return n;
}