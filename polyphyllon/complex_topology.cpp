#include "complex.h"

using namespace DIAPLEXIS;

void Complex::compute_simplicial_dimension()
{
  unsigned int i;
  const unsigned int nv = events.size();

  for(i=0; i<nv; ++i) {
    if (!events[i].active()) continue;
    events[i].topological_dimension = vertex_dimension(i,-1);
  }
}

void Complex::get_edge_topology(std::vector<std::set<int> >& vx) const
{
  unsigned int i;
  const unsigned int nv = events.size();

  vx.clear();
  for(i=0; i<nv; ++i) {
    if (!events[i].active()) continue;
    vx.push_back(events[i].neighbours);
  }
}

void Complex::compute_degree_distribution(bool logarithmic,int sheet) const
{
  std::vector<double> histogram;
  SYNARMOSMA::Graph* G = new SYNARMOSMA::Graph;

  compute_graph(G,sheet);

  G->degree_distribution(logarithmic,histogram);

  std::cout << "Vertex degree histogram (sheet = " << sheet << "):" << std::endl;
  for(int i=0; i<(signed) histogram.size(); ++i) {
    if (histogram[i] > 0) std::cout << i+1 << "  " << 100.0*histogram[i] << std::endl;
  }
  delete G;
}

void Complex::compute_connectivity_distribution(bool direct,int sheet) const
{
  int i,j,m = 0,l = cardinality(0,sheet);
  std::vector<int> pcount;
  const int nv = (signed) events.size();
  const double ntotal = double(l*(l-1)/2);

  for(i=0; i<nv; ++i) {
    pcount.push_back(0);
  }

  // If we know memory is abundant relative to the spacetime 
  // size, we can compute the distances all at once
  if (direct) {
    int offset[nv];
    SYNARMOSMA::Graph G;
    SYNARMOSMA::pair_index distances;
    SYNARMOSMA::pair_index::const_iterator qt;

    compute_graph(&G,offset,sheet);
    G.compute_distances(distances);
    if (sheet == -1) {
      for(i=0; i<nv; ++i) {
        if (!events[i].active()) continue;
        for(j=1+i; j<nv; ++j) {
          if (!events[j].active()) continue;
          qt = distances.find(std::pair<int,int>(offset[i],offset[j]));
          l = qt->second;
          if (l > m) m = l;
          pcount[l] += 1;
        }
      }
    }
    else {
      for(i=0; i<nv; ++i) {
        if (!events[i].active(sheet)) continue;
        for(j=1+i; j<nv; ++j) {
          if (!events[j].active(sheet)) continue;
          qt = distances.find(std::pair<int,int>(offset[i],offset[j]));
          l = qt->second;
          if (l > m) m = l;
          pcount[l] += 1;
        }
      }
    }
  }
  else {
    if (sheet == -1) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i,j,l) reduction(max:m) schedule(dynamic,1)
#endif
      for(i=0; i<nv; ++i) {
        if (!events[i].active()) continue;
        for(j=1+i; j<nv; ++j) {
          if (!events[j].active()) continue;
          l = combinatorial_distance(i,j,-1);
          if (l > m) m = l;
#ifdef _OPENMP
#pragma omp critical
          {
#endif
          pcount[l] += 1;
#ifdef _OPENMP
          }
#endif
        }
      }
    }
    else {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i,j,l) reduction(max:m) schedule(dynamic,1)
#endif
      for(i=0; i<nv; ++i) {
        if (!events[i].active(sheet)) continue;
        for(j=1+i; j<nv; ++j) {
          if (!events[j].active(sheet)) continue;
          l = combinatorial_distance(i,j,sheet);
          if (l > m) m = l;
#ifdef _OPENMP
#pragma omp critical
          {
#endif
          pcount[l] += 1;
#ifdef _OPENMP
          }
#endif
        }
      }
    }
  }
  std::cout << "Combinatorial distance histogram (sheet = " << sheet << "):" << std::endl;  
  for(i=1; i<=m; ++i) {
    std::cout << i << "  " << 100.0*double(pcount[i])/ntotal << std::endl;
  }
}

std::pair<double,double> Complex::random_walk(int sheet) const
{
  std::pair<double,double> output;
  SYNARMOSMA::Graph G;

  compute_graph(&G,sheet);
  int L = G.size()/4;
  // Do a random walk with a length equal to a quarter of the graph's size and 
  // a fifth of its vertices as starting point.  
  output = G.random_walk(L,0.2);

  return output;
}

void Complex::compute_dependent_simplices(const std::set<int>& modified_vertices)
{
  // A method that takes a set of vertices whose coordinates have
  // been modified and determines which simplices will have their
  // modified property set to true...
  if (modified_vertices.empty()) return;
  int i,j,m,n;
  std::set<int> current,next;
  std::set<int>::const_iterator it,jt,kt;

#ifdef VERBOSE
  std::cout << "There are " << modified_vertices.size() << " vertices directly implicated." << std::endl;
#endif

  for(i=1; i<=Complex::ND; ++i) {
    n = (signed) simplices[i].size();
    for(j=0; j<n; ++j) {
      simplices[i][j].modified = false;
    }
  }
  // We assume that the number of modified vertices is small relative
  // to the total number of vertices in the complex...
  for(it=modified_vertices.begin(); it!=modified_vertices.end(); ++it) {
    n = *it;
    current = events[n].entourage;
    for(i=1; i<=Complex::ND; ++i) {
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

double Complex::dimensional_stress(int v,int sheet) const
{
  // This method should sum the dimensional stress associated with
  // each d-simplex (d > 0) that contains the vertex v and exists on
  // the specified sheet.
  if (sheet == -1) {
    if (!events[v].active()) return 0.0;
  }
  else {
    if (!events[v].active(sheet)) return 0.0;
  }

  int i,j,ds;
  const int nd = dimension(-1);
  double sum = 0.0;

  for(i=nd; i>0; i--) {
    ds = (signed) simplices[i].size();
    for(j=0; j<ds; ++j) {
      if (!simplices[i][j].active()) continue;
      if (simplices[i][j].contains(v)) sum += dimensional_stress(i,j,sheet);
    }
  }
  return sum;
}

double Complex::dimensional_stress(int d,int n,int sheet) const
{
  // This method measures the standard deviation of the
  // simplicial dimensions of the vertices of a given
  // d-simplex
  if (sheet == -1) {
    if (!simplices[d][n].active()) return 0.0;
  }
  else {
    if (!simplices[d][n].active(sheet)) return 0.0;
  }

  int i,vx[1+d];
  double alpha,v[1+d],sigma = 0.0,mu = 0.0;

  simplices[d][n].get_vertices(vx);

  for(i=0; i<1+d; ++i) {
    alpha = double(vertex_dimension(vx[i],sheet));
    v[i] = alpha;
    mu += alpha;
  }
  mu = mu/double(1+d);
  for(i=0; i<1+d; ++i) {
    sigma += (v[i] - mu)*(v[i] - mu);
  }
  sigma = std::sqrt(sigma/double(1+d));
  return sigma;
}

void Complex::recompute_parity(int n) 
{
  // The edge n has had its parity changed, so now we need to recompute 
  // the parity of the higher-dimensional simplices that depend on it.
  std::set<int> temp;
  
  temp.insert(n);
  recompute_parity(temp);
}

void Complex::recompute_parity(const std::set<int>& edges) 
{
  // The 1-simplices in the argument have had their parity changed, so now 
  // we need to recompute the parity of the higher-dimensional simplices that 
  // depend on them.
  if (edges.empty()) return;
  int i,d = 1;
  std::set<int> current,next;
  std::set<int>::const_iterator it,jt;
  
  current = edges;

  do {
    for(it=current.begin(); it!=current.end(); ++it) {
      i = *it;
      compute_simplex_parity(d,i);
      for(jt=simplices[d][i].entourage.begin(); jt!=simplices[d][i].entourage.end(); ++jt) {
        next.insert(*jt);
      }
    }
    if (next.empty()) break;
    current = next;
    next.clear();
    d += 1;
  } while(true);
}

void Complex::compute_parity() 
{
  int i,j;
  const int d = dimension(-1);

  for(i=2; i<=d; ++i) {
    for(j=0; j<(signed) simplices[i].size(); ++j) {
      if (!simplices[i][j].modified) continue;
      compute_simplex_parity(i,j);
    }
  }
}

void Complex::compute_fvector(std::vector<int>& f,std::vector<int>& fstar,int sheet) const
{
  int i,j,s = 0;
  const int D = dimension(sheet);

  f.clear();
  fstar.clear();

  // First the f-vector
  for(i=0; i<=D; ++i) {
    f.push_back(cardinality(i,sheet));
  }

  // Now the dual vector, f*
  for(i=0; i<=D; ++i) {
    s += SYNARMOSMA::ipow(-1,i)*f[i];
  }
  fstar.push_back(s);
  for(i=1; i<D; ++i) {
    s = 0;
    for(j=i; j<=D; ++j) {
      s += SYNARMOSMA::ipow(-1,j-i)*int(SYNARMOSMA::binomial(j,i))*f[j]; 
    }
    fstar.push_back(s);
  }
  fstar.push_back(f[D]);
}

void Complex::compute_hvector(std::vector<int>& h,int sheet) const
{
  int i,j,sum;
  std::vector<int> f,fstar;
  const int D = dimension(sheet);

  compute_fvector(f,fstar,sheet);
  for(i=0; i<=D; ++i) {
    sum = SYNARMOSMA::ipow(-1,i)*int(SYNARMOSMA::binomial(D+1,i));
    for(j=1; j<=i; ++j) {
      sum += SYNARMOSMA::ipow(-1,i-j)*int(SYNARMOSMA::binomial(D+1-j,i-j))*f[j-1];
    }
    h.push_back(sum);
  }
}

void Complex::vertex_degree_statistics(double* output,int sheet) const
{
  SYNARMOSMA::Graph G;  
  compute_graph(&G,sheet);
  output[0] = double(G.max_degree());
  output[1] = double(G.min_degree());
  output[2] = G.average_degree();
}

int Complex::cyclicity(int sheet) const
{
  // This method should calculate the number of cyclic edges in the complex
  // associated with the sheet, assuming the complex is connected.
  SYNARMOSMA::Graph G;

  compute_graph(&G,sheet);
  if (!G.connected()) return 0;
  return (G.size() - G.bridge_count());
}

void Complex::compute_simplex_parity(int d,int n)
{
  if (d < 2) return;
  int i,j,vx[1+d];
  std::set<int> S;
  SYNARMOSMA::hash_map::const_iterator qt;

  simplices[d][n].get_vertices(vx);

  simplices[d][n].parity = 1;
  // The parity of the simplex is simply the product of the parity of its edges...
  for(i=0; i<=d; ++i) {
    for(j=1+i; j<=d; ++j) {
      S.insert(vx[i]); S.insert(vx[j]);
      qt = index_table[1].find(S);
      simplices[d][n].parity *= simplices[1][qt->second].parity;
      S.clear();     
    }
  }
}

void Complex::simplex_membership(int v,std::vector<int>& output) const
{
  int i,k = 0;
  unsigned int j,m;

  output.clear();
  for(i=1; i<=Complex::ND; ++i) {
    m = simplices[i].size();
    for(j=0; j<m; ++j) {
      if (!simplices[i][j].active()) continue;
      if (simplices[i][j].contains(v)) ++k;
    }
    output.push_back(k);
    k = 0;
  }
}

void Complex::compute_graph(SYNARMOSMA::Graph* G,int* offset,int sheet) const
{
  int i,vx[2];
  const int nv = (signed) events.size();
  const int ne = (signed) simplices[1].size();

  G->clear();

  if (sheet == -1) {
    for(i=0; i<nv; ++i) {
      offset[i] = -1;
      if (!events[i].active()) continue;
      offset[i] = G->add_vertex();
    }
    for(i=0; i<ne; ++i) {
      if (!simplices[1][i].active()) continue;
      simplices[1][i].get_vertices(vx);
#ifdef DEBUG
      assert(offset[vx[0]] >= 0 && offset[vx[1]] >= 0);
#endif
      G->add_edge(offset[vx[0]],offset[vx[1]]);
    }
  }
  else {
    for(i=0; i<nv; ++i) {
      offset[i] = -1;
      if (!events[i].active(sheet)) continue;
      offset[i] = G->add_vertex();
    }
    for(i=0; i<ne; ++i) {
      if (!simplices[1][i].active(sheet)) continue;
      simplices[1][i].get_vertices(vx);
#ifdef DEBUG
      assert(offset[vx[0]] >= 0 && offset[vx[1]] >= 0);
#endif
      G->add_edge(offset[vx[0]],offset[vx[1]]);
    }
  }
}

void Complex::compute_graph(SYNARMOSMA::Graph* G,int base,int steps,int sheet) const
{
  int i,v,w,hop = 1;
  const int nv = (signed) events.size();
  SYNARMOSMA::hash_map::const_iterator qt;
  std::set<int>::const_iterator it,jt;
  std::set<int> current,next,S;
  std::vector<int> offset;

  G->clear();
  for(i=0; i<nv; ++i) {
    offset.push_back(-1);
  }

  current.insert(base);
  offset[base] = G->add_vertex();

  do {
    for(it=current.begin(); it!=current.end(); ++it) {
      v = *it;
      for(jt=events[v].neighbours.begin(); jt!=events[v].neighbours.end(); ++jt) {
        w = *jt;
        S.clear();
        S.insert(v);
        S.insert(w);
        qt = index_table[1].find(S);
        if (!simplices[1][qt->second].active(sheet)) continue;
        if (offset[w] == -1) {
          offset[w] = G->add_vertex();
          next.insert(w);
        }
        G->add_edge(offset[v],offset[w]);
      }
    }
    if (next.empty() || hop == steps) break;
    hop++;
    current = next;
    next.clear();
  } while(true);
}

void Complex::compute_global_nexus(SYNARMOSMA::Nexus* NX,int sheet) const
{
  int i,j;
  std::vector<int> offset;
  std::set<int> vx;
  std::set<int>::const_iterator it;
  const int nv = (signed) events.size();
  const int n = dimension(sheet);

  NX->initialize(n);
  for(i=0; i<nv; ++i) {
    offset.push_back(-1);
  }
  if (sheet == -1) {
    for(i=0; i<nv; ++i) {
      if (!events[i].active()) continue;
      offset[i] = NX->add_vertex();
    }
    for(i=n; i>=1; --i) {
      for(j=0; j<(signed) simplices[i].size(); ++j) {
        if (!simplices[i][j].active()) continue;
        for(it=simplices[i][j].vertices.begin(); it!=simplices[i][j].vertices.end(); ++it) {
          vx.insert(offset[*it]);
        }
        NX->paste(vx);
        vx.clear();
      }
    }
  }
  else {
    for(i=0; i<nv; ++i) {
      if (!events[i].active(sheet)) continue;
      offset[i] = NX->add_vertex();
    }
    for(i=n; i>=1; --i) {
      for(j=0; j<(signed) simplices[i].size(); ++j) {
        if (!simplices[i][j].active(sheet)) continue;
        for(it=simplices[i][j].vertices.begin(); it!=simplices[i][j].vertices.end(); ++it) {
          vx.insert(offset[*it]);
        }
        NX->paste(vx);
        vx.clear();
      }
    }
  }
  NX->regularization();
}

void Complex::compute_local_nexus(SYNARMOSMA::Nexus* NX,int base,int sheet) const
{
  int i,j;
  std::vector<int> offset;
  std::set<int>::const_iterator it;
  std::set<int> vx;
  const int n = vertex_dimension(base,sheet);
  const int nv = (signed) events.size();

  for(i=0; i<nv; ++i) {
    offset.push_back(-1);
  }
  NX->initialize(n);
  offset[base] = NX->add_vertex();
  if (sheet == -1) {
    for(i=1; i<=n; ++i) {
      for(j=0; j<(signed) simplices[i].size(); ++j) {
        if (!simplices[i][j].active()) continue;
        if (simplices[i][j].contains(base)) {
          for(it=simplices[i][j].vertices.begin(); it!=simplices[i][j].vertices.end(); ++it) {
            if (offset[*it] == -1) offset[*it] = NX->add_vertex();
            vx.insert(offset[*it]);
          }
          NX->paste(vx);
          vx.clear();
        }
      }
    }
  }
  else {
    for(i=1; i<=n; ++i) {
      for(j=0; j<(signed) simplices[i].size(); ++j) {
        if (!simplices[i][j].active(sheet)) continue;
        if (simplices[i][j].contains(base)) {
          for(it=simplices[i][j].vertices.begin(); it!=simplices[i][j].vertices.end(); ++it) {
            if (offset[*it] == -1) offset[*it] = NX->add_vertex();
            vx.insert(offset[*it]);
          }
          NX->paste(vx);
          vx.clear();
        }
      }
    }
  }
  NX->regularization();
}

int Complex::chromatic_number(int sheet) const
{
  // Computes the chromatic number chi of the graph associated with the
  // colour "p"; we already know that 1 <= chi <= max_degree+1.  
  SYNARMOSMA::Graph G;

  compute_graph(&G,sheet);
  // It makes no sense to try to compute a single chromatic number for
  // a disconnected graph.
  if (!G.connected()) return 0;
  return G.chromatic_number();
}

bool Complex::edge_parity_mutation(int base,int sheet)
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

bool Complex::edge_parity_mutation(int u,int v,int sheet)
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

void Complex::simplex_deletion(int d,int n,int sheet)
{
  std::set<int>::const_iterator it;
  std::set<int> parents;
  int i,dp1 = d + 1;

  if (sheet == -1) {
    if (simplices[d][n].active()) return;
    simplices[d][n].deactivate();
    simplices[d][n].modified = true;
  } 
  else {
    if (!simplices[d][n].active(sheet)) return;
    simplices[d][n].deactivate(sheet);
    simplices[d][n].modified = true;
    for(it=simplices[d][n].vertices.begin(); it!=simplices[d][n].vertices.end(); ++it) {
      events[*it].topology_modified = true;
    }
  }
  parents = simplices[d][n].entourage;
  for(it=parents.begin(); it!=parents.end(); ++it) {
    i = *it;
    simplex_deletion(dp1,i,sheet);
  }
}

bool Complex::simplex_addition(int u,int v,const std::set<int>& locus,int n)
{
  bool output = false;
  std::set<int> S;
  std::set<int>::const_iterator it;
  SYNARMOSMA::hash_map::const_iterator qt;

  S.insert(u); S.insert(v);
  qt = index_table[1].find(S);
  if (qt == index_table[1].end()) {
    output = true;
    Simplex s(S,locus);
    s.incept = n;
    simplices[1].push_back(s);
    index_table[1][S] = simplices[1].size() - 1;
    events[v].neighbours.insert(u);
    events[u].neighbours.insert(v);
  }
  else {
    for(it=locus.begin(); it!=locus.end(); ++it) {
      if (!simplices[1][qt->second].active(*it)) {
        simplices[1][qt->second].activate(*it);
        output = true;
      }
    }
  }
  for(it=locus.begin(); it!=locus.end(); ++it) {
    if (!events[u].active(*it)) {
      events[u].activate(*it);
      output = true;
    }
    if (!events[v].active(*it)) {
      events[v].activate(*it);
      output = true;
    }
  }
  if (output) {
    events[u].topology_modified = true;
    events[v].topology_modified = true;
  }
  return output;
}

bool Complex::simplex_addition(const std::set<int>& S,const std::set<int>& locus,int n)
{
  int i,j;
  std::set<int> fc;
  std::set<int>::const_iterator it,jt;
  std::vector<int> vec,vx;
  SYNARMOSMA::hash_map::const_iterator qt;
  const int d = (signed) S.size() - 1;

  Simplex s(S,locus);
  s.incept = n;

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
    bool output = false;

    for(it=locus.begin(); it!=locus.end(); ++it) {
      if (!simplices[d][qt->second].active(*it)) {
        simplices[d][qt->second].activate(*it);
        output = true;
      }
    }
    if (!output) return false;
  }

#ifdef VERBOSE
  std::cout << "Adding a " << d << "-simplex to the spacetime complex..." << std::endl;
#endif

  for(it=S.begin(); it!=S.end(); ++it) {
    events[*it].topology_modified = true;
    for(jt=locus.begin(); jt!=locus.end(); ++jt) {
      events[*it].activate(*jt);
    }
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
      s.initialize(fc,locus);
      s.incept = n;
      simplices[i].push_back(s);
      index_table[i][fc] = simplices[i].size() - 1;
    }
    else {
      for(it=locus.begin(); it!=locus.end(); ++it) {
        if (!simplices[i][qt->second].active(*it)) simplices[i][qt->second].activate(*it);
      }
    }
    fc.clear();
    while(SYNARMOSMA::next_combination(vec,1+d)) {
      for(j=0; j<=i; ++j) {
        fc.insert(vx[vec[j]]);
      }
      qt = index_table[i].find(fc);
      if (qt == index_table[i].end()) {
        s.initialize(fc,locus);
        s.incept = n;
        simplices[i].push_back(s);
        index_table[i][fc] = simplices[i].size() - 1;
      }
      else {
        for(it=locus.begin(); it!=locus.end(); ++it) {
          if (!simplices[i][qt->second].active(*it)) simplices[i][qt->second].activate(*it);
        }
      }
      fc.clear();
    }
    vec.clear();
  }

  for(it=locus.begin(); it!=locus.end(); ++it) {
    simplicial_implication(*it);
  }
  return true;
}

void Complex::compute_modified_events()
{
  int i,j,n,m,l,nhop;
  std::set<int> S,vx,current,next;
  std::set<int>::const_iterator it,jt,kt;
  const int nv = (signed) events.size();
  bool done[nv];

  for(i=0; i<nv; ++i) {
    done[i] = false;
    if (!events[i].active()) continue;
    if (events[i].topology_modified) S.insert(i);
  }

  for(i=1; i<=Complex::ND; ++i) {
    n = (signed) simplices[i].size();
    for(j=0; j<n; ++j) {
      if (!simplices[i][j].modified) continue;
      simplices[i][j].get_vertices(vx);
      for(it=vx.begin(); it!=vx.end(); ++it) {
        events[*it].topology_modified = true;
        S.insert(*it);
      }
    }
  }

  // Now with this know we can calculate which events need to have their
  // entwinement and/or dimensional stress recalculated
#ifdef VERBOSE
  std::cout << "There are " << S.size() << " vertices directly modified in this relaxation step." << std::endl;
#endif
  for(it=S.begin(); it!=S.end(); ++it) {
    n = *it;
    nhop = 0;
    // Every vertex within topological_radius hops of n is labelled as
    // modified
    current.insert(n);
    done[n] = true;
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
    } while(nhop < Complex::topological_radius);
    current.clear();
    next.clear();
    for(i=0; i<nv; ++i) {
      if (done[i]) events[i].topology_modified = true;
      done[i] = false;
    }
  }
  int nmod = 0;
  for(i=0; i<nv; ++i) {
    if (events[i].topology_modified) nmod++;
  }
#ifdef VERBOSE
  std::cout << "There are " << nmod << " modified vertices out of " << nv << " in this relaxation step." << std::endl;
#endif
}

void Complex::simplicial_implication(int sheet)
{
  int i,j,k,n,m,vx[2];
  Simplex S;
  std::string sx;
  std::set<int>::const_iterator it;
  SYNARMOSMA::hash_map::const_iterator qt;
  const int ulimit = dimension(sheet);

  if (sheet == -1) {
    std::set<int> ubiquity;

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
              simplices[i-1][qt->second].activate(*it);
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
        events[vx[0]].activate(*it);
        events[vx[1]].activate(*it);
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
              simplices[i-1][qt->second].activate(sheet);
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
        events[vx[0]].activate(sheet);
      }
      if (!events[vx[1]].active(sheet)) {
#ifdef VERBOSE
        std::cout << "Restoring vertex " << vx[1] << " to regularize the complex" << std::endl;
#endif
        events[vx[1]].activate(sheet);
      }
    }
  }
}

void Complex::simplicial_implication(int base,int sheet) const
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

int Complex::entourage(int base,int sheet) const
{
  if (!events[base].active(sheet)) return -1;
  // Calculates a measure of this event's integration/implication in its
  // spacetime neighbourhood
  int i,j,n,m,output = 0;

  for(i=1; i<=Complex::ND; ++i) {
    n = 0;
    m = (signed) simplices[i].size();
    for(j=0; j<m; ++j) {
      if (!simplices[i][j].active(sheet)) continue;
      if (simplices[i][j].contains(base)) n++;
    }
    output += i*n;
  }
  return output;
}

int Complex::max_degree() const
{
  int i,n,output = 0;
  for(i=0; i<(signed) events.size(); ++i) {
    if (!events[i].active()) continue;
    n = (signed) events[i].neighbours.size();
    if (output < n) output = n;
  }
  return output;
}

double Complex::entwinement(int sheet) const
{
  // This method produces a real number between 0 and 1 that measures the
  // degree of "labyrinthicity" of the graph
  SYNARMOSMA::Graph G;
  compute_graph(&G,sheet);
  return G.entwinement();
}

double Complex::cyclic_resistance(int sheet) const
{
  SYNARMOSMA::Graph G;
  compute_graph(&G,sheet);
  return G.cyclic_resistance();
}

int Complex::combinatorial_distance(int v1,int v2,int sheet) const
{
  // A method to calculate the topological distance
  // between the two vertices v1 and v2
  if (v1 == v2) return 0;

  // A sanity check...
#ifdef DEBUG
  if (sheet == -1) {
    assert(events[v1].active() && events[v2].active());
  }
  else {
    assert(events[v1].active(sheet) && events[v2].active(sheet));
  }
#endif

  int d,offset[events.size()];
  SYNARMOSMA::Graph G;
  compute_graph(&G,offset,sheet);
  d = G.distance(offset[v1],offset[v2]);
  return d;
}

int Complex::cardinality(int d,int sheet) const
{
  if (d < 0) return 0;
  unsigned int i;
  int n = 0;

  if (sheet == -1) {
    if (d == 0) {
      const unsigned int M = events.size();
      for(i=0; i<M; ++i) {
        if (events[i].active()) n++;
      }
    }
    else {
      const unsigned int M = simplices[d].size();
      for(i=0; i<M; ++i) {
        if (simplices[d][i].active()) n++;
      }
    }
  }
  else {
    if (d == 0) {
      const unsigned int M = events.size();
      for(i=0; i<M; ++i) {
        if (events[i].active(sheet)) n++;
      }
    }
    else {
      const unsigned int M = simplices[d].size();
      for(i=0; i<M; ++i) {
        if (simplices[d][i].active(sheet)) n++;
      }
    }
  }
  return n;
}

int Complex::weighted_entourage(int n1,int n2) const
{
  int i,nfound,output = 0;
  unsigned int n,j;

  for(i=2; i<=Complex::ND; ++i) {
    nfound = 0;
    n = simplices[i].size();
    for(j=0; j<n; ++j) {
      if (!simplices[i][j].active()) continue;
      if (simplices[i][j].vertices.count(n1) > 0 && simplices[i][j].vertices.count(n2) > 0) nfound++;
    }
    output += i*nfound;
  }
  return output;
}

int Complex::vertex_valence(int v,int sheet) const
{
  int nd = 0;
  std::set<int>::const_iterator it;
  for(it=events[v].entourage.begin(); it!=events[v].entourage.end(); ++it) {
    if (simplices[1][*it].active(sheet)) nd++;
  }
  return nd;
}

int Complex::vertex_dimension(int v,int sheet) const
{
  int i,j,n;
  if (v < 0 || v >= (signed) events.size()) return -1;
  if (sheet == -1) {
    for(i=Complex::ND; i>=1; i--) {
      n = (signed) simplices[i].size();
      for(j=0; j<n; ++j) {
        if (!simplices[i][j].active()) continue;
        if (simplices[i][j].contains(v)) return i;
      }
    }
    if (!events[v].active()) return -1;
    return 0;
  }
  else {
    for(i=Complex::ND; i>=1; i--) {
      n = (signed) simplices[i].size();
      for(j=0; j<n; ++j) {
        if (!simplices[i][j].active(sheet)) continue;
        if (simplices[i][j].contains(v)) return i;
      }
    }
    if (events[v].active(sheet)) return 0;
  }
  return -1;
}

double Complex::dimensional_frontier(int D,int sheet) const
{
  unsigned int i,s = 0,ne = 0;
  int vx[2];
  const unsigned int n = simplices[1].size();

  for(i=0; i<n; ++i) {
    if (!simplices[1][i].active(sheet)) continue;
    simplices[1][i].get_vertices(vx);
    vx[0] = std::max(vertex_dimension(vx[0],sheet),D);
    vx[1] = std::max(vertex_dimension(vx[1],sheet),D);
    if (vx[0] != vx[1]) s++;
    ne++;
  }
  if (ne == 0) return 0.0;

  return double(s)/double(ne);
}

double Complex::dimensional_uniformity(int geometric_dimension,int sheet) const
{
  int i,d,nv = 0,sdimension = dimension(sheet),sum = 0;
  if (sdimension < geometric_dimension) sdimension = geometric_dimension;
  const int n = (signed) events.size();

  if (sheet == -1) {
    for(i=0; i<n; ++i) {
      if (!events[i].active()) continue;
      d = vertex_dimension(i,sheet);
      d = std::max(d,geometric_dimension);
      sum += d - sdimension;
      nv++;
    }
  }
  else {
    for(i=0; i<n; ++i) {
      if (!events[i].active(sheet)) continue;
      d = vertex_dimension(i,sheet);
      d = std::max(d,geometric_dimension);
      sum += d - sdimension;
      nv++;
    }
  }
  if (nv == 0) return 0.0;

  return double(sum)/double(nv);
}

bool Complex::connected(int sheet) const
{
  // For this method we calculate the 1-skeleton of the simplicial complex and then,
  // as a graph, determine its connectedness.
  SYNARMOSMA::Graph G; 

  compute_graph(&G,sheet);

  return G.connected();
}

