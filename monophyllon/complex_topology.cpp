#include "complex.h"

using namespace DIAPLEXIS;

void Complex::compute_simplicial_dimension()
{
  int i;
  const int nv = (signed) events.size();

  for(i=0; i<nv; ++i) {
    if (!events[i].active) continue;
    events[i].topological_dimension = vertex_dimension(i);
  }
}

void Complex::get_edge_topology(std::vector<std::set<int> >& vx) const
{
  int i;
  const int nv = cardinality(0);

  vx.clear();
  for(i=0; i<nv; ++i) {
    if (!events[i].active) continue;
    vx.push_back(events[i].neighbours);
  }
}

void Complex::compute_degree_distribution(bool logarithmic) const
{
  std::vector<double> histogram;
  SYNARMOSMA::Graph* G = new SYNARMOSMA::Graph;

  compute_graph(G);

  G->degree_distribution(logarithmic,histogram);

  std::cout << "Vertex degree histogram:" << std::endl;
  for(int i=0; i<(signed) histogram.size(); ++i) {
    if (histogram[i] > 0) std::cout << i+1 << "  " << 100.0*histogram[i] << std::endl;
  }
  delete G;
}

void Complex::compute_connectivity_distribution(bool direct) const
{
  int i,j,m = 0,l = cardinality(0);
  std::vector<int> pcount;
  const int nv = (signed) events.size();
  const double ntotal = double(l*(l-1)/2);

  for(i=0; i<nv; ++i) {
    pcount.push_back(0);
  }

  // If we know memory is abundant relative to the spacetime 
  // size, we can compute the distances all at once
  if (direct) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i,j,l) reduction(max:m) schedule(dynamic,1)
#endif
    for(i=0; i<nv; ++i) {
      if (!events[i].active) continue;
      for(j=1+i; j<nv; ++j) {
        if (!events[j].active) continue;
        l = combinatorial_distance(i,j);
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
    int offset[nv];
    SYNARMOSMA::Graph G;
    SYNARMOSMA::pair_index distances;
    SYNARMOSMA::pair_index::const_iterator qt;

    compute_graph(&G,offset);
    G.compute_distances(distances);
    for(i=0; i<nv; ++i) {
      if (!events[i].active) continue;
      for(j=1+i; j<nv; ++j) {
        if (!events[j].active) continue;
        qt = distances.find(std::pair<int,int>(offset[i],offset[j]));
        l = qt->second;
        if (l > m) m = l;
        pcount[l] += 1;
      }
    }
  }
  std::cout << "Combinatorial distance histogram:" << std::endl;  
  for(i=1; i<=m; ++i) {
    std::cout << i << "  " << 100.0*double(pcount[i])/ntotal << std::endl;
  }
}

bool Complex::active_element(int d,int n) const
{
  if (d == 0) {
    return events[n].active;
  }
  else {
    return simplices[d][n].active;
  }
}

std::pair<double,double> Complex::random_walk() const
{
  std::pair<double,double> output;
  SYNARMOSMA::Graph G;

  compute_graph(&G);
  int L = G.size()/4;
  // Do a random walk with a length equal to a quarter of the graph's size and 
  // a fifth of its vertices as starting point.  
  output = G.random_walk(L,0.2);

  return output;
}

bool Complex::reduction(int base)
{
  const int d = vertex_dimension(base);
  if (d < 2) return false;
  int i,n,m,vx[1+d];
  std::set<int> candidates,s1;
  SYNARMOSMA::hash_map::const_iterator qt;
  const int N = (signed) simplices[d].size();

  // Find which edges this vertex possesses are used in n-simplices (n > 1) and
  // eliminate one of them...
  for(i=0; i<N; ++i) {
    if (!simplices[d][i].active) continue;
    if (simplices[d][i].contains(base)) candidates.insert(i);
  }
  m = RND->irandom(candidates);
  simplices[d][m].get_vertices(vx);
  do {
    n = RND->irandom(1+d);
    if (vx[n] != base) break;
  } while(true);
  // Delete the edge v:vx[n]...
#ifdef VERBOSE
  std::cout << "Deleting edge with key " << base << ":" << vx[n] << std::endl;
#endif
  s1.insert(base); s1.insert(vx[n]);
  qt = index_table[1].find(s1);
  simplex_deletion(1,qt->second);
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

void Spacetime::compute_topological_dependency(const std::set<int>& vx)
{
  int i,n,m,l,nhop,nmod = 0;
  std::set<int> next,current;
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
    } while(nhop < Complex::topological_radius);
    current.clear();
    next.clear();
    for(i=0; i<nv; ++i) {
      if (done[i] == 1) events[i].topology_modified = true;
      done[i] = 0;
    }
  }
  for(i=0; i<nv; ++i) {
    if (events[i].topology_modified) nmod++;
  }
#ifdef VERBOSE
  std::cout << "There are " << nmod << " modified vertices out of " << nv << std::endl;
#endif
}

double Complex::dimensional_stress(int d,int n) const
{
  // This method measures the standard deviation of the
  // simplicial dimensions of the vertices of a given
  // d-simplex
  int i,vx[1+d];
  double alpha,v[1+d],sigma = 0.0,mu = 0.0;

  simplices[d][n].get_vertices(vx);

  for(i=0; i<1+d; ++i) {
    alpha = double(vertex_dimension(vx[i]));
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
  // The edge n has had its orientation changed, so now we need to recompute 
  // the orientation of the higher-dimensional simplices that depend on it.
  std::set<int> temp;
  
  temp.insert(n);
  recompute_parity(temp);
}

void Complex::recompute_parity(const std::set<int>& edges) 
{
  // The 1-simplices in the argument have had their orientation changed, so now 
  // we need to recompute the orientation of the higher-dimensional simplices that 
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
  const int d = dimension();

  for(i=2; i<=d; ++i) {
    for(j=0; j<(signed) simplices[i].size(); ++j) {
      if (!simplices[i][j].modified) continue;
      compute_simplex_parity(i,j);
    }
  }
}

void Complex::compute_fvector(std::vector<int>& f,std::vector<int>& fstar) const
{
  int i,j,s = 0;
  const int D = dimension();

  f.clear();
  fstar.clear();

  // First the f-vector
  for(i=0; i<=D; ++i) {
    f.push_back(cardinality(i));
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

void Complex::compute_hvector(std::vector<int>& h) const
{
  int i,j,sum;
  std::vector<int> f,fstar;
  const int D = dimension();

  compute_fvector(f,fstar);
  for(i=0; i<=D; ++i) {
    sum = SYNARMOSMA::ipow(-1,i)*int(SYNARMOSMA::binomial(D+1,i));
    for(j=1; j<=i; ++j) {
      sum += SYNARMOSMA::ipow(-1,i-j)*int(SYNARMOSMA::binomial(D+1-j,i-j))*f[j-1];
    }
    h.push_back(sum);
  }
}

void Complex::vertex_degree_statistics(double* output) const
{
  SYNARMOSMA::Graph G;  
  compute_graph(&G);
  output[0] = double(G.max_degree());
  output[1] = double(G.min_degree());
  output[2] = G.average_degree();
}

int Complex::cyclicity() const
{
  // This method should calculate the number of cyclic edges in the complex
  // associated with the sheet, assuming the complex is connected.
  SYNARMOSMA::Graph G;

  compute_graph(&G);
  if (!G.connected()) return 0;
  return (G.size() - G.bridge_count());
}

double Complex::set_logical_atoms(int n)
{ 
#ifdef DEBUG
  assert(n > 0);
#endif
  int i,j,natoms;
  double sigma,output = 0.0;
  std::set<int> cset;
  const int nvertex = (signed) events.size();

  for(int i=0; i<nvertex; ++i) {
    events[i].theorem.clear();
  }
 
  // Set the logical atoms in a purely random manner, the argument 
  // "n" represents the total number of propositional atoms in the 
  // entire spacetime
  for(i=0; i<nvertex; ++i) {
    if (!events[i].active) continue;
    cset.clear();
    // The more energetic and the higher the topological dimension of a 
    // vertex, the greater the number of atomic propositions in its theorem 
    // property.
    sigma = (1.0 + events[i].get_energy())*double(4 + events[i].topological_dimension);
    sigma *= RND->drandom(1.0,1.5);
    natoms = int(sigma);
    if (natoms >= n) {
      for(j=0; j<n; ++j) {
        cset.insert(j); 
      }
    }
    else {
      do {
        cset.insert(RND->irandom(n));
        if ((signed) cset.size() == natoms) break;
      } while(true);
    }
    events[i].theorem.set_atoms(cset);
    output += double(cset.size());   
  }
  output = output/double(cardinality(0));
  return output;
}

double Complex::logical_energy(int v) const
{
  if (events[v].neighbours.empty()) return 0.0;
  double sum = 0.0;
  std::set<int>::const_iterator it;
  SYNARMOSMA::Proposition q,p = events[v].theorem;

  for(it=events[v].neighbours.begin(); it!=events[v].neighbours.end(); ++it) {
    q = p & events[*it].theorem;
    sum += double(q.satisfiable());
  }
  sum = sum/double(events[v].neighbours.size());
  return sum;
}

bool Complex::logical_conformity(int v) const
{
  std::set<int>::const_iterator it;
  SYNARMOSMA::Proposition Q = events[v].theorem;

  for(it=events[v].neighbours.begin(); it!=events[v].neighbours.end(); ++it) {
    Q = Q & events[*it].theorem;
  }
  return Q.satisfiable();
}

void Complex::compute_simplex_energy(int d,int n)
{
  int i,vx[1+d];
  double alpha = 0.0;

  simplices[d][n].get_vertices(vx);

  for(i=0; i<1+d; ++i) {
    alpha += events[vx[i]].get_energy();
  }
  simplices[d][n].energy = alpha/double(1+d);
}

void Complex::compute_simplex_parity(int d,int n)
{
  if (d < 2) return;
  int i,j,vx[1+d];
  std::set<int> S;
  SYNARMOSMA::hash_map::const_iterator qt;

  simplices[d][n].get_vertices(vx);

  simplices[d][n].parity = 1;
  // A single undirected edge means the whole simplex is undirected...
  for(i=0; i<=d; ++i) {
    for(j=1+i; j<=d; ++j) {
      S.insert(vx[i]); S.insert(vx[j]);
      qt = index_table[1].find(S);
      simplices[d][n].parity *= simplices[1][qt->second].parity;
      S.clear();     
    }
  }
}

double Complex::dimensional_stress(int v) const
{
  // This method should sum the dimensional stress associated with
  // each d-simplex (d > 0) that contains the vertex v and exists on
  // the specified sheet.
  int i,j,ds;
  const int nd = dimension();
  double sum = 0.0;

  for(i=nd; i>0; i--) {
    ds = (signed) simplices[i].size();
    for(j=0; j<ds; ++j) {
      if (!simplices[i][j].active) continue;
      if (simplices[i][j].contains(v)) sum += dimensional_stress(i,j);
    }
  }
  return sum;
}

void Complex::simplex_membership(int v,std::vector<int>& output) const
{
  int i,j,m,k = 0;

  output.clear();
  for(i=1; i<=Complex::ND; ++i) {
    m = (signed) simplices[i].size();
    for(j=0; j<m; ++j) {
      if (!simplices[i][j].active) continue;
      if (simplices[i][j].contains(v)) ++k;
    }
    output.push_back(k);
    k = 0;
  }
}

void Complex::compute_graph(SYNARMOSMA::Graph* G) const
{
  int offset[events.size()];
  compute_graph(G,offset);
}

void Complex::compute_graph(SYNARMOSMA::Graph* G,int* offset) const
{
  int i,vx[2];
  const int nv = (signed) events.size();
  const int ne = (signed) simplices[1].size();

  G->clear();

  for(i=0; i<nv; ++i) {
    offset[i] = -1;
    if (!events[i].active) continue;
    offset[i] = G->add_vertex();
  }
  for(i=0; i<ne; ++i) {
    if (!simplices[1][i].active) continue;
    simplices[1][i].get_vertices(vx);
    G->add_edge(offset[vx[0]],offset[vx[1]]);
  }
}

void Complex::compute_graph(SYNARMOSMA::Graph* G,int base,int steps) const
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
        if (!simplices[1][qt->second].active) continue;
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

void Complex::compute_global_nexus(SYNARMOSMA::Nexus* NX) const
{
  int i,j;
  std::vector<int> offset;
  std::set<int> vx;
  std::set<int>::const_iterator it;
  const int nv = (signed) events.size();
  const int n = dimension();

  NX->initialize(n);
  for(i=0; i<nv; ++i) {
    offset.push_back(-1);
  }

  for(i=0; i<nv; ++i) {
    if (!events[i].active) continue;
    offset[i] = NX->add_vertex();
  }
  for(i=n; i>=1; --i) {
    for(j=0; j<(signed) simplices[i].size(); ++j) {
      if (!simplices[i][j].active) continue;
      for(it=simplices[i][j].vertices.begin(); it!=simplices[i][j].vertices.end(); ++it) {
        vx.insert(offset[*it]);
      }
      NX->paste(vx);
      vx.clear();
    }
  }
  NX->regularization();
}

void Complex::compute_local_nexus(SYNARMOSMA::Nexus* NX,int base) const
{
  int i,j;
  std::vector<int> offset;
  std::set<int>::const_iterator it;
  std::set<int> vx;
  const int n = vertex_dimension(base);
  const int nv = (signed) events.size();

  for(i=0; i<nv; ++i) {
    offset.push_back(-1);
  }
  NX->initialize(n);
  offset[base] = NX->add_vertex();

  for(i=1; i<=n; ++i) {
    for(j=0; j<(signed) simplices[i].size(); ++j) {
      if (!simplices[i][j].active) continue;
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
  NX->regularization();
}

int Complex::chromatic_number() const
{
  // Computes the chromatic number chi of the graph associated with the
  // colour "p"; we already know that 1 <= chi <= max_degree+1.  
  SYNARMOSMA::Graph G;

  compute_graph(&G);
  // It makes no sense to try to compute a single chromatic number for
  // a disconnected graph.
  if (!G.connected()) return 0;
  return G.chromatic_number();
}

bool Complex::edge_parity_mutation(int base)
{
  int n;
  std::set<int> candidates;
  std::set<int>::const_iterator it;

  for(it=events[base].entourage.begin(); it!=events[base].entourage.end(); ++it) {
    if (simplices[1][*it].active) {
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

bool Complex::edge_parity_mutation(int u,int v)
{
  // This is the method used for dynamic hyphansis where there is a single call 
  // to recompute the parity of all the higher-dimensional simplices, so no 
  // need to include one at the method's end. 
  int n;
  std::set<int> S;
  SYNARMOSMA::hash_map::const_iterator qt;

  S.insert(u); S.insert(v);
  qt = index_table[1].find(S);
  if (qt == index_table[1].end()) return false;
  if (!simplices[1][qt->second].active) return false;
  n = qt->second;
  if (simplices[1][n].parity == 0) {
    simplices[1][n].parity = (RND->irandom(2) == 0) ? 1 : -1;
  }
  else {
    simplices[1][n].parity *= -1;
  }
  return true;
}

void Complex::compute_entourages()
{
  int i,j,k,ns,vx[2];
  std::set<int> s,v;
  std::set<int>::const_iterator it;
  const int ulimit = dimension();
  SYNARMOSMA::hash_map::const_iterator qt;

  // What about removing items from the entourage of a d-simplex, when this item has
  // changed its ubiquity?
  for(i=1; i<=Complex::ND; ++i) {
    for(j=0; j<(signed) simplices[i].size(); ++j) {
      if (!simplices[i][j].active) continue;
      for(it=simplices[i][j].entourage.begin(); it!=simplices[i][j].entourage.end(); ++it) {
        if (!simplices[i+1][*it].active) continue;
        s.insert(*it);
      }
      simplices[i][j].entourage = s;
      s.clear();
    }
  }
  for(i=0; i<(signed) events.size(); ++i) {
    if (!events[i].active) continue;
    for(it=events[i].entourage.begin(); it!=events[i].entourage.end(); ++it) {
      if (!simplices[1][*it].active) continue;
      s.insert(*it);
    }
    events[i].entourage = s;
    s.clear();
  }

  for(i=ulimit; i>=2; i--) {
    ns = (signed) simplices[i].size();
    for(j=0; j<ns; ++j) {
      if (!simplices[i][j].active) continue;
      for(k=0; k<1+i; ++k) {
        v = simplices[i][j].faces[k];
        qt = index_table[i-1].find(v);
        if (qt == index_table[i-1].end()) throw std::runtime_error("Missing entourage element!");
        simplices[i-1][qt->second].active = true;
        simplices[i-1][qt->second].entourage.insert(j);
      }
    }
  }
  // Now the edges...
  ns = (signed) simplices[1].size();
  for(i=0; i<ns; ++i) {
    if (!simplices[1][i].active) continue;
    simplices[1][i].get_vertices(vx);

    events[vx[0]].active = true;
    events[vx[0]].entourage.insert(i);

    events[vx[1]].active = true;
    events[vx[1]].entourage.insert(i);
  }
}

void Complex::compute_neighbours()
{
  int i,vx[2];

  for(i=0; i<(signed) events.size(); ++i) {
    events[i].neighbours.clear();
  }
  for(i=0; i<(signed) simplices[1].size(); ++i) {
    if (!simplices[1][i].active) continue;
    simplices[1][i].get_vertices(vx);
    events[vx[0]].neighbours.insert(vx[1]);
    events[vx[1]].neighbours.insert(vx[0]);
  }
}

void Complex::simplex_deletion(int d,int n)
{
  std::set<int>::const_iterator it;
  std::set<int> parents;
  int i,dp1 = d + 1;
  
  simplices[d][n].active = false;
  parents = simplices[d][n].entourage;
  for(it=parents.begin(); it!=parents.end(); ++it) {
    i = *it;
    simplex_deletion(dp1,i);
  }
}

bool Complex::simplex_addition(const std::set<int>& S,std::set<int>& vx_delta)
{
  int i,j;
  std::set<int> fc;
  std::set<int>::const_iterator it;
  std::vector<int> vec,vx;
  SYNARMOSMA::hash_map::const_iterator qt;
  const int d = (signed) S.size() - 1;
  Simplex s(S);

  qt = index_table[d].find(S);
  if (qt == index_table[d].end()) {
    simplices[d].push_back(s);
    index_table[d][S] = simplices[d].size() - 1;
  }
  else {
    if (simplices[d][qt->second].active) {
      return false;
    }
    else {
      simplices[d][qt->second].active = true;
    }
  }

#ifdef VERBOSE
  std::cout << "Adding a " << d << "-simplex to the spacetime complex..." << std::endl;
#endif

  for(it=S.begin(); it!=S.end(); ++it) {
    vx_delta.insert(*it);
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
      simplices[i].push_back(Simplex(fc));
      index_table[i][fc] = simplices[i].size() - 1;
    }
    else {
      simplices[i][qt->second].active = true;
    }
    fc.clear();
    while(SYNARMOSMA::next_combination(vec,1+d)) {
      for(j=0; j<=i; ++j) {
        fc.insert(vx[vec[j]]);
      }
      qt = index_table[i].find(fc);
      if (qt == index_table[i].end()) {
        simplices[i].push_back(Simplex(fc));
        index_table[i][fc] = simplices[i].size() - 1;
      }
      else {
        simplices[i][qt->second].active = true;
      }
      fc.clear();
    }
    vec.clear();
  }
  simplicial_implication();
  return true;
}

void Complex::inversion()
{
  int i,j;
  std::set<int> hold;
  std::set<int>::const_iterator it;
  Simplex S;
  const int nv = (signed) events.size();

  for(i=0; i<nv; ++i) {
    // hold = V / events[i].neighbours
    for(j=0; j<nv; ++j) {
      if (!events[j].active) continue;
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
        S.initialize(i,j);
        simplices[1].push_back(S);
        index_table[1][S.vertices] = (signed) simplices[1].size() - 1;
      }
    }
  }
}

void Complex::simplicial_implication()
{
  int i,j,k,n,m,vx[2];
  const int ulimit = dimension();
  Simplex S;
  SYNARMOSMA::hash_map::const_iterator qt;
  std::set<int>::const_iterator itz;

  for(i=ulimit; i>=2; i--) {
    n = (signed) simplices[i].size();
    m = (signed) simplices[i-1].size();
    for(j=0; j<n; ++j) {
      if (!simplices[i][j].active) continue;
      for(k=0; k<1+i; ++k) {
        qt = index_table[i-1].find(simplices[i][j].faces[k]);
        if (qt == index_table[i-1].end()) {
          S.initialize(simplices[i][j].faces[k]);
          simplices[i-1].push_back(S);
          index_table[i-1][S.vertices] = m;
          m++;
        }
        else {
          simplices[i-1][qt->second].active = true;
        }
      }
    }
  }
  n = (signed) simplices[1].size();
  for(i=0; i<n; ++i) {
    if (!simplices[1][i].active) continue;
    simplices[1][i].get_vertices(vx);
    events[vx[0]].active = true;
    events[vx[1]].active = true;
  }
}

void Complex::simplicial_implication(int base) const
{
  // This method will calculate all of the n-simplices (n > 1) that are "implied" by
  // the base vertex and its neighbours (via their mutual edges) and list them, as well
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
            if (!simplices[1][qt->second].active) {
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
  // So how many of these implied simplices actually exist? And having found one or more which officially don't
  // exist, what to do with them? If deleting a single edge can eliminate an entire tower of dependent n-simplices
  // (n > 1) then shouldn't recreating this same edge cause the tower of n-simplices to be restored?
  for(i=M; i>2; --i) {
    nsimp = 0;
    nfound = 0;
    for(it=implied_simplex[i].begin(); it!=implied_simplex[i].end(); ++it) {
      qt = index_table[i-1].find(*it);
      if (qt != index_table[i-1].end()) {
        if (simplices[i-1][qt->second].active) nfound++;
      }
      nsimp++;
    }
#ifdef VERBOSE
    std::cout << "There are " << nsimp << " implied " << i-1 << "-simplices of which " << nfound << " already exist in the spacetime complex." << std::endl;
#endif
  }
  delete[] implied_simplex;
}

int Complex::entourage(int base) const
{
  // Calculates a measure of this event's integration/implication in its
  // spacetime neighbourhood
  int i,j,n,m,output = 0;

  for(i=1; i<=Complex::ND; ++i) {
    n = 0;
    m = (signed) simplices[i].size();
    for(j=0; j<m; ++j) {
      if (!simplices[i][j].active) continue;
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
    if (!events[i].active) continue;
    n = (signed) events[i].neighbours.size();
    if (output < n) output = n;
  }
  return output;
}

double Complex::entwinement() const
{
  // This method produces a real number between 0 and 1 that measures the
  // degree of "labyrinthicity" of the graph
  SYNARMOSMA::Graph G;
  compute_graph(&G);
  return G.entwinement();
}

double Complex::cyclic_resistance() const
{
  SYNARMOSMA::Graph G;
  compute_graph(&G);
  return G.cyclic_resistance();
}

int Complex::combinatorial_distance(int v1,int v2) const
{
  // A method to calculate the topological distance
  // between the two vertices v1 and v2
  if (v1 == v2) return 0;

  // A sanity check...
#ifdef DEBUG
  assert(events[v1].active && events[v2].active);
#endif

  int d,offset[events.size()];
  SYNARMOSMA::Graph G;
  compute_graph(&G,offset);
  d = G.distance(offset[v1],offset[v2]);
  return d;
}

int Complex::cardinality_safe(int d) const
{
  int i,n = 0;
  if (d == 0) {
    const int M = (signed) events.size();
    for(i=0; i<M; ++i) {
      if (!events[i].active) continue;
      n++;
    }
  }
  else {
    const int M = (signed) simplices[d].size();
    for(i=0; i<M; ++i) {
      if (!simplices[d][i].active) continue;
      n++;
    }
  }
  return n;
}

int Complex::circuit_rank() const
{
  int output = cardinality(1) - cardinality(0);
  if (connected()) {
    return 1 + output;
  }
  std::vector<int> components;
  int n = component_analysis(components);
  return n + output;
}

int Complex::euler_characteristic() const
{
  int i,pf = 1,chi = 0;
  const int D = dimension();
  for(i=0; i<=D; ++i) {
    chi += pf*cardinality(i);
    pf *= -1;
  }
  return chi;
}

int Complex::dimension() const
{
  int i,j;

  for(i=Complex::ND; i>0; i--) {
    for(j=0; j<(signed) simplices[i].size(); ++j) {
      if (simplices[i][j].active) return i;
    }
  }
  for(i=0; i<(signed) events.size(); ++i) {
    if (events[i].active) return 0;
  }
  return -1;
}

int Complex::total_dimension() const
{
  int i,sum = 0;

  for(i=0; i<(signed) events.size(); ++i) {
    if (!events[i].active) continue;
    sum += vertex_dimension(i);
  }
  return sum;
}

int Complex::structural_index() const
{
  int i,j,l,n,sum = 0,d = dimension();


  for(i=0; i<(signed) events.size(); ++i) {
    if (!events[i].active) continue;
    sum++;
  }
  for(i=1; i<d; ++i) {
    l = 0;
    n = (signed) simplices[i].size();
    for(j=0; j<n; ++j) {
      if (!simplices[i][j].active) continue;
      l++;
    }
    sum += (1+i)*l;
  }
  return sum;
}

int Complex::weighted_entourage(int n1,int n2) const
{
  int i,j,nfound,output = 0;
  bool f1,f2;
  std::set<int>::const_iterator it;

  for(i=2; i<=Complex::ND; ++i) {
    nfound = 0;
    for(j=0; j<(signed) simplices[i].size(); ++j) {
      it = std::find(simplices[i][j].vertices.begin(),simplices[i][j].vertices.end(),n1);
      f1 = (it == simplices[i][j].vertices.end()) ? false : true;

      it = std::find(simplices[i][j].vertices.begin(),simplices[i][j].vertices.end(),n2);
      f2 = (it == simplices[i][j].vertices.end()) ? false : true;

      if (f1 && f2) nfound++;
    }
    output += i*nfound;
  }
  return output;
}

int Complex::vertex_valence(int v) const
{
  int nd = 0;
  std::set<int>::const_iterator it;
  for(it=events[v].entourage.begin(); it!=events[v].entourage.end(); ++it) {
    if (simplices[1][*it].active) nd++;
  }
  return nd;
}

int Complex::vertex_dimension(int v) const
{
  int i,j,n;
  if (v < 0 || v >= (signed) events.size()) return -1;

  if (!events[v].active) return -1;

  for(i=Complex::ND; i>=1; i--) {
    n = (signed) simplices[i].size();
    for(j=0; j<n; ++j) {
      if (!simplices[i][j].active) continue;
      if (simplices[i][j].contains(v)) return i;
    }
  }
  return 0;
}

double Complex::dimensional_frontier(int D) const
{
  int i,d[2],s = 0;
  for(i=0; i<(signed) simplices[1].size(); ++i) {
    if (!simplices[1][i].active) continue;
    simplices[1][i].get_vertices(d);
    d[0] = (d[0] < D) ? D : d[0];
    d[1] = (d[1] < D) ? D : d[1];
    if (d[0] != d[1]) s++;
  }
  return double(s)/double(simplices[1].size());
}

double Complex::dimensional_uniformity(int geometric_dimension) const
{
  int i,n,nv,sdimension = dimension(),sum = 0;
  if (sdimension < geometric_dimension) sdimension = geometric_dimension;

  nv = 0;
  for(i=0; i<(signed) events.size(); ++i) {
    if (!events[i].active) continue;
    n = vertex_dimension(i);
    n = (n < geometric_dimension) ? geometric_dimension : n;
    sum += sdimension - n;
    nv++;
  }
  return double(sum)/double(nv);
}

bool Complex::consistent() const
{
  int i,j,k,l,n,m,vx[2];
  bool found;
  std::set<int>::const_iterator it;
  std::set<int> S;
  SYNARMOSMA::hash_map::const_iterator qt;
  const int nv = (signed) events.size();
  const int ulimit = dimension();

#ifdef DEBUG
  assert(simplices[0].empty());
  assert(index_table[0].empty());
#endif

  for(i=Complex::ND; i>=2; i--) {
    n = (signed) simplices[i].size();
    for(j=0; j<n; ++j) {
      if (!simplices[i][j].active) continue;
      for(it=simplices[i][j].entourage.begin(); it!=simplices[i][j].entourage.end(); ++it) {
        if (!simplices[i+1][*it].active) {
          std::cout << "Error with entourage: " << i << "  " << j << "  " << *it << "  " << ulimit << std::endl;
          return false;
        }
      }
      if (simplices[i][j].dimension() != i) {
        std::cout << i << "-simplex  " << j << " has dimension " << simplices[i][j].dimension() << std::endl;
        return false;
      }
      for(k=0; k<1+i; ++k) {
        S = simplices[i][j].faces[k];
        found = false;
        m = (signed) simplices[i-1].size();
        for(l=0; l<m; ++l) {
          if (!simplices[i-1][l].active) continue;
          if (S == simplices[i-1][l].vertices) {
            found = true;
            break;
          }
        }
        if (!found) {
          std::cout << "Can't find face " << SYNARMOSMA::make_key(S) << " of simplex " << SYNARMOSMA::make_key(simplices[i][j].vertices) << std::endl;
          for(l=0; l<(signed) simplices[i-1].size(); ++l) {
            if (!simplices[i-1][l].active) continue;
            std::cout << SYNARMOSMA::make_key(simplices[i-1][l].vertices) << "  " << simplices[i-1][l].active << std::endl;
          }
          return false;
        }
        qt = index_table[i-1].find(S);
        if (qt == index_table[i-1].end()) {
          std::cout << "Problem with the index tables for " << SYNARMOSMA::make_key(simplices[i][j].vertices) << " and " << SYNARMOSMA::make_key(S) << std::endl;
          return false;
        }
      }
      for(it=simplices[i][j].vertices.begin(); it!=simplices[i][j].vertices.end(); ++it) {
        l = *it;
        if (l < 0 || l >= nv) {
          std::cout << "Nonexistent vertex: " << l << std::endl;
          return false;
        }
        if (!events[l].active) {
          std::cout << "Inactive vertex: " << l << std::endl;
          return false;
        }
      }
    }
  }
  for(i=0; i<(signed) simplices[1].size(); ++i) {
    if (!simplices[1][i].active) continue;
    simplices[1][i].get_vertices(vx);
    if (vx[0] < 0 || vx[0] >= nv) return false;
    if (vx[1] < 0 || vx[0] >= nv) return false;
    if (!events[vx[0]].active) {
      std::cout << "Inactive vertex: " << vx[0] << std::endl;
      return false;
    }
    if (!events[vx[1]].active) {
      std::cout << "Inactive vertex: " << vx[1] << std::endl;
      return false;
    }
    if (events[vx[0]].neighbours.count(vx[1]) == 0) {
      std::cout << "Edge in simplices array but not neighbour set: " << vx[0] << "  " << vx[1] << std::endl;
      std::cout << SYNARMOSMA::make_key(simplices[1][i].vertices) << "  " << simplices[1][i].incept << std::endl;
      return false;
    }
    if (events[vx[1]].neighbours.count(vx[0]) == 0) {
      std::cout << "Edge in simplices array but not neighbour set: " << vx[1] << "  " << vx[0] << std::endl;
      std::cout << SYNARMOSMA::make_key(simplices[1][i].vertices) << "  " << simplices[1][i].incept << std::endl;
      return false;
    }
  }
  for(i=0; i<nv; ++i) {
    if (!events[i].active) continue;
    if (events[i].neighbours.empty()) 
    for(it=events[i].entourage.begin(); it!=events[i].entourage.end(); ++it) {
      if (!simplices[1][*it].active) {
        std::cout << "Error with entourage: " << i << "  " << *it << std::endl;
        return false;
      }
    }
    for(it=events[i].neighbours.begin(); it!=events[i].neighbours.end(); ++it) {
      n = *it;
      if (n == i) {
        std::cout << "Illegal link: " << i << "  " << i << std::endl;
        return false;
      }
      else if (n < 0) {
        std::cout << "Illegal link: " << i << "  " << n << std::endl;
        return false;
      }
      else if (n >= nv) {
        std::cout << "Illegal link: " << i << "  " << n << std::endl;
        return false;
      }
      S.clear();
      S.insert(i);
      S.insert(n);
      qt = index_table[1].find(S);
      if (qt == index_table[1].end()) {
        std::cout << "Edge in neighbour set but not in index table!" << std::endl;
        return false;
      }
    }
  }
  // Make sure that each n-simplex only exists once...
  for(i=1; i<=Complex::ND; ++i) {
    n = (signed) simplices[i].size();
    for(j=0; j<n; ++j) {
      S = simplices[i][j].vertices;
      for(k=0; k<n; ++k) {
        if (k == j) continue;
        if (S == simplices[i][k].vertices) {
          std::cout << "Illegal " << i << "-simplex duplication with " << SYNARMOSMA::make_key(S) << std::endl;
          return false;
        }
      }
    }
  }
  return true;
}

bool Complex::connected() const
{
  // For this method we calculate the 1-skeleton of the simplicial complex and then,
  // as a graph, determine its connectedness.
  SYNARMOSMA::Graph G;

  compute_graph(&G);

  return G.connected();
}

int Complex::component_analysis(std::vector<int>& component) const
{
  int i,n,ct = 0;
  std::vector<int> cvalue;
  SYNARMOSMA::Graph G;

  compute_graph(&G);
  n = G.component_analysis(cvalue);

  // We now need to include the unused vertices in the output component vector
  component.clear();

  for(i=0; i<(signed) events.size(); ++i) {
    if (!events[i].active) {
      component.push_back(-1);
      continue;
    }
    component.push_back(cvalue[ct]);
    ct++;
  }

  return n;
}

