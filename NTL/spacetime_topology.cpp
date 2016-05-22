#include "spacetime.h"

extern SYNARMOSMA::Random RND;

using namespace DIAPLEXIS;

void Spacetime::compute_simplicial_dimension()
{
  int i;
  const int nv = (signed) events.size();

  for(i=0; i<nv; ++i) {
    if (events[i].ubiquity == 1) continue;
    events[i].topological_dimension = vertex_dimension(i,-1);
  }
}

void Spacetime::get_edge_topology(std::vector<std::set<int> >& vx) const
{
  int i;
  const int nv = cardinality(0,-1);

  vx.clear();
  for(i=0; i<nv; ++i) {
    if (events[i].ubiquity == 1) continue;
    vx.push_back(events[i].neighbours);
  }
}

void Spacetime::get_ubiquity(int d,int n,std::string& output) const
{
  const int nt = (signed) codex.size();

  if (d == 0) {
    output = "{";
    for(int i=0; i<nt-1; ++i) {
      if (NTL::divide(events[n].ubiquity,codex[i].colour) == 1) {
        output += "1,";
      }
      else {
        output += "0,";
      }
    }
    if (NTL::divide(events[n].ubiquity,codex[nt-1].colour) == 1) {
      output += "1}";
    }
    else {
      output += "0}";
    }
  }
  else {
    for(int i=0; i<nt-1; ++i) {
      if (NTL::divide(simplices[d][n].ubiquity,codex[i].colour) == 1) {
        output += "1,";
      }
      else {
        output += "0,";
      }
    }
    if (NTL::divide(simplices[d][n].ubiquity,codex[nt-1].colour) == 1) {
      output += "1}";
    }
    else {
      output += "0}";
    }
  }
}

void Spacetime::compute_degree_distribution(bool logarithmic,int sheet) const
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

void Spacetime::compute_connectivity_distribution(int sheet) const
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
  if (geometry->get_memory_type()) {
    int offset[nv];
    SYNARMOSMA::Graph G;
    SYNARMOSMA::edge_hash distances;
    SYNARMOSMA::edge_hash::const_iterator qt;

    compute_graph(&G,offset,sheet);
    G.compute_distances(distances);
    if (sheet == -1) {
      for(i=0; i<nv; ++i) {
        if (events[i].ubiquity == 1) continue;
        for(j=1+i; j<nv; ++j) {
          if (events[j].ubiquity == 1) continue;
          qt = distances.find(std::pair<int,int>(offset[i],offset[j]));
          l = qt->second;
          if (l > m) m = l;
          pcount[l] += 1;
        }
      }
    }
    else {
      for(i=0; i<nv; ++i) {
        if (NTL::divide(events[i].ubiquity,codex[sheet].colour) == 0) continue;
        for(j=1+i; j<nv; ++j) {
          if (NTL::divide(events[j].ubiquity,codex[sheet].colour) == 0) continue;
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
        if (events[i].ubiquity == 1) continue;
        for(j=1+i; j<nv; ++j) {
          if (events[j].ubiquity == 1) continue;
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
        if (NTL::divide(events[i].ubiquity,codex[sheet].colour) == 0) continue;
        for(j=1+i; j<nv; ++j) {
          if (NTL::divide(events[j].ubiquity,codex[sheet].colour) == 0) continue;
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

bool Spacetime::active_element(int d,int n) const
{
  bool output = true;

  if (d == 0) {
    if (events[n].ubiquity == 1) output = false;
  }
  else {
    if (simplices[d][n].ubiquity == 1) output = false;
  }
  return output;
}

void Spacetime::random_walk(double* mean,double* sdeviation,int sheet) const
{
  SYNARMOSMA::Graph G;
  compute_graph(&G,sheet);
  G.random_walk(mean,sdeviation,geometry->dimension());
}

double Spacetime::dimensional_stress(int d,int n,int sheet) const
{
  // This method measures the standard deviation of the
  // simplicial dimensions of the vertices of a given
  // d-simplex
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

void Spacetime::compute_fvector(std::vector<int>& f,std::vector<int>& fstar,int sheet) const
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
      s += SYNARMOSMA::ipow(-1,j-i)*int(SYNARMOSMA::factorial(j)/(SYNARMOSMA::factorial(i)*SYNARMOSMA::factorial(j-i)))*f[j];
    }
    fstar.push_back(s);
  }
  fstar.push_back(f[D]);
}

void Spacetime::vertex_degree_statistics(double* output,int sheet) const
{
  SYNARMOSMA::Graph G;  
  compute_graph(&G,sheet);
  output[0] = double(G.max_degree());
  output[1] = double(G.min_degree());
  output[2] = G.average_degree();
}

int Spacetime::cyclicity(int sheet) const
{
  // This method should calculate the number of cyclic edges in the complex
  // associated with the sheet, assuming the complex is connected.
  SYNARMOSMA::Graph G;

  compute_graph(&G,sheet);
  if (!G.connected()) return 0;
  return (G.size() - G.bridge_count());
}

bool Spacetime::logical_conformity(int v) const
{
  std::set<int>::const_iterator it;
  SYNARMOSMA::Proposition Q = events[v].theorem;
  for(it=events[v].neighbours.begin(); it!=events[v].neighbours.end(); ++it) {
    Q = Q*events[*it].theorem;
  }
  return Q.satisfiable();
}

void Spacetime::compute_simplex_energy(int d,int n)
{
  int i,vx[1+d];
  double alpha = 0.0;

  simplices[d][n].get_vertices(vx);

  for(i=0; i<1+d; ++i) {
    alpha += events[vx[i]].get_energy();
  }
  simplices[d][n].energy = alpha/double(1+d);
}

double Spacetime::dimensional_stress(int v,int sheet) const
{
  // This method should sum the dimensional stress associated with
  // each d-simplex (d > 0) that contains the vertex v and exists on
  // the specified sheet.
  int i,j,ds;
  const int nd = dimension(-1);
  double sum = 0.0;

  for(i=nd; i>0; i--) {
    ds = (signed) simplices[i].size();
    for(j=0; j<ds; ++j) {
      if (simplices[i][j].ubiquity == 1) continue;
      if (simplices[i][j].contains(v)) sum += dimensional_stress(i,j,sheet);
    }
  }
  return sum;
}

void Spacetime::simplex_membership(int v,std::vector<int>& output) const
{
  int i,j,m,k = 0;

  output.clear();
  for(i=1; i<=Spacetime::ND; ++i) {
    m = (signed) simplices[i].size();
    for(j=0; j<m; ++j) {
      if (simplices[i][j].ubiquity == 1) continue;
      if (simplices[i][j].contains(v)) ++k;
    }
    output.push_back(k);
    k = 0;
  }
}

void Spacetime::compute_graph(SYNARMOSMA::Graph* G,int sheet) const
{
  int offset[events.size()];

  compute_graph(G,offset,sheet);
}

void Spacetime::compute_graph(SYNARMOSMA::Graph* G,int* offset,int sheet) const
{
  int i,vx[2];
  const int nv = (signed) events.size();
  const int ne = (signed) simplices[1].size();

  G->clear();

  if (sheet == -1) {
    for(i=0; i<nv; ++i) {
      offset[i] = -1;
      if (events[i].ubiquity == 1) continue;
      offset[i] = G->add_vertex();
    }
    for(i=0; i<ne; ++i) {
      if (simplices[1][i].ubiquity == 1) continue;
      simplices[1][i].get_vertices(vx);
      G->add_edge(offset[vx[0]],offset[vx[1]]);
    }
  }
  else {
    for(i=0; i<nv; ++i) {
      offset[i] = -1;
      if (NTL::divide(events[i].ubiquity,codex[sheet].colour) == 0) continue;
      offset[i] = G->add_vertex();
    }
    for(i=0; i<ne; ++i) {
      if (NTL::divide(simplices[1][i].ubiquity,codex[sheet].colour) == 0) continue;
      simplices[1][i].get_vertices(vx);
      G->add_edge(offset[vx[0]],offset[vx[1]]);
    }
  }
}

void Spacetime::compute_graph(SYNARMOSMA::Graph* G,int base,int steps,int sheet) const
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
        if (NTL::divide(simplices[1][qt->second].ubiquity,codex[sheet].colour) == 0) continue;
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

void Spacetime::compute_causal_graph(SYNARMOSMA::Directed_Graph* G,int base,int sheet) const
{
  int i,j,l,v;
  SYNARMOSMA::RELATION rho;
  SYNARMOSMA::hash_map::const_iterator qt;
  std::set<int> S;
  std::set<int>::const_iterator it;
  std::vector<int> offset,current,next;
  std::vector<int>::const_iterator v_it;
  const int nv = (signed) events.size();

  G->clear();
  current.push_back(base);
  for(i=0; i<nv; ++i) {
    offset.push_back(-1);
  }
  offset[base] = G->add_vertex();

  if (sheet == -1) {
    do {
      for(v_it=current.begin(); v_it!=current.end(); ++v_it) {
        v = *v_it;
        for(it=events[v].neighbours.begin(); it!=events[v].neighbours.end(); ++it) {
          j = *it;
          S.clear();
          S.insert(v);
          S.insert(j);
          qt = index_table[1].find(S);
          if (simplices[1][qt->second].ubiquity == 1) continue;
          rho = simplices[1][qt->second].orientation;
          if (rho == SYNARMOSMA::DISPARATE) continue;
          if (offset[j] == -1) {
            offset[j] = G->add_vertex();
            next.push_back(j);
          }
          G->add_edge(offset[v],offset[j],rho);
        }
      }
      if (next.empty()) break;
      current = next;
      next.clear();
    } while(true);
  }
  else {
    do {
      for(v_it=current.begin(); v_it!=current.end(); ++v_it) {
        v = *v_it;
        for(it=events[v].neighbours.begin(); it!=events[v].neighbours.end(); ++it) {
          j = *it;
          S.clear();
          S.insert(v);
          S.insert(j);
          qt = index_table[1].find(S);
          l = qt->second;
          if (NTL::divide(simplices[1][l].ubiquity,codex[sheet].colour) == 0) continue;
          rho = simplices[1][qt->second].orientation;
          if (rho == SYNARMOSMA::DISPARATE) continue;
          if (offset[j] == -1) {
            offset[j] = G->add_vertex();
            next.push_back(j);
          }
          G->add_edge(offset[v],offset[j],rho);
        }
      }
      if (next.empty()) break;
      current = next;
      next.clear();
    } while(true);
  }
}

void Spacetime::compute_global_nexus(SYNARMOSMA::Nexus* NX,int sheet) const
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
      if (events[i].ubiquity == 1) continue;
      offset[i] = NX->add_vertex();
    }
    for(i=n; i>=1; --i) {
      for(j=0; j<(signed) simplices[i].size(); ++j) {
        if (simplices[i][j].ubiquity == 1) continue;
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
      if (NTL::divide(events[i].ubiquity,codex[sheet].colour) == 0) continue;
      offset[i] = NX->add_vertex();
    }
    for(i=n; i>=1; --i) {
      for(j=0; j<(signed) simplices[i].size(); ++j) {
        if (NTL::divide(simplices[i][j].ubiquity,codex[sheet].colour) == 0) continue;
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

void Spacetime::compute_local_nexus(SYNARMOSMA::Nexus* NX,int base,int sheet) const
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
        if (simplices[i][j].ubiquity == 1) continue;
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
        if (NTL::divide(simplices[i][j].ubiquity,codex[sheet].colour) == 0) continue;
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

void Spacetime::compute_lightcones()
{
  // This method only makes sense for the entire polycosmic spacetime complex, i.e. 
  // sheet = -1, with a Lorentzian metric
  if (geometry->get_euclidean()) return;
  int i,j,k;
  SYNARMOSMA::RELATION rho;
  SYNARMOSMA::hash_map::const_iterator qt;
  std::set<int>::const_iterator it,jt;
  std::set<int> pcurrent,fcurrent,old,v,null,S;
  const int nv = (signed) events.size();

  // First we compute the past of each vertex
  for(i=0; i<nv; ++i) {
    if (events[i].ubiquity == 1) continue;
    pcurrent.clear();
    fcurrent.clear();
    for(it=events[i].neighbours.begin(); it!=events[i].neighbours.end(); ++it) {
      j = *it;
      S.clear();
      S.insert(i);
      S.insert(j);
      qt = index_table[1].find(S);
      rho = simplices[1][qt->second].orientation;
      if (rho == SYNARMOSMA::DISPARATE) continue;
      if (rho == SYNARMOSMA::BEFORE) {
        pcurrent.insert(j);
      }
      else {
        fcurrent.insert(j);
      }
    }
    if (pcurrent.empty()) {
      events[i].anterior = null;
    }
    else {
      do {
        for(it=pcurrent.begin(); it!=pcurrent.end(); ++it) {
          j = *it;
          if (j == i) continue;
          for(jt=events[j].neighbours.begin(); jt!=events[j].neighbours.end(); ++jt) {
            k = *jt;
            S.clear();
            S.insert(j);
            S.insert(k);
            qt = index_table[1].find(S);
            rho = simplices[1][qt->second].orientation;
            if (rho == SYNARMOSMA::DISPARATE) continue;
            if (rho == SYNARMOSMA::BEFORE) {
              // Check to see if k doesn't already exist in the set "old"
              if (old.count(k) == 0) v.insert(k);
            }
          }
        }
        for(it=pcurrent.begin(); it!=pcurrent.end(); ++it) {
          old.insert(*it);
        }
        if (v.empty()) break;
        pcurrent = v;
        v.clear();
      } while(true);
      events[i].anterior = old;
      old.clear();
    }
    if (fcurrent.empty()) {
      events[i].posterior = null;
      continue;
    }
    do {
      for(it=fcurrent.begin(); it!=fcurrent.end(); ++it) {
        j = *it;
        if (j == i) continue;
        for(jt=events[j].neighbours.begin(); jt!=events[j].neighbours.end(); ++jt) {
          k = *jt;
          S.clear();
          S.insert(j);
          S.insert(k);
          qt = index_table[1].find(S);
          rho = simplices[1][qt->second].orientation;
          if (rho == SYNARMOSMA::DISPARATE) continue;
          if (rho == SYNARMOSMA::AFTER) {
            if (old.count(k) == 0) v.insert(k);
          }
        }
      }
      for(it=fcurrent.begin(); it!=fcurrent.end(); ++it) {
        old.insert(*it);
      }
      if (v.empty()) break;
      fcurrent = v;
      v.clear();
    } while(true);
    events[i].posterior = old;
    old.clear();
  }
}

double Spacetime::compute_temporal_nonlinearity(int sheet) const
{
  if (geometry->get_euclidean()) return 0.0;
  int i,j,k,nsink = 0,nsource = 0,causal_loop = 0;
  double output,nlinearity = 0.0;
  SYNARMOSMA::RELATION rho;
  std::set<int> past,future,fcurrent,pcurrent,v,null,old,S;
  std::set<int>::const_iterator it,jt;
  SYNARMOSMA::hash_map::const_iterator qt;
  SYNARMOSMA::Directed_Graph G;
  const int nv = (signed) events.size();
  const double na = double(cardinality(0,sheet));

  if (sheet == -1) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i,j,k,it,jt,qt,rho,G,old,v,pcurrent,fcurrent,past,future,S) reduction(+:nlinearity,nsource,nsink,causal_loop)
#endif
    for(i=0; i<nv; ++i) {
      if (events[i].ubiquity == 1) continue;
      // Now calculate the future and past lightcones for this vertex on this sheet...
      pcurrent.clear();
      fcurrent.clear();
      for(it=events[i].neighbours.begin(); it!=events[i].neighbours.end(); ++it) {
        j = *it;
        S.clear();
        S.insert(i);
        S.insert(j);
        qt = index_table[1].find(S);
        rho = simplices[1][qt->second].orientation;
        if (rho == SYNARMOSMA::DISPARATE) continue;
        if (rho == SYNARMOSMA::BEFORE) {
          pcurrent.insert(j);
        }
        else {
          fcurrent.insert(j);
        }
      }
      if (pcurrent.empty()) {
        past = null;
      }
      else {
        do {
          for(it=pcurrent.begin(); it!=pcurrent.end(); ++it) {
            j = *it;
            if (j == i) continue;
            for(jt=events[j].neighbours.begin(); jt!=events[j].neighbours.end(); ++jt) {
              k = *jt;
              S.clear();
              S.insert(j);
              S.insert(k);
              qt = index_table[1].find(S);
              rho = simplices[1][qt->second].orientation;
              if (rho == SYNARMOSMA::DISPARATE) continue;
              if (rho == SYNARMOSMA::BEFORE) {
                // Check to see if k doesn't already exist in the set "old"
                if (old.count(k) == 0) v.insert(k);
              }
            }
          }
          for(it=pcurrent.begin(); it!=pcurrent.end(); ++it) {
            old.insert(*it);
          }
          if (v.empty()) break;
          pcurrent = v;
          v.clear();
        } while(true);
        past = old;
        old.clear();
      }
      if (fcurrent.empty()) {
        future = null;
        continue;
      }
      do {
        for(it=fcurrent.begin(); it!=fcurrent.end(); ++it) {
          j = *it;
          if (j == i) continue;
          for(jt=events[j].neighbours.begin(); jt!=events[j].neighbours.end(); ++jt) {
            k = *jt;
            S.clear();
            S.insert(j);
            S.insert(k);
            qt = index_table[1].find(S);
            rho = simplices[1][qt->second].orientation;
            if (rho == SYNARMOSMA::DISPARATE) continue;
            if (rho == SYNARMOSMA::AFTER) {
              if (old.count(k) == 0) v.insert(k);
            }
          }
        }
        for(it=fcurrent.begin(); it!=fcurrent.end(); ++it) {
          old.insert(*it);
        }
        if (v.empty()) break;
        fcurrent = v;
        v.clear();
      } while(true);
      future = old;
      old.clear();
      if (past.count(i) == 1 || future.count(i) == 1) causal_loop++;
      // If it's a sink, that means all of its edges have orientation equal to -1;
      // for a source, the edges must all have an orientation equal to +1.
      if (past.empty() && !future.empty()) {
        compute_causal_graph(&G,i,-1);
        nlinearity += G.cyclicity();
        nsource++;
        continue;
      }
      else if (!past.empty() && future.empty()) {
        compute_causal_graph(&G,i,-1);
        nlinearity += G.cyclicity();
        nsink++;
        continue;
      }
    }
  }
  else {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i,j,k,it,jt,qt,rho,G,old,v,pcurrent,fcurrent,past,future,S) reduction(+:nlinearity,nsource,nsink,causal_loop)
#endif
    for(i=0; i<nv; ++i) {
      if (NTL::divide(events[i].ubiquity,codex[sheet].colour) == 0) continue;
      // Now calculate the future and past lightcones for this vertex on this sheet...
      pcurrent.clear();
      fcurrent.clear();
      for(it=events[i].neighbours.begin(); it!=events[i].neighbours.end(); ++it) {
        j = *it;
        S.clear();
        S.insert(i);
        S.insert(j);
        qt = index_table[1].find(S);
        if (NTL::divide(simplices[1][qt->second].ubiquity,codex[sheet].colour) == 0) continue;
        rho = simplices[1][qt->second].orientation;
        if (rho == SYNARMOSMA::DISPARATE) continue;
        if (rho == SYNARMOSMA::BEFORE) {
          pcurrent.insert(j);
        }
        else {
          fcurrent.insert(j);
        }
      }
      if (pcurrent.empty()) {
        past = null;
      }
      else {
        do {
          for(it=pcurrent.begin(); it!=pcurrent.end(); ++it) {
            j = *it;
            if (j == i) continue;
            for(jt=events[j].neighbours.begin(); jt!=events[j].neighbours.end(); ++jt) {
              k = *jt;
              S.clear();
              S.insert(j);
              S.insert(k);
              qt = index_table[1].find(S);
              if (NTL::divide(simplices[1][qt->second].ubiquity,codex[sheet].colour) == 0) continue;
              rho = simplices[1][qt->second].orientation;
              if (rho == SYNARMOSMA::DISPARATE) continue;
              if (rho == SYNARMOSMA::BEFORE) {
                // Check to see if k doesn't already exist in the set "old"
                if (old.count(k) == 0) v.insert(k);
              }
            }
          }
          for(it=pcurrent.begin(); it!=pcurrent.end(); ++it) {
            old.insert(*it);
          }
          if (v.empty()) break;
          pcurrent = v;
          v.clear();
        } while(true);
        past = old;
        old.clear();
      }
      if (fcurrent.empty()) {
        future = null;
        continue;
      }
      do {
        for(it=fcurrent.begin(); it!=fcurrent.end(); ++it) {
          j = *it;
          if (j == i) continue;
          for(jt=events[j].neighbours.begin(); jt!=events[j].neighbours.end(); ++jt) {
            k = *jt;
            S.clear();
            S.insert(j);
            S.insert(k);
            qt = index_table[1].find(S);
            if (NTL::divide(simplices[1][qt->second].ubiquity,codex[sheet].colour) == 0) continue;
            rho = simplices[1][qt->second].orientation;
            if (rho == SYNARMOSMA::DISPARATE) continue;
            if (rho == SYNARMOSMA::AFTER) {
              if (old.count(k) == 0) v.insert(k);
            }
          }
        }
        for(it=fcurrent.begin(); it!=fcurrent.end(); ++it) {
          old.insert(*it);
        }
        if (v.empty()) break;
        fcurrent = v;
        v.clear();
      } while(true);
      future = old;
      old.clear();
      if (past.count(i) == 1 || future.count(i) == 1) causal_loop++;
      // If it's a sink, that means all of its edges have orientation equal to -1;
      // for a source, the edges must all have an orientation equal to +1.
      if (past.empty() && !future.empty()) {
        compute_causal_graph(&G,i,sheet);
        nlinearity += G.cyclicity();
        nsource++;
        continue;
      }
      else if (!past.empty() && future.empty()) {
        compute_causal_graph(&G,i,sheet);
        nlinearity += G.cyclicity();
        nsink++;
        continue;
      }
    }
  }
  output = nlinearity/double(nsink + nsource) + double(causal_loop)/na;
  return output;
}

int Spacetime::chromatic_number(int sheet) const
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

int Spacetime::entourage(int base,int sheet) const
{
  // Calculates a measure of this event's integration/implication in its
  // spacetime neighbourhood
  int i,j,n,m,output = 0;

  for(i=1; i<=Spacetime::ND; ++i) {
    n = 0;
    m = (signed) simplices[i].size();
    for(j=0; j<m; ++j) {
      if (NTL::divide(simplices[i][j].ubiquity,codex[sheet].colour) == 0) continue;
      if (simplices[i][j].contains(base)) n++;
    }
    output += i*n;
  }
  return output;
}

int Spacetime::max_degree() const
{
  int i,n,output = 0;
  for(i=0; i<(signed) events.size(); ++i) {
    if (events[i].ubiquity == 1) continue;
    n = (signed) events[i].neighbours.size();
    if (output < n) output = n;
  }
  return output;
}

double Spacetime::entwinement(int sheet) const
{
  // This method produces a real number between 0 and 1 that measures the
  // degree of "labyrinthicity" of the graph
  SYNARMOSMA::Graph G;
  compute_graph(&G,sheet);
  return G.entwinement();
}

double Spacetime::cyclic_resistance(int sheet) const
{
  SYNARMOSMA::Graph G;
  compute_graph(&G,sheet);
  return G.cyclic_resistance();
}

int Spacetime::combinatorial_distance(int v1,int v2,int sheet) const
{
  // A method to calculate the topological distance
  // between the two vertices v1 and v2
  if (v1 == v2) return 0;

  // A sanity check...
#ifdef DEBUG
  if (sheet == -1) {
    assert(events[v1].ubiquity > 1 && events[v2].ubiquity > 1);
  }
  else {
    assert(NTL::divide(events[v1].ubiquity,codex[sheet].colour) == 1 && NTL::divide(events[v2].ubiquity,codex[sheet].colour) == 1);
  }
#endif

  int d,offset[events.size()];
  SYNARMOSMA::Graph G;
  compute_graph(&G,offset,sheet);
  d = G.distance(offset[v1],offset[v2]);
  return d;
}

int Spacetime::cardinality_safe(int d,int sheet) const
{
  int i,n = 0;
  if (sheet == -1) {
    if (d == 0) {
      const int M = (signed) events.size();
      for(i=0; i<M; ++i) {
        if (events[i].ubiquity == 1) continue;
        n++;
      }
    }
    else {
      const int M = (signed) simplices[d].size();
      for(i=0; i<M; ++i) {
        if (simplices[d][i].ubiquity == 1) continue;
        n++;
      }
    }
  }
  else {
    if (d == 0) {
      const int M = (signed) events.size();
      for(i=0; i<M; ++i) {
        if (NTL::divide(events[i].ubiquity,codex[sheet].colour) == 0) continue;
        n++;
      }
    }
    else {
      const int M = (signed) simplices[d].size();
      for(i=0; i<M; ++i) {
        if (NTL::divide(simplices[d][i].ubiquity,codex[sheet].colour) == 0) continue;
        n++;
      }
    }
  }
  return n;
}

int Spacetime::circuit_rank(int sheet) const
{
  int output = cardinality(1,sheet) - cardinality(0,sheet);
  if (connected(sheet)) {
    return 1 + output;
  }
  std::vector<int> components;
  int n = component_analysis(components,sheet);
  return n + output;
}

int Spacetime::euler_characteristic(int sheet) const
{
  int i,pf = 1,chi = 0;
  const int D = dimension(sheet);
  for(i=0; i<=D; ++i) {
    chi += pf*cardinality(i,sheet);
    pf *= -1;
  }
  return chi;
}

bool Spacetime::active_simplex(int d,int i,int sheet) const
{
  bool output = (NTL::divide(simplices[d][i].ubiquity,codex[sheet].colour) == 1) ? true : false;
  return output;
}

int Spacetime::dimension(int sheet) const
{
  int i,j;

  if (sheet == -1) {
    for(i=Spacetime::ND; i>0; i--) {
      for(j=0; j<(signed) simplices[i].size(); ++j) {
        if (simplices[i][j].ubiquity > 1) return i;
      }
    }
    for(i=0; i<(signed) events.size(); ++i) {
      if (events[i].ubiquity > 1) return 0;
    }
  }
  else {
    for(i=Spacetime::ND; i>0; i--) {
      for(j=0; j<(signed) simplices[i].size(); ++j) {
        if (NTL::divide(simplices[i][j].ubiquity,codex[sheet].colour) == 1) return i;
      }
    }
    for(i=0; i<(signed) events.size(); ++i) {
      if (NTL::divide(events[i].ubiquity,codex[sheet].colour) == 1) return 0;
    }
  }
  return -1;
}

int Spacetime::total_dimension(int sheet) const
{
  int i,sum = 0;

  if (sheet == -1) {
    for(i=0; i<(signed) events.size(); ++i) {
      if (events[i].ubiquity == 1) continue;
      sum += dimension(i);
    }
  }
  else {
    for(i=0; i<(signed) events.size(); ++i) {
      if (NTL::divide(events[i].ubiquity,codex[sheet].colour) == 1) sum += vertex_dimension(i,sheet);
    }
  }
  return sum;
}

int Spacetime::structural_index(int sheet) const
{
  int i,j,l,n,sum = 0,d = dimension(sheet);

  if (sheet == -1) {
    for(i=0; i<(signed) events.size(); ++i) {
      if (events[i].ubiquity == 1) continue;
      sum++;
    }
    for(i=1; i<d; ++i) {
      l = 0;
      n = (signed) simplices[i].size();
      for(j=0; j<n; ++j) {
        if (simplices[i][j].ubiquity == 1) continue;
        l++;
      }
      sum += (1+i)*l;
    }
  }
  else {
    for(i=0; i<(signed) events.size(); ++i) {
      if (NTL::divide(events[i].ubiquity,codex[sheet].colour) == 1) sum++;
    }
    for(i=1; i<d; ++i) {
      l = 0;
      n = (signed) simplices[i].size();
      for(j=0; j<n; ++j) {
        if (NTL::divide(simplices[i][j].ubiquity,codex[sheet].colour) == 1) l++;
      }
      sum += (1+i)*l;
    }
  }
  return sum;
}

int Spacetime::weighted_entourage(int n1,int n2) const
{
  int i,j,nfound,output = 0;
  bool f1,f2;
  std::set<int>::const_iterator it;

  for(i=2; i<=Spacetime::ND; ++i) {
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

int Spacetime::vertex_valence(int v,int sheet) const
{
  int nd = 0;
  std::set<int>::const_iterator it;
  for(it=events[v].entourage.begin(); it!=events[v].entourage.end(); ++it) {
    if (NTL::divide(simplices[1][*it].ubiquity,codex[sheet].colour) == 1) nd++;
  }
  return nd;
}

int Spacetime::vertex_dimension(int v,int sheet) const
{
  int i,j,n;
  if (v < 0 || v >= (signed) events.size()) return -1;
  if (sheet == -1) {
    for(i=Spacetime::ND; i>=1; i--) {
      n = (signed) simplices[i].size();
      for(j=0; j<n; ++j) {
        if (simplices[i][j].ubiquity == 1) continue;
        if (simplices[i][j].contains(v)) return i;
      }
    }
    if (events[v].ubiquity == 1) return -1;
    return 0;
  }
  else {
    for(i=Spacetime::ND; i>=1; i--) {
      n = (signed) simplices[i].size();
      for(j=0; j<n; ++j) {
        if (NTL::divide(simplices[i][j].ubiquity,codex[sheet].colour) == 0) continue;
        if (simplices[i][j].contains(v)) return i;
      }
    }
    if (NTL::divide(events[v].ubiquity,codex[sheet].colour) == 1) return 0;
  }
  return -1;
}

double Spacetime::dimensional_frontier(int D,int sheet) const
{
  int i,d[2],s = 0;
  for(i=0; i<(signed) simplices[1].size(); ++i) {
    if (NTL::divide(simplices[1][i].ubiquity,codex[sheet].colour) == 0) continue;
    simplices[1][i].get_vertices(d);
    d[0] = (d[0] < D) ? D : d[0];
    d[1] = (d[1] < D) ? D : d[1];
    if (d[0] != d[1]) s++;
  }
  return double(s)/double(simplices[1].size());
}

double Spacetime::dimensional_uniformity(int sheet) const
{
  int i,n,nv,sdimension = dimension(sheet),sum = 0;
  const int D = geometry->dimension();
  if (sdimension < D) sdimension = D;

  nv = 0;
  if (sheet == -1) {
    for(i=0; i<(signed) events.size(); ++i) {
      if (events[i].ubiquity == 1) continue;
      n = vertex_dimension(i,sheet);
      n = (n < D) ? D : n;
      sum += sdimension - n;
      nv++;
    }
  }
  else {
    for(i=0; i<(signed) events.size(); ++i) {
      if (NTL::divide(events[i].ubiquity,codex[sheet].colour) == 0) continue;
      n = vertex_dimension(i,sheet);
      n = (n < D) ? D : n;
      sum += sdimension - n;
      nv++;
    }
  }
  return double(sum)/double(nv);
}

bool Spacetime::consistent(int sheet) const
{
  int i,j,k,l,n,m,vx[2];
  bool found;
  std::set<int>::const_iterator it;
  std::set<int> S;
  SYNARMOSMA::hash_map::const_iterator qt;
  const int nv = (signed) events.size();
  const int ulimit = dimension(sheet);

#ifdef DEBUG
  assert(simplices[0].empty());
  assert(index_table[0].empty());
#endif

  if (sheet == -1) {
    for(i=Spacetime::ND; i>=2; i--) {
      n = (signed) simplices[i].size();
      for(j=0; j<n; ++j) {
        if (simplices[i][j].ubiquity == 1) continue;
        for(it=simplices[i][j].entourage.begin(); it!=simplices[i][j].entourage.end(); ++it) {
          if (simplices[i+1][*it].ubiquity == 1) {
            std::cout << "Error with entourage ubiquity: " << i << "  " << j << "  " << *it << "  " << ulimit << std::endl;
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
            if (simplices[i-1][l].ubiquity == 1) continue;
            if (S == simplices[i-1][l].vertices) {
              found = true;
              break;
            }
          }
          if (!found) {
            std::cout << "Can't find face " << SYNARMOSMA::make_key(S) << " of simplex " << SYNARMOSMA::make_key(simplices[i][j].vertices) << std::endl;
            for(l=0; l<(signed) simplices[i-1].size(); ++l) {
              if (simplices[i-1][l].ubiquity == 1) continue;
              std::cout << SYNARMOSMA::make_key(simplices[i-1][l].vertices) << "  " << simplices[i-1][l].ubiquity << std::endl;
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
          if (events[l].ubiquity == 1) {
            std::cout << "Inactive vertex: " << l << std::endl;
            return false;
          }
        }
      }
    }
    for(i=0; i<(signed) simplices[1].size(); ++i) {
      if (simplices[1][i].ubiquity == 1) continue;
      simplices[1][i].get_vertices(vx);
      if (vx[0] < 0 || vx[0] >= nv) return false;
      if (vx[1] < 0 || vx[0] >= nv) return false;
      if (events[vx[0]].ubiquity == 1) {
        std::cout << "Inactive vertex: " << vx[0] << std::endl;
        return false;
      }
      if (events[vx[1]].ubiquity == 1) {
        std::cout << "Inactive vertex: " << vx[1] << std::endl;
        return false;
      }
    }
    for(i=0; i<nv; ++i) {
      if (events[i].ubiquity == 1) continue;
      for(it=events[i].entourage.begin(); it!=events[i].entourage.end(); ++it) {
        if (simplices[1][*it].ubiquity == 1) {
          std::cout << "Error with entourage ubiquity: " << i << "  " << *it << std::endl;
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
      }
    }
  }
  else {
    for(i=Spacetime::ND; i>=2; i--) {
      n = (signed) simplices[i].size();
      for(j=0; j<n; ++j) {
        if (NTL::divide(simplices[i][j].ubiquity,codex[sheet].colour) == 0) continue;
        if (simplices[i][j].dimension() != i) {
          std::cout << i << "-simplex  " << j << " has dimension " << simplices[i][j].dimension() << std::endl;
          return false;
        }
        for(k=0; k<1+i; ++k) {
          S = simplices[i][j].faces[k];
          found = false;
          for(l=0; l<(signed) simplices[i-1].size(); ++l) {
            if (NTL::divide(simplices[i-1][l].ubiquity,codex[sheet].colour) == 0) continue;
            if (S == simplices[i-1][l].vertices) {
              found = true;
              break;
            }
          }
          if (!found) {
            std::cout << "For sheet " << sheet << ", can't find face " << SYNARMOSMA::make_key(S) << " of simplex " << SYNARMOSMA::make_key(simplices[i][j].vertices) << std::endl;
            for(l=0; l<(signed) simplices[i-1].size(); ++l) {
              if (NTL::divide(simplices[i-1][l].ubiquity,codex[sheet].colour) == 0) continue;
              std::cout << SYNARMOSMA::make_key(simplices[i-1][l].vertices) << "  " << simplices[i-1][l].ubiquity << std::endl;
            }
            return false;
          }
        }
        for(it=simplices[i][j].vertices.begin(); it!=simplices[i][j].vertices.end(); ++it) {
          l = *it;
          if (l < 0 || l >= nv) {
            std::cout << "Nonexistent vertex: " << l << std::endl;
            return false;
          }
          if (NTL::divide(events[l].ubiquity,codex[sheet].colour) == 0) {
            std::cout << "Inactive vertex: " << l << std::endl;
            return false;
          }
        }
      }
    }
    for(i=0; i<(signed) simplices[1].size(); ++i) {
      if (NTL::divide(simplices[1][i].ubiquity,codex[sheet].colour) == 0) continue;
      simplices[1][i].get_vertices(vx);
      if (vx[0] < 0 || vx[0] >= nv) return false;
      if (vx[1] < 0 || vx[0] >= nv) return false;
      if (NTL::divide(events[vx[0]].ubiquity,codex[sheet].colour) == 0) {
        std::cout << "Inactive vertex: " << vx[0] << std::endl;
        return false;
      }
      if (NTL::divide(events[vx[1]].ubiquity,codex[sheet].colour) == 0) {
        std::cout << "Inactive vertex: " << vx[1] << std::endl;
        return false;
      }
    }
    for(i=0; i<nv; ++i) {
      if (NTL::divide(events[i].ubiquity,codex[sheet].colour) == 0) continue;
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
      }
    }
  }
  // Make sure that each n-simplex only exists once...
  for(i=1; i<=Spacetime::ND; ++i) {
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

bool Spacetime::connected(int sheet) const
{
  // For this method we calculate the 1-skeleton of the simplicial complex and then,
  // as a graph, determine its connectedness.
  SYNARMOSMA::Graph G;

  compute_graph(&G,sheet);

  return G.connected();
}

int Spacetime::component_analysis(std::vector<int>& component,int sheet) const
{
  int i,n,ct = 0;
  std::vector<int> cvalue;
  SYNARMOSMA::Graph G;

  compute_graph(&G,sheet);
  n = G.component_analysis(cvalue);
  // We now need to include the unused vertices in the output component vector
  component.clear();
  if (sheet == -1) {
    for(i=0; i<(signed) events.size(); ++i) {
      if (events[i].ubiquity == 1) {
        component.push_back(-1);
        continue;
      }
      component.push_back(cvalue[ct]);
      ct++;
    }
  }
  else {
    for(i=0; i<(signed) events.size(); ++i) {
      if (NTL::divide(events[i].ubiquity,codex[sheet].colour) == 0) {
        component.push_back(-1);
        continue;
      }
      component.push_back(cvalue[ct]);
      ct++;
    }
  }
  return n;
}

