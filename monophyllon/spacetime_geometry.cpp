#include "spacetime.h"

using namespace DIAPLEXIS;

void Spacetime::get_coordinates(int v,std::vector<double>& x) const
{
  geometry->get_coordinates(v,x);
}

void Spacetime::get_coordinates(std::vector<double>& x) const
{
  int i;
  unsigned int j;
  std::vector<double> vx;
  const int nv = skeleton->cardinality(0);

  x.clear();
  for(i=0; i<nv; ++i) {
    if (!skeleton->active_event(i)) continue;
    get_coordinates(i,vx);
    for(j=0; j<geometry->dimension(); ++j) {
      x.push_back(vx[j]);
    }
  }
}

void Spacetime::arclength_statistics(double* output) const
{
  int i;
  double wm,avg_length = 0.0,max_length = 0.0,min_length = 1000.0;
  const int Ne = (signed) skeleton->simplices[1].size();
  const double ne = double(skeleton->cardinality(1));

  for(i=0; i<Ne; ++i) {
    if (!skeleton->active_simplex(1,i)) continue;
    wm = skeleton->simplices[1][i].get_volume();
    avg_length += wm;
    if (wm > max_length) max_length = wm;
    if (wm < min_length) min_length = wm;
  }
  if (ne > 0) avg_length /= ne;
  output[0] = max_length;
  output[1] = min_length;
  output[2] = avg_length;
}

bool Spacetime::adjust_dimension()
{
  // This method handles the geometry and energy changes to
  // minimize the structural deficiency...
  int i,n;
  std::vector<int> vdimension;
  const int nv = (signed) skeleton->events.size();
  const int D = (signed) geometry->dimension();
  const bool uniform = geometry->get_uniform();

  system_size = 0;

  if (uniform) {
    for(i=0; i<nv; ++i) {
      if (!skeleton->active_event(i)) {
        vdimension.push_back(-1);
        continue;
      }
      n = skeleton->events[i].get_topological_dimension();
      system_size += D;
      vdimension.push_back(n);
    }
  }
  else {
    for(i=0; i<nv; ++i) {
      if (!skeleton->active_event(i)) {
        vdimension.push_back(-1);
        continue;
      }
      n = skeleton->events[i].get_topological_dimension();
      system_size += (n <= D) ? D : n;
      vdimension.push_back(n);
    }
  }
  return geometry->adjust_dimension(vdimension);
}

void Spacetime::compute_causal_graph(SYNARMOSMA::Directed_Graph* G,int base) const
{
  int i,j,v;
  std::set<int> S,N;
  std::set<int>::const_iterator it;
  std::vector<int> offset,current,next;
  std::vector<int>::const_iterator v_it;
  SYNARMOSMA::hash_map::const_iterator qt;
  const int nv = (signed) skeleton->events.size();

  G->clear();
  current.push_back(base);
  for(i=0; i<nv; ++i) {
    offset.push_back(-1);
  }
  offset[base] = G->add_vertex();

  do {
    for(v_it=current.begin(); v_it!=current.end(); ++v_it) {
      v = *v_it;
      skeleton->events[v].get_neighbours(N);
      for(it=N.begin(); it!=N.end(); ++it) {
        j = *it;
        S.clear();
        S.insert(v);
        S.insert(j);
        qt = skeleton->index_table[1].find(S);
        if (!skeleton->active_simplex(1,qt->second)) continue;
        if (!skeleton->simplices[1][qt->second].timelike()) continue;
        if (offset[j] == -1) {
          offset[j] = G->add_vertex();
          next.push_back(j);
        }
        G->add_edge(offset[v],offset[j],geometry->get_temporal_order(v,j));
      }
    }
    if (next.empty()) break;
    current = next;
    next.clear();
  } while(true);
}

void Spacetime::compute_total_lightcone(int v,std::set<int>& past_cone,std::set<int>& future_cone) const
{
  // This method assumes that the compute_lightcones method has already been called to fill 
  // the anterior and posterior properties of the events!
  int i,j;
  std::set<int> S,current,next;
  std::set<int>::const_iterator it,jt;

  // First the past cone...
  current.insert(v);
  past_cone.clear();
  do {
    for(it=current.begin(); it!=current.end(); ++it) {
      i = *it;
      skeleton->events[i].get_anterior(S);
      for(jt=S.begin(); jt!=S.end(); ++jt) {
        j = *jt;
        if (past_cone.count(j) > 0) continue;
        next.insert(j);
      }
    }
    if (next.empty()) break;
    for(it=next.begin(); it!=next.end(); ++it) {
      past_cone.insert(*it);
    }
    current = next;
    next.clear();
  } while(true);

  // Now the future cone...
  current.clear();
  current.insert(v);
  future_cone.clear();
  do {
    for(it=current.begin(); it!=current.end(); ++it) {
      i = *it;
      skeleton->events[i].get_posterior(S);
      for(jt=S.begin(); jt!=S.end(); ++jt) {
        j = *jt;
        if (future_cone.count(j) > 0) continue;
        next.insert(j);
      }
    }
    if (next.empty()) break;
    for(it=next.begin(); it!=next.end(); ++it) {
      future_cone.insert(*it);
    }
    current = next;
    next.clear();
  } while(true);
}

double Spacetime::compute_temporal_nonlinearity() const
{
  int i,nsink = 0,nsource = 0,causal_loop = 0;
  double output,nlinearity = 0.0;
  std::set<int> past,future;
  SYNARMOSMA::Directed_Graph G;
  const int nv = (signed) skeleton->events.size();
  const double na = double(skeleton->cardinality(0));

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i,G,past,future) reduction(+:nlinearity,nsource,nsink,causal_loop)
#endif
  for(i=0; i<nv; ++i) {
    if (!skeleton->active_event(i)) continue;
    // Now calculate the future and past lightcones for this vertex on this sheet...
    compute_total_lightcone(i,past,future);
    if (past.count(i) == 1 || future.count(i) == 1) causal_loop++;
    // If it's a sink, that means all of its edges have orientation equal to -1;
    // for a source, the edges must all have an orientation equal to +1.
    if (past.empty() && !future.empty()) {
      compute_causal_graph(&G,i);
      nlinearity += G.cyclicity();
      nsource++;
    }
    else if (!past.empty() && future.empty()) {
      compute_causal_graph(&G,i);
      nlinearity += G.cyclicity();
      nsink++;
    }
  }
  output = nlinearity/double(nsink + nsource) + double(causal_loop)/na;
  return output;
}

void Spacetime::compute_lightcones()
{
  int i,vx[2];
  const int n = (signed) skeleton->events.size();
  const int m = (signed) skeleton->simplices[1].size();

  for(i=0; i<n; ++i) {
    skeleton->events[i].clear_anterior();
    skeleton->events[i].clear_posterior();
  }
  for(i=0; i<m; ++i) {
    if (!skeleton->active_simplex(1,i)) continue;
    if (!skeleton->simplices[1][i].timelike()) continue;
    skeleton->simplices[1][i].get_vertices(vx);
    if (geometry->get_temporal_order(vx[0],vx[1]) == SYNARMOSMA::Relation::before) {
      skeleton->events[vx[0]].add_posterior(vx[1]);
      skeleton->events[vx[1]].add_anterior(vx[0]);
    }
    else {
      skeleton->events[vx[1]].add_posterior(vx[0]);
      skeleton->events[vx[0]].add_anterior(vx[1]);
    }
  }
}

void Spacetime::chorogenesis(int nsteps)
{
  assert(solver == Geometry_Solver::mechanical);
  assert(edge_probability > 0.3);
  assert(geometry->get_euclidean());
  assert(skeleton->connected());
  int dpopulation[1+geometry->dimension()];
  unsigned int i,j;
  std::vector<std::pair<long,int> > factors;
  const unsigned int D = geometry->dimension();
  const unsigned int nv = skeleton->events.size();

  SYNARMOSMA::factorize(nv,factors);
  j = 1;
  for(i=0; i<factors.size(); ++i) {
    assert(factors[i].second%D == 0);
    j *= SYNARMOSMA::ipow(factors[i].first,factors[i].second/D);
  }
  const int n = j;
  for(i=0; i<=D; ++i) {
    dpopulation[i] = SYNARMOSMA::ipow(2,i)*int(SYNARMOSMA::binomial(D,i))*SYNARMOSMA::ipow(n - 2,D - i);
  }
  // Zero out the spacetime energy and make sure the vertex geometry is dimensionally homogeneous...
  std::vector<double> x,y;
  system_size = 0;
  for(i=0; i<nv; ++i) {
    skeleton->events[i].nullify_energy();
    geometry->get_coordinates(i,x);
    system_size += D;
    if (D == x.size()) continue;
    for(j=0; j<D; ++j) {
      y.push_back(x[j]);
    }
    geometry->set_coordinates(i,y);
    y.clear();
  }
  // We have now verified all the necessary conditions and can begin  
  skeleton->compute_degree_distribution(false);
  
  skeleton->compute_connectivity_distribution(!geometry->get_memory_type());

  int vx[2];
  unsigned int d,cfactor,derror,csize;
  double temperature = 1.0;
  std::vector<int> candidates,reorder;
  std::vector<double> hgram;
  SYNARMOSMA::Graph* G = new SYNARMOSMA::Graph;
  const unsigned int ne = skeleton->simplices[1].size();

#ifdef VERBOSE
  std::cout << "Beginning chorogenesis with " << ne << " edges..." << std::endl;
#endif  

  do {
    iterations += 1;
    for(i=0; i<ne; ++i) {
      if (!skeleton->active_simplex(1,i)) continue;
      skeleton->simplices[1][i].get_vertices(vx);
      d = skeleton->vertex_valence(vx[0]);
      if (d <= 2*geometry->dimension()) continue;
      d = skeleton->vertex_valence(vx[1]);
      if (d <= 2*geometry->dimension()) continue;
      candidates.push_back(i);
    }
    if (candidates.empty()) {
      std::cout << "No viable candidates for topological operations, exiting loop..." << std::endl;
      break;
    }
    csize = candidates.size();
    cfactor = int(skeleton->RND->drandom(0.1*temperature,temperature)*csize);
    if (cfactor == 0) {
      std::cout << "Temperature too low for topological operations, exiting loop..." << std::endl;
      break;
    }
    cfactor = std::max(cfactor,2*nv);
#ifdef VERBOSE
    std::cout << "There are " << csize << " and " << cfactor << " edges to delete at temperature " << temperature << std::endl;
#endif
    d = 0;
    skeleton->RND->shuffle(reorder,csize);
    for(i=0; i<csize; ++i) {
      j = candidates[reorder[i]];
      skeleton->simplex_deletion(1,j);
      if (!skeleton->connected()) {
        skeleton->simplices[1][j].activate();
        continue;
      }
      d++;
      if (d == cfactor) break;
    }
#ifdef VERBOSE
    std::cout << "Deleted " << d << " edges from the spacetime complex." << std::endl;
#endif
    regularization(false);
    skeleton->compute_graph(G);
    G->degree_distribution(false,hgram);
    derror = 0;
    for(i=0; i<hgram.size(); ++i) {
      d = int(nv*hgram[i]);
      if (i > D) {
        derror += d;
      }
      else {
        derror += std::abs(d - dpopulation[i]);
      }
    }
#ifdef VERBOSE
    std::cout << "Vertex degree error is " << derror << std::endl;
#endif
    skeleton->compute_connectivity_distribution(!geometry->get_memory_type());
    // Now the geometry...
    optimize();
    // Prepare for the next iteration...
    temperature *= 0.95; 
    candidates.clear();    
  } while(true);
  if (iterations == 1) optimize();
  compute_volume();
  compute_obliquity();
  structural_deficiency();
  write_state();
  delete G;
}

double Spacetime::compute_abnormality(const std::vector<int>& flexible_edge) const
{
  int i,vx[2];
  double d,ell,output = 0.0;
  const int ne = (signed) skeleton->simplices[1].size();
  const double pfactor = (2.0/M_PI)*5.0;
  const double sq_tolerance = edge_flexibility_threshold*edge_flexibility_threshold;

  for(i=0; i<ne; ++i) {
    if (!skeleton->active_simplex(1,i)) continue;
    skeleton->simplices[1][i].get_vertices(vx);
    d = geometry->get_squared_distance(vx[0],vx[1],false);
    if (flexible_edge[i] == 1) {
      if (d > sq_tolerance) {
        d = std::sqrt(d);
        output += (d - edge_flexibility_threshold)*(d - edge_flexibility_threshold);
      }
      continue;
    }
    d = std::sqrt(d);
    ell = 1.0/(1.0 + pfactor*std::atan(0.5*(skeleton->events[vx[0]].get_energy() + skeleton->events[vx[1]].get_energy())));
    output += (d - ell)*(d - ell);
  }
  return output;
}

double Spacetime::compute_abnormality(const std::vector<double>& x,const std::vector<int>& flexible_edge) const
{
  int i,vx[2];
  double d,ell,output = 0.0;
  const int ne = (signed) skeleton->simplices[1].size();
  const double pfactor = (2.0/M_PI)*5.0;
  const double sq_tolerance = edge_flexibility_threshold*edge_flexibility_threshold;

  for(i=0; i<system_size; ++i) {
    geometry->set_element(i,x[i]);
  }
  geometry->compute_squared_distances();

  for(i=0; i<ne; ++i) {
    if (!skeleton->active_simplex(1,i)) continue;
    skeleton->simplices[1][i].get_vertices(vx);
    d = geometry->get_squared_distance(vx[0],vx[1],false);
    if (flexible_edge[i] == 1) {
      if (d > sq_tolerance) {
        d = std::sqrt(d);
        output += (d - edge_flexibility_threshold)*(d - edge_flexibility_threshold);
      }
      continue;
    }
    d = std::sqrt(d);
    ell = 1.0/(1.0 + pfactor*std::atan(0.5*(skeleton->events[vx[0]].get_energy() + skeleton->events[vx[1]].get_energy())));
    output += (d - ell)*(d - ell);
  }
  return output;
}

void Spacetime::compute_geometric_gradient(std::vector<double>& df,bool negate,const std::vector<int>& flexible_edge)
{
  int i;
  std::set<int>::const_iterator it;

  df.clear();
  for(i=0; i<system_size; ++i) {
    df.push_back(0.0);
  }
  if (geometry->get_relational()) {
    std::set<int> modified_vertices,last,current;
    double E,alpha;
    const double B = compute_abnormality(flexible_edge);

    for(i=0; i<system_size; ++i) {
      geometry->add(i,geometry_tolerance);
      geometry->get_implied_vertices(i,current);
      skeleton->compute_dependent_simplices(current);
      modified_vertices = last;
      for(it=current.begin(); it!=current.end(); ++it) {
        modified_vertices.insert(*it);
      }
      geometry->compute_squared_distances(modified_vertices);
      compute_volume();
      E = compute_abnormality(flexible_edge);
      alpha = (E - B)/geometry_tolerance;
      df[i] = alpha;
      geometry->add(i,-geometry_tolerance);
      last = current;
      current.clear();
    }
    geometry->compute_squared_distances(last);
  }
  else {
    int j,k,na = 0;
    double l,ell;
    std::set<int> S,N;
    std::vector<double> x1,x2;
    SYNARMOSMA::hash_map::const_iterator qt;
    const int nv = (signed) skeleton->events.size();
    const int D = geometry->dimension();
    const double pfactor = (2.0/M_PI)*5.0;
    const double sq_tolerance = edge_flexibility_threshold*edge_flexibility_threshold;
    double alpha[D];

    for(i=0; i<nv; ++i) {
      if (!skeleton->active_event(i)) continue;
      na++;
    }
#ifdef DEBUG
    assert(system_size == D*na);
#endif
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i,j,k,l,S,N,x1,x2,ell,alpha,it,qt)
#endif
    for(i=0; i<nv; ++i) {
      if (!skeleton->active_event(i)) continue;
      geometry->get_coordinates(i,x1);
      for(j=0; j<D; ++j) {
        alpha[j] = 0.0;
      }
      skeleton->events[i].get_neighbours(N);
      for(it=N.begin(); it!=N.end(); ++it) {
        k = *it;
        geometry->get_coordinates(k,x2);
        l = geometry->get_squared_distance(i,k,false);
        S.clear();
        S.insert(i);
        S.insert(k);
        qt = skeleton->index_table[1].find(S);
        if (flexible_edge[qt->second] == 1) {
          if (l > sq_tolerance) {
            l = std::sqrt(l);
            for(j=0; j<D; ++j) {
              alpha[j] += 2.0*(l - edge_flexibility_threshold)*(x1[j] - x2[j])/l;
            }
          }
        }
        else {
          ell = 1.0/(1.0 + pfactor*std::atan(0.5*(skeleton->events[i].get_energy() + skeleton->events[k].get_energy())));
          l = std::sqrt(l);
          for(j=0; j<D; ++j) {
            alpha[j] += 2.0*(l - ell)*(x1[j] - x2[j])/l;
          }
        }
      }
      for(j=0; j<D; ++j) {
        df[D*i+j] = alpha[j];
      }
    }
  }
  if (negate) {
    for(i=0; i<system_size; ++i) {
      df[i] = -df[i];
    }
  }
}

double Spacetime::minimize_lengths(const std::vector<int>& S1,const std::vector<int>& S2,int* vx) const
{
  int v1,v2;
  unsigned int i,j;
  double delta,mdelta = std::numeric_limits<double>::infinity();
  const unsigned int n1 = S1.size();
  const unsigned int n2 = S2.size();

  for(i=0; i<n1; ++i) {
    v1 = S1[i];
    for(j=0; j<n2; ++j) {
      v2 = S2[j];
      delta = std::abs(geometry->get_squared_distance(v1,v2,false));
      if (delta < mdelta) {
        mdelta = delta;
        vx[0] = v1;
        vx[1] = v2;
      }
    }
  }
  return mdelta;
}

double Spacetime::compute_temporal_vorticity(int v) const
{
  // A value near zero means that the current topology is in accord with the chronogeometry
  // whereas a large positive value means much more topological entwinement is needed and a
  // large negative value means the topology must be altered to resemble that of the Cartesian
  // initial state.
  int u,w,tcount;
  double l,tipsy,vorticity;
  std::set<int> S,N;
  std::vector<int> jset;
  std::set<int>::const_iterator it,jt;
  std::vector<int>::const_iterator vit;
  SYNARMOSMA::Relation d1,d2;
  SYNARMOSMA::Directed_Graph G;
  SYNARMOSMA::hash_map::const_iterator qt;

  compute_causal_graph(&G,v);
  vorticity = G.cyclicity();

  tipsy = 0.0;
  skeleton->events[v].get_neighbours(N);
  for(it=N.begin(); it!=N.end(); ++it) {
    u = *it;
    S.clear();
    S.insert(v);
    S.insert(u);
    qt = skeleton->index_table[1].find(S);
    l = skeleton->simplices[1][qt->second].get_volume();
    if (!skeleton->simplices[1][qt->second].spacelike()) continue;
    // This edge is spacelike, so it will contribute to the temporal vorticity
    jset.clear();
    for(jt=N.begin(); jt!=N.end(); ++jt) {
      w = *jt;
      if (skeleton->events[u].is_neighbour(w)) jset.push_back(w);
    }
    if (jset.empty()) continue;
    tcount = 0;
    for(vit=jset.begin(); vit!=jset.end(); ++vit) {
      w = *vit;
      if (w == u || w == v) continue;
      d1 = geometry->get_temporal_order(v,w);
      d2 = geometry->get_temporal_order(u,w);
      if (d1 != d2) tcount++;
    }
    tipsy += double(tcount)/(1.0 + l*double(jset.size()));
  }
  vorticity += tipsy;
 
  return vorticity;
}

void Spacetime::compute_obliquity()
{
  const int nv = (signed) skeleton->events.size();
  if (geometry->get_relational()) {
    for(int i=0; i<nv; ++i) {
      skeleton->events[i].set_obliquity(0.0);
    }
    return;
  }
  int i,j;
  double rho,theta,alpha;
  bool first;
  std::vector<double> vx,vy;
  std::set<int> N;
  std::set<int>::const_iterator it;

  const double A = 2.5;

  for(i=0; i<nv; ++i) {
    if (!skeleton->active_event(i)) continue;

    if (skeleton->vertex_valence(i) < 2 || !skeleton->events[i].get_geometry_modified()) continue;

    skeleton->events[i].get_neighbours(N);
    j = *(N.begin());

    geometry->vertex_difference(i,j,vx);

    rho = 0.0;
    first = true;
    for(it=N.begin(); it!=N.end(); ++it) {
      if (first) {
        first = false;
        continue;
      }
      j = *it;
      geometry->vertex_difference(i,j,vy);
      alpha = geometry->get_argument(vx,vy);
      theta = std::acos(alpha);
      rho += A*std::sin(2.0*theta)*std::sin(2.0*theta);
    }
    rho = rho/double(skeleton->vertex_valence(i) - 1);
    if (rho < std::numeric_limits<double>::epsilon()) rho = 0.0;
    skeleton->events[i].set_obliquity(rho);
  }
}

double Spacetime::representational_energy(bool weighted) const
{
  // A routine that follows the algorithm outlined in C. Godsil and G. Royle, "Algebraic
  // Graph Theory" (Springer, 2001), Section 13.3 (pp. 284--286)
  if (skeleton->dimension() < 1) return 0.0;
  int i,j,k,l,nelements,n = 0;
  double w,E = 0.0;
  std::vector<double> xc,Lx;
  std::vector<unsigned int> Lc;
  std::vector<int> offset;
  SYNARMOSMA::Graph G;
  const int nv = (signed) skeleton->events.size();
  const int N = skeleton->cardinality(0);
  const int D = geometry->dimension();  
  double coordinates[N][D],result[N][D],sum[D];
  SYNARMOSMA::Matrix<double> L(N);

  skeleton->compute_graph(&G);
  
  for(i=0; i<nv; ++i) {
    if (!skeleton->active_event(i)) {
      offset.push_back(-1);
      continue;
    }
    offset.push_back(n);
    geometry->get_coordinates(i,xc);
    for(j=0; j<D; ++j) {
      coordinates[n][j] = xc[j];
    }
    n++;
  }

  if (weighted) {
    int vx[2];
    const int ne = (signed) skeleton->simplices[1].size();
  
    for(i=0; i<ne; ++i) {
      if (!skeleton->active_simplex(1,i)) continue;
      w = 2.0/std::abs(skeleton->simplices[1][i].get_volume());
      skeleton->simplices[1][i].get_vertices(vx);
      j = offset[vx[0]];
      k = offset[vx[1]];
      // The diagonal elements need to be incremented...
      L.set(j,j,w,true);
      L.set(k,k,w,true);
      // but not the off-diagonal elements...
      L.set(j,k,-w);
      L.set(k,j,-w);
    }
  }
  else {
    G.compute_laplacian(&L);
  }

  for(i=0; i<N; ++i) {
    for(j=0; j<D; ++j) {
      sum[j] = 0.0;
    }
    L.get_row(Lx,Lc,i);
    nelements = (signed) Lx.size();
    for(j=0; j<nelements; ++j) {
      for(k=0; k<D; ++k) {
        sum[k] += coordinates[Lc[j]][k]*Lx[j];
      }
    }
    for(j=0; j<D; ++j) {
      result[i][j] = sum[j];
    }
  }

  for(i=0; i<D; ++i) {
    w = 0.0;
    for(l=0; l<N; ++l) {
      w += result[l][i]*coordinates[l][i];
    }
    E += w;
  }
  return E;
}

bool Spacetime::realizable(int d,int n) const
{
  // Can the d*(d+1)/2 edge lengths lead to a geometrically
  // realizable d-simplex?
  // An edge or vertex is always geometrically realizable...
  if (d < 2) return true;
  int i,j,info,dp1 = d + 1;
  SYNARMOSMA::hash_map::const_iterator qt;
  double alpha;
  bool output = true;
  char jtype = 'N';
  char uplo = 'U';
  int nwork = 3*dp1 - 1;
  int vx[dp1];
  std::set<int> S;
  double* A = new double[dp1*dp1];
  double* w = new double[dp1];
  double* work = new double[nwork];

  skeleton->simplices[d][n].get_vertices(vx);

  for(i=0; i<dp1; ++i) {
    alpha = 0.0;
    if (i != 0) {
      S.clear();
      S.insert(vx[i]);
      S.insert(vx[0]);
      qt = skeleton->index_table[1].find(S);
      alpha = skeleton->simplices[1][qt->second].get_volume();
    }
    for(j=0; j<dp1; ++j) {
      if (j != 0) {
        S.clear();
        S.insert(vx[j]);
        S.insert(vx[0]);
        qt = skeleton->index_table[1].find(S);
        alpha += skeleton->simplices[1][qt->second].get_volume();
      }
      S.clear();
      S.insert(vx[i]);
      S.insert(vx[j]);
      qt = skeleton->index_table[1].find(S);
      alpha -= skeleton->simplices[1][qt->second].get_volume();
      A[dp1*i+j] = alpha;
    }
  }
  dsyev_(&jtype,&uplo,&dp1,A,&dp1,w,work,&nwork,&info);
#ifdef DEBUG
  assert(info == 0);
#endif
  for(i=0; i<dp1; ++i) {
    if (w[i] < 0.0) {
      output = false;
      break;
    }
  }
  delete[] A;
  delete[] w;
  delete[] work;
  return output;
}

void Spacetime::compute_volume()
{
  int i,j,k,l,n,m,vx[Complex::ND+3];
  double prefactor,V,l1,l2,l3;
  std::set<int> S;
  std::vector<std::set<int> > F;
  SYNARMOSMA::UINT64 q,p = 8;
  SYNARMOSMA::Matrix<double> A;
  SYNARMOSMA::hash_map::const_iterator qt;

  compute_lengths();

  for(i=0; i<(signed) skeleton->simplices[2].size(); ++i) {
    // The triangles are a very simple case we can handle
    // without LAPACK thanks to Heron's formula
    if (!skeleton->active_simplex(2,i) || !skeleton->simplices[2][i].get_modified()) continue;
    skeleton->simplices[2][i].get_faces(F);
    qt = skeleton->index_table[1].find(F[0]);
    l1 = skeleton->simplices[1][qt->second].get_squared_volume();
    qt = skeleton->index_table[1].find(F[1]);
    l2 = skeleton->simplices[1][qt->second].get_squared_volume();
    qt = skeleton->index_table[1].find(F[2]);
    l3 = skeleton->simplices[1][qt->second].get_squared_volume();
    V = -(l3*l3 - 2.0*l3*(l1 + l2) + (l2 - l1)*(l2 - l1))/16.0;
    skeleton->simplices[2][i].set_squared_volume(V);
    skeleton->simplices[2][i].set_volume(std::sqrt(std::abs(V)));
    skeleton->simplices[2][i].set_modified(false);
  }

  for(i=3; i<=Complex::ND; ++i) {
    if (skeleton->simplices[i].empty()) continue;
    m = i + 2;
    n = (signed) skeleton->simplices[i].size();
    q = SYNARMOSMA::factorial(i);
    q *= q;
    prefactor = 1.0/(double(p)*double(q));
    prefactor *= ((i+1)%2 == 0) ? 1.0 : -1.0;
    A.initialize(m,m);
    for(j=1; j<m; ++j) {
      A.set(0,j,1.0);
      A.set(j,0,1.0);
    }
    for(j=0; j<n; ++j) {
      if (!skeleton->active_simplex(i,j) || !skeleton->simplices[i][j].get_modified()) continue;
      skeleton->simplices[i][j].get_vertices(vx);
      for(k=0; k<1+i; ++k) {
        for(l=k+1; l<1+i; ++l) {
          S.clear();
          S.insert(vx[k]);
          S.insert(vx[l]);
          qt = skeleton->index_table[1].find(S);
          V = skeleton->simplices[1][qt->second].get_squared_volume();
          A.set(1+k,1+l,V);
          A.set(1+l,1+k,V);
        }
      }
      V = prefactor*A.determinant();
      skeleton->simplices[i][j].set_squared_volume(V);
      skeleton->simplices[i][j].set_volume(std::sqrt(std::abs(V)));
      skeleton->simplices[i][j].set_modified(false);
    }
    p *= 2;
  }
}

void Spacetime::compute_lengths()
{
  int i,vx[2];
  double delta;
  const int ne = (signed) skeleton->simplices[1].size();

  for(i=0; i<ne; ++i) {
    if (!skeleton->active_simplex(1,i) || !skeleton->simplices[1][i].get_modified()) continue;
    skeleton->simplices[1][i].get_vertices(vx);
    delta = geometry->get_squared_distance(vx[0],vx[1],false);
    skeleton->simplices[1][i].set_squared_volume(delta);
    skeleton->simplices[1][i].set_volume(std::sqrt(std::abs(delta)));
    skeleton->simplices[1][i].set_modified(false);
  }
}

