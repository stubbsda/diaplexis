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
  const int nv = cardinality(0);

  x.clear();
  for(i=0; i<nv; ++i) {
    if (!events[i].active) continue;
    get_coordinates(i,vx);
    for(j=0; j<geometry->dimension(); ++j) {
      x.push_back(vx[j]);
    }
  }
}

double Spacetime::get_geometric_distance(int n,int m) const
{
  return geometry->get_squared_distance(n,m,false);
}

double Spacetime::total_energy() const
{
  const int nv = (signed) events.size();
  double sum = 0.0;

  for(int i=0; i<nv; ++i) {
    if (events[i].active) sum += events[i].get_energy();
  }
  return sum;
}

void Spacetime::arclength_statistics(double* output) const
{
  int i;
  double wm,avg_length = 0.0,max_length = 0.0,min_length = 1000.0;
  const int Ne = (signed) simplices[1].size();
  const double ne = double(cardinality(1));

  for(i=0; i<Ne; ++i) {
    if (!simplices[1][i].active) continue;
    wm = simplices[1][i].volume;
    avg_length += wm;
    if (wm > max_length) max_length = wm;
    if (wm < min_length) min_length = wm;
  }
  if (ne > 0) avg_length /= ne;
  output[0] = max_length;
  output[1] = min_length;
  output[2] = avg_length;
}

void Spacetime::compute_causal_graph(SYNARMOSMA::Directed_Graph* G,int base) const
{
  int i,j,v;
  std::set<int> S;
  std::set<int>::const_iterator it;
  std::vector<int> offset,current,next;
  std::vector<int>::const_iterator v_it;
  SYNARMOSMA::hash_map::const_iterator qt;
  const int nv = (signed) events.size();

  G->clear();
  current.push_back(base);
  for(i=0; i<nv; ++i) {
    offset.push_back(-1);
  }
  offset[base] = G->add_vertex();

  do {
    for(v_it=current.begin(); v_it!=current.end(); ++v_it) {
      v = *v_it;
      for(it=events[v].neighbours.begin(); it!=events[v].neighbours.end(); ++it) {
        j = *it;
        S.clear();
        S.insert(v);
        S.insert(j);
        qt = index_table[1].find(S);
        if (!simplices[1][qt->second].active) continue;
        if (!simplices[1][qt->second].timelike()) continue;
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
  std::set<int> current,next;
  std::set<int>::const_iterator it,jt;

  // First the past cone...
  current.insert(v);
  past_cone.clear();
  do {
    for(it=current.begin(); it!=current.end(); ++it) {
      i = *it;
      for(jt=events[i].anterior.begin(); jt!=events[i].anterior.end(); ++jt) {
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
      for(jt=events[i].posterior.begin(); jt!=events[i].posterior.end(); ++jt) {
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
  const int nv = (signed) events.size();
  const double na = double(cardinality(0));

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i,G,past,future) reduction(+:nlinearity,nsource,nsink,causal_loop)
#endif
  for(i=0; i<nv; ++i) {
    if (!events[i].active) continue;
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
  const int n = (signed) events.size();
  const int m = (signed) simplices[1].size();

  for(i=0; i<n; ++i) {
    events[i].anterior.clear();
    events[i].posterior.clear();
  }
  for(i=0; i<m; ++i) {
    if (!simplices[1][i].active) continue;
    if (!simplices[1][i].timelike()) continue;
    simplices[1][i].get_vertices(vx);
    if (geometry->get_temporal_order(vx[0],vx[1]) == SYNARMOSMA::Relation::before) {
      events[vx[0]].posterior.insert(vx[1]);
      events[vx[1]].anterior.insert(vx[0]);
    }
    else {
      events[vx[1]].posterior.insert(vx[0]);
      events[vx[0]].anterior.insert(vx[1]);
    }
  }
}

void Spacetime::chorogenesis(int nsteps)
{
  assert(solver == Geometry_Solver::mechanical);
  assert(edge_probability > 0.3);
  assert(geometry->get_euclidean());
  assert(connected());
  int dpopulation[1+geometry->dimension()];
  unsigned int i,j;
  std::vector<std::pair<long,int> > factors;
  const unsigned int D = geometry->dimension();
  const unsigned int nv = events.size();

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
    events[i].nullify_energy();
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
  compute_degree_distribution(false);
  
  compute_connectivity_distribution();

  int vx[2];
  unsigned int d,cfactor,derror,csize;
  double temperature = 1.0;
  std::vector<int> candidates,reorder;
  std::vector<double> hgram;
  SYNARMOSMA::Graph* G = new SYNARMOSMA::Graph;
  const unsigned int ne = simplices[1].size();

#ifdef VERBOSE
  std::cout << "Beginning chorogenesis with " << ne << " edges..." << std::endl;
#endif  

  do {
    iterations += 1;
    for(i=0; i<ne; ++i) {
      if (!simplices[1][i].active) continue;
      simplices[1][i].get_vertices(vx);
      d = events[vx[0]].neighbours.size();
      if (d <= 2*geometry->dimension()) continue;
      d = events[vx[1]].neighbours.size();
      if (d <= 2*geometry->dimension()) continue;
      candidates.push_back(i);
    }
    if (candidates.empty()) {
      std::cout << "No viable candidates for topological operations, exiting loop..." << std::endl;
      break;
    }
    csize = candidates.size();
    cfactor = int(RND->drandom(0.1*temperature,temperature)*csize);
    if (cfactor == 0) {
      std::cout << "Temperature too low for topological operations, exiting loop..." << std::endl;
      break;
    }
    cfactor = std::max(cfactor,2*nv);
#ifdef VERBOSE
    std::cout << "There are " << csize << " and " << cfactor << " edges to delete at temperature " << temperature << std::endl;
#endif
    d = 0;
    RND->shuffle(reorder,csize);
    for(i=0; i<csize; ++i) {
      j = candidates[reorder[i]];
      simplex_deletion(1,j);
      if (!connected()) {
        simplices[1][j].active = true;
        continue;
      }
      d++;
      if (d == cfactor) break;
    }
#ifdef VERBOSE
    std::cout << "Deleted " << d << " edges from the spacetime complex." << std::endl;
#endif
    regularization(false);
    compute_graph(G);
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
    compute_connectivity_distribution();
    // Now the geometry...
    optimize();
    // Prepare for the next iteration...
    temperature *= 0.95; 
    candidates.clear();    
  } while(true);
  if (iterations == 1) optimize();
  compute_volume();
  compute_curvature();
  compute_obliquity();
  structural_deficiency();
  write_state();
  delete G;
}

void Spacetime::determine_flexible_edges()
{
  int i,n,vx[2];
  bool e_neighbour;
  std::set<int>::const_iterator it;
  const int ne = (signed) simplices[1].size();

  flexible_edge.clear();
  for(i=0; i<ne; ++i) {
    flexible_edge.push_back(0);
  }

  for(i=0; i<ne; ++i) {
    if (!simplices[1][i].active) continue;
    simplices[1][i].get_vertices(vx);
    if (!events[vx[0]].zero_energy() || !events[vx[1]].zero_energy()) continue;
    // Now check to see if one of these vertices has a neighbour with energy > 0
    e_neighbour = false;
    n = vx[0];
    for(it=events[n].neighbours.begin(); it!=events[n].neighbours.end(); ++it) {
      if (!events[*it].zero_energy()) e_neighbour = true;
    }
    if (e_neighbour) {
      flexible_edge[i] = 1;
      continue;
    }
    n = vx[1];
    for(it=events[n].neighbours.begin(); it!=events[n].neighbours.end(); ++it) {
      if (!events[*it].zero_energy()) e_neighbour = true;
    }
    if (e_neighbour) flexible_edge[i] = 1;
  }
}

double Spacetime::compute_abnormality() const
{
  int i,vx[2];
  double d,ell,output = 0.0;
  const int ne = (signed) simplices[1].size();
  const double pfactor = (2.0/M_PI)*5.0;
  const double sq_tolerance = edge_flexibility_threshold*edge_flexibility_threshold;

  for(i=0; i<ne; ++i) {
    if (!simplices[1][i].active) continue;
    simplices[1][i].get_vertices(vx);
    d = geometry->get_squared_distance(vx[0],vx[1],false);
    if (flexible_edge[i] == 1) {
      if (d > sq_tolerance) {
        d = std::sqrt(d);
        output += (d - edge_flexibility_threshold)*(d - edge_flexibility_threshold);
      }
      continue;
    }
    d = std::sqrt(d);
    ell = 1.0/(1.0 + pfactor*std::atan(0.5*(events[vx[0]].get_energy() + events[vx[1]].get_energy())));
    output += (d - ell)*(d - ell);
  }
  return output;
}

double Spacetime::compute_abnormality(const std::vector<double>& x) const
{
  int i,vx[2];
  double d,ell,output = 0.0;
  const int ne = (signed) simplices[1].size();
  const double pfactor = (2.0/M_PI)*5.0;
  const double sq_tolerance = edge_flexibility_threshold*edge_flexibility_threshold;

  for(i=0; i<system_size; ++i) {
    geometry->set_element(i,x[i]);
  }
  geometry->compute_squared_distances();

  for(i=0; i<ne; ++i) {
    if (!simplices[1][i].active) continue;
    simplices[1][i].get_vertices(vx);
    d = geometry->get_squared_distance(vx[0],vx[1],false);
    if (flexible_edge[i] == 1) {
      if (d > sq_tolerance) {
        d = std::sqrt(d);
        output += (d - edge_flexibility_threshold)*(d - edge_flexibility_threshold);
      }
      continue;
    }
    d = std::sqrt(d);
    ell = 1.0/(1.0 + pfactor*std::atan(0.5*(events[vx[0]].get_energy() + events[vx[1]].get_energy())));
    output += (d - ell)*(d - ell);
  }
  return output;
}

void Spacetime::compute_geometric_gradient(std::vector<double>& df,bool negate)
{
  int i;
  std::set<int>::const_iterator it;

  df.clear();
  for(i=0; i<system_size; ++i) {
    df.push_back(0.0);
  }
  if (geometry->get_relational()) {
    std::set<int> vmodified,last,current;
    double E,alpha;
    const double B = compute_abnormality();

    for(i=0; i<system_size; ++i) {
      geometry->add(i,geometry_tolerance);
      geometry->get_implied_vertices(i,current);
      compute_geometric_dependency(current);
      vmodified = last;
      for(it=current.begin(); it!=current.end(); ++it) {
        vmodified.insert(*it);
      }
      geometry->compute_squared_distances(vmodified);
      compute_volume();
      E = compute_abnormality();
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
    std::set<int> S;
    std::vector<double> x1,x2;
    SYNARMOSMA::hash_map::const_iterator qt;
    const int nv = (signed) events.size();
    const int D = geometry->dimension();
    const double pfactor = (2.0/M_PI)*5.0;
    const double sq_tolerance = edge_flexibility_threshold*edge_flexibility_threshold;
    double alpha[D];

    for(i=0; i<nv; ++i) {
      if (!events[i].active) continue;
      na++;
    }
#ifdef DEBUG
    assert(system_size == D*na);
#endif
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i,j,k,l,S,x1,x2,ell,alpha,it,qt)
#endif
    for(i=0; i<nv; ++i) {
      if (!events[i].active) continue;
      geometry->get_coordinates(i,x1);
      for(j=0; j<D; ++j) {
        alpha[j] = 0.0;
      }
      for(it=events[i].neighbours.begin(); it!=events[i].neighbours.end(); ++it) {
        k = *it;
        geometry->get_coordinates(k,x2);
        l = geometry->get_squared_distance(i,k,false);
        S.clear();
        S.insert(i);
        S.insert(k);
        qt = index_table[1].find(S);
        if (flexible_edge[qt->second] == 1) {
          if (l > sq_tolerance) {
            l = std::sqrt(l);
            for(j=0; j<D; ++j) {
              alpha[j] += 2.0*(l - edge_flexibility_threshold)*(x1[j] - x2[j])/l;
            }
          }
        }
        else {
          ell = 1.0/(1.0 + pfactor*std::atan(0.5*(events[i].get_energy() + events[k].get_energy())));
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

bool Spacetime::delaunay() const
{
  // A method that checks if the spacetime, viewed as a simplicial complex,
  // satisfies the Delaunay property, that is the circumsphere of each d-simplex
  // contains no vertices, where d is the dimension of the complex.
  return false;
}

double Spacetime::compute_temporal_vorticity(int v) const
{
  // A value near zero means that the current topology is in accord with the chronogeometry
  // whereas a large positive value means much more topological entwinement is needed and a
  // large negative value means the topology must be altered to resemble that of the Cartesian
  // initial state.
  int u,w,tcount;
  double l,tipsy,vorticity;
  std::set<int> S;
  std::vector<int> jset;
  std::set<int>::const_iterator it,jt;
  std::vector<int>::const_iterator vit;
  SYNARMOSMA::Relation d1,d2;
  SYNARMOSMA::Directed_Graph G;
  SYNARMOSMA::hash_map::const_iterator qt;

  compute_causal_graph(&G,v);
  vorticity = G.cyclicity();

  tipsy = 0.0;
  for(it=events[v].neighbours.begin(); it!=events[v].neighbours.end(); ++it) {
    u = *it;
    S.clear();
    S.insert(v);
    S.insert(u);
    qt = index_table[1].find(S);
    l = simplices[1][qt->second].volume;
    if (!simplices[1][qt->second].spacelike()) continue;
    // This edge is spacelike, so it will contribute to the temporal vorticity
    jset.clear();
    for(jt=events[v].neighbours.begin(); jt!=events[v].neighbours.end(); ++jt) {
      w = *jt;
      if (events[u].neighbours.count(w) > 0) jset.push_back(w);
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

int Spacetime::simplex_embedding(int d,int n) const
{
  int i,j,in1,vx[1+d],ns = 0,nt = 0,p = 0,nwork = 3*n - 1;
  double delta,sum,s1,s2,c,dmatrix[n*n],A[n*n],a[n],b[n],w[n],work[nwork];
  char jtype = 'V',uplo = 'U';
  std::set<int> S;
  std::set<int>::const_iterator it;
  SYNARMOSMA::hash_map::const_iterator qt;

  simplices[d][n].get_vertices(vx);

  for(i=0; i<(1+d)*(1+d); ++i) {
    dmatrix[i] = 0.0;
  }

  // First we need to check if all the edges of this simplex are
  // of the same type
  for(i=0; i<1+d; ++i) {
    for(j=i+1; j<1+d; ++j) {
      S.clear();
      S.insert(vx[i]);
      S.insert(vx[j]);
      qt = index_table[1].find(S);
      delta = std::abs(simplices[1][qt->second].volume);
      dmatrix[n*i+j] = delta;
      dmatrix[n*j+i] = delta;
      if (simplices[1][qt->second].spacelike()) {
        ns++;
      }
      else if (simplices[1][qt->second].timelike()) {
        nt++;
      }
    }
  }
  if (nt > 0 && ns > 0) return -1;
  if (ns == 0) {
    // If the simplex's edges are all timelike, we should change the
    // sign from negative to positive before proceeding with the
    // calculation
    for(i=0; i<n*n; ++i) {
      dmatrix[i] = -dmatrix[i];
    }
  }

  for(i=0; i<n*n; ++i) {
    dmatrix[i] = -0.5*dmatrix[i]*dmatrix[i];
  }

  sum = 0.0;
  for(i=0; i<n; ++i) {
    s1 = 0.0;
    s2 = 0.0;
    for(j=0; j<n; ++j) {
      sum += dmatrix[n*i+j];
      s1 += dmatrix[n*i+j];
      s2 += dmatrix[n*j+i];
    }
    a[i] = s1/double(n);
    b[i] = s2/double(n);
  }
  c = sum/double(n*n);
  for(i=0; i<n; ++i) {
    for(j=0; j<n; ++j) {
      A[n*i+j] = dmatrix[n*i+j] - a[i] - b[j] + c;
    }
  }
  dsyev_(&jtype,&uplo,&n,A,&n,w,work,&nwork,&in1);
  for(i=0; i<n; ++i) {
    if (std::abs(w[i]) > 0.0) p++;
  }
  return p;
}

void Spacetime::compute_obliquity()
{
  const int nv = (signed) events.size();
  if (geometry->get_relational()) {
    for(int i=0; i<nv; ++i) {
      events[i].obliquity = 0.0;
    }
    return;
  }
  int i,j;
  double rho,theta,alpha;
  bool first;
  std::vector<double> vx,vy;
  std::set<int>::const_iterator it;

  const double A = 2.5;

  for(i=0; i<nv; ++i) {
    if (!events[i].active || events[i].neighbours.size() < 2 || !events[i].geometry_modified) continue;

    j = *(events[i].neighbours.begin());

    geometry->vertex_difference(i,j,vx);

    rho = 0.0;
    first = true;
    for(it=events[i].neighbours.begin(); it!=events[i].neighbours.end(); ++it) {
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
    rho = rho/double(events[i].neighbours.size() - 1);
    if (rho < std::numeric_limits<double>::epsilon()) rho = 0.0;
    events[i].obliquity = rho;
  }
}

void Spacetime::compute_curvature()
{
  int i;
  double alpha;
  const int nv = (signed) events.size();

  for(i=0; i<nv; ++i) {
    events[i].curvature = 0.0;
  }
  for(i=0; i<nv; ++i) {
    // To calculate the curvature of v, we need to find the greatest
    // simplicial dimension of v subject to the requirement that its
    // entourage at this dimensionality be greater than one.
    if (!events[i].active) continue;
    alpha = 0.0;
    events[i].curvature = alpha;
  }
}

double Spacetime::representational_energy(bool weighted) const
{
  // A routine that follows the algorithm outlined in C. Godsil and G. Royle, "Algebraic
  // Graph Theory" (Springer, 2001), Section 13.3 (pp. 284--286)
  if (dimension() < 1) return 0.0;
  int i,n = 0;
  std::vector<int> offset;
  SYNARMOSMA::Graph G;
  const int nv = (signed) events.size();
  
  compute_graph(&G);

  SYNARMOSMA::Matrix<double> L(G.order());
  
  for(i=0; i<nv; ++i) {
    if (!events[i].active) {
      offset.push_back(-1);
      continue;
    }
    offset.push_back(n);
    n++;
  }

  if (weighted) {
    int j,k,vx[2];
    double w;
    const int ne = (signed) simplices[1].size();
  
    for(i=0; i<ne; ++i) {
      if (!simplices[1][i].active) continue;
      w = 2.0/std::abs(simplices[1][i].volume);
      simplices[1][i].get_vertices(vx);
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

  return geometry->inner_product(L,offset);
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

  simplices[d][n].get_vertices(vx);

  for(i=0; i<dp1; ++i) {
    alpha = 0.0;
    if (i != 0) {
      S.clear();
      S.insert(vx[i]);
      S.insert(vx[0]);
      qt = index_table[1].find(S);
      alpha = simplices[1][qt->second].volume;
    }
    for(j=0; j<dp1; ++j) {
      if (j != 0) {
        S.clear();
        S.insert(vx[j]);
        S.insert(vx[0]);
        qt = index_table[1].find(S);
        alpha += simplices[1][qt->second].volume;
      }
      S.clear();
      S.insert(vx[i]);
      S.insert(vx[j]);
      qt = index_table[1].find(S);
      alpha -= simplices[1][qt->second].volume;
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
  int i,j,k,l,n,m,vx[Spacetime::ND+3];
  double prefactor,V,l1,l2,l3;
  std::set<int> S;
  SYNARMOSMA::UINT64 q,p = 8;
  SYNARMOSMA::Matrix<double> A;
  SYNARMOSMA::hash_map::const_iterator qt;

  compute_lengths();

  for(i=0; i<(signed) simplices[2].size(); ++i) {
    // The triangles are a very simple case we can handle
    // without LAPACK
    if (!simplices[2][i].active || !simplices[2][i].modified) continue;
    qt = index_table[1].find(simplices[2][i].faces[0]);
    l1 = simplices[1][qt->second].sq_volume;
    qt = index_table[1].find(simplices[2][i].faces[1]);
    l2 = simplices[1][qt->second].sq_volume;
    qt = index_table[1].find(simplices[2][i].faces[2]);
    l3 = simplices[1][qt->second].sq_volume;
    V = -(l3*l3 - 2.0*l3*(l1 + l2) + (l2 - l1)*(l2 - l1))/16.0;
    simplices[2][i].volume = std::sqrt(std::abs(V));
    simplices[2][i].sq_volume = V;
    simplices[2][i].modified = false;
  }

  for(i=3; i<=Spacetime::ND; ++i) {
    if (simplices[i].empty()) continue;
    m = i + 2;
    n = (signed) simplices[i].size();
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
      if (!simplices[i][j].active || !simplices[i][j].modified) continue;
      simplices[i][j].get_vertices(vx);
      for(k=0; k<1+i; ++k) {
        for(l=k+1; l<1+i; ++l) {
          S.clear();
          S.insert(vx[k]);
          S.insert(vx[l]);
          qt = index_table[1].find(S);
          V = simplices[1][qt->second].sq_volume;
          A.set(1+k,1+l,V);
          A.set(1+l,1+k,V);
        }
      }
      V = prefactor*A.determinant();
      simplices[i][j].volume = std::sqrt(std::abs(V));
      simplices[i][j].sq_volume = V;
      simplices[i][j].modified = false;
    }
    p *= 2;
  }
}

void Spacetime::compute_lengths()
{
  int i,vx[2];
  double delta;
  const int ne = (signed) simplices[1].size();

  for(i=0; i<ne; ++i) {
    if (!simplices[1][i].active || !simplices[1][i].modified) continue;
    simplices[1][i].get_vertices(vx);
    delta = geometry->get_squared_distance(vx[0],vx[1],false);
    simplices[1][i].sq_volume = delta;
    simplices[1][i].volume = std::sqrt(std::abs(delta));
    simplices[1][i].modified = false;
  }
}

