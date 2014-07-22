#include "spacetime.h"

extern Random RND;

void Spacetime::get_coordinates(int v,std::vector<double>& x) const
{
  geometry->get_coordinates(v,x);
}

void Spacetime::get_coordinates(std::vector<double>& x) const
{
  int i,j;
  std::vector<double> vx;
  const int nv = cardinality(0,-1);

  x.clear();
  for(i=0; i<nv; ++i) {
    if (events[i].ubiquity == 1) continue;
    get_coordinates(i,vx);
    for(j=0; j<geometry->dimension(); ++j) {
      x.push_back(vx[j]);
    }
  }
}

double Spacetime::get_geometric_distance(int n,int m,bool t) const
{
  return geometry->get_distance(n,m,t);
}

double Spacetime::total_energy(int sheet) const
{
  const int nv = (signed) events.size();
  double sum = 0.0;

  if (sheet == -1) {
    for(int i=0; i<nv; ++i) {
      if (events[i].ubiquity > 1) sum += events[i].energy;
    }
  }
  else {
    for(int i=0; i<nv; ++i) {
      if (NTL::divide(events[i].ubiquity,codex[sheet].colour) == 1) sum += events[i].energy;
    }
  }
  return sum;
}

void Spacetime::chorogenesis(int nsteps)
{
  //assert(initial_state == RANDOM);
  assert(solver == MECHANICAL);
  assert(edge_probability > 0.3);
  assert(geometry->get_euclidean());
  assert(connected(-1));
  int i,j,dpopulation[1+geometry->dimension()];
  std::vector<std::pair<long,int> > factors;
  const int D = geometry->dimension();
  const int nv = (signed) events.size();

  factorize(nv,factors);
  j = 1;
  for(i=0; i<(signed) factors.size(); ++i) {
    assert(factors[i].second%D == 0);
    j *= ipow(factors[i].first,factors[i].second/D);
  }
  const int n = j;
  for(i=0; i<=D; ++i) {
    dpopulation[i] = ipow(2,i)*factorial(D)/(factorial(i)*factorial(D - i))*ipow(n - 2,D - i);
  }
  // Zero out the spacetime energy and make sure the vertex geometry is dimensionally homogeneous...
  std::vector<double> x,y;
  system_size = 0;
  for(i=0; i<nv; ++i) {
    events[i].energy = 0.0;
    geometry->get_coordinates(i,x);
    system_size += D;
    if (D == (signed) x.size()) continue;
    for(j=0; j<D; ++j) {
      y.push_back(x[j]);
    }
    geometry->set_coordinates(i,y);
    y.clear();
  }
  // We have now verified all the necessary conditions and can begin  
  compute_degree_distribution(false,-1);
  
  compute_connectivity_distribution(-1);

  int d,cfactor,derror,csize,vx[2];
  double temperature = 1.0;
  std::vector<int> candidates,reorder;
  std::vector<double> hgram;
  Graph* G = new Graph(cardinality(0,-1));
  const int ne = (signed) simplices[1].size();

#ifdef VERBOSE
  std::cout << "Beginning chorogenesis with " << ne << " edges..." << std::endl;
#endif  

  do {
    iterations += 1;
    for(i=0; i<ne; ++i) {
      if (simplices[1][i].ubiquity == 1) continue;
      simplices[1][i].get_vertices(vx);
      d = (signed) events[vx[0]].neighbours.size();
      if (d <= 2*geometry->dimension()) continue;
      d = (signed) events[vx[1]].neighbours.size();
      if (d <= 2*geometry->dimension()) continue;
      candidates.push_back(i);
    }
    if (candidates.empty()) {
      std::cout << "No viable candidates for topological operations, exiting loop..." << std::endl;
      break;
    }
    csize = (signed) candidates.size();
    cfactor = int(RND.drandom(0.1*temperature,temperature)*csize);
    if (cfactor == 0) {
      std::cout << "Temperature too low for topological operations, exiting loop..." << std::endl;
      break;
    }
    cfactor = std::max(cfactor,2*nv);
    std::cout << "There are " << csize << " and " << cfactor << " edges to delete at temperature " << temperature << std::endl;
    d = 0;
    RND.shuffle(reorder,csize);
    for(i=0; i<csize; ++i) {
      j = candidates[reorder[i]];
      simplex_deletion(1,j,-1);
      if (!connected(-1)) {
        simplices[1][j].ubiquity = 2;
        continue;
      }
      d++;
      if (d == cfactor) break;
    }
    std::cout << "Deleted " << d << " edges from the spacetime complex." << std::endl;
    regularization(false,-1);
    compute_graph(G,-1);
    G->degree_distribution(false,hgram);
    derror = 0;
    for(i=0; i<(signed) hgram.size(); ++i) {
      d = int(nv*hgram[i]);
      if (i > D) {
        derror += d;
      }
      else {
        derror += std::abs(d - dpopulation[i]);
      }
    }
    std::cout << "Vertex degree error is " << derror << std::endl;
    compute_connectivity_distribution(-1);
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
    if (simplices[1][i].ubiquity == 1) continue;
    simplices[1][i].get_vertices(vx);
    if (events[vx[0]].energy > 0.0 || events[vx[1]].energy > 0.0) continue;
    // Now check to see if one of these vertices has a neighbour with energy > 0
    e_neighbour = false;
    n = vx[0];
    for(it=events[n].neighbours.begin(); it!=events[n].neighbours.end(); ++it) {
      if (events[*it].energy > 0.0) e_neighbour = true;
    }
    if (e_neighbour) {
      flexible_edge[i] = 1;
      continue;
    }
    n = vx[1];
    for(it=events[n].neighbours.begin(); it!=events[n].neighbours.end(); ++it) {
      if (events[*it].energy > 0.0) e_neighbour = true;
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
  const double sq_cutoff = edge_flexibility_threshold*edge_flexibility_threshold;

  for(i=0; i<ne; ++i) {
    if (simplices[1][i].ubiquity == 1) continue;
    simplices[1][i].get_vertices(vx);
    d = geometry->get_distance(vx[0],vx[1],false);
    if (flexible_edge[i] == 1) {
      if (d > sq_cutoff) {
        d = std::sqrt(d);
        output += (d - edge_flexibility_threshold)*(d - edge_flexibility_threshold);
      }
      continue;
    }
    d = std::sqrt(d);
    ell = 1.0/(1.0 + pfactor*std::atan(0.5*(events[vx[0]].energy + events[vx[1]].energy)));
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
  const double sq_cutoff = edge_flexibility_threshold*edge_flexibility_threshold;

  for(i=0; i<system_size; ++i) {
    geometry->set_element(i,x[i]);
  }
  geometry->compute_distances();

  for(i=0; i<ne; ++i) {
    if (simplices[1][i].ubiquity == 1) continue;
    simplices[1][i].get_vertices(vx);
    d = geometry->get_distance(vx[0],vx[1],false);
    if (flexible_edge[i] == 1) {
      if (d > sq_cutoff) {
        d = std::sqrt(d);
        output += (d - edge_flexibility_threshold)*(d - edge_flexibility_threshold);
      }
      continue;
    }
    d = std::sqrt(d);
    ell = 1.0/(1.0 + pfactor*std::atan(0.5*(events[vx[0]].energy + events[vx[1]].energy)));
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
      geometry->add(i,geometry_cutoff);
      geometry->get_implied_vertices(i,current);
      compute_geometric_dependency(current);
      vmodified = last;
      for(it=current.begin(); it!=current.end(); it++) {
        vmodified.insert(*it);
      }
      geometry->compute_distances(vmodified);
      compute_volume();
      E = compute_abnormality();
      alpha = (E - B)/geometry_cutoff;
      df[i] = alpha;
      geometry->add(i,-geometry_cutoff);
      last = current;
      current.clear();
    }
    geometry->compute_distances(last);
  }
  else {
    int j,k,na = 0;
    double l,ell;
    std::vector<double> x1,x2;
    hash_map::const_iterator qt;
    const int nvertex = (signed) events.size();
    const int D = geometry->dimension();
    const double pfactor = (2.0/M_PI)*5.0;
    const double sq_cutoff = edge_flexibility_threshold*edge_flexibility_threshold;
    double alpha[D];

    for(i=0; i<nvertex; ++i) {
      if (events[i].ubiquity == 1) continue;
      na++;
    }

    assert(system_size == D*na);

#ifdef PARALLEL
#pragma omp parallel for default(shared) private(i,j,k,l,x1,x2,ell,alpha,it,qt)
#endif
    for(i=0; i<nvertex; ++i) {
      if (events[i].ubiquity == 1) continue;
      geometry->get_coordinates(i,x1);
      for(j=0; j<D; ++j) {
        alpha[j] = 0.0;
      }
      for(it=events[i].neighbours.begin(); it!=events[i].neighbours.end(); ++it) {
        k = *it;
        geometry->get_coordinates(k,x2);
        l = geometry->get_distance(i,k,false);
        qt = index_table[1].find(make_key(i,k));
        if (flexible_edge[qt->second] == 1) {
          if (l > sq_cutoff) {
            l = std::sqrt(l);
            for(j=0; j<D; ++j) {
              alpha[j] += 2.0*(l - edge_flexibility_threshold)*(x1[j] - x2[j])/l;
            }
          }
        }
        else {
          ell = 1.0/(1.0 + pfactor*std::atan(0.5*(events[i].energy + events[k].energy)));
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
  int i,j,v1,v2;
  double delta,mdelta = 1000.0;
  const int n1 = (signed) S1.size();
  const int n2 = (signed) S2.size();

  for(i=0; i<n1; ++i) {
    v1 = S1[i];
    for(j=0; j<n2; ++j) {
      v2 = S2[j];
      delta = geometry->get_distance(v1,v2,true);
      if (delta < 0.0) delta = -delta;
      if (delta < mdelta) {
        mdelta = delta;
        vx[0] = v1;
        vx[1] = v2;
      }
    }
  }
  return mdelta;
}

void Spacetime::reciprocate()
{
  geometry->reciprocate();

  compute_volume();
  compute_curvature();
  compute_obliquity();
}

bool Spacetime::delaunay() const
{
  // A method that checks if the spacetime, viewed as a simplicial complex,
  // satisfies the Delaunay property, that is the circumsphere of each d-simplex
  // contains no vertices, where d is the dimension of the complex.
  return false;
}

void Spacetime::compute_obliquity()
{
  const int nvertex = (signed) events.size();
  if (geometry->get_relational()) {
    for(int i=0; i<nvertex; ++i) {
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

  for(i=0; i<nvertex; ++i) {
    if (events[i].ubiquity == 1 || events[i].neighbours.size() < 2 || !events[i].geometry_modified) continue;

    j = *(events[i].neighbours.begin());

    geometry->vertex_difference(i,j,vx);

    rho = 0.0;
    first = true;
    for(it=events[i].neighbours.begin(); it!=events[i].neighbours.end(); it++) {
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
    if (rho < Spacetime::epsilon) rho = 0.0;
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
    if (events[i].ubiquity == 1) continue;
    alpha = 0.0;
    events[i].curvature = alpha;
  }
}

double Spacetime::representational_energy(bool weighted) const
{
  // A routine that follows the logic of C. Godsil and G. Royle, "Algebraic
  // Graph Theory" (Springer, 2001), Section 13.3 (pp. 284--286)
  double value,energy = 0.0;
  int i,j,k,l,vx[2],nv = 0;
  std::vector<int> offset;
  std::set<int>::const_iterator it;
  const int nvertex = (signed) events.size();
  const int nedge = (signed) simplices[1].size();

  if (dimension(-1) < 1) return energy;

  for(i=0; i<nvertex; ++i) {
    if (events[i].ubiquity == 1) {
      offset.push_back(-1);
      continue;
    }
    offset.push_back(nv);
    nv++;
  }

  std::vector<double>* laplacian = new std::vector<double>[nv];
  double diagonal[nv];
  for(i=0; i<nv; ++i) {
    diagonal[i] = 0.0;
  }

  if (weighted) {
    // What is the appropriate value for l?
    double w;
    l = 0;
    for(i=0; i<nedge; ++i) {
      if (simplices[1][i].ubiquity == 1) continue;
      value = std::abs(simplices[1][i].volume);
      w = double(1+l)/value;
      simplices[1][i].get_vertices(vx);
      j = offset[vx[0]];
      k = offset[vx[1]];
      diagonal[j] += w;
      diagonal[k] += w;
      laplacian[j].push_back(-w);
      laplacian[j].push_back(double(j));
      laplacian[k].push_back(-w);
      laplacian[k].push_back(double(i));
    }
  }
  else {
    // We can construct the Laplacian of the spacetime graph directly in this
    // case
    for(i=0; i<nedge; ++i) {
      if (simplices[1][i].ubiquity == 1) continue;
      simplices[1][i].get_vertices(vx);
      j = offset[vx[0]];
      k = offset[vx[1]];
      diagonal[j] += 1.0;
      diagonal[k] += 1.0;
      laplacian[j].push_back(-1.0);
      laplacian[j].push_back(double(vx[0]));
      laplacian[k].push_back(-1.0);
      laplacian[k].push_back(double(vx[1]));
    }
  }
  for(i=0; i<nvertex; ++i) {
    if (offset[i] == -1) continue;
    laplacian[offset[i]].push_back(diagonal[offset[i]]);
    laplacian[offset[i]].push_back(double(i));
  }

  energy = geometry->inner_product(laplacian,offset,nv);
  delete[] laplacian;

  return energy;
}

bool Spacetime::realizable(int d,int n) const
{
  // Can the d*(d+1)/2 edge lengths lead to a geometrically
  // realizable d-simplex?
  // An edge or vertex is always geometrically realizable...
  if (d < 2) return true;
  int i,j,info,dp1 = d + 1;
  hash_map::const_iterator qt;
  double alpha;
  bool output = true;
  char jtype = 'N';
  char uplo = 'U';
  int nwork = 3*dp1 - 1;
  int vx[dp1];
  double* A = new double[dp1*dp1];
  double* w = new double[dp1];
  double* work = new double[nwork];

  simplices[d][n].get_vertices(vx);

  for(i=0; i<dp1; ++i) {
    alpha = 0.0;
    if (i != 0) {
      qt = index_table[1].find(make_key(vx[i],vx[0]));
      alpha = simplices[1][qt->second].volume;
    }
    for(j=0; j<dp1; ++j) {
      if (j != 0) {
        qt = index_table[1].find(make_key(vx[j],vx[0]));
        alpha += simplices[1][qt->second].volume;
      }
      qt = index_table[1].find(make_key(vx[i],vx[j]));
      alpha -= simplices[1][qt->second].volume;
      A[dp1*i+j] = alpha;
    }
  }
  dsyev_(&jtype,&uplo,&dp1,A,&dp1,w,work,&nwork,&info);
  assert(info == 0);
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
  int i,j,k,l,n,m,parity,info,vx[3],pivots[ND+3];
  UINT64 q,p = 8;
  double prefactor,V,l1,l2,l3,A[(ND+3)*(ND+3)];
  hash_map::const_iterator qt;

  compute_lengths();

  for(i=0; i<(signed) simplices[2].size(); ++i) {
    // The triangles are a very simple case we can handle
    // without LAPACK
    if (simplices[2][i].ubiquity == 1 || !simplices[2][i].modified) continue;
    simplices[2][i].get_vertices(vx);
    qt = index_table[1].find(make_key(vx[0],vx[1]));
    l1 = simplices[1][qt->second].sq_volume;
    qt = index_table[1].find(make_key(vx[0],vx[2]));
    l2 = simplices[1][qt->second].sq_volume;
    qt = index_table[1].find(make_key(vx[1],vx[2]));
    l3 = simplices[1][qt->second].sq_volume;
    V = -(l3*l3 - 2.0*l3*(l1 + l2) + (l2 - l1)*(l2 - l1))/16.0;
    simplices[2][i].volume = std::sqrt(std::abs(V));
    simplices[2][i].sq_volume = V;
    simplices[2][i].modified = false;
  }

  for(i=3; i<=ND; ++i) {
    if (simplices[i].empty()) continue;
    m = i + 2;
    n = (signed) simplices[i].size();
    q = factorial(i);
    q *= q;
    prefactor = 1.0/(double(p)*double(q));
    prefactor *= ((i+1)%2 == 0) ? 1.0 : -1.0;
    A[0] = 0.0;
    for(j=1; j<m; ++j) {
      A[j*(1+m)] = 0.0;
      A[j] = 1.0;
      A[m*j] = 1.0;
    }
    for(j=0; j<n; ++j) {
      if (simplices[i][j].ubiquity == 1 || !simplices[i][j].modified) continue;
      simplices[i][j].get_vertices(pivots);
      for(k=0; k<1+i; ++k) {
        for(l=k+1; l<1+i; ++l) {
          qt = index_table[1].find(make_key(pivots[k],pivots[l]));
          V = simplices[1][qt->second].sq_volume;
          A[m*(1+k)+1+l] = V;
          A[m*(1+l)+1+k] = V;
        }
      }
      dgetrf_(&m,&m,A,&m,pivots,&info);
      assert(info >= 0);
      V = prefactor;
      parity = 1;
      for(k=0; k<m; ++k) {
        V *= A[k*(1+m)];
        if (pivots[i] < i) parity *= -1;
      }
      V = double(parity)*V;
      /*
      if (info > 0) {
        std::set<int>::const_iterator it;
        std::cout << "Null simplex: " << i << "  " << V << std::endl;
        for(it=simplices[i][j].vertices.begin(); it!=simplices[i][j].vertices.end(); it++) {
          std::cout << "[ ";
          for(k=0; k<(signed) events[*it].x.size(); ++k) {
            std::cout << events[*it].x[k] << " ";
          }
          std::cout << "]" << std::endl;
        }
        std::exit(1);
      }
      */
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
    if (simplices[1][i].ubiquity == 1 || !simplices[1][i].modified) continue;
    simplices[1][i].get_vertices(vx);
    delta = geometry->get_distance(vx[0],vx[1],true);
    simplices[1][i].sq_volume = delta;
    simplices[1][i].volume = std::sqrt(std::abs(delta));
    simplices[1][i].modified = false;
  }
}

