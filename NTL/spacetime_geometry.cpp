#include "spacetime.h"

extern SYNARMOSMA::Random RND;

using namespace DIAPLEXIS;

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

void Spacetime::arclength_statistics(double* output,int sheet) const
{
  int i;
  double wm,avg_length = 0.0,max_length = 0.0,min_length = 1000.0;
  const int Ne = (signed) simplices[1].size();
  const double ne = double(cardinality(1,sheet));

  if (sheet == -1) {
    for(i=0; i<Ne; ++i) {
      if (simplices[1][i].ubiquity == 1) continue;
      wm = simplices[1][i].volume;
      avg_length += wm;
      if (wm > max_length) max_length = wm;
      if (wm < min_length) min_length = wm;
    }
    if (ne > 0) avg_length /= ne;
  }
  else {
    for(i=0; i<Ne; ++i) {
      if (NTL::divide(simplices[1][i].ubiquity,codex[sheet].colour) == 0) continue;
      wm = simplices[1][i].volume;
      avg_length += wm;
      if (wm > max_length) max_length = wm;
      if (wm < min_length) min_length = wm;
    }
    if (ne > 0) avg_length /= ne;
  }
  output[0] = max_length;
  output[1] = min_length;
  output[2] = avg_length;
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

  SYNARMOSMA::factorize(nv,factors);
  j = 1;
  for(i=0; i<(signed) factors.size(); ++i) {
    assert(factors[i].second%D == 0);
    j *= SYNARMOSMA::ipow(factors[i].first,factors[i].second/D);
  }
  const int n = j;
  for(i=0; i<=D; ++i) {
    dpopulation[i] = SYNARMOSMA::ipow(2,i)*SYNARMOSMA::factorial(D)/(SYNARMOSMA::factorial(i)*SYNARMOSMA::factorial(D - i))*SYNARMOSMA::ipow(n - 2,D - i);
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
  SYNARMOSMA::Graph* G = new SYNARMOSMA::Graph(cardinality(0,-1));
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
#ifdef VERBOSE
    std::cout << "There are " << csize << " and " << cfactor << " edges to delete at temperature " << temperature << std::endl;
#endif
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
#ifdef VERBOSE
    std::cout << "Deleted " << d << " edges from the spacetime complex." << std::endl;
#endif
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
#ifdef VERBOSE
    std::cout << "Vertex degree error is " << derror << std::endl;
#endif
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
      for(it=current.begin(); it!=current.end(); ++it) {
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
    std::set<int> S;
    std::vector<double> x1,x2;
    SYNARMOSMA::hash_map::const_iterator qt;
    const int nv = (signed) events.size();
    const int D = geometry->dimension();
    const double pfactor = (2.0/M_PI)*5.0;
    const double sq_cutoff = edge_flexibility_threshold*edge_flexibility_threshold;
    double alpha[D];

    for(i=0; i<nv; ++i) {
      if (events[i].ubiquity == 1) continue;
      na++;
    }

    assert(system_size == D*na);

#ifdef PARALLEL
#pragma omp parallel for default(shared) private(i,j,k,l,x1,x2,ell,alpha,it,S,qt)
#endif
    for(i=0; i<nv; ++i) {
      if (events[i].ubiquity == 1) continue;
      geometry->get_coordinates(i,x1);
      for(j=0; j<D; ++j) {
        alpha[j] = 0.0;
      }
      for(it=events[i].neighbours.begin(); it!=events[i].neighbours.end(); ++it) {
        k = *it;
        geometry->get_coordinates(k,x2);
        l = geometry->get_distance(i,k,false);
        S.clear();
        S.insert(i);
        S.insert(k);
        qt = index_table[1].find(S);
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

double Spacetime::compute_temporal_vorticity(int v,int sheet) const
{
  // A value near zero means that the current topology is in accord with the chronogeometry
  // whereas a large positive value means much more topological entwinement is needed and a
  // large negative value means the topology must be altered to resemble that of the Cartesian
  // initial state.
  if (geometry->get_euclidean()) return 0.0;
  int in1,in2,tcount;
  double l,w1,w2,tipsy,vorticity;
  CAUSALITY d1,d2;
  std::set<int> n1,n2,S;
  std::vector<int> jset;
  std::set<int>::const_iterator it;
  std::vector<int>::const_iterator vit;
  SYNARMOSMA::hash_map::const_iterator qt;
  SYNARMOSMA::Graph* G = new SYNARMOSMA::Graph;

  // This makes a graph out of the entire forward light cone of this vertex
  compute_causal_graph(G,v,FUTURE,sheet);
  w1 = G->cyclicity();
  // And this makes one out of the entire backward light cone
  compute_causal_graph(G,v,PAST,sheet);
  w2 = G->cyclicity();
  // Now assign the arithmetic average of the cyclicity of these graphs to be the 
  // vertex's temporal vorticity
  vorticity = 0.5*(w1 + w2);
  delete G;
  n1 = events[v].neighbours;
  tipsy = 0.0;
  if (sheet == -1) {
    for(it=n1.begin(); it!=n1.end(); ++it) {
      in1 = *it;
      S.clear();
      S.insert(v);
      S.insert(in1);
      qt = index_table[1].find(S);
      l = simplices[1][qt->second].volume;
      if (!(simplices[1][qt->second].sq_volume > 0.0)) continue;
      // This edge is spacelike, so it will contribute to the temporal vorticity
      n2 = events[in1].neighbours;
      std::set_intersection(n1.begin(),n1.end(),n2.begin(),n2.end(),jset.begin());
      if (jset.empty()) continue;
      tcount = 0;
      for(vit=jset.begin(); vit!=jset.end(); ++vit) {
        in2 = *vit;

        S.clear();
        S.insert(v);
        S.insert(in2);
        qt = index_table[1].find(S);
        d1 = simplices[1][qt->second].orientation;

        S.clear();
        S.insert(in1);
        S.insert(in2);
        qt = index_table[1].find(S);
        d2 = simplices[1][qt->second].orientation;
        if (d1 != d2) tcount++;
      }
      tipsy += double(tcount)/(l*double(jset.size()));
    }
  }
  else {
    for(it=n1.begin(); it!=n1.end(); ++it) {
      in1 = *it;
      S.clear();
      S.insert(v);
      S.insert(in1);
      qt = index_table[1].find(S);
      if (NTL::divide(simplices[1][qt->second].ubiquity,codex[sheet].colour) == 0) continue;
      l = simplices[1][qt->second].volume;
      if (!(simplices[1][qt->second].sq_volume > 0.0)) continue;
      // This edge is spacelike, so it will contribute to the temporal vorticity
      n2 = events[in1].neighbours;
      std::set_intersection(n1.begin(),n1.end(),n2.begin(),n2.end(),jset.begin());
      if (jset.empty()) continue;
      tcount = 0;
      for(vit=jset.begin(); vit!=jset.end(); ++vit) {
        in2 = *vit;

        S.clear();
        S.insert(v);
        S.insert(in2);
        qt = index_table[1].find(S);
        if (NTL::divide(simplices[1][qt->second].ubiquity,codex[sheet].colour) == 0) continue;
        d1 = simplices[1][qt->second].orientation;

        S.clear();
        S.insert(in1);
        S.insert(in2);
        qt = index_table[1].find(S);
        if (NTL::divide(simplices[1][qt->second].ubiquity,codex[sheet].colour) == 0) continue;
        d2 = simplices[1][qt->second].orientation;
        if (d1 != d2) tcount++;
      }
      tipsy += double(tcount)/(l*double(jset.size()));
    }
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
      if (simplices[1][qt->second].sq_volume < 0.0) {
        nt++;
      }
      else {
        ns++;
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
    if (events[i].ubiquity == 1 || events[i].neighbours.size() < 2 || !events[i].geometry_modified) continue;

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
  int i,j,k,l,vx[2],n = 0;
  std::vector<int> offset;
  std::set<int>::const_iterator it;
  const int nv = (signed) events.size();
  const int ne = (signed) simplices[1].size();

  if (dimension(-1) < 1) return energy;

  for(i=0; i<nv; ++i) {
    if (events[i].ubiquity == 1) {
      offset.push_back(-1);
      continue;
    }
    offset.push_back(n);
    n++;
  }

  std::vector<double>* laplacian = new std::vector<double>[n];
  double diagonal[n];
  for(i=0; i<n; ++i) {
    diagonal[i] = 0.0;
  }

  if (weighted) {
    // What is the appropriate value for l?
    double w;
    l = 0;
    for(i=0; i<ne; ++i) {
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
    for(i=0; i<ne; ++i) {
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
  for(i=0; i<nv; ++i) {
    if (offset[i] == -1) continue;
    laplacian[offset[i]].push_back(diagonal[offset[i]]);
    laplacian[offset[i]].push_back(double(i));
  }

  energy = geometry->inner_product(laplacian,offset,n);
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
  std::set<int> S;
  SYNARMOSMA::hash_map::const_iterator qt;
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
  int i,j,k,l,n,m,vx[Spacetime::ND+3];
  SYNARMOSMA::UINT64 q,p = 8;
  double prefactor,V,l1,l2,l3,A[(Spacetime::ND+3)*(Spacetime::ND+3)];
  std::set<int> S;
  SYNARMOSMA::hash_map::const_iterator qt;

  compute_lengths();

  for(i=0; i<(signed) simplices[2].size(); ++i) {
    // The triangles are a very simple case we can handle
    // without LAPACK
    if (simplices[2][i].ubiquity == 1 || !simplices[2][i].modified) continue;    
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
    A[0] = 0.0;
    for(j=1; j<m; ++j) {
      A[j*(1+m)] = 0.0;
      A[j] = 1.0;
      A[m*j] = 1.0;
    }
    for(j=0; j<n; ++j) {
      if (simplices[i][j].ubiquity == 1 || !simplices[i][j].modified) continue;
      simplices[i][j].get_vertices(vx);
      for(k=0; k<1+i; ++k) {
        for(l=k+1; l<1+i; ++l) {
          S.clear();
          S.insert(vx[k]);
          S.insert(vx[l]);
          qt = index_table[1].find(S);
          V = simplices[1][qt->second].sq_volume;
          A[m*(1+k)+1+l] = V;
          A[m*(1+l)+1+k] = V;
        }
      }
      V = prefactor*SYNARMOSMA::determinant(A,m);
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

