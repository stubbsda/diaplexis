#include "spacetime.h"

using namespace DIAPLEXIS;

bool Spacetime::interplication(int centre,double size,int D,int sheet)
{
  // Method to construct a highly-entwined knot after carving out a hole for 
  // in the spacetime complex - it really only makes sense if D is large 
  // enough...  
  if (D <= 3) return false;

  const double dm = double(D);

  int i,j,k,d,p,q,m,nc,nbridge,its;
  double l,alpha,width = size/dm;
  std::vector<double> xc,cvertex;
  std::vector<std::pair<int,double> > ambient;
  std::set<int> S,base,nvertex,bvertex,kvertex;
  std::set<int>::const_iterator it;
  const int g_dim = (signed) geometry->dimension();

  // We begin by eliminating the central vertex along with any vertices on this 
  // sheet that are within the a sphere of radius "size"
  geometry->get_coordinates(centre,cvertex); 
  assert(vertex_deletion(centre,sheet));
  for(i=0; i<(signed) events.size(); ++i) {
    if (!events[i].active(sheet)) continue;
    l = geometry->get_squared_distance(centre,i,false);
    if (l < size) {
      assert(vertex_deletion(i,sheet));
      continue;
    }
    ambient.push_back(std::pair<int,double>(i,l - size));
  }
  std::sort(ambient.begin(),ambient.end(),SYNARMOSMA::pair_predicate_dbl);

  // Ideally the number of boundary vertices should be the surface area of a d-sphere
  d = g_dim;
  alpha = double(d)/2.0;
  l = std::pow(M_PI,alpha)*std::pow(size,double(d - 1))/boost::math::tgamma(1.0 + alpha);  
  d *= int(l);
  q = (signed) ambient.size();
  if (d > q) d = q;
  for(i=0; i<d; ++i) {
    // If we can't get enough boundary vertices that are close enough, exit
    if (ambient[i].second > 1.5*size) break;
    bvertex.insert(ambient[i].first);
  }
  assert(bvertex.size() > 1);

  // First add the central element of the knot, a D-simplex
  for(i=0; i<1+D; ++i) {
    for(j=0; j<g_dim; ++j) {
      xc.push_back(cvertex[j] + width*(RND->drandom() - 0.5));
    }
    if (!geometry->get_uniform()) {
      for(j=g_dim; j<D; ++j) {
        xc.push_back(-dm + dm/2.0*RND->drandom()); 
      }
    }
    S.insert(vertex_addition(xc,sheet));
    xc.clear();
  }
  simplex_addition(S,sheet);

  // Now add a set of lower-dimensional simplices to this knot and also tie it 
  // in to the existing spacetime complex...
  base = S;
  kvertex = S;
  S.clear();
  for(i=D-1; i>=3; --i) {
    // As the dimension shrinks, the number of simplices should grow
    l = double(i);
    m = int(RND->nrandom(dm - l,2.0));
#ifdef VERBOSE
    std::cout << "Adding " << m << " " << i << "-simplices to this complex" << std::endl;
#endif
    if (m <= 0) continue;
    nc = 0;
    width = 2.0*size/l;
    if (width > size) width = size;
    do {    
      for(j=0; j<1+i; ++j) {
        q = -1;
        if (RND->drandom() < (0.5*l/dm)) {
          for(k=0; k<g_dim; ++k) {
            xc.push_back(cvertex[k] + width*(RND->drandom() - 0.5));
          }
          if (!geometry->get_uniform()) {
            for(k=g_dim; k<i; ++k) {
              alpha = RND->nrandom(-l,dm);
              if (alpha < -dm) alpha = -dm + 0.1 + 1.5*RND->drandom();
              if (alpha > -0.1) alpha = -0.1 - RND->drandom();
              xc.push_back(alpha); 
            }
          }
          q = vertex_addition(xc,sheet);
          xc.clear();
          nvertex.insert(q);
          kvertex.insert(q);
        }
        else {
          its = 0;
          do {
            its++;
            q = RND->irandom(base);
            // Make sure it isn't already in S,
            if (S.count(q) == 0) break;
          } while(its < 2*((signed) base.size()));
        }
        if (q == -1) {
          // Create a new vertex from scratch
          for(k=0; k<g_dim; ++k) {
            xc.push_back(cvertex[k] + width*(RND->drandom() - 0.5));
          }
          if (!geometry->get_uniform()) {
            for(k=g_dim; k<i; ++k) {
              alpha = RND->nrandom(-l,dm);
              if (alpha < -dm) alpha = -dm + 0.1 + 1.5*RND->drandom();
              if (alpha > -0.1) alpha = -0.1 - RND->drandom();
              xc.push_back(alpha); 
            }
          }
          q = vertex_addition(xc,sheet);
          xc.clear();
          nvertex.insert(q);
          kvertex.insert(q);          
        }
        S.insert(q);
      }
#ifdef VERBOSE
      std::cout << "Created " << i << "-simplex with " << nvertex.size() << " new vertices" << std::endl;
#endif
      simplex_addition(S,sheet);
      S.clear();
      nc++;
    } while(nc < m);
    base = nvertex;
    nvertex.clear();
  }
  // Finally the 2-simplexes to bridge the knot with the ambient spacetime complex
  d = bvertex.size()*(bvertex.size() - 1)/2;
  nbridge = int((0.1 + 0.15*RND->drandom())*d);
#ifdef VERBOSE
  std::cout << "Adding " << nbridge << " 2-simplices" << std::endl;
#endif
  d = 0;
  do {
    q = RND->irandom(bvertex);
    if (RND->drandom() < 0.33) {
      do { 
        p = RND->irandom(bvertex);
        if (p != q) break;
      } while(true);
    }
    else {
      alpha = 10.0*size;
      for(it=kvertex.begin(); it!=kvertex.end(); ++it) {
        l = geometry->get_squared_distance(q,*it,false);
        if (l < alpha) {
          p = *it;
          alpha = l;
        }
      }
    }
    S.insert(q);
    S.insert(p);
    // Find the closest new vertex...
    ambient.clear();
    for(it=kvertex.begin(); it!=kvertex.end(); ++it) {
      if (p == *it) continue;
      l = geometry->get_squared_distance(q,*it,false);
      alpha = geometry->get_squared_distance(p,*it,false);
      ambient.push_back(std::pair<int,double>(*it,l + alpha));
    }
    std::sort(ambient.begin(),ambient.end(),SYNARMOSMA::pair_predicate_dbl);
    if (RND->drandom() < 0.5) {
      S.insert(ambient[0].first);
    }
    else {
      S.insert(ambient[1].first);
    }
    if (simplex_addition(S,sheet)) d++;
    S.clear();
  } while(d < nbridge);

  regularization(true,sheet);
  return true;
}

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

bool Spacetime::stellar_addition(int base,int sheet)
{
  if (vertex_dimension(base,sheet) < 2) return false;
  // Take one of these 2-simplices and eliminate it in favour of a 
  // new vertex
  int i,m = (signed) simplices[2].size();
  unsigned int l;
  std::vector<int> sx;
  std::vector<double> xc,xtemp;
  std::set<int> candidates,S;
  SYNARMOSMA::hash_map::const_iterator qt;

  for(i=0; i<m; ++i) {
    if (!simplices[2][i].active(sheet)) continue;
    if (simplices[2][i].contains(base)) candidates.insert(i);
  }
  // Choose a 2-simplex at random, get its three vertices and then 
  // eliminate each of the three edges in succession...
  m = RND->irandom(candidates);
  simplices[2][m].get_vertices(sx);

  S.clear();
  S.insert(sx[0]); S.insert(sx[1]);
  qt = index_table[1].find(S);
  simplex_deletion(1,qt->second,sheet);

  S.clear();
  S.insert(sx[0]); S.insert(sx[2]);
  qt = index_table[1].find(S);
  simplex_deletion(1,qt->second,sheet);

  S.clear();
  S.insert(sx[1]); S.insert(sx[2]);
  qt = index_table[1].find(S);
  simplex_deletion(1,qt->second,sheet);    
  // Now add the new vertex and the three edges...
  for(l=0; l<geometry->dimension(); ++l) {
    xc.push_back(0.0);
  }
  for(i=0; i<3; ++i) {
    geometry->get_coordinates(sx[i],xtemp);
    for(l=0; l<geometry->dimension(); ++l) {
      xc[l] += xtemp[l]; 
    }
  }
  for(l=0; l<geometry->dimension(); ++l) {
    xc[l] = xc[l]/3.0;
  }
  m = vertex_addition(xc,sheet);
  for(i=0; i<3; ++i) {
    S.clear();
    S.insert(m);
    S.insert(sx[i]);
    simplex_addition(S,sheet);
  }
  regularization(true,sheet);
  return true;    
}

bool Spacetime::germination(int base,int sheet)
{
  // This method constructs new neighbour vertices w_i for the vertex base which are
  // unit distance from v and orthogonal to v's existing edges, if possible.
  if (events[base].boundary) return false;

  bool good,modified = false;
  double a,b,d_min,delta;
  std::set<int> N,Dm2,Dm1,locus,current,tset,free_dims;
  std::set<int>::const_iterator it,jt;
  std::vector<double> x,y,z,xc,bvector;
  SYNARMOSMA::hash_map::const_iterator qt;
  Simplex S;
  int i,j,vx[2],D1 = -1,D2 = -1,m,in1,n = -1,mi = 0;
  const int nv = (signed) events.size();
  const int ne = (signed) simplices[1].size();
  const int D = (signed) geometry->dimension();

  locus.insert(sheet);

  for(i=0; i<ne; ++i) {
    if (!simplices[1][i].active(sheet)) continue;
    simplices[1][i].get_vertices(vx);
    if (base != vx[0] && base != vx[1]) continue;
    j = (vx[0] == base) ? vx[1] : vx[0];
    N.insert(j);
  }
  for(it=N.begin(); it!=N.end(); ++it) {
    n = *it;
    if (events[n].active()) break;
  }
#ifdef VERBOSE
  std::cout << "Germination with " << base << "  " << n << std::endl;
#endif
  Dm2.insert(n);

  // The problem here is that this assumes that a^2 + b^2 != 0, which we are
  // by no means guaranteed if the background dimension > 2
  geometry->vertex_difference(n,base,z);
  geometry->get_coordinates(base,bvector);
  xc = bvector;
  for(i=0; i<D; ++i) {
    x.push_back(0.0);
  }
  for(i=0; i<D; ++i) {
    for(j=1+i; j<D; ++j) {
      a = z[i];
      b = z[j];
      delta = a*a + b*b;
      if (delta < std::numeric_limits<double>::epsilon()) continue;
      D1 = i;
      D2 = j;
      x[D1] = a/std::sqrt(a*a+b*b);
      x[D2] = b/std::sqrt(a*a+b*b);
      break;
    }
  }
  if (D1 == -1) return false;

  // Check the three possible locations for a vertex orthogonal to this edge:
  // First, a 90 degree rotation...
  d_min = 0.5;
  in1 = -1;
  xc[D1] = -x[D2] + bvector[D1];
  xc[D2] =  x[D1] + bvector[D2];
  for(i=0; i<nv; ++i) {
    if (i == base) continue;
    delta = geometry->get_squared_distance(i,xc);
    delta = std::sqrt(delta);
    if (delta < d_min) {
      in1 = i;
      d_min = delta;
    }
  }
  if (in1 >= 0) {
    if (!events[in1].active(sheet)) {
      modified = true;
      events[in1].set_active(sheet);
    }
    S.initialize(base,in1,locus);
    qt = index_table[1].find(S.vertices);
    if (qt == index_table[1].end()) {
      modified = true;
      simplices[1].push_back(S);
      index_table[1][S.vertices] = (signed) simplices[1].size() - 1;
    }
    else {
      if (!simplices[1][qt->second].active(sheet)) {
        modified = true;
        simplices[1][qt->second].set_active(sheet);
      }
    }
  }
  else {
    modified = true;
    m = vertex_addition(base,sheet);
    Dm1.insert(m);
    geometry->set_coordinates(m,xc);
    S.initialize(base,m,locus);
    simplices[1].push_back(S);
    index_table[1][S.vertices] = (signed) simplices[1].size() - 1;
  }

  // Second, a 180 degree rotation...
  d_min = 0.5;
  in1 = -1;
  xc[D1] = -x[D1] + bvector[D1];
  xc[D2] = -x[D2] + bvector[D2];
  for(i=0; i<nv; ++i) {
    if (i == base) continue;
    delta = geometry->get_squared_distance(i,xc);
    delta = std::sqrt(delta);
    if (delta < d_min) {
      in1 = i;
      d_min = delta;
    }
  }
  if (in1 >= 0) {
    if (!events[in1].active(sheet)) {
      modified = true;
      events[in1].set_active(sheet);
    }
    S.initialize(base,in1,locus);
    qt = index_table[1].find(S.vertices);
    if (qt == index_table[1].end()) {
      modified = true;
      simplices[1].push_back(S);
      index_table[1][S.vertices] = (signed) simplices[1].size() - 1;
    }
    else {
      if (!simplices[1][qt->second].active(sheet)) {
        modified = true;
        simplices[1][qt->second].set_active(sheet);
      }
    }
  }
  else {
    modified = true;
    m = vertex_addition(base,sheet);
    geometry->set_coordinates(m,xc);
    S.initialize(base,m,locus);
    simplices[1].push_back(S);
    index_table[1][S.vertices] = (signed) simplices[1].size() - 1;
  }

  // And finally the 270 degree rotation...
  d_min = 0.5;
  in1 = -1;
  xc[D1] =  x[D2] + bvector[D1];
  xc[D2] = -x[D1] + bvector[D2];
  for(i=0; i<nv; ++i) {
    if (i == base) continue;
    delta = geometry->get_squared_distance(i,xc);
    delta = std::sqrt(delta);
    if (delta < d_min) {
      in1 = i;
      d_min = delta;
    }
  }
  if (in1 >= 0) {
    if (!events[in1].active(sheet)) {
      modified = true;
      events[in1].set_active(sheet);
    }
    S.initialize(base,in1,locus);
    qt = index_table[1].find(S.vertices);
    if (qt == index_table[1].end()) {
      modified = true;
      simplices[1].push_back(S);
      index_table[1][S.vertices] = (signed) simplices[1].size() - 1;
    }
    else {
      if (!simplices[1][qt->second].active(sheet)) {
        simplices[1][qt->second].set_active(sheet);
      }
    }
  }
  else {
    modified = true;
    m = vertex_addition(base,sheet);
    Dm1.insert(m);
    geometry->set_coordinates(m,xc);
    S.initialize(base,m,locus);
    simplices[1].push_back(S);
    index_table[1][S.vertices] = (signed) simplices[1].size() - 1;
  }

  if (modified) {
    regularization(false,sheet);
    regularization(false,-1);
  }

  // Now we use rolling cross products to branch out into higher
  // dimensions...
  for(i=0; i<D; ++i) {
    if (i == D1 || i == D2) continue;
    free_dims.insert(i);
  }

  do {
    if (Dm1.empty() || Dm2.empty()) break;

    good = false;
    for(it=Dm2.begin(); it!=Dm2.end(); ++it) {
      n = *it;
      geometry->vertex_difference(n,base,y);
      x.clear();
      x.push_back(y[D1]);
      x.push_back(y[D2]);
      x.push_back(0.0);
      for(jt=free_dims.begin(); jt!=free_dims.end(); ++jt) {
        mi = *jt;
        x[2] = y[mi];
        delta = SYNARMOSMA::norm(x);
        if (delta > std::numeric_limits<double>::epsilon()) {
          good = true;
          break;
        }
      }
    }
    if (good) {
      for(i=0; i<3; ++i) {
        x[i] = x[i]/delta;
      }
    }
    else {
      break;
    }

    // What is m equal to?
    good = false;
    for(it=Dm1.begin(); it!=Dm1.end(); ++it) {
      m = *it;
      geometry->vertex_difference(m,base,z);
      y.clear();
      y.push_back(z[D1]);
      y.push_back(z[D2]);
      y.push_back(z[mi]);
      delta = SYNARMOSMA::norm(y);
      if (delta > std::numeric_limits<double>::epsilon()) {
        good = true;
        break;
      }
    }
    if (good) {
      for(i=0; i<3; ++i) {
        y[i] = y[i]/delta;
      }
    }
    else {
      break;
    }

    SYNARMOSMA::cross_product(x,y,z);

    // Now look to see if there are any existing vertices near xc
    // and -xc...
    xc.clear();
    for(i=0; i<D; ++i) {
      xc.push_back(bvector[i]);
    }
    xc[D1] = z[0] + bvector[D1];
    xc[D2] = z[1] + bvector[D2];
    xc[mi] = z[2] + bvector[mi];

    d_min = 0.5;
    in1 = -1;
    for(i=0; i<nv; ++i) {
      if (i == base) continue;
      delta = geometry->get_squared_distance(i,xc);
      delta = std::sqrt(delta);
      if (delta < d_min) {
        in1 = i;
        d_min = delta;
      }
    }
    if (in1 >= 0) {
      if (!events[in1].active(sheet)) {
        modified = true;
        events[in1].set_active(sheet);
      }
      S.initialize(base,in1,locus);
      qt = index_table[1].find(S.vertices);
      if (qt == index_table[1].end()) {
        modified = true;
        simplices[1].push_back(S);
        index_table[1][S.vertices] = (signed) simplices[1].size() - 1;
      }
      else {
        if (!simplices[1][qt->second].active(sheet)) {
          simplices[1][qt->second].set_active(sheet);
        }
      }
    }
    else {
      modified = true;
      m = vertex_addition(base,sheet);
      current.insert(m);
      geometry->set_coordinates(m,xc);
      S.initialize(base,m,locus);
      simplices[1].push_back(S);
      index_table[1][S.vertices] = (signed) simplices[1].size() - 1;
    }

    // Now the mirror image in a three-dimensional subspace...
    xc[D1] = -z[0] + bvector[D1];
    xc[D2] = -z[1] + bvector[D2];
    xc[mi] = -z[2] + bvector[mi];
    d_min = 0.5;
    in1 = -1;
    for(i=0; i<nv; ++i) {
      if (i == base) continue;
      delta = geometry->get_squared_distance(i,xc);
      delta = std::sqrt(delta);
      if (delta < d_min) {
        in1 = i;
        d_min = delta;
      }
    }
    if (in1 >= 0) {
      if (!events[in1].active(sheet)) {
        modified = true;
        events[in1].set_active(sheet);
      }
      S.initialize(base,in1,locus);
      qt = index_table[1].find(S.vertices);
      if (qt == index_table[1].end()) {
        modified = true;
        simplices[1].push_back(S);
        index_table[1][S.vertices] = (signed) simplices[1].size() - 1;
      }
      else {
        if (!simplices[1][qt->second].active(sheet)) {
          simplices[1][qt->second].set_active(sheet);
        }
      }
    }
    else {
      modified = true;
      m = vertex_addition(base,sheet);
      current.insert(m);
      geometry->set_coordinates(m,xc);
      S.initialize(base,m,locus);
      simplices[1].push_back(S);
      index_table[1][S.vertices] = (signed) simplices[1].size() - 1;
    }

    for(it=free_dims.begin(); it!=free_dims.end(); ++it) {
      if (*it == mi) continue;
      tset.insert(*it);
    }
    free_dims = tset;
    tset.clear();

    Dm2 = Dm1;
    Dm1 = current;
    current.clear();

    D1 = D2;
    D2 = mi;
    mi = -1;
  } while(!free_dims.empty());

  return modified;
}

bool Spacetime::correction(int base,int sheet)
{
  // This method needs to loop over all vertices in a given sheet and find those which
  // are capable of adding another vertex at a distance of (roughly) one and which is
  // orthogonal to the vertex's current set of edges.
  int i,j,n,m,in1;
  unsigned int r;
  bool active,modified = false;
  double l,d1,d2,dbest;
  std::set<int> locus,candidates;
  Simplex s1;
  SYNARMOSMA::hash_map::const_iterator qt;
  const int nv = (signed) events.size();

  locus.insert(sheet);

  for(i=0; i<nv; ++i) {
    if (i == base) continue;
    if (!events[i].active(sheet)) continue;
    if (std::abs(events[i].deficiency) < std::numeric_limits<double>::epsilon()) continue;
    l = geometry->get_squared_distance(base,i,true);
    if (l < 3.8025 || l > 4.2025) continue;
    // See if there is a third vertex that lies between these two...
    in1 = -1;
    active = false;
    dbest = 100.0;
    for(j=0; j<nv; ++j) {
      if (j == base || j == i) continue;
      d1 = geometry->get_squared_distance(base,j,true);
      d2 = geometry->get_squared_distance(i,j,true);
      if ((d1 > 0.81 && d1 < 1.21) && (d2 > 0.81 && d2 < 1.21)) {
        l = (d1 - 1.0)*(d1 - 1.0) + (d2 - 1.0)*(d2 - 1.0);
        if (l < dbest) {
          dbest = l;
          in1 = j;
        }
        if (events[j].active()) {
          active = true;
          break;
        }
      }
      if (active) continue;
      if (in1 >= 0) {
        if (!events[in1].active(sheet)) {
#ifdef VERBOSE
          std::cout << "Restoring vertex " << in1 << " in correction between " << base << "  " << i << std::endl;
#endif
          events[in1].set_active(sheet);
          modified = true;
        }
        // Connect this vertex to i and j if necessary
        s1.initialize(base,in1,locus);
        qt = index_table[1].find(s1.vertices);
        if (qt == index_table[1].end()) {
#ifdef VERBOSE
          std::cout << "Adding edge connecting " << base << " and " << in1 << " in correction" << std::endl;
#endif
          simplices[1].push_back(s1);
          index_table[1][s1.vertices] = (signed) simplices[1].size() - 1;
          modified = true;
        }
        else {
          if (!simplices[1][qt->second].active(sheet)) {
#ifdef VERBOSE
            std::cout << "Restoring edge with key " << SYNARMOSMA::make_key(simplices[1][qt->second].vertices) << " in correction" << std::endl;
#endif
            simplices[1][qt->second].set_active(sheet);
            modified = true;
          }
        }
        s1.initialize(i,in1,locus);
        qt = index_table[1].find(s1.vertices);
        if (qt == index_table[1].end()) {
#ifdef VERBOSE
          std::cout << "Adding edge connecting " << i << " and " << in1 << " in correction" << std::endl;
#endif
          simplices[1].push_back(s1);
          index_table[1][s1.vertices] = (signed) simplices[1].size() - 1;
          modified = true;
        }
        else {
          if (!simplices[1][qt->second].active(sheet)) {
#ifdef VERBOSE
            std::cout << "Restoring edge with key " << SYNARMOSMA::make_key(simplices[1][qt->second].vertices) << " in correction" << std::endl;
#endif
            simplices[1][qt->second].set_active(sheet);
            modified = true;
          }
        }
        continue;
      }
      candidates.insert(i);
    }
  }
  if (modified) regularization(true,sheet);
  if (candidates.empty()) return false;
  n = RND->irandom(candidates);
  std::vector<double> xc,x1,x2;
  // Perform the vertex fission on the n'th vertex pair...
  geometry->get_coordinates(base,x1);
  geometry->get_coordinates(n,x2);
  for(r=0; r<geometry->dimension(); ++r) {
    l = 0.5*(x1[r] + x2[r]);
    xc.push_back(l);
  }
  // Add something to check here if this new vertex is within 0.5 of an existing vertex...
  for(i=0; i<nv; ++i) {
    if (!events[i].active()) continue;
    if (geometry->get_squared_distance(i,xc) <= 0.5) return false;
  }
#ifdef VERBOSE
  std::cout << "Adding vertex between " << base << " and " << n << std::endl;
#endif
  m = vertex_addition(xc,sheet);
  s1.initialize(base,m,locus);
  simplices[1].push_back(s1);
  index_table[1][s1.vertices] = (signed) simplices[1].size() - 1;
  s1.initialize(n,m,locus);
  simplices[1].push_back(s1);
  index_table[1][s1.vertices] = (signed) simplices[1].size() - 1;

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

bool Spacetime::reduction(int base,int sheet)
{
  const int d = vertex_dimension(base,sheet);
  if (d < 2) return false;
  int i,n,m,vx[1+d];
  std::set<int> candidates,s1;
  SYNARMOSMA::hash_map::const_iterator qt;
  const int N = (signed) simplices[d].size();

  // Find which edges this vertex possesses are used in n-simplices (n > 1) and
  // eliminate one of them...
  for(i=0; i<N; ++i) {
    if (!simplices[d][i].active(sheet)) continue;
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
  simplex_deletion(1,qt->second,sheet);
  return true;
}

bool Spacetime::contraction(int base,double l,int sheet)
{
  int i,vx[2];
  std::set<int> pool;
  const double L2 = l*l;
  const int ne = (signed) simplices[1].size();

  if (sheet == -1) {
    for(i=0; i<ne; ++i) {
      if (!simplices[1][i].active()) continue;
      simplices[1][i].get_vertices(vx);
      if (vx[0] != base && vx[1] != base) continue;
      if (simplices[1][i].sq_volume > L2) pool.insert(i);
    }
  }
  else {
    for(i=0; i<ne; ++i) {
      if (!simplices[1][i].active(sheet)) continue;
      simplices[1][i].get_vertices(vx);
      if (vx[0] != base && vx[1] != base) continue;
      if (simplices[1][i].sq_volume >= L2) pool.insert(i);
    }
  }
  if (pool.empty()) return false;
  i = RND->irandom(pool);
#ifdef VERBOSE
  std::cout << "Deleting edge with key " << SYNARMOSMA::make_key(simplices[1][i].vertices) << " in contraction" << std::endl;
#endif
  simplex_deletion(1,i,sheet);
  return true;
}

bool Spacetime::compensation_m(int base,int sheet)
{
  int i,j,vx[2];
  double l;
  std::set<int> candidates,s1;
  SYNARMOSMA::hash_map::const_iterator qt;
  const int ne = (signed) simplices[1].size();

  if (vertex_dimension(base,sheet) < 2) return false;
#ifdef VERBOSE
  std::cout << "Compensation with " << dimension(sheet) << std::endl;
#endif

  // Remove an edge from this vertex: ideally one that connects with a vertex whose degree
  // is also excessive *and* which is relatively far away (edge length >> 1)
  for(i=0; i<ne; ++i) {
    if (!simplices[1][i].active(sheet)) continue;
    simplices[1][i].get_vertices(vx);
    if (vx[0] != base && vx[1] != base) continue;
    j = (vx[0] == base) ? vx[1] : vx[0];
    if (vertex_dimension(j,sheet) < 2) continue;
    if (events[j].deficiency < -std::numeric_limits<double>::epsilon()) continue;
    l = geometry->get_squared_distance(base,j,true);
    if (l >= 0.81 && l <= 1.21) continue;
    s1.clear();
    s1.insert(base);
    s1.insert(j);
    qt = index_table[1].find(s1);
    candidates.insert(qt->second);
  }

  if (candidates.empty()) {
#ifdef VERBOSE
    std::cout << "Unable to perform dimensional trimming, no viable candidates..." << std::endl;
#endif
    return false;
  }
  j = RND->irandom(candidates);
#ifdef VERBOSE
  std::cout << "Deleting edge with key " << SYNARMOSMA::make_key(simplices[1][j].vertices) << " to reduce the complex's dimensionality" << std::endl;
#endif
  simplex_deletion(1,j,sheet);
  return true;
}

bool Spacetime::compensation_g(int base,int sheet)
{
  int i,j,vx[2];
  double l;
  std::set<int> candidates,s1;
  SYNARMOSMA::hash_map::const_iterator qt;
  const int ne = (signed) simplices[1].size();

  unsigned int sdegree = 0;
  for(i=0; i<ne; ++i) {
    if (!simplices[1][i].active(sheet)) continue;
    simplices[1][i].get_vertices(vx);
    if (base != vx[0] && base != vx[1]) continue;
    sdegree++;
  }
  if (sdegree == 2*geometry->dimension()) return false;
  bool up;
  std::set<int> locus;
  const unsigned int idegree = 2*geometry->dimension();
  const int nv = (signed) events.size();

  locus.insert(sheet);

  up = (sdegree < idegree) ? true : false;
#ifdef VERBOSE
  std::cout << "Degree trim with " << base << "  " << up << "  " << sdegree - idegree << std::endl;
#endif
  if (up) {
    // Add an edge to the vertex base: ideally one that connects with a vertex whose degree is
    // also too low and which is relatively close at hand. If the only candidates are distant
    // then do nothing and return false!
    for(i=0; i<nv; ++i) {
      if (i == base) continue;
      if (!events[i].active(sheet)) continue;
      if (edge_exists(base,i,sheet)) continue;
      l = geometry->get_squared_distance(base,i,true);
      if (l >= 0.81 && l <= 1.21) candidates.insert(i);
    }

    if (candidates.empty()) {
#ifdef VERBOSE
      std::cout << "Failure in positive compensation ..." << std::endl;
#endif
      return false;
    }
    j = RND->irandom(candidates);
    codex[sheet].vx_delta.insert(base);
    codex[sheet].vx_delta.insert(j);
    s1.clear();
    s1.insert(base);
    s1.insert(j);
    qt = index_table[1].find(s1);
    if (qt != index_table[1].end()) {
#ifdef DEBUG
      assert(!simplices[1][qt->second].active(sheet));
#endif
#ifdef VERBOSE
      std::cout << "Restoring edge with key " << SYNARMOSMA::make_key(simplices[1][qt->second].vertices) << " in positive compensation" << std::endl;
#endif
      simplices[1][qt->second].set_active(sheet);
    }
    else {
#ifdef VERBOSE
      std::cout << "Adding edge connecting " << base << " and " << j << " in positive compensation" << std::endl;
#endif
      Simplex S(base,j,locus);
      // Now add this simplex to the spacetime complex...
      simplices[1].push_back(S);
      index_table[1][S.vertices] = (signed) simplices[1].size() - 1;
    }
  }
  else {
    // Remove an edge from this vertex: ideally one that connects with a vertex whose degree
    // is also excessive *and* which is relatively far away (edge length >> 1)
    for(i=0; i<ne; ++i) {
      if (!simplices[1][i].active(sheet)) continue;
      simplices[1][i].get_vertices(vx);
      if (vx[0] != base && vx[1] != base) continue;
      j = (vx[0] == base) ? vx[1] : vx[0];
      l = geometry->get_squared_distance(base,j,true);
      if (l >= 0.81 && l <= 1.21) continue;
      s1.clear();
      s1.insert(base);
      s1.insert(j);
      qt = index_table[1].find(s1);
      candidates.insert(qt->second);
    }

    if (candidates.empty()) {
#ifdef VERBOSE
      std::cout << "Failure in negative compensation..." << std::endl;
#endif
      return false;
    }
    i = RND->irandom(candidates);
#ifdef VERBOSE
    std::cout << "Deleting edge with key " << SYNARMOSMA::make_key(simplices[1][i].vertices) << " in negative compensation" << std::endl;
#endif
    simplex_deletion(1,i,sheet);
  }
  return true;
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

int Spacetime::compression(double threshold,std::set<int>& vmodified)
{
  int i,n,nc,v[2];
  std::set<int> locus;
  std::set<int>::const_iterator it;
  double s,alpha;
  std::vector<int> candidates;
  const int nv = (signed) events.size();
  const int ne = (signed) simplices[1].size();
  const int nt = (signed) codex.size();

  for(i=0; i<ne; ++i) {
    if (!simplices[1][i].active()) continue;
    alpha = std::abs(simplices[1][i].volume);
    if (alpha > threshold) {
      s = std::exp(-(alpha - threshold));
      if (RND->drandom() > s) candidates.push_back(i);
    }
  }
  nc = (signed) candidates.size();
  for(i=0; i<nc; ++i) {
    n = candidates[i];
    simplex_deletion(1,n,-1);
    // Remove the references in events[v/w].neighbours and
    // events[v/w].entourage
    simplices[1][n].get_vertices(v);

    it = std::find(events[v[0]].neighbours.begin(),events[v[0]].neighbours.end(),v[1]);
    if (it != events[v[0]].neighbours.end()) events[v[0]].neighbours.erase(it);
    it = std::find(events[v[1]].neighbours.begin(),events[v[1]].neighbours.end(),v[0]);
    if (it != events[v[1]].neighbours.end()) events[v[1]].neighbours.erase(it);

    it = std::find(events[v[0]].entourage.begin(),events[v[0]].entourage.end(),n);
    if (it != events[v[0]].entourage.end()) events[v[0]].entourage.erase(it);
    it = std::find(events[v[1]].entourage.begin(),events[v[1]].entourage.end(),n);
    if (it != events[v[1]].entourage.end()) events[v[1]].entourage.erase(it);

    vmodified.insert(v[0]);
    vmodified.insert(v[1]);
  }
  if (!(connected(-1))) {
    // We need to add enough (short!) edges to keep the spacetime connected...
    int j,k,nsedge,ncomp;
    double l,tlength;
    std::set<int>::const_iterator jt,kt;
    std::set<int> linked;
    std::vector<int> component,sedge;
    std::vector<std::tuple<int,int,double> > connect;
    Simplex S;

    ncomp = component_analysis(component,-1);
#ifdef VERBOSE
    std::cout << "There are " << ncomp << " components in this spacetime." << std::endl;
#endif
    std::vector<int>* cvertex = new std::vector<int>[ncomp];
    for(i=0; i<nv; ++i) {
      if (!events[i].active()) continue;
      cvertex[component[i]].push_back(i);
    }

    // Now we want to find the shortest possible edge that connects all of these
    // components together, first the brute force solution
    for(i=0; i<ncomp; ++i) {
      for(j=1+i; j<ncomp; ++j) {
        l = minimize_lengths(cvertex[i],cvertex[j],v);
        connect.push_back(std::tuple<int,int,double>(v[0],v[1],l));
      }
    }
    n = (signed) connect.size();
    bool used[n];
    for(i=0; i<n; ++i) {
      used[i] = false;
    }
    std::sort(connect.begin(),connect.end(),SYNARMOSMA::tuple_predicate);
    // Assuming the array "connect" has been sorted in ascending order for the last
    // element...
    j = std::get<0>(connect[0]);
    k = std::get<1>(connect[0]);
    sedge.push_back(j);
    sedge.push_back(k);
    linked.insert(component[j]);
    linked.insert(component[k]);
    tlength = std::get<2>(connect[0]);
    used[0] = true;
    do {
      for(i=1; i<n; ++i) {
        if (used[i]) continue;
        j = component[std::get<0>(connect[i])];
        k = component[std::get<1>(connect[i])];
        jt = std::find(linked.begin(),linked.end(),j);
        kt = std::find(linked.begin(),linked.end(),k);
        if (jt != linked.end() && kt != linked.end()) continue;
        if (jt == linked.end() && kt == linked.end()) continue;
        sedge.push_back(std::get<0>(connect[i]));
        sedge.push_back(std::get<1>(connect[i]));
        tlength += std::get<2>(connect[i]);
        linked.insert(j);
        linked.insert(k);
        break;
      }
    } while((signed) linked.size() < ncomp);
    nsedge = (signed) sedge.size()/2;
#ifdef DEBUG
    assert((nsedge+1) == ncomp);
#endif
#ifdef VERBOSE
    std::cout << "Adding " << nsedge << " edge(s) to reconnect the spacetime complex..." << std::endl;
#endif
    // Now add the edges in sedge...
    for(i=0; i<nsedge; ++i) {
      j = sedge[2*i];
      k = sedge[2*i+1];
      l = RND->irandom(nt);
      locus.insert(l);
      S.initialize(j,k,locus);
      simplices[1].push_back(S);
      index_table[1][S.vertices] = (signed) simplices[1].size() - 1;
      events[j].neighbours.insert(k);
      events[k].neighbours.insert(j);
      vmodified.insert(k);
      vmodified.insert(j);
      locus.erase(l);
    }
    nc -= nsedge;
#ifdef DEBUG
    assert(connected(-1));
#endif
    delete[] cvertex;
  }
  return nc;
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

bool Spacetime::foliation_x(int base,int sheet)
{
  int i,p,n1,n2,vx[2];
  std::set<int> candidates,S;
  SYNARMOSMA::hash_map::iterator qt;
  const int ne = (signed) simplices[1].size();

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
  S.insert(n1); S.insert(n2);
  qt = index_table[1].find(S);
  if (qt == index_table[1].end()) return false;
  if (!simplices[1][qt->second].active(sheet)) return false;
  simplex_deletion(1,qt->second,sheet);
  return true;
}

bool Spacetime::amputation(int base,double tolerance,int sheet)
{
  int i,j,n,p,vx[2];
  std::set<int> candidates;
  const int ne = (signed) simplices[1].size();
  const int ulimit = dimension(sheet);

  if (events[base].deficiency > tolerance) candidates.insert(base);
  for(i=0; i<ne; ++i) {
    if (!simplices[1][i].active(sheet)) continue;
    simplices[1][i].get_vertices(vx);
    if (base != vx[0] && base != vx[1]) continue;
    n = (vx[0] == base) ? vx[1] : vx[0];
    if (events[n].deficiency > tolerance) candidates.insert(n);
  }
  
  if (candidates.empty()) return false;
  n = RND->irandom(candidates);

#ifdef VERBOSE
  std::cout << "Amputating vertex " << n << " and all its dependent simplices" << std::endl;
#endif
  // Delete this vertex and all its edges...
  events[n].set_inactive(sheet);
  codex[sheet].vx_delta.insert(n);
  for(i=ulimit; i>=1; --i) {
    p = (signed) simplices[i].size();
    for(j=0; j<p; ++j) {
      if (!simplices[i][j].active(sheet)) continue;
      if (simplices[i][j].contains(n)) simplex_deletion(i,j,sheet);
    }
  }
  if (!events[n].active()) {
    if (!events[n].zero_energy()) {
      events[base].increment_energy(events[n].get_energy());
      events[n].nullify_energy();
    }
  }
  return true;
}

bool Spacetime::fusion_x(int base,double tolerance,int sheet)
{
  int i,j;
  std::vector<std::pair<int,int> > candidates;
  const int nv = (signed) events.size();

  for(i=0; i<nv; ++i) {
    if (i == base) continue;
    if (!events[i].active(sheet)) continue;
    if (events[i].deficiency < std::numeric_limits<double>::epsilon()) continue;
    if (geometry->get_squared_distance(base,i,true) > tolerance) continue;
    candidates.push_back(std::pair<int,int>(base,i));
  }
  if (candidates.empty()) return false;
  i = (signed) candidates.size();
  j = RND->irandom(i);
#ifdef VERBOSE
  std::cout << "Fusing deficient vertices: " << candidates[j].second << " => " << candidates[j].first << std::endl;
#endif
  vertex_fusion(candidates[j].first,candidates[j].second,sheet);
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

void Spacetime::regularization(bool minimal,int sheet)
{
  if (!minimal) simplicial_implication(sheet);

  compute_neighbours();
  if (connected(sheet)) {
    compute_entourages(sheet);
    return;
  }

  int i,j,k,nc,v1,v2,n1,n2,loc = -1;
  double l,mdelta;
  unsigned int max = 0;
  Simplex S;
  std::set<int> locus,colours;
  std::vector<int> component;
  SYNARMOSMA::hash_map::const_iterator qt;
  const int nv = (signed) events.size();
  const int nt = (signed) codex.size();

  for(i=0; i<nt; ++i) {
    if (codex[i].active) colours.insert(i);
  }
  nc = component_analysis(component,sheet);
#ifdef VERBOSE
  std::cout << "There are " << nc << " components in this spacetime." << std::endl;
#endif
  // So there are nc (>1) components in this graph, we will unite them in the
  // most minimal fashion, using (nc-1) edges
  // November 2, 2012: We need to consider ensuring that the minimalist linking of these
  // distinct components is also minimalist in a geometric sense (i.e. the shortest possible
  // edge) and if a "component" is nothing more than an isolated vertex with zero energy
  // then drop the vertex (set ubiquity = 1).
  std::vector<int>* cvertex = new std::vector<int>[nc];
  if (sheet == -1) {
    for(i=0; i<nv; ++i) {
      if (!events[i].active()) continue;
      cvertex[component[i]].push_back(i);
    }
  }
  else {
    for(i=0; i<nv; ++i) {
      if (!events[i].active(sheet)) continue;
      cvertex[component[i]].push_back(i);
    }
  }

  for(i=0; i<nc; ++i) {
    if (cvertex[i].size() > max) {
      max = cvertex[i].size();
      loc = i;
    }
  }
  for(i=0; i<nc; ++i) {
    if (cvertex[i].size() == 1) {
      // Isolated vertex, let's see what its energy is...
      if (events[cvertex[i][0]].zero_energy()) {
        // Eliminate this vertex altogether...
#ifdef VERBOSE
        std::cout << "Eliminating isolated vertex " << cvertex[i][0] << std::endl;
#endif
        for(j=0; j<nt; ++j) {
          events[cvertex[i][0]].set_inactive(j);
        }
        cvertex[i].erase(cvertex[i].begin());
      }
    }
  }
  n1 = (signed) cvertex[loc].size();
  if (sheet == -1) {
    for(i=0; i<nc; ++i) {
      if (i == loc) continue;
      if (cvertex[i].empty()) continue;
      // We need to choose the closest pair of vertices, not just a random pair
      n2 = (signed) cvertex[i].size();
      mdelta = std::numeric_limits<double>::infinity();
      v1 = -1;
      v2 = -1;
      for(j=0; j<n1; ++j) {
        for(k=0; k<n2; ++k) {
          l = geometry->get_squared_distance(cvertex[loc][j],cvertex[i][k],true);
          if (l < mdelta) {
            mdelta = l;
            v1 = cvertex[loc][j];
            v2 = cvertex[i][k];
          }
        }
      }
#ifdef DEBUG
      assert(v1 > -1);
      assert(v2 > -1);
#endif
      j = RND->irandom(colours);
      locus.insert(j);
      S = Simplex(v1,v2,locus);
      locus.erase(j);
      qt = index_table[1].find(S.vertices);
      if (qt == index_table[1].end()) {
#ifdef VERBOSE
        std::cout << "Adding edge connecting " << v1 << " and " << v2 << " to maintain the complex's connectedness" << std::endl;
#endif
        simplices[1].push_back(S);
        index_table[1][S.vertices] = (signed) simplices[1].size() - 1;
      }
      else {
#ifdef DEBUG
        assert(!simplices[1][qt->second].active());
#endif
#ifdef VERBOSE
        std::cout << "Restoring the simplex with key " << SYNARMOSMA::make_key(simplices[1][qt->second].vertices) << " to maintain the complex's connectedness" << std::endl;
#endif
        simplices[1][qt->second].set_active(j);
      }
      // Need to ensure that v1 and v2 are also coloured by sheet...
      if (!events[v1].active(j)) events[v1].set_active(j); 
      if (!events[v2].active(j)) events[v2].set_active(j); 
    }
  }
  else {
    locus.insert(sheet);
    for(i=0; i<nc; ++i) {
      if (i == loc) continue;
      if (cvertex[i].empty()) continue;
      n2 = (signed) cvertex[i].size();
      mdelta = std::numeric_limits<double>::infinity();
      v1 = -1;
      v2 = -1;
      for(j=0; j<n1; ++j) {
        for(k=0; k<n2; ++k) {
          l = geometry->get_squared_distance(cvertex[loc][j],cvertex[i][k],true);
          if (l < mdelta) {
            mdelta = l;
            v1 = cvertex[loc][j];
            v2 = cvertex[i][k];
          }
        }
      }
#ifdef VERBOSE
      std::cout << "Linking together " << v1 << " and " << v2 << " at a distance of " << mdelta << std::endl;
#endif
#ifdef DEBUG
      assert(v1 > -1);
      assert(v2 > -1);
#endif
      S = Simplex(v1,v2,locus);
      qt = index_table[1].find(S.vertices);
      if (qt == index_table[1].end()) {
        simplices[1].push_back(S);
        index_table[1][S.vertices] = (signed) simplices[1].size() - 1;
      }
      else {
        if (!simplices[1][qt->second].active(sheet)) simplices[1][qt->second].set_active(sheet);
      }
    }
  }
  compute_neighbours();
#ifdef DEBUG
  assert(connected(sheet));
#endif
  compute_entourages(sheet);
  delete[] cvertex;
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

void Spacetime::superposition_fusion(std::set<int>& vmodified)
{
  int i,j,m,nf,v1,v2,nfail = 0,nfused = 0;
  double delta,pfusion = 0.0,alpha = 0.1;
  std::vector<std::pair<int,int> > candidates;
  std::vector<int> modified;
  const double na = double(cardinality(0,-1));
  const int nv = (signed) events.size();
  const int ulimit = dimension(-1);

  // Vertex fusion - if two vertices are close enough together, they
  // should coalesce.
  for(i=0; i<nv; ++i) {
    modified.push_back(0);
    if (!events[i].active()) continue;
    for(j=1+i; j<nv; ++j) {
      if (!events[j].active()) continue;
      delta = geometry->get_squared_distance(i,j,false);
      if (delta < alpha) candidates.push_back(std::pair<int,int>(i,j));
    }
  }
  if (candidates.empty()) return;

  // Now we need to carry out the fusions...
  nf = (signed) candidates.size();
  do {
    i = RND->irandom(nf);
    v1 = candidates[i].first;
    v2 = candidates[i].second;
    if (modified[v1] || modified[v2]) {
      nfail++;
      continue;
    }
    modified[v1] = 1;
    modified[v2] = 1;
#ifdef VERBOSE
    std::cout << "Fusing vertices: " << v2 << " => " << v1 << " via superposition" << std::endl;
#endif
    vertex_fusion(v1,v2,-1);
    vmodified.insert(v1);
    nfused++;
    pfusion = double(nfused)/na;
  } while(pfusion < 0.05 && nfused < nf && nfail < 20);

  // Then recalculate the index_table hash map..
  for(i=1; i<=ulimit; ++i) {
    index_table[i].clear();
    m = (signed) simplices[i].size();
    for(j=0; j<m; ++j) {
      simplices[i][j].entourage.clear();
      index_table[i][simplices[i][j].vertices] = j;
    }
  }
  for(i=0; i<nv; ++i) {
    events[i].entourage.clear();
  }

  // Then recalculate the entourages...
  compute_entourages(-1);
  compute_neighbours();
#ifdef DEBUG
  assert(consistent(-1));
#endif
}

void Spacetime::superposition_fission(std::set<int>& vmodified)
{
  int i,j,k,nc;
  bool change;
  std::set<int> locus;
  std::set<int>::const_iterator it;
  Event vx;
  Simplex S;
  const int nv = (signed) events.size();
  const int nt = (signed) codex.size();

  // Finally, the opposite possibility - that a given vertex might undergo
  // spontaneous fission, creating a new vertex in its immediate vicinity...
  nc = nv;
  do {
    change = false;
    for(i=0; i<nc; ++i) {
      if (!events[i].active()) continue;
      if (RND->drandom() > 0.01) continue;
      vx = events[i];      
      for(j=0; j<nt; ++j) {
        if (RND->drandom() > 0.9) {
          if (vx.active(j)) {
            vx.set_inactive(j);
          }
          else {
            vx.set_active(j);
          }
        }
      }
      // Check that vx.ubiquity > 1
      if (!vx.active()) vx.set_active(RND->irandom(nt));
      for(it=vx.neighbours.begin(); it!=vx.neighbours.end(); ++it) {
        j = *it;
        events[j].neighbours.insert(nc);
        vmodified.insert(j);
        for(k=0; k<nt; ++k) {
          if (RND->drandom() > 0.33) locus.insert(k);
        }
        if (locus.empty()) locus.insert(RND->irandom(nt));
        S.initialize(nc,j,locus);
        simplices[1].push_back(S);
        index_table[1][S.vertices] = (signed) simplices[1].size() - 1;
        S.clear();
      }
      vmodified.insert(nc);
      geometry->vertex_addition(i);
      events.push_back(vx);
      change = true;
      break;
    }
    if (change) nc++;
#ifdef VERBOSE
    std::cout << "Created vertex " << nc - 1 << " from " << i << std::endl;
#endif
  } while((nc - nv) < 4);
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

bool Spacetime::deflation(int base,int sheet)
{
  int d = vertex_dimension(base,sheet);
  if (d < 2) return false;
  int i,n,dw = 1;
  std::set<int> candidates;
  if (d > 2) dw = RND->irandom(1,d);
  const int m = (signed) simplices[dw].size();
  for(i=0; i<m; ++i) {
    if (!simplices[dw][i].active(sheet)) continue;
    if (simplices[dw][i].contains(base)) candidates.insert(i);
  }
  n = RND->irandom(candidates);
#ifdef VERBOSE
  std::cout << "Deflating a " << d << "-simplex with base " << base << " to a " << dw << "-simplex." << std::endl;
#endif
  simplex_deletion(dw,n,sheet);
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

