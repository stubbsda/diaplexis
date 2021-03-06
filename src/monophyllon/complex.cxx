#include "complex.h"

using namespace DIAPLEXIS;

template<class kind>
Complex<kind>::Complex()
{
  allocate();
}

template<class kind>
Complex<kind>::~Complex()
{
  delete[] simplices;
  delete[] index_table;
  delete H;
  delete pi1;
  delete RND;
}

template<class kind>
void Complex<kind>::allocate()
{
  H = new SYNARMOSMA::Homology(SYNARMOSMA::Homology::Field::mod2,SYNARMOSMA::Homology::Method::native);
  pi1 = new SYNARMOSMA::Homotopy;
  simplices = new std::vector<Simplex>[1+Complex<kind>::ND];
  index_table = new SYNARMOSMA::hash_map[1+Complex<kind>::ND];
  RND = new SYNARMOSMA::Random;
}

template<class kind>
void Complex<kind>::clear()
{
  H->clear();
  pi1->clear();
  events.clear();
  for(int i=1; i<=Complex<kind>::ND; ++i) {
    simplices[i].clear();
    index_table[i].clear();
  }
}

template<class kind>
void Complex<kind>::compute_entourages()
{
  int i,j,k,ns,vx[2];
  std::set<int> s,v;
  std::set<int>::const_iterator it;
  const int ulimit = dimension();
  SYNARMOSMA::hash_map::const_iterator qt;

  // What about removing items from the entourage of a d-simplex, when this item has
  // changed its ubiquity?
  for(i=1; i<=Complex<kind>::ND; ++i) {
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

template<class kind>
void Complex<kind>::compute_neighbours()
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

template<class kind>
bool Complex<kind>::consistent() const
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

  for(i=Complex<kind>::ND; i>=2; i--) {
    n = (signed) simplices[i].size();
    for(j=0; j<n; ++j) {
      if (!simplices[i][j].active) continue;
        for(it=simplices[i][j].entourage.begin(); it!=simplices[i][j].entourage.end(); ++it) {
          if (!simplices[i+1][*it].active) {
            std::cout << "Error with entourage ubiquity: " << i << "  " << j << "  " << *it << "  " << ulimit << std::endl;
          }
        }
      for(k=0; k<1+i; ++k) {
        qt = index_table[i-1].find(simplices[i][j].faces[k]);
        if (!simplices[i-1][qt->second].active) {
          std::cout << "Error with entailment: " << i << "  " << j << "  " << SYNARMOSMA::make_key(simplices[i][j].faces[k]) << "  " << ulimit << std::endl;
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
          std::cout << "Error with entailment for vertex: " << l << std::endl;
          return false;
        }
      }
    }
  }
  for(i=0; i<(signed) simplices[1].size(); ++i) {
    if (!simplices[1][i].active) continue;
    for(it=simplices[1][i].entourage.begin(); it!=simplices[1][i].entourage.end(); ++it) {
      if (!simplices[2][*it].active) {
        std::cout << "Error with entourage ubiquity: 1 " << i << "  " << *it << "  " << ulimit << std::endl;
      }
    }
    simplices[1][i].get_vertices(vx);
    if (vx[0] < 0 || vx[0] >= nv) return false;
    if (vx[1] < 0 || vx[0] >= nv) return false;
    if (!events[vx[0]].active) {
      std::cout << "Error with entailment for vertex: " << vx[0] << std::endl;
      return false;
    }
    if (!events[vx[1]].active) {
      std::cout << "Error with entailment for vertex: " << vx[1] << std::endl;
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
    for(it=events[i].entourage.begin(); it!=events[i].entourage.end(); ++it) {
      if (!simplices[1][*it].active) {
        std::cout << "Error with entourage ubiquity: 0 " << i << "  " << *it << std::endl;
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
  for(i=1; i<=Complex<kind>::ND; ++i) {
    n = (signed) simplices[i].size();
    for(j=0; j<n; ++j) {
      S = simplices[i][j].vertices;
      for(k=0; k<n; ++k) {
        if (k == j) continue;
        if (S == simplices[i][k].vertices) {
          std::cout << "Illegal " << i << "-simplex duplication with " << SYNARMOSMA::make_key(S) << " for indices " << k << " and " << j << std::endl;
          return false;
        }
      }
    }
  }
  return true;
}

template<class kind>
int Complex<kind>::circuit_rank() const
{
  int output = cardinality(1) - cardinality(0);
  if (connected()) {
    return 1 + output;
  }
  std::vector<int> components;
  int n = component_analysis(components);
  return n + output;
}

template<class kind>
int Complex<kind>::euler_characteristic() const
{
  int i,pf = 1,chi = 0;
  const int D = dimension();
  for(i=0; i<=D; ++i) {
    chi += pf*cardinality(i);
    pf *= -1;
  }
  return chi;
}

template<class kind>
int Complex<kind>::dimension() const
{
  int i,j;

  for(i=Complex<kind>::ND; i>0; i--) {
    for(j=0; j<(signed) simplices[i].size(); ++j) {
      if (simplices[i][j].active) return i;
    }
  }
  for(i=0; i<(signed) events.size(); ++i) {
    if (events[i].active) return 0;
  }
  return -1;
}

template<class kind>
int Complex<kind>::total_dimension() const
{
  int i,d,sum = 0;

  for(i=0; i<(signed) events.size(); ++i) {
    d = vertex_dimension(i);
    if (d > 0) sum += d;
  }
  return sum;
}

template<class kind>
int Complex<kind>::structural_index() const
{
  int i,l,sum = 0;
  unsigned int j,n = events.size();
  const int d = dimension();

  for(j=0; j<n; ++j) {
    if (events[j].active) sum++;
  }
  for(i=1; i<d; ++i) {
    l = 0;
    n = simplices[i].size();
    for(j=0; j<n; ++j) {
      if (simplices[i][j].active) l++;
    }
    sum += (1+i)*l;
  }
  return sum;
}

template<class kind>
int Complex<kind>::component_analysis(std::vector<int>& component) const
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

template<class kind>
void Complex<kind>::compute_global_topology(bool high_memory)
{
  // To calculate the global deficiency, we need to compute the Betti numbers and
  // the fundamental group, for the total spacetime, operations that are serial...
  if (cardinality(0) < 2) return;

  SYNARMOSMA::Nexus* NX = new SYNARMOSMA::Nexus;

  compute_global_nexus(NX);

  H->compute(NX);
  if (high_memory) pi1->compute(NX);
  // Finally, the pseudomanifold and orientability properties
  pseudomanifold = NX->pseudomanifold(&boundary);
  if (pseudomanifold) orientable = NX->orientable();

  delete NX;
}

template<class kind>
void Complex<kind>::write_topology() const
{
  int i,j,n;
  std::vector<int> vx;
  SYNARMOSMA::UINT64 q;
  const int ns = cardinality(0);
  const int ulimit = dimension();

  for(i=0; i<(signed) events.size(); ++i) {
    if (!events[i].active) continue;
    simplex_membership(i,vx);
    std::cout << i+1 << ": [";
    for(j=1; j<ulimit; ++j) {
      std::cout << vx[j-1] << ",";
    }
    std::cout << vx[ulimit-1] << "]" << std::endl;
  }
  for(i=Complex<kind>::ND; i>=1; i--) {
    n = cardinality(i);
    if (n > 0) {
      q = int(SYNARMOSMA::binomial(ns,i+1));
      std::cout << "There are " << n << " (" << q << ") " << i << "-simplices in this complex." << std::endl;
    }
  }

  std::cout << "There are " << ns << " (" << ns << ") 0-simplices in this complex." << std::endl;
  std::cout << "The Euler characteristic is " << euler_characteristic() << std::endl;
}

template<class kind>
void Complex<kind>::write_incastrature(const std::string& filename) const
{
  // A method that generates the Hasse diagram corresponding to the
  // complex's simplicial structure, with the diagram stored in the
  // PDF "hasse.pdf"; the method assumes that the Graphviz library
  // has been installed on the system.
  int i,j,k;

  std::ofstream s(filename,std::ios::trunc);

  s << "digraph G {" << std::endl;

  for(i=Complex<kind>::ND; i>0; i--) {
    for(j=0; j<(signed) simplices[i].size(); ++j) {
      if (!simplices[i][j].active) continue;
      for(k=0; k<1+i; ++k) {
        s << "  \"" << SYNARMOSMA::make_key(simplices[i][j].vertices) << "\" -> \"" << SYNARMOSMA::make_key(simplices[i][j].faces[k]) << "\";" << std::endl;
      }
    }
  }
  for(i=0; i<(signed) events.size(); ++i) {
    if (!events[i].active) continue;
    s << "  \"" << i << "\" -> \"NULL\";" << std::endl;
  }
  s << "}" << std::endl;
  s.close();
}

template<class kind>
bool Complex<kind>::energy_check() const
{
  bool output = true;
  const int nv = (signed) events.size();

  for(int i=0; i<nv; ++i) {
    if (!events[i].zero_energy()) {
      if (!events[i].active) {
#ifdef VERBOSE
        std::cout << "An inactive event has non-zero energy! " << i << "  " << events[i].get_energy() << std::endl;
#endif
        output = false;
      }
    }
  }
  return output;
}

template<class kind>
void Complex<kind>::get_energy_values(std::vector<double>& output) const
{
  int i;
  const int nv = (signed) events.size();

  output.clear();

  for(i=0; i<nv; ++i) {
    if (!events[i].active) continue;
    output.push_back(events[i].get_energy());
  }
}

template<class kind>
void Complex<kind>::get_deficiency_values(std::vector<double>& output) const
{
  int i;
  const int nv = (signed) events.size();

  output.clear();

  for(i=0; i<nv; ++i) {
    if (!events[i].active) continue;
    output.push_back(events[i].deficiency);
  }
}

template<class kind>
void Complex<kind>::write_event_data(int v) const
{
  if (v < 0 || v >= (signed) events.size()) return;
  std::cout << "For event " << v << " we have:" << std::endl;
  std::cout << "    Incept = " << events[v].incept << std::endl;
  std::cout << "    Structural deficiency = " << events[v].deficiency << std::endl;
  std::cout << "    Energy = " << events[v].get_energy() << std::endl;
  std::cout << "    Obliquity = " << events[v].obliquity << std::endl;
  std::cout << std::endl;
}

template<class kind>
void Complex<kind>::write_graph(const std::string& filename) const
{
  int i,j,vx[2],N1 = 0;
  double E;
  std::vector<int> offset;
  const int nv = (signed) events.size();
  const int ne = (signed) simplices[1].size();

  for(i=0; i<nv; ++i) {
    offset.push_back(-1);
  }

  std::ofstream s(filename,std::ios::out | std::ios::trunc | std::ios::binary);
  i = cardinality(0);
  s.write((char*)(&i),sizeof(int));
  i = cardinality(1);
  s.write((char*)(&i),sizeof(int));
  // First calculate the appropriate vertex offset and write out the energy...
  for(i=0; i<nv; ++i) {
    if (!events[i].active) continue;
    offset[i] = N1;
    N1++;
    E = events[i].get_energy();
    s.write((char*)(&E),sizeof(double));
  }
  // Finally the edges...
  for(i=0; i<ne; ++i) {
    if (!simplices[1][i].active) continue;
    simplices[1][i].get_vertices(vx);
    j = offset[vx[0]];
    s.write((char*)(&j),sizeof(int));
    j = offset[vx[1]];
    s.write((char*)(&j),sizeof(int));
  }
  s.close();
}

template<class kind>
int Complex<kind>::serialize(std::ofstream& s) const
{
  int i,j,n,output = 0;
  unsigned long q;

  s.write((char*)(&pseudomanifold),sizeof(bool)); output += sizeof(bool);
  s.write((char*)(&boundary),sizeof(bool)); output += sizeof(bool);
  s.write((char*)(&orientable),sizeof(bool)); output += sizeof(bool);

  // Write the initial random number seed to disk...
  q = RND->get_seed();
  s.write((char*)(&q),sizeof(long)); output += sizeof(long);

  n = (signed) events.size();
  s.write((char*)(&n),sizeof(int)); output += sizeof(int);
  for(i=0; i<n; ++i) {
    output += events[i].serialize(s);
  }

  for(i=1; i<=Complex<kind>::ND; ++i) {
    n = (signed) simplices[i].size();
    s.write((char*)(&n),sizeof(int)); output += sizeof(int);
    for(j=0; j<n; ++j) {
      output += simplices[i][j].serialize(s);
    }
  }
  // Now the algebraic properties...
  output += H->serialize(s);
  output += pi1->serialize(s);

  return output;
}

template<class kind>
int Complex<kind>::deserialize(std::ifstream& s)
{
  int i,j,n,output = 0;
  unsigned long q;
  Event<kind> v;
  Simplex S;

  clear();

  s.read((char*)(&pseudomanifold),sizeof(bool)); output += sizeof(bool);
  s.read((char*)(&boundary),sizeof(bool)); output += sizeof(bool);
  s.read((char*)(&orientable),sizeof(bool)); output += sizeof(bool);

  s.read((char*)(&q),sizeof(long)); output += sizeof(long);
  RND->set_seed(q);

  s.read((char*)(&n),sizeof(int)); output += sizeof(int);
  for(i=0; i<n; ++i) {
    output += v.deserialize(s);
    events.push_back(v);
  }

  for(i=1; i<=Complex<kind>::ND; ++i) {
    s.read((char*)(&n),sizeof(int)); output += sizeof(int);
    for(j=0; j<n; ++j) {
      output += S.deserialize(s);
      simplices[i].push_back(S);
    }
  }

  for(i=1; i<=Complex<kind>::ND; ++i) {
    for(j=0; j<(signed) simplices[i].size(); ++j) {
      index_table[i][simplices[i][j].vertices] = j;
    }
  }
  compute_entourages();

  // Now the algebraic properties...
  output += H->deserialize(s);
  output += pi1->deserialize(s);

  return output;
}

template<class kind>
double Complex<kind>::total_energy() const
{
  unsigned int i,nv = events.size();
  double sum = 0.0;

  for(i=0; i<nv; ++i) {
    if (events[i].active) sum += events[i].get_energy();
  }
  return sum;
}

template<class kind>
int Complex<kind>::simplex_embedding(int d,int n) const
{
#ifdef DEBUG
  assert(d > 1);
#endif
  int N = 1 + d;  
  int i,j,info,vx[N],ns = 0,nt = 0,p = 0,nwork = 3*N - 1;
  double delta,sum,s1,s2,c,dmatrix[N*N],A[N*N],a[N],b[N],w[N],work[nwork];
  char jtype = 'V',uplo = 'U';
  std::set<int> S;
  std::set<int>::const_iterator it;
  SYNARMOSMA::hash_map::const_iterator qt;

  simplices[d][n].get_vertices(vx);

  for(i=0; i<N*N; ++i) {
    dmatrix[i] = 0.0;
  }

  // First we need to check if all the edges of this simplex are
  // of the same type
  for(i=0; i<N; ++i) {
    for(j=i+1; j<N; ++j) {
      S.clear();
      S.insert(vx[i]);
      S.insert(vx[j]);
      qt = index_table[1].find(S);
      delta = simplices[1][qt->second].volume;
      dmatrix[N*i+j] = delta;
      dmatrix[N*j+i] = delta;
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
    for(i=0; i<N*N; ++i) {
      dmatrix[i] = -dmatrix[i];
    }
  }

  for(i=0; i<N*N; ++i) {
    dmatrix[i] = -0.5*dmatrix[i]*dmatrix[i];
  }

  sum = 0.0;
  for(i=0; i<N; ++i) {
    s1 = 0.0;
    s2 = 0.0;
    for(j=0; j<N; ++j) {
      sum += dmatrix[N*i+j];
      s1 += dmatrix[N*i+j];
      s2 += dmatrix[N*j+i];
    }
    a[i] = s1/double(N);
    b[i] = s2/double(N);
  }
  c = sum/double(N*N);
  for(i=0; i<N; ++i) {
    for(j=0; j<N; ++j) {
      A[N*i+j] = dmatrix[N*i+j] - a[i] - b[j] + c;
    }
  }
  dsyev_(&jtype,&uplo,&N,A,&N,w,work,&nwork,&info);
#ifdef DEBUG
  assert(info == 0);
#endif
  for(i=0; i<N; ++i) {
    if (std::abs(w[i]) > std::numeric_limits<double>::epsilon()) p++;
  }
  return p;
}

template<class kind>
void Complex<kind>::determine_flexible_edges(std::vector<int>& flexible_edge)
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
