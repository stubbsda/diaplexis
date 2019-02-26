#include "complex.h"

using namespace DIAPLEXIS;

const int Complex::topological_radius;
const int Complex::ND;

Complex::Complex()
{
  allocate();
}

Complex::Complex(const Complex& source)
{

}

Complex::~Complex()
{
  delete[] simplices;
  delete[] index_table;
  delete H;
  delete pi1;
  delete RND;
}

void Complex::allocate()
{
  H = new SYNARMOSMA::Homology(SYNARMOSMA::Homology::Field::mod2,SYNARMOSMA::Homology::Method::native);
  pi1 = new SYNARMOSMA::Homotopy;
  simplices = new std::vector<Simplex>[1+Complex::ND];
  index_table = new SYNARMOSMA::hash_map[1+Complex::ND];
  RND = new SYNARMOSMA::Random;
}

void Complex::clear()
{
  H->clear();
  pi1->clear();
  events.clear();
  for(int i=1; i<=Complex::ND; ++i) {
    simplices[i].clear();
    index_table[i].clear();
  }
}

void Complex::compute_global_topology(bool high_memory)
{
  // To calculate the global deficiency, we need to compute the Betti numbers and
  // the fundamental group, for the total spacetime, operations that are serial...
  SYNARMOSMA::Nexus* NX = new SYNARMOSMA::Nexus;

  compute_global_nexus(NX);

  H->compute(NX);
  if (high_memory) pi1->compute(NX);
  // Finally, the pseudomanifold and orientability properties
  pseudomanifold = NX->pseudomanifold(&boundary);
  if (pseudomanifold) orientable = NX->orientable();

  delete NX;
}

void Complex::write_topology() const
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
  for(i=Complex::ND; i>=1; i--) {
    n = cardinality(i);
    if (n > 0) {
      q = int(SYNARMOSMA::binomial(ns,i+1));
      std::cout << "There are " << n << " (" << q << ") " << i << "-simplices in this complex." << std::endl;
    }
  }

  std::cout << "There are " << ns << " (" << ns << ") 0-simplices in this complex." << std::endl;
  std::cout << "The Euler characteristic is " << euler_characteristic() << std::endl;
}

void Complex::write_incastrature(const std::string& filename) const
{
  // A method that generates the Hasse diagram corresponding to the
  // complex's simplicial structure, with the diagram stored in the
  // PDF "hasse.pdf"; the method assumes that the Graphviz library
  // has been installed on the system.
  int i,j,k;

  std::ofstream s(filename,std::ios::trunc);

  s << "digraph G {" << std::endl;

  for(i=Complex::ND; i>0; i--) {
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

bool Complex::energy_check() const
{
  bool output = true;
  const int nv = (signed) events.size();

  for(int i=0; i<nv; ++i) {
    if (!events[i].zero_energy()) {
      if (!events[i].active) {
        std::cout << "Potential problem here: " << i << "  " << events[i].get_energy() << std::endl;
        output = false;
      }
    }
  }
  return output;
}

void Complex::get_energy_values(std::vector<double>& output) const
{
  int i;
  const int nv = (signed) events.size();

  output.clear();

  for(i=0; i<nv; ++i) {
    if (!events[i].active) continue;
    output.push_back(events[i].get_energy());
  }
}

void Complex::get_deficiency_values(std::vector<double>& output) const
{
  int i;
  const int nv = (signed) events.size();

  output.clear();

  for(i=0; i<nv; ++i) {
    if (!events[i].active) continue;
    output.push_back(events[i].deficiency);
  }
}

void Complex::write_vertex_data(int v) const
{
  if (v < 0 || v >= (signed) events.size()) return;
  std::cout << "For vertex " << v << " we have:" << std::endl;
  std::cout << "    Incept = " << events[v].incept << std::endl;
  std::cout << "    Structural deficiency = " << events[v].deficiency << std::endl;
  std::cout << "    Energy = " << events[v].get_energy() << std::endl;
  std::cout << "    Orthogonality = " << events[v].obliquity << std::endl;
  std::cout << std::endl;
}

void Complex::write_graph(const std::string& filename) const
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

int Complex::serialize(std::ofstream& s) const
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

  for(i=1; i<=Complex::ND; ++i) {
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

int Complex::deserialize(std::ifstream& s)
{
  int i,j,n,output = 0;
  unsigned long q;
  Event v;
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

  for(i=1; i<=Complex::ND; ++i) {
    s.read((char*)(&n),sizeof(int)); output += sizeof(int);
    for(j=0; j<n; ++j) {
      output += S.deserialize(s);
      simplices[i].push_back(S);
    }
  }

  for(i=1; i<=Complex::ND; ++i) {
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

double Complex::total_energy() const
{
  unsigned int i,nv = events.size();
  double sum = 0.0;

  for(i=0; i<nv; ++i) {
    if (events[i].active) sum += events[i].get_energy();
  }
  return sum;
}

int Complex::simplex_embedding(int d,int n) const
{
  int i,j,info,vx[1+d],ns = 0,nt = 0,p = 0,nwork = 3*n - 1;
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
  dsyev_(&jtype,&uplo,&n,A,&n,w,work,&nwork,&info);
#ifdef DEBUG
  assert(info == 0);
#endif
  for(i=0; i<n; ++i) {
    if (std::abs(w[i]) > 0.0) p++;
  }
  return p;
}

void Complex::determine_flexible_edges(std::vector<int>& flexible_edge)
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
