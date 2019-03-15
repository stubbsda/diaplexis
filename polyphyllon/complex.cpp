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

void Complex::write_topology(int sheet) const
{
  int i,j,n;
  std::vector<int> vx;
  SYNARMOSMA::UINT64 q;
  const int ns = cardinality(0,sheet);
  const int ulimit = dimension(sheet);

  if (sheet == -1) {
    for(i=0; i<(signed) events.size(); ++i) {
      if (!events[i].active()) continue;
      simplex_membership(i,vx);
      std::cout << i+1 << ": [";
      for(j=1; j<ulimit; ++j) {
        std::cout << vx[j-1] << ",";
      }
      std::cout << vx[ulimit-1] << "]" << std::endl;
    }
    for(i=Complex::ND; i>=1; i--) {
      n = cardinality(i,sheet);
      if (n > 0) {
        q = int(SYNARMOSMA::binomial(ns,i+1));
        std::cout << "There are " << n << " (" << q << ") " << i << "-simplices in this complex." << std::endl;
      }
    }
  }
  else {
    for(i=0; i<(signed) events.size(); ++i) {
      if (!events[i].active(sheet)) continue;
      simplex_membership(i,vx);
      std::cout << i+1 << ": [";
      for(j=1; j<ulimit; ++j) {
        std::cout << vx[j-1] << ",";
      }
      std::cout << vx[ulimit-1] << "]" << std::endl;
    }
    for(i=Complex::ND; i>=1; i--) {
      n = cardinality(i,sheet);
      if (n > 0) {
        q = int(SYNARMOSMA::binomial(ns,i+1));
        std::cout << "There are " << n << " (" << q << ") " << i << "-simplices in this complex." << std::endl;
      }
    }
  }
  std::cout << "There are " << ns << " (" << ns << ") 0-simplices in this complex." << std::endl;
  std::cout << "The Euler characteristic is " << euler_characteristic(sheet) << std::endl;
}

void Complex::write_incastrature(const std::string& filename,int sheet) const
{
  // A method that generates the Hasse diagram corresponding to the
  // complex's simplicial structure, with the diagram stored in the
  // PDF "hasse.pdf"; the method assumes that the Graphviz library
  // has been installed on the system.
  int i,j,k;

  std::ofstream s(filename,std::ios::trunc);

  s << "digraph G {" << std::endl;
  if (sheet == -1) {
    for(i=Complex::ND; i>0; i--) {
      for(j=0; j<(signed) simplices[i].size(); ++j) {
        if (!simplices[i][j].active()) continue;
        for(k=0; k<1+i; ++k) {
          s << "  \"" << SYNARMOSMA::make_key(simplices[i][j].vertices) << "\" -> \"";
          s << SYNARMOSMA::make_key(simplices[i][j].faces[k]) << "\";" << std::endl;
        }
      }
    }
    for(i=0; i<(signed) events.size(); ++i) {
      if (!events[i].active()) continue;
      s << "  \"" << i << "\" -> \"NULL\";" << std::endl;
    }
  }
  else {
    for(i=Complex::ND; i>0; i--) {
      for(j=0; j<(signed) simplices[i].size(); ++j) {
        if (!simplices[i][j].active(sheet)) continue;
        for(k=0; k<1+i; ++k) {
          s << "  \"" << SYNARMOSMA::make_key(simplices[i][j].vertices) << "\" -> \"";
          s << SYNARMOSMA::make_key(simplices[i][j].faces[k]) << "\";" << std::endl;
        }
      }
    }
    for(i=0; i<(signed) events.size(); ++i) {
      if (!events[i].active(sheet)) continue;
      s << "  \"" << i << "\" -> \"NULL\";" << std::endl;
    }
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
      if (!events[i].active()) {
        std::cout << "Potential problem here: " << i << "  " << events[i].get_energy() << std::endl;
        output = false;
      }
    }
  }
  return output;
}