#include "complex.h"

using namespace DIAPLEXIS;

void Complex::inversion()
{
  int i,j;
  std::set<int> hold;
  std::set<int>::const_iterator it;
  Simplex S;
  const int nv = (signed) events.size();

  for(i=0; i<nv; ++i) {
    if (!events[i].active) {
      // In principle this should be unnecessary...
      events[i].neighbours.clear();
      continue;
    }
    // hold = V / events[i].neighbours
    for(j=0; j<nv; ++j) {
      if (!events[j].active) continue;
      if (j == i) continue;
      if (events[i].neighbours.count(j) == 0) hold.insert(j);
    }
    events[i].neighbours = hold;
    hold.clear();
  }
  for(i=1; i<=Complex::ND; ++i) {
    simplices[i].clear();
    index_table[i].clear();
  }
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
  compute_entourages();
}

double Complex::distribution_fitness(int* volume,const std::vector<int>& affinity,int nprocs) const
{
  int i,sum = 0,bcount = 0;
  double mu,sigma = 0.0;
  std::set<int>::const_iterator it;
  const int nv = (signed) events.size();

  for(i=0; i<nprocs; ++i) {
    sum += volume[i];
  }
  mu = double(sum)/double(nprocs);
  for(i=0; i<nprocs; ++i) {
    sigma += (volume[i] - mu)*(volume[i] - mu);
  }
  sigma = std::sqrt(sigma/double(nprocs));
  for(i=0; i<nv; ++i) {
    if (!events[i].active) continue;
    for(it=events[i].neighbours.begin(); it!=events[i].neighbours.end(); ++it) {
      if (*it < i) continue;
      if (affinity[*it] != affinity[i]) bcount++;
    }
  }
  return sigma + double(bcount);
}

void Complex::distribute(int nprocs) const
{
  if (nprocs < 1) throw std::invalid_argument("The number of processors must be greater than zero!");
  int i,j,k,n,p,p_old,ecount,bcount,its = 0,volume[nprocs],max_dim = 0,current = -1,cproc = 0,nreal = 0;
  int cneighbour;
  bool done,bdry;
  double cost,current_cost;
  std::vector<int> affinity;
  std::vector<std::pair<int,int> > candidates;
  std::set<int> vx,neg,next;
  std::set<int>::const_iterator it;
  const int nv = (signed) events.size();
  const int D = dimension();

  // The algorithm should begin by making current equal to the 
  // vertex with the highest simplicial dimension...
  for(i=0; i<nv; ++i) {
    affinity.push_back(-1);
    if (!events[i].active) continue;
    nreal++;
    if (events[i].topological_dimension > max_dim) {
      current = i;
      max_dim = events[i].topological_dimension;
    }
  }
  for(i=0; i<nprocs; ++i) {
    volume[i] = 0;
  }
#ifdef VERBOSE
  std::cout << "Distributing " << nreal << " vertices across " << nprocs << " processor elements, with maximum simplicial dimension of " << max_dim << "." << std::endl;
#endif
  // First do the higher-dimensional vertices...
  do {
    // So take all of the n-skeleton->simplices (n > 1) that are active and which 
    // contain the current vertex and assign all of their vertices to the 
    // current processor...
    if (events[current].topological_dimension > 1) {
      affinity[current] = cproc;
      volume[cproc] += 1;
      // Need to loop over all d-skeleton->simplices, d > 1
      for(i=2; i<=events[current].topological_dimension; ++i) {
        for(j=0; j<(signed) simplices[i].size(); ++j) {
          if (!simplices[i][j].active) continue;
          if (simplices[i][j].contains(current)) {
            simplices[i][j].get_vertices(vx);
            for(it=vx.begin(); it!=vx.end(); ++it) {
              if (affinity[*it] == -1) {
                affinity[*it] = cproc;
                volume[cproc] += 1;
              }
#ifdef DEBUG
              else {
                assert(affinity[*it] == cproc);
              }
#endif
            }
          }
        }
      }
      // Now handle the neighbouring d-skeleton->simplices (d > 1)...
      do {
        for(i=2; i<=D; ++i) {
          for(j=0; j<(signed) simplices[i].size(); ++j) {
            if (!simplices[i][j].active) continue;
            simplices[i][j].get_vertices(vx);
            bdry = false;
            for(it=vx.begin(); it!=vx.end(); ++it) {
              if (affinity[*it] == -1) {
                neg.insert(*it);
              }
              else {
#ifdef DEBUG
                assert(affinity[*it] == cproc);
#endif
                bdry = true;
              }
            }
            if (!neg.empty() && bdry) {
              for(it=neg.begin(); it!=neg.end(); ++it) {
                next.insert(*it);
              }
            }
            neg.clear(); 
          }
        }
        if (next.empty()) break;
        for(it=next.begin(); it!=next.end(); ++it) {
          affinity[*it] = cproc;
        }  
        next.clear();
      } while(true);
    }    
    done = true;
    for(i=0; i<nv; ++i) {
      if (!events[i].active) continue;
      if (affinity[i] == -1 && events[i].topological_dimension > 1) {
        done = false;
        current = i;
        break;
      }
    }
    cproc = (cproc + 1)%nprocs;
  } while(!done);
  // Now everything else, adding vertices to equilibrate the population counts 
  // and minimize the number of boundary edges...
#ifdef VERBOSE
  std::cout << "Done handling higher-dimensional vertices, now doing edges and vertices..." << std::endl;
#endif  
  for(i=0; i<nprocs; ++i) {
    if (volume[i] == 0) {
      // Find an initial vertex, ideally far from any existing vertices that have 
      // been assigned to a processor...
      candidates.clear();
      for(j=0; j<nv; ++j) {
        if (!events[j].active) continue;
        if (affinity[j] > -1) continue;
        cneighbour = 0;
        for(it=events[j].neighbours.begin(); it!=events[j].neighbours.end(); ++it) {
          if (affinity[*it] > -1) cneighbour++; 
        }
        candidates.push_back(std::pair<int,int>(j,cneighbour));         
      }
      if (candidates.empty()) break;
      std::sort(candidates.begin(),candidates.end(),SYNARMOSMA::pair_predicate_int);
      n = candidates[0].first;
      affinity[n] = i;
      volume[i] += 1;
    }
    done = false;
    do {
      next.clear();
      for(j=0; j<nv; ++j) {
        if (!events[j].active) continue;
        if (affinity[j] > -1) continue;
        for(it=events[j].neighbours.begin(); it!=events[j].neighbours.end(); ++it) {
          if (affinity[*it] == i) next.insert(j);
        }
      }
      if (next.empty()) break;
      for(it=next.begin(); it!=next.end(); ++it) {
#ifdef DEBUG
        assert(affinity[*it] == -1);
#endif
        affinity[*it] = i;
        volume[i] += 1;
        if (volume[i] >= nreal/nprocs) {
          done = true;
          break;
        }
      }
    } while(!done);
  }
  for(i=0; i<nv; ++i) {
    if (!events[i].active) continue;
    if (affinity[i] == -1) {
      for(it=events[i].neighbours.begin(); it!=events[i].neighbours.end(); ++it) {
        p = affinity[*it];
        if (p > -1) {
          affinity[i] = p;
          volume[p] += 1;
          break;
        }
      }
    }
  }
  current_cost = distribution_fitness(volume,affinity,nprocs);
#ifdef VERBOSE
  std::cout << "At iteration 0 the cost is " << current_cost << std::endl;
#endif
  do {
    do {
      n = RND->irandom(nv);
      if (!events[n].active) continue;
      if (events[n].topological_dimension > 1) continue;
      break;
    } while(true);
    p_old = affinity[n];
    do {
      p = RND->irandom(nprocs);
      if (p != p_old) break;
    } while(true);
    affinity[n] = p;
    volume[p] += 1;
    volume[p_old] -= 1;
    cost = distribution_fitness(volume,affinity,nprocs);
    if (cost < current_cost) {
      current_cost = cost;
    }
    else {
      affinity[n] = p_old;
      volume[p] -= 1;
      volume[p_old] += 1;
    }
    its++;
#ifdef VERBOSE
    if (its % 250 == 0) std::cout << "At iteration " << its << " the cost is " << current_cost << std::endl; 
#endif
  } while(its < 10000);
  // Some basic sanity checks: every vertex has a processor...
#ifdef DEBUG
  n = 0;
  for(i=0; i<nv; ++i) {
    if (!events[i].active) continue;
    if (affinity[i] == -1) n++;
  }
  assert(n == 0);
  // And the processors haven't overcounted the vertices.
  n = 0;
  for(i=0; i<nprocs; ++i) {
    n += volume[i];
  }
  assert(n == nreal);
#endif
  // Analysis of the distribution of vertices and edges among the processors...
  for(i=0; i<nprocs; ++i) {
    ecount = 0; bcount = 0;
    for(j=0; j<nv; ++j) {
      if (!events[j].active) continue;
      if (affinity[j] != i) continue;
      for(it=events[j].neighbours.begin(); it!=events[j].neighbours.end(); ++it) {
        k = *it;
        if (affinity[k] == i) {
          ecount++;
        }
        else {
          bcount++;
        }
      }
    }
#ifdef VERBOSE 
    std::cout << "Processor " << i << " owns " << volume[i] << " vertices and " << ecount << " edges, with " << bcount << " boundary edges." << std::endl;
#endif
  }
}

double Complex::set_logical_atoms(int n)
{ 
  if (n < 1) throw std::invalid_argument("The number of atoms must be greater than zero!");

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
  if (events[v].neighbours.empty() || !events[v].active) return 0.0;
  int n = 0;
  double sum = 0.0;
  std::set<int>::const_iterator it;
  SYNARMOSMA::Proposition q,p = events[v].theorem;

  for(it=events[v].neighbours.begin(); it!=events[v].neighbours.end(); ++it) {
    if (!events[*it].active) continue;
    q = p & events[*it].theorem;
    sum += double(q.satisfiable());
    n++;
  }
  sum = sum/double(n);
  return sum;
}

bool Complex::logical_conformity(int v) const
{
  if (!events[v].active) return false;

  std::set<int>::const_iterator it;
  SYNARMOSMA::Proposition p = events[v].theorem;

  for(it=events[v].neighbours.begin(); it!=events[v].neighbours.end(); ++it) {
    if (!events[*it].active) continue;
    p = p & events[*it].theorem;
  }
  return p.satisfiable();
}

void Complex::compute_simplex_energy(int d,int n)
{
  int i,vx[1+d];
  double alpha = 0.0;

  simplices[d][n].get_vertices(vx);

  for(i=0; i<1+d; ++i) {
    alpha += events[vx[i]].get_energy();
  }
  simplices[d][n].energy = alpha/double(1 + d);
}

double Complex::parity_hamiltonian(double J,bool ferromagnetic) const 
{
  // An Ising-like model based on the simplex parity...
  int i;
  unsigned int j,n;
  SYNARMOSMA::INT64 H = 0;
  const int D = dimension();

  for(i=D; i>=1; --i) {
    n = simplices[i].size();
    for(j=0; j<n; ++j) {
      if (!simplices[i][j].active) continue;
      H += i*simplices[i][j].parity;
    }
  }
  if (!ferromagnetic) H = -H;
  return -J*double(H);
}

void Complex::energy_diffusion(int nchip)
{
  // This algorithm, based on the parallel chip-firing game for graphs
  // (cf. T-Y Jiang et al., SIAM J. Disc. Math., 29:615-630, (2015)) is 
  // simple but suffers from the flaw that it does not conserve energy 
  // in the spacetime complex unless nchip > E_total/epsilon, where the 
  // epsilon value is used to test the difference of E_total before and 
  // after the diffusion process. For a fixed graph topology the chip-
  // firing game always settles into a limit cycle but with the dynamic 
  // topology employed here there is no such guarantee.
  const unsigned int nv = events.size();
  const double dE = total_energy()/double(nchip);
  int j;
  unsigned int i,d,chip_count[nv];
  bool fired[nv];
  double E;
  std::set<int>::const_iterator it;

  for(i=0; i<nv; ++i) {
    chip_count[i] = 0;
    fired[i] = false;
    if (!events[i].active) continue;
    E = events[i].get_energy();
    chip_count[i] = int(E/dE);
  }
  for(i=0; i<nv; ++i) {
    d = (unsigned) vertex_valence(i);
    if (chip_count[i] < d) continue;
    fired[i] = true;
    chip_count[i] -= d;
  }
  for(i=0; i<nv; ++i) {
    for(it=events[i].neighbours.begin(); it!=events[i].neighbours.end(); ++it) {
      j = *it;
      if (fired[j]) chip_count[i] += 1;
    }
  }
  for(i=0; i<nv; ++i) {
    if (chip_count[i] == 0) {
      events[i].nullify_energy();
      continue;
    }
    events[i].set_energy(dE*double(chip_count[i]));
  } 
}

void Complex::energy_diffusion(double Lambda,double cutoff)
{
  const int nv = (signed) events.size();
  int i,j,v,n,m,nc; 
  double Enew[nv],E,En,l,d,E_tx,residue; 
  std::set<int>::const_iterator it;
  std::vector<std::tuple<int,double,double> > tvertex;
  std::vector<std::pair<int,double> > candidates,ivertex;

  // I want to try a new method here that begins by considering 
  // those vertices that have non-zero energy. If such a vertex 
  // has a negative deficiency, that means I want to lose energy, 
  // so I next look for a neighbour with positive deficiency (so 
  // it needs energy), with a transfer that is dependent on both 
  // the geometric distance and the relative energy delta. 
#ifdef DEBUG
  double Esum1 = 0.0;
  for(i=0; i<nv; ++i) {
    Esum1 += events[i].get_energy();
  }
#endif
  for(i=0; i<nv; ++i) {
    Enew[i] = -1.0;
    // Inactive vertex...
    if (!events[i].active) continue;
    l = std::abs(events[i].deficiency);
    if (l > std::numeric_limits<double>::epsilon()) candidates.push_back(std::pair<int,double>(i,l));
  }
  if (candidates.empty()) return;  
  std::sort(candidates.begin(),candidates.end(),SYNARMOSMA::pair_predicate_dbl);
  nc = (signed) candidates.size();
  for(i=nc-1; i>=0; --i) {    
    v = candidates[i].first;
    if (Enew[v] > -1.0) continue;
    E = events[v].get_energy();
    // Look for neighbours with an energy value less than 
    // mine that don't already have a new energy value...
    tvertex.clear();
    for(it=events[v].neighbours.begin(); it!=events[v].neighbours.end(); ++it) {
      n = *it;
      if (Enew[n] > -1.0) continue;
      tvertex.push_back(std::tuple<int,double,double>(n,events[n].get_energy(),events[n].deficiency));
    }
    if (tvertex.empty()) continue;
    m = (signed) tvertex.size();
    // If the deficiency > 0, this vertex needs to absorb energy, while 
    // if the deficiency < 0 it wants to donate energy
    if (events[v].deficiency < -std::numeric_limits<double>::epsilon()) {
      // The greater my deficiency the more energy I want to transfer and ideally it needs to be 
      // transferred to the vertices with the largest deficiency
      if (E < std::numeric_limits<double>::epsilon()) continue;
      ivertex.clear();
      for(j=0; j<m; ++j) {
       n = std::get<0>(tvertex[j]);
       En = std::get<1>(tvertex[j]);
       d = std::get<2>(tvertex[j]);
       if (d > std::numeric_limits<double>::epsilon()) ivertex.push_back(std::pair<int,double>(n,d));
      }
      d = -events[v].deficiency/Lambda;
      E_tx = (E < d) ? E : d;
      if (ivertex.empty()) {
        // Just do this the mindless way - a fraction of the vertex's energy is divided amongst 
        // itself and its neighbours, with the fraction decreasing as the energy 
        E_tx = E_tx/double(1 + m);
        Enew[v] = E - double(m)*E_tx;
        for(j=0; j<m; ++j) {
          n = std::get<0>(tvertex[j]);
          Enew[n] = events[n].get_energy() + E_tx;
        }
      }
      else {
        // We have some likely candidates
        m = (signed) ivertex.size();
        // Sum up the outstanding energy need among the neighbours...
        d = 0.0;
        for(j=0; j<m; ++j) {
          d += ivertex[j].second/Lambda;
        }
        // If it's less than the total energy this vertex can donate...
        if (d < E_tx) {
          // We will eliminate each neighbour's need and then spread the remainder equally around
          residue = (E_tx - d)/double(1 + m);
          for(j=0; j<m; ++j) {
            n = ivertex[j].first;
            Enew[n] = events[n].get_energy() + ivertex[j].second/Lambda + residue;
          }
          Enew[v] = E - d - double(m)*residue;
        }
        else {
          // Here there isn't enough to go around, so we will choose random neighbours to 
          // transfer energy to until we exhaust this vertex's spare energy
          Enew[v] = E - E_tx;
          residue = E_tx;
          do {
            j = RND->irandom(m);
            if (Enew[ivertex[j].first] > -1.0) continue;
            l = ivertex[j].second/Lambda;
            n = ivertex[j].first;
            if (l < residue) {
              Enew[n] = events[n].get_energy() + l;
              residue -= l;
            }
            else {
              Enew[n] = events[n].get_energy() + residue;
              break;
            }
          } while(true);
        } 
      }
    }
    else {
      // Look for a candidate which has a positive deficiency
      ivertex.clear();
      for(j=0; j<m; ++j) {
       n = std::get<0>(tvertex[j]);
       En = std::get<1>(tvertex[j]);
       d = std::get<2>(tvertex[j]);
       if (En > std::numeric_limits<double>::epsilon()) ivertex.push_back(std::pair<int,double>(n,En));
      }
      // If none of my neighbours have any energy there is nothing to do but 
      // skip to another vertex
      if (ivertex.empty()) continue;
      // This vertex needs to grab as much energy as it can from its neighbours, 
      // giving priority to those neighbours with a negative deficiency
      m = (signed) ivertex.size();
      d = 0.0;
      for(j=0; j<m; ++j) {
        n = ivertex[j].first;
        if (!(events[n].deficiency < -std::numeric_limits<double>::epsilon())) continue;
        En = events[n].get_energy();
        E_tx = -events[n].deficiency/Lambda;
        d += (En < E_tx) ? En : E_tx;
      }
      if (d < (events[v].deficiency/Lambda)) {
        for(j=0; j<m; ++j) {
          n = ivertex[j].first;
          if (!(events[n].deficiency < -std::numeric_limits<double>::epsilon())) continue;
          En = events[n].get_energy();
          E_tx = -events[n].deficiency/Lambda;
          if (En < E_tx) {
            Enew[n] = 0.0;
          }
          else {
            Enew[n] = events[n].get_energy() - E_tx;
          }
        }
        Enew[v] = E + d;
      }
      else {
        Enew[v] = E + events[v].deficiency/Lambda;
        residue = events[v].deficiency/Lambda;
        do {
          j = RND->irandom(m);
          n = ivertex[j].first;
          if (Enew[n] > -1.0 || !(events[n].deficiency < -std::numeric_limits<double>::epsilon())) continue;
          En = events[n].get_energy();
          E_tx = -events[n].deficiency/Lambda;
          l = (En < E_tx) ? En : E_tx;
          if (l < residue) {
            Enew[n] = (En < E_tx) ? 0.0 : events[n].get_energy() - l;
            residue -= l;
          }
          else {
            Enew[n] = events[n].get_energy() - residue;
            break;
          }
        } while(true);
      }
    }
  }
#ifdef DEBUG
  // A check to verify that the energy has been conserved during its diffusion through the spacetime
  // network.
  double Esum2 = 0.0;
  for(i=0; i<nv; ++i) {
    if (Enew[i] > std::numeric_limits<double>::epsilon()) {
      Esum2 += Enew[i];
    }
    else {
      Esum2 += events[i].get_energy();
    }
  }
  if (std::abs(Esum2 - Esum1) > cutoff) throw std::runtime_error("Energy conservation error!");
#endif
  for(i=0; i<nv; ++i) {
    if (Enew[i] > std::numeric_limits<double>::epsilon()) events[i].set_energy(Enew[i]);
  }
}
