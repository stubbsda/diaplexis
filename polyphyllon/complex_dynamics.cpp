#include "complex.h"

using namespace DIAPLEXIS;

void Complex::compute_delta(std::set<int>& modified_vertices)
{
  int i,j,m,n,l,nhop;
  std::set<int> vx,current,next;
  std::set<int>::const_iterator it,jt,kt;
  const int nv = (signed) events.size();
  const int nt = (signed) codex.size();
  int done[nv];

  for(i=0; i<nv; ++i) {
    done[i] = 0;
    events[i].topology_modified = false;
  }

  // All of the new vertices...
  for(i=0; i<nv; ++i) {
    if (events[i].incept == -1) vx.insert(i);
  }
  for(i=1; i<=Spacetime::ND; ++i) {
    n = (signed) simplices[i].size();
    // And all the new d-simplices (d >= 1)...
    for(j=0; j<n; ++j) {
      if (simplices[i][j].incept >= 0) continue;
      for(it=simplices[i][j].vertices.begin(); it!=simplices[i][j].vertices.end(); ++it) {
        vx.insert(*it);
      }
    }
  }

  for(i=0; i<nt; ++i) {
    if (!codex[i].active) continue;
    for(it=codex[i].vx_delta.begin(); it!=codex[i].vx_delta.end(); ++it) {
      vx.insert(*it);
    }
  }
  for(i=0; i<nt; ++i) {
    codex[i].vx_delta.clear();
  }
  // Now with this know we can calculate which events need to have their
  // entwinement and/or dimensional stress recalculated
#ifdef VERBOSE
  std::cout << "There are " << vx.size() << " vertices directly implicated" << std::endl;
#endif
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
#ifdef VERBOSE
  std::cout << "There are " << nmod << " modified vertices out of " << nv << std::endl;
#endif
}

double Complex::parity_hamiltonian(double J,bool ferromagnetic,int sheet) const 
{
  // An Ising-like model based on the edge parity...
  int i,ND = dimension(sheet);
  unsigned int j,n;
  SYNARMOSMA::INT64 H = 0;

  if (sheet == -1) {
    for(i=ND; i>=1; --i) {
      n = simplices[i].size();
      for(j=0; j<n; ++j) {
        if (!simplices[i][j].active()) continue;
        H += i*simplices[i][j].parity;
      }
    }
  }
  else {
    for(i=ND; i>=1; --i) {
      n = simplices[i].size();
      for(j=0; j<n; ++j) {
        if (!simplices[i][j].active(sheet)) continue;
        H += i*simplices[i][j].parity;
      }
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
  const int nv = (signed) events.size();
  const double dE = total_energy(-1)/double(nchip);
  int i,j;
  unsigned int d,chip_count[nv];
  bool fired[nv];
  double E;
  std::set<int>::const_iterator it;

  for(i=0; i<nv; ++i) {
    chip_count[i] = 0;
    fired[i] = false;
    if (!events[i].active()) continue;
    E = events[i].get_energy();
    chip_count[i] = int(E/dE);
  }
  for(i=0; i<nv; ++i) {
    d = (unsigned) vertex_valence(i,-1);
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

void Complex::energy_diffusion(double Lambda)
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
    if (!events[i].active()) continue;
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
  // Use the float epsilon because the double epsilon is much too sensitive given the likelihood of round-off
  // error from the various multiplications and divisions carried out in the diffusion algorithm.
  if (std::abs(Esum2 - Esum1) > std::numeric_limits<float>::epsilon()) throw std::runtime_error("Energy conservation error!");
#endif
  for(i=0; i<nv; ++i) {
    if (Enew[i] > std::numeric_limits<double>::epsilon()) events[i].set_energy(Enew[i]);
  }
}