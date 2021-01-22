#include "spacetime.h"

using namespace DIAPLEXIS;

const int Spacetime::N_EXP;
const int Spacetime::N_IMP;
const std::string Spacetime::EXP_OP[] = {"D","Ux","Ox","R","C","N","A","G","Sg","Sm","Y"};
const std::string Spacetime::IMP_OP[] = {"I","Um","Om","E","F","P","V","Δ"};

std::string Spacetime::implicative_scale(int key,std::vector<double>& parameters) const
{
  // This consists of twelve "notes", eight of which belong to the scale itself
  // (diatonic notes) and four chromatic notes
  // The implicative scale is the treble clef in the F major scale, so the following
  // twelve piano keys in ascending pitch,
  // F4, G4, G4 sharp, A4, A4 sharp, C5, C5 sharp, D5, E5, F5, F5 sharp, G5
  // The four chromatic notes are G4 sharp, C5 sharp, F5 sharp and G5
  // The twelve implicative operations are Um, Om, V, P, I1, I2, E1, E2, E3, F1, F2 and F3
  // V <=> G5*
  // P <=> F5 sharp*
  // E3 <=> F5
  // F3 <=> E5
  // Om <=> D5
  // F2 <=> C5 sharp*
  // F1 <=> C5
  // E2 <=> A4 sharp
  // Um <=> A4
  // I2 <=> G4 sharp*
  // E1 <=> G4
  // I1 <=> F4
  std::string output = "";
  parameters.clear();
  switch (key) {
    case 45:
      output = "I";
      parameters.push_back(0.25);
      break;
    case 47:
      output = "E";
      parameters.push_back(0.25);
      break;
    case 48:
      output = "I";
      parameters.push_back(0.75);
      break;
    case 49:
      output = "Um";
      break;
    case 50:
      output = "E";
      parameters.push_back(0.5);
      break;
    case 52:
      output = "F";
      parameters.push_back(0.25);
      break;
    case 53:
      output = "F";
      parameters.push_back(0.5);
      break;
    case 54:
      output = "Om";
      break;
    case 56:
      output = "F";
      parameters.push_back(0.75);
      break;
    case 57:
      output = "E";
      parameters.push_back(0.75);
      break;
    case 58:
      output = "P";
      parameters.push_back(double(3));
      break;
    case 59:
      output = "V";
      break;
    default:
      throw std::invalid_argument("Illegal key value in implicative scale!");
  }
  return output;
}

std::string Spacetime::explicative_scale(int key,std::vector<double>& parameters) const
{
  // This consists of twelve "notes", eight of which belong to the scale itself
  // (diatonic notes) and four chromatic notes
  // The explicative scale is the bass clef in the F major scale, so the following
  // twelve piano keys in descending pitch,
  // A3, G3, F3 sharp, F3, E3, D3, C3 sharp, C3, A2 sharp, A2, G2, F2
  // The four chromatic notes are A3, G3, F3 sharp and C3 sharp
  // The twelve explicative operations are D, Ox, R, C, G, Sg, Sm, A, N1, N2, Ux1 and Ux2
  // G <=> F2
  // A <=> G2
  // D <=> A2
  // C <=> A2 sharp
  // Sg <=> C3
  // N2 <=> C3 sharp*
  // Ux2 <=> D3
  // Sm <=> E3
  // R <=> F3
  // Ox <=> F3 sharp*
  // N1 <=> G3*
  // Ux1 <=> A3*
  std::string output = "";
  parameters.clear();
  switch (key) {
    case 37:
      output = "Ux";
      parameters.push_back(0.2);
      break;
    case 35:
      output = "N";
      parameters.push_back(2.5);
      break;
    case 34:
      output = "Ox";
      break;
    case 33:
      output = "R";
      break;
    case 32:
      output = "Sm";
      break;
    case 30:
      output = "Ux";
      parameters.push_back(0.5);
      break;
    case 29:
      output = "N";
      parameters.push_back(1.2);
      break;
    case 28:
      output = "Sg";
      break;
    case 26:
      output = "C";
      break;
    case 25:
      output = "D";
      break;
    case 23:
      output = "A";
      break;
    case 21:
      output = "G";
      break;
    default:
      throw std::invalid_argument("Illegal key value in explicative scale!");
  }
  return output;
}

void Spacetime::implication(std::string& output) const
{
  // Should return one of the implicative operators: {F,Um,Om,E,I,P,V,Δ}
  double alpha;
  if (iterations < 50) {
    alpha = skeleton->RND->drandom();
    if (alpha < 0.3) {
      output = "F";
    }
    else if (alpha < 0.6) {
      if (skeleton->RND->drandom() < 0.5) {
        output = "E";
      }
      else {
        output = "I";
      }
    }
    else if (alpha < 0.75) {
      output = "Om";
    }
    else if (alpha < 0.9) {
      output = "Um";
    }
    else {
      if (skeleton->RND->drandom() < 0.33) {
        output = "P";
      }
      else {
        output = "V";
      }
    }
  }
  else {
    if (skeleton->RND->drandom() < 0.5) {
      if (skeleton->RND->drandom() < 0.5) {
        output = "P";
      }
      else {
        output = "V";
      }
    }
    else {
      alpha = skeleton->RND->drandom();
      if (alpha < 0.4) {
        output = "Om";
      }
      else if (alpha < 0.6) {
        output = "Um";
      }
      else if (alpha < 0.8) {
        output = "F";
      }
      else {
        if (skeleton->RND->drandom() < 0.67) {
          output = "E";
        }
        else {
          output = "I";
        }
      }
    }
  }
}

void Spacetime::explication(std::string& output) const
{
  // Should return one of the explicative operators: {G,C,A,Sg,Sm,D,N,Y,R,Ox,Ux}
  if (skeleton->RND->drandom() < 0.1) {
    output = "G";
    return;
  }
  double alpha;
  if (iterations < 50) {
    if (skeleton->RND->drandom() < (0.35 + 1.0/(1+iterations/2))) {
      output = "Sg";
    }
    else {
      if (iterations <= 10) {
        output = "C";
      }
      else {
        if (skeleton->RND->drandom() < 0.5) {
          output = "C";
        }
        else {
          output = "A";
        }
      }
    }
  }
  else {
    alpha = skeleton->RND->drandom();
    if (alpha < 0.25) {
      output = "C";
    }
    else if (alpha < 0.5) {
      output = "Sg";
    }
    else if (alpha < 0.75) {
      output = "N";
    }
    else {
      output = "Ux";
    }
  }
}

void Spacetime::pitch_mapping(std::set<int>* candidates,int sheet) const
{
  int i;
  double sigma,drange,min_value = 1e10,max_value=-1e10;
  const int nv = (signed) skeleton->events.size();

  // We use the negation of the structural deficiency because deficiency < 0 => implication (higher pitch)
  // and deficiency > 0 => explication (lower pitch).
  for(i=0; i<nv; ++i) {
    if (!skeleton->active_event(i,sheet)) continue;
    sigma = -skeleton->events[i].get_deficiency();
    if (sigma < min_value) min_value = sigma;
    if (sigma > max_value) max_value = sigma;
  }
  drange = max_value - min_value;
  for(i=0; i<nv; ++i) {
    if (!skeleton->active_event(i,sheet)) continue;
    sigma = 20.0 + 40.0/drange*(-skeleton->events[i].get_deficiency() - min_value);
    if (sigma < 22.0) {
      candidates[0].insert(i);
    }
    else if (sigma < 24.0) {
      candidates[1].insert(i);
    }
    else if (sigma < 25.5) {
      candidates[2].insert(i);
    }
    else if (sigma < 27.0) {
      candidates[3].insert(i);
    }
    else if (sigma < 28.5) {
      candidates[4].insert(i);
    }
    else if (sigma < 29.5) {
      candidates[5].insert(i);
    }
    else if (sigma < 31.0) {
      candidates[6].insert(i);
    }
    else if (sigma < 32.5) {
      candidates[7].insert(i);
    }
    else if (sigma < 33.5) {
      candidates[8].insert(i);
    }
    else if (sigma < 34.5) {
      candidates[9].insert(i);
    }
    else if (sigma < 36.0) {
      candidates[10].insert(i);
    }
    else if (sigma < 38.0) {
      candidates[11].insert(i);
    }
    else if (sigma < 42.0) {
      candidates[12].insert(i);
    }
    else if (sigma < 46.0) {
      candidates[13].insert(i);
    }
    else if (sigma < 47.5) {
      candidates[14].insert(i);
    }
    else if (sigma < 48.5) {
      candidates[15].insert(i);
    }
    else if (sigma < 49.5) {
      candidates[16].insert(i);
    }
    else if (sigma < 51.0) {
      candidates[17].insert(i);
    }
    else if (sigma < 52.5) {
      candidates[18].insert(i);
    }
    else if (sigma < 53.5) {
      candidates[19].insert(i);
    }
    else if (sigma < 55.0) {
      candidates[20].insert(i);
    }
    else if (sigma < 56.5) {
      candidates[21].insert(i);
    }
    else if (sigma < 57.5) {
      candidates[22].insert(i);
    }
    else if (sigma < 58.5) {
      candidates[23].insert(i);
    }
    else {
      candidates[24].insert(i);
    }
  }
}

int Spacetime::musical_hyphansis(int sheet)
{
  int i,j,v,w,nsuccess = 0;
  bool success = false;
  std::string op,tag,opstring;
  std::vector<int> key_list = codex[sheet].hyphantic_notes[iterations];
  std::vector<double> pvalues;
  std::set<int> candidates[25];
  const int opcount = (signed) key_list.size();

  // Open the hyphantic log file
  std::ofstream s(hyphansis_file,std::ios::app);

  if (opcount == 0) {
    // We're done!
    s << "  </Sheet>" << std::endl;
    s.close();
    return 0;
  }

  pitch_mapping(candidates,sheet);

  // Start "playing" the notes for this voice - our instrument is the topology of spacetime...
  for(i=0; i<opcount; ++i) {
    j = key_list[i];
    if (candidates[key_mapping[j]].empty()) continue;
    if (j == 40) {
      // Edge reorientation for this neutral key, which requires 
      // two distinct vertices since it operates on an edge...
      if (candidates[key_mapping[40]].size() < 2) continue;
      v = skeleton->RND->irandom(candidates[key_mapping[40]]);
      do {
        w = skeleton->RND->irandom(candidates[key_mapping[40]]);
        if (w != v) break; 
      } while(true);
      if (!skeleton->active_event(v,sheet) || !skeleton->active_event(w,sheet)) continue;
      success = skeleton->edge_parity_mutation(v,w,sheet);
      if (success) {
        s << "  <Operation>T," << v << ":" << w << "</Operation>" << std::endl;
        regularization(false,sheet);
        codex[sheet].hyphantic_ops += "T";
        nsuccess++;
        candidates[key_mapping[40]].erase(v);
        candidates[key_mapping[40]].erase(w);
      }
#ifdef DEBUG
      // A very costly assertion...
      assert(skeleton->consistent(sheet));
#endif
      continue;
    }
    v = skeleton->RND->irandom(candidates[key_mapping[j]]);
    // An earlier hyphantic operation may have rendered this
    // event inactive...
    if (!skeleton->active_event(v,sheet)) continue;
    // Now we have the base event v, next we need to get the operator and 
    // parameters for this piano key
    op = (j > 40) ? implicative_scale(j,pvalues) : explicative_scale(j,pvalues);
    opstring = op + "," + std::to_string(v);
    if (op == "F") {
      success = fission(v,pvalues[0],sheet);
      opstring += "," + std::to_string(pvalues[0]); 
    }
    else if (op == "Um") {
      success = fusion_m(v,sheet);
    }
    else if (op == "Om") {
      success = foliation_m(v,sheet);
    }
    else if (op == "E") {
      success = expansion(v,pvalues[0],sheet);
      opstring += "," + std::to_string(pvalues[0]);
    }
    else if (op == "I") {
      success = inflation(v,pvalues[0],sheet);
      opstring += "," + std::to_string(pvalues[0]); 
    }
    else if (op == "P") {
      success = perforation(v,0,sheet);
      opstring += ",0"; 
    }
    else if (op == "V") {
      success = circumvolution(v,sheet);
    }
    else if (op == "D") {
      success = deflation(v,sheet);
    }
    else if (op == "Ux") {
      success = fusion_x(v,pvalues[0],sheet);
      opstring += "," + std::to_string(pvalues[0]); 
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
      opstring += "," + std::to_string(pvalues[0]); 
    }
    else if (op == "A") {
      success = amputation(v,10.0,sheet);
      opstring += ",10.0"; 
    }
    else if (op == "G") {
      success = germination(v,sheet);
    }
    if (success) {
      s << "    <Operation>" << opstring << "</Operation>" << std::endl;
      regularization(false,sheet);
      codex[sheet].hyphantic_ops += op;
      nsuccess++;
    }
#ifdef DEBUG
    // This is an extremely costly assertion...
    assert(skeleton->consistent(sheet));
#endif
  }

  // We're done, so close the hyphantic log file and return
  s << "  </Sheet>" << std::endl;
  s.close(); 

  return nsuccess;
}

int Spacetime::dynamic_hyphansis(int sheet)
{
  int i,v,n,vx[2],nsuccess = 0;
  double alpha;
  std::string op,opstring;
  std::set<int> r_edges;
  bool success = false;
  std::vector<std::pair<int,double> > candidates;
  const int nv = (signed) skeleton->events.size();

  for(i=0; i<nv; ++i) {
    if (!skeleton->events[i].active(sheet)) continue;
    alpha = skeleton->events[i].get_deficiency();
    // Eliminate events whose structural deficiency is close to zero...
    if (std::abs(alpha) < std::numeric_limits<double>::epsilon()) continue;
    candidates.push_back(std::pair<int,double>(i,alpha));
  }

  std::ofstream s(hyphansis_file,std::ios::app);

  if (candidates.size() == 1) {
    v = candidates[0].first;
    if (skeleton->events[v].get_deficiency() < -std::numeric_limits<double>::epsilon()) {
      if (!expansion(v,sheet)) throw std::runtime_error("Unable to expand singleton spacetime complex!");
      codex[sheet].hyphantic_ops += 'E';
      regularization(false,sheet);
      s << "    <Operation>E," << v << "</Operation>" << std::endl;
      nsuccess++;
    }
  }
  else if (candidates.size() > 1) {
    std::sort(candidates.begin(),candidates.end(),SYNARMOSMA::pair_predicate_dbl);
    const int nc = (signed) candidates.size();
    for(i=nc-1; i>0; --i) {
      v = candidates[i].first;
      // An earlier hyphantic operation may have rendered this
      // event inactive...
      if (!skeleton->events[v].active(sheet)) continue;
      alpha = skeleton->events[v].get_deficiency();
      if (alpha < -std::numeric_limits<double>::epsilon()) {
        implication(op);
        opstring = op + "," + std::to_string(v);
        if (op == "F") {
          success = fission(v,0.4,sheet);
          opstring += ",0.4";
        }
        else if (op == "Um") {
          success = fusion_m(v,sheet);
        }
        else if (op == "Om") {
          success = foliation_m(v,sheet);
        }
        else if (op == "E") {
          success = expansion(v,0.15,sheet);
          opstring += ",0.15";
        }
        else if (op == "I") {
          success = inflation(v,0.25,sheet);
          opstring += ",0.25";
        }
        else if (op == "P") {
          success = perforation(v,0,sheet);
          opstring += ",0";
        }
        else if (op == "V") {
          success = circumvolution(v,sheet);
        }
        else if (op == "Δ") {
          success = stellar_deletion(v,sheet);
        }
      }
      else if (alpha > std::numeric_limits<double>::epsilon()) {
        if (skeleton->vertex_dimension(v,sheet) > 1) {
          if (skeleton->RND->drandom() < alpha/10.0) {
            op = "D";
            success = deflation(v,sheet);
            opstring = "D," + std::to_string(v);
          }
          else {
            op = "R";
            success = reduction(v,sheet);
            opstring = "R," + std::to_string(v);
          }
        }
        else {
          explication(op);
          opstring = op + "," + std::to_string(v);
          if (op == "C") {
            success = correction(v,sheet);
          }
          else if (op == "N") {
            success = contraction(v,1.2,sheet);
            opstring += ",1.2";
          }
          else if (op == "Ux") {
            success = fusion_x(v,0.5,sheet);
            opstring += "0.5";
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
            opstring += ",10.0";
          }
          if (!success && op == "R") {
            if (skeleton->RND->irandom(2) == 0) {
              op = "N";
              success = contraction(v,1.2,sheet);
              opstring = "N," + std::to_string(v) + ",1.2"; 
            }
            else {
              op = "Ux";
              success = fusion_x(v,0.5,sheet);
              opstring = "Ux," + std::to_string(v) + ",0.5";
            }
          }
          if (!success) {
            op = "Sg";
            success = compensation_g(v,sheet);
            opstring = "Sg," + std::to_string(v); 
          }
        }
      }
      if (success) {
        s << "    <Operation>" << opstring << "</Operation>" << std::endl;
        regularization(false,sheet);
        codex[sheet].hyphantic_ops += op;
        nsuccess++;
      }
#ifdef DEBUG
      // This is an extremely costly assertion - in one test with an optimized 
      // version it increased the topological hyphansis time from 18.1 seconds 
      // to 823.5 seconds (June 12, 2020). 
      assert(skeleton->consistent(sheet));
#endif
      // If more than 10% of the initial candidate vertices have been successfully used in hyphantic 
      // operations, it's time to exit - we don't want to modify the topology too profoundly in any 
      // given relaxation step...
      if (double(nsuccess)/nc > 0.1) break;
    }
  }
  n = (signed) skeleton->simplices[1].size();
  for(i=0; i<n; ++i) {
    if (!skeleton->simplices[1][i].active(sheet)) continue;
    skeleton->simplices[1][i].get_vertices(vx);
    alpha = 0.5*(std::abs(skeleton->events[vx[0]].get_deficiency()) + std::abs(skeleton->events[vx[1]].get_deficiency()));
    if (alpha > std::numeric_limits<double>::epsilon()) continue;
    if (skeleton->RND->drandom() < parity_mutation) {
      skeleton->simplices[1][i].get_vertices(vx);
      if (skeleton->edge_parity_mutation(vx[0],vx[1],sheet)) {
        s << "    <Operation>T," << vx[0] << ":" << vx[1] << "</Operation>" << std::endl;
        codex[sheet].hyphantic_ops += "T";
        nsuccess++;
      }
    }
  }

  // We're done, so close the hyphantic log file and return
  s << "  </Sheet>" << std::endl;
  s.close();

  return nsuccess;
}

void Spacetime::hyphansis(int sheet)
{
  codex[sheet].hyphantic_ops = "";

  std::ofstream s(hyphansis_file,std::ios::app);
  s << "  <Sheet>" << std::endl;
  s << "    <Index>" << sheet << "</Index>" << std::endl;
  s.close();

#ifdef VERBOSE
  int i,npos = 0,nneg = 0,nze = 0,nv = (signed) skeleton->events.size();
  double alpha;
  for(i=0; i<nv; ++i) {
    if (!skeleton->events[i].active(sheet)) continue;
    alpha = skeleton->events[i].get_deficiency();
    if (!skeleton->events[i].zero_energy()) nze++;
    if (alpha > std::numeric_limits<double>::epsilon()) npos++;
    if (alpha < -std::numeric_limits<double>::epsilon()) nneg++;
  }
  std::cout << "There are " << nze << " events with positive energy or " << 100.0*double(nze)/double(skeleton->cardinality(0,sheet)) << " percent of the total." << std::endl;
  std::cout << "There are " << npos << " positive events and " << nneg << " negative events in the spacetime complex." << std::endl;
#endif

  if (weaving == Hyphansis::dynamic) {
    dynamic_hyphansis(sheet);
  }
  else if (weaving == Hyphansis::musical) {
    musical_hyphansis(sheet);
  }
}

bool Spacetime::event_fusion(int n1,int n2,int sheet)
{
  int i,j,l;
  std::set<int> vx,locus;
  const int ulimit = skeleton->dimension(sheet);

  if (sheet == -1) {
    // A method that converts all occurences of n2 into n1...
    int k,m,n,im1,in1;
    std::set<int> S,duplicate;
    std::set<int>::const_iterator jt;
    std::set<int>::reverse_iterator it;
    std::vector<int> mindex;
    std::vector<int>* mutation = new std::vector<int>[ulimit+1];

    // Do the actual event swap...
    for(i=1; i<=ulimit; ++i) {
      m = (signed) skeleton->simplices[i].size();
      for(j=0; j<m; ++j) {
        if (skeleton->simplices[i][j].exchange(n1,n2)) mutation[i].push_back(j);
      }
    }

    // Place skeleton->simplices in the right dimension...
    for(i=ulimit; i>1; i--) {
      im1 = i - 1;
      n = (signed) skeleton->simplices[i-1].size();
      for(j=0; j<(signed) mutation[i].size(); ++j) {
        in1 = mutation[i][j];
        if (skeleton->simplices[i][in1].dimension() == im1) duplicate.insert(in1);
      }
      for(it=duplicate.rbegin(); it!=duplicate.rend(); ++it) {
        in1 = *it;
        skeleton->simplices[im1].push_back(skeleton->simplices[i][in1]);
        mutation[im1].push_back(n);
        n++;
        skeleton->simplices[i].erase(skeleton->simplices[i].begin() + in1);
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
      if (skeleton->simplices[1][in1].dimension() == 0) duplicate.insert(in1);
    }
    for(it=duplicate.rbegin(); it!=duplicate.rend(); ++it) {
      in1 = *it;
      skeleton->simplices[1].erase(skeleton->simplices[1].begin() + in1);
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
      n = (signed) skeleton->simplices[i].size();
      for(j=0; j<(signed) mutation[i].size(); ++j) {
        in1 = mutation[i][j];
#ifdef DEBUG
        assert(in1 < n);
#endif
        skeleton->simplices[i][in1].get_vertices(vx);
        for(l=0; l<n; ++l) {
          if (l == in1) continue;
          skeleton->simplices[i][l].get_vertices(S);
          if (S == vx) {
            duplicate.insert(in1);
            skeleton->simplices[i][in1].get_ubiquity(locus);
            for(jt=locus.begin(); jt!=locus.end(); ++jt) {
              skeleton->simplices[i][l].activate(*jt);
            }
          }
        }
      }
      for(it=duplicate.rbegin(); it!=duplicate.rend(); ++it) {
        in1 = *it;
        skeleton->simplices[i].erase(skeleton->simplices[i].begin() + in1);
      }
      duplicate.clear();
    }
    skeleton->events[n2].get_ubiquity(locus);
    for(jt=locus.begin(); jt!=locus.end(); ++jt) {
      skeleton->events[n1].activate(*jt);
    }
    skeleton->events[n2].deactivate();
    delete[] mutation;

    // Then recalculate the index_table hash map..
    for(i=1; i<=ulimit; ++i) {
      skeleton->index_table[i].clear();
      m = (signed) skeleton->simplices[i].size();
      for(j=0; j<m; ++j) {
        skeleton->simplices[i][j].clear_entourage();
        skeleton->simplices[i][j].get_vertices(S);
        skeleton->index_table[i][S] = j;
      }
    }
    for(i=0; i<(signed) skeleton->events.size(); ++i) {
      skeleton->events[i].clear_entourage();
    }

    // Then recalculate the entourages...
    skeleton->compute_entourages(-1);
    skeleton->compute_neighbours();
  }
  else {
    bool found;
    std::set<int> v;
    std::set<int>::const_iterator it;
    SYNARMOSMA::hash_map::const_iterator qt;

    locus.insert(sheet);

    for(i=1; i<=Complex::ND; ++i) {
      for(j=0; j<(signed) skeleton->simplices[i].size(); ++j) {
        if (!skeleton->simplices[i][j].active(sheet)) continue;
        // Check of this simplex contains the event "n2", if so
        // divide its colour by sheet and create a new simplex with the
        // right events. Check to see if it already exists, if so
        // multiply its colour by sheet otherwise add the new simplex.
        skeleton->simplices[i][j].get_vertices(v);
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
        skeleton->simplices[i][j].deactivate(sheet);
        for(it=v.begin(); it!=v.end(); ++it) {
          skeleton->events[*it].set_topology_modified(true);
        }
        l = (signed) vx.size() - 1;
        if (l >= 1) {
          qt = skeleton->index_table[l].find(vx);
          if (qt != skeleton->index_table[l].end()) {
            skeleton->simplices[l][qt->second].activate(sheet); 
          }
          else {
            skeleton->simplices[l].push_back(Simplex(vx,locus));
            skeleton->index_table[l][vx] = (signed) skeleton->simplices[l].size() - 1;
          }
        }
        else {
          l = SYNARMOSMA::element(vx);
          skeleton->events[l].activate(sheet); 
        }
        vx.clear();
      }
    }
    skeleton->events[n2].deactivate(sheet); 
  }
  if (!skeleton->events[n2].active()) {
    if (!skeleton->events[n2].zero_energy()) {
      skeleton->events[n1].increment_energy(skeleton->events[n2].get_energy());
      skeleton->events[n2].nullify_energy();
    }
  }
  return true;
}

bool Spacetime::event_twist(int sheet)
{
  // This method will fuse two 0-skeleton->simplices with each other, so as to twist the
  // complex's topology, creating non-orientability and torsion groups in the
  // simplicial homology.
  int i,j,n1,n2,nc,delta,tolerance,n = (signed) skeleton->events.size();
  std::vector<std::string> ekeys;
  SYNARMOSMA::hash_map::const_iterator qt;
  std::set<int> parents;
  std::set<int>::const_iterator it;
  std::vector<int> candidates,used;
  std::vector<std::pair<int,int> > slist;
  static int first = 1;

#ifdef DEBUG
  assert(skeleton->consistent(sheet));
#endif

  tolerance = (first == 1) ? 2 : 1;

  for(i=0; i<n; ++i) {
    if (skeleton->events[i].active(sheet)) candidates.push_back(i);
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
  // to be inactive and then go through all d-skeleton->simplices, d > 0, setting every
  // occurrence of v1 to v2.
  for(i=0; i<nc; ++i) {
    for(j=1+i; j<nc; ++j) {
      delta = skeleton->combinatorial_distance(candidates[i],candidates[j],-1);
      if (delta > tolerance) {
        slist.push_back(std::pair<int,int>(candidates[i],candidates[j]));
      }
    }
  } 
  if (slist.empty()) return false;
  n1 = skeleton->RND->irandom(slist.size());
  n1 = slist[i].first;
  n2 = slist[i].second;
#ifdef VERBOSE
  std::cout << "Twist is " << n1 << " => " << n2 << std::endl;
#endif
  skeleton->events[n2].deactivate(sheet);
  if (!skeleton->events[n2].active()) {
    if (!skeleton->events[n2].zero_energy()) {
      skeleton->events[n1].increment_energy(skeleton->events[n2].get_energy());
      skeleton->events[n2].nullify_energy();
    }
  }
  event_fusion(n1,n2,sheet);
  if (first == 1) first = 0;
  return true;
}

bool Spacetime::event_deletion(int n,int sheet)
{
  if (!skeleton->events[n].active(sheet)) return false;
  int i,vx[2],ne = (signed) skeleton->simplices[1].size();
  skeleton->events[n].deactivate(sheet);
  for(i=0; i<ne; ++i) {
    if (!skeleton->simplices[1][i].active(sheet)) continue;
    skeleton->simplices[1][i].get_vertices(vx);
    if (vx[0] == n || vx[1] == n) skeleton->simplex_deletion(1,i,sheet);
  }
  return true;  
}

int Spacetime::event_addition(const std::vector<double>& xc,int sheet)
{
  int n = (signed) skeleton->events.size();
  Event vt;

  geometry->vertex_addition(xc);
  vt.activate(sheet); 
  skeleton->events.push_back(vt);

  return n;
}

int Spacetime::event_addition(const std::set<int>& antecedents,int sheet)
{
  int n = (signed) skeleton->events.size();
  Event vt;

  geometry->vertex_addition(antecedents);
  vt.activate(sheet);
  skeleton->events.push_back(vt);

  return n;
}

int Spacetime::event_addition(int base,int sheet)
{
  int n = (signed) skeleton->events.size();
  Event vt;

  geometry->vertex_addition(base);
  vt.activate(sheet);
  skeleton->events.push_back(vt);

  return n;
}
