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

int Spacetime::select_event(const std::vector<int>& candidates,double intensity) const
{
  if (candidates.empty()) return -1;
  // The closer the intensity is to unity, the more we should try to choose an element 
  // of the candidates vector which is close to the beginning
  int i,output,n = (signed) candidates.size();
  double cdeficit,tdeficit;
  std::vector<int> vcandidates;

  for(i=0; i<n; ++i) {
    if (!skeleton->active_event(candidates[i])) continue;
    vcandidates.push_back(candidates[i]);
  }
  if (vcandidates.empty()) return -1;
  if (vcandidates.size() == 1) return vcandidates[0];
  n = (signed) vcandidates.size();
  tdeficit = std::abs(skeleton->events[vcandidates[0]].get_deficiency() - skeleton->events[vcandidates[n-1]].get_deficiency());
  for(i=0; i<n; ++i) {
    output = vcandidates[i];
    cdeficit = std::abs(skeleton->events[output].get_deficiency() - skeleton->events[vcandidates[0]].get_deficiency());
    if (intensity <= cdeficit/tdeficit) break;
  }
  return output;
}

int Spacetime::musical_hyphansis(const std::vector<std::pair<int,double> >& candidates)
{
  int i,j,v,its,opcount,nsuccess = 0;
  double m_width = 1.0,x_width = 1.0;
  bool success = false;
  std::string line,temp,op,opstring;
  std::stringstream sstream;
  std::vector<int> key_list,m_events,x_events,neutral_events,m_keys,x_keys;
  std::vector<std::string> elements;
  std::vector<double> pvalues;
  const int nc = (signed) candidates.size();

  // Open the file containing the hyphantic score 
  std::ifstream mscore;
  mscore.exceptions(std::ifstream::badbit);
  try {
    mscore.open(hyphansis_score);

    // Now read the measure that corresponds to this iteration and sheet...
    while(mscore.good()) {
      getline(mscore,line);
      // If the line is empty or doesn't contain a forward slash, ignore it...
      if (line.empty()) continue;
      if (line.find('/') == std::string::npos) continue;
      // Tokenize the line at the forward slash...
      elements.clear();
      sstream.str(line);
      while(getline(sstream,temp,'/')) { 
        elements.push_back(temp); 
      }
#ifdef DEBUG
      assert(elements.size() == 3);
#endif
      its = std::stoi(elements[0]) - 1;
      if (its < iterations) continue;
      if (its > iterations) break;
      // So this is a line for this relaxation step, check if it is the right sheet/voice...
      v = std::stoi(elements[1]);
      // So, grab the piano key...
      key_list.push_back(std::stoi(elements[2]));
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
    return 0; 
  }

  // Start "playing" the notes for this voice - our instrument is the topology of spacetime...
  opcount = (signed) key_list.size();

  // This is a fairly complicated operation - we need to play the notes the in the right order, 
  // while at the same time highest pitched key above 44 is assigned to the event with the most 
  // negative deficiency, the next highest pitched key above 44 acts upon the event with the  
  // second most negative deficiency etc.
  // We need to create a list of the distinct explicative and implicative piano keys in this measure, 
  // paired to the appropriate event
  for(i=0; i<opcount; ++i) {
    if (key_list[i] > 40) {
      if (std::count(m_keys.begin(),m_keys.end(),key_list[i]) == 0) m_keys.push_back(key_list[i]);
    }
    else if (key_list[i] < 40) {
      if (std::count(x_keys.begin(),x_keys.end(),key_list[i]) == 0) x_keys.push_back(key_list[i]);
    }
  }
  // Now sort these piano key values in the correct order, meaning ascending 
  // for explicative (1 to 44) and descending for implicative (88 to 45)
  std::sort(m_keys.begin(),m_keys.end(),std::greater<int>());
  std::sort(x_keys.begin(),x_keys.end());

  for(i=nc-1; i>0; --i) {
    v = candidates[i].first;
    if (!skeleton->active_event(v)) continue;
    neutral_events.push_back(v);
    if (skeleton->events[v].get_deficiency() < -std::numeric_limits<double>::epsilon()) {
      m_events.push_back(v);
    }
    else if (skeleton->events[v].get_deficiency() > std::numeric_limits<double>::epsilon()) {
      x_events.push_back(v);
    }
  }

  if (!m_keys.empty()) m_width = double(m_keys[0] - m_keys.back());
  if (!x_keys.empty()) x_width = double(x_keys.back() - x_keys[0]);
  for(i=0; i<opcount; ++i) {
    j = key_list[i];
    if (j > 40) {
      v = select_event(m_events,double(j - m_keys.back())/m_width);
    }
    else if (j < 40) {
      v = select_event(x_events,double(j - x_keys[0])/x_width);
    }
    else {
      // Any active event will do for a neutral operator
      v = skeleton->RND->irandom(neutral_events);
    }
    if (v == -1) continue;
    // Now we have the base event v, next we need to get the operator and 
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
    opstring = op + "," + std::to_string(v);
    if (op == "F") {
      success = fission(v,pvalues[0]);
      opstring += "," + std::to_string(pvalues[0]); 
    }
    else if (op == "Um") {
      success = fusion_m(v);
    }
    else if (op == "Om") {
      success = foliation_m(v);
    }
    else if (op == "E") {
      success = expansion(v,pvalues[0]);
      opstring += "," + std::to_string(pvalues[0]); 
    }
    else if (op == "I") {
      success = inflation(v,pvalues[0]);
      opstring += "," + std::to_string(pvalues[0]); 
    }
    else if (op == "P") {
      success = perforation(v,0);
      opstring += ",0"; 
    }
    else if (op == "V") {
      success = circumvolution(v);
    }
    else if (op == "D") {
      success = deflation(v);
    }
    else if (op == "Ux") {
      success = fusion_x(v,pvalues[0]);
      opstring += "," + std::to_string(pvalues[0]); 
    }
    else if (op == "Ox") {
      success = foliation_x(v);
    }
    else if (op == "Sg") {
      success = compensation_g(v);
    }
    else if (op == "Sm") {
      success = compensation_m(v);
    }
    else if (op == "R") {
      success = reduction(v);
    }
    else if (op == "C") {
      success = correction(v);
    }
    else if (op == "Δ") {
      success = stellar_deletion(v);
    }
    else if (op == "Y") {
      success = stellar_addition(v);
    }
    else if (op == "N") {
      success = contraction(v,pvalues[0]);
      opstring += "," + std::to_string(pvalues[0]); 
    }
    else if (op == "A") {
      success = amputation(v,10.0);
      opstring += ",10.0"; 
    }
    else if (op == "G") {
      success = germination(v);
    }
    else if (op == "T") {
      success = skeleton->edge_parity_mutation(v);
    }
    if (success) {
      s << "  <Operation>" << opstring << "</Operation>" << std::endl;
      regularization(false);
      hyphantic_ops += op;
      nsuccess++;
    }
#ifdef DEBUG
    // This is an extremely costly assertion...
    assert(skeleton->consistent());
#endif
  }

  // We're done, so close the hyphantic log file and return
  s.close(); 

  return nsuccess;
}

int Spacetime::dynamic_hyphansis(const std::vector<std::pair<int,double> >& candidates)
{
  int i,v,n,vx[2],nsuccess = 0;
  double alpha;
  std::string op,opstring;
  std::set<int> r_edges;
  bool success = false;
  const int nc = (signed) candidates.size();

  std::ofstream s(hyphansis_file,std::ios::app);

  for(i=nc-1; i>0; --i) {
    v = candidates[i].first;
    // An earlier hyphantic operation may have rendered this
    // event inactive...
    if (!skeleton->active_event(v)) continue;
    alpha = skeleton->events[v].get_deficiency();
    if (alpha < -std::numeric_limits<double>::epsilon()) {
      implication(op);
      opstring = op + "," + std::to_string(v);
      if (op == "F") {
        success = fission(v,0.4);
        opstring += ",0.4";
      }
      else if (op == "Um") {
        success = fusion_m(v);
      }
      else if (op == "Om") {
        success = foliation_m(v);
      }
      else if (op == "E") {
        success = expansion(v,0.15);
        opstring += ",0.15";
      }
      else if (op == "I") {
        success = inflation(v,0.25);
        opstring += ",0.25";
      }
      else if (op == "P") {
        success = perforation(v,0);
        opstring += ",0";
      }
      else if (op == "V") {
        success = circumvolution(v);
      }
      else if (op == "Δ") {
        success = stellar_deletion(v);
      }
    }
    else if (alpha > std::numeric_limits<double>::epsilon()) {
      if (skeleton->vertex_dimension(v) > 1) {
        if (skeleton->RND->drandom() < alpha/10.0) {
          op = "D";
          success = deflation(v);
          opstring = "D," + std::to_string(v);
        }
        else {
          op = "R";
          success = reduction(v);
          opstring = "R," + std::to_string(v);
        }
      }
      else {
        explication(op);
        opstring = op + "," + std::to_string(v);
        if (op == "C") {
          success = correction(v);
        }
        else if (op == "N") {
          success = contraction(v,1.2);
          opstring += ",1.2";
        }
        else if (op == "Ux") {
          success = fusion_x(v,0.5);
          opstring += "0.5";
        }
        else if (op == "Sg") {
          success = compensation_g(v);
        }
        else if (op == "Sm") {
          success = compensation_m(v);
        }
        else if (op == "G") {
          success = germination(v);
        }
        else if (op == "Y") {
          success = stellar_addition(v);
        }
        else if (op == "A") {
          success = amputation(v,10.0);
          opstring += ",10.0";
        }
        if (!success && op == "R") {
          if (skeleton->RND->irandom(2) == 0) {
            op = "N";
            success = contraction(v,1.2);
            opstring = "N," + std::to_string(v) + ",1.2"; 
          }
          else {
            op = "Ux";
            success = fusion_x(v,0.5);
            opstring = "Ux," + std::to_string(v) + ",0.5";
          }
        }
        if (!success) {
          op = "Sg";
          success = compensation_g(v);
          opstring = "Sg," + std::to_string(v); 
        }
      }
    }
    if (success) {
      s << "  <Operation>" << opstring << "</Operation>" << std::endl;
      regularization(false);
      hyphantic_ops += op;
      nsuccess++;
    }
#ifdef DEBUG
    // This is an extremely costly assertion...
    assert(skeleton->consistent());
#endif
    // If more than 10% of the initial candidate vertices have been successfully used in hyphantic
    // operations, it's time to exit - we don't want to modify the topology too profoundly in any
    // given relaxation step...
    if (double(nsuccess)/nc > 0.1) break;
  }

  n = (signed) skeleton->simplices[1].size();
  for(i=0; i<n; ++i) {
    if (!skeleton->active_simplex(1,i)) continue;
    if (skeleton->RND->drandom() < parity_mutation) {
      skeleton->simplices[1][i].get_vertices(vx);
      if (skeleton->edge_parity_mutation(vx[0],vx[1])) {
        s << "  <Operation>T," << vx[0] << "</Operation>" << std::endl;
        r_edges.insert(i);
      }
    }
  }
  skeleton->recompute_parity(r_edges);

  // We're done, so close the hyphantic log file and return
  s.close();

  return nsuccess;
}

void Spacetime::hyphansis()
{
  int v;
  double alpha;
  std::vector<std::pair<int,double> > candidates;
  const int nv = (signed) skeleton->events.size();

  hyphantic_ops = "";

  std::ofstream s(hyphansis_file,std::ios::app);
#ifdef VERBOSE
  int npos = 0,nneg = 0,nze = 0;
#endif
  for(int i=0; i<nv; ++i) {
    if (!skeleton->active_event(i)) continue;
    alpha = skeleton->events[i].get_deficiency();
#ifdef VERBOSE
    if (!skeleton->events[i].zero_energy()) nze++;
    if (alpha > std::numeric_limits<double>::epsilon()) npos++;
    if (alpha < -std::numeric_limits<double>::epsilon()) nneg++;
#endif
    // Eliminate events whose structural deficiency is close to zero...
    alpha = std::abs(alpha);
    if (alpha < std::numeric_limits<double>::epsilon()) continue;
    candidates.push_back(std::pair<int,double>(i,alpha));
  }
#ifdef VERBOSE
  std::cout << "There are " << nze << " events with positive energy or " << 100.0*double(nze)/double(skeleton->cardinality(0)) << " percent of the total." << std::endl;
  std::cout << "There are " << npos << " positive events and " << nneg << " negative events in the spacetime complex." << std::endl;
#endif
  if (candidates.empty()) {
    s.close();
    return;
  }

  if (candidates.size() == 1) {
    v = candidates[0].first;
    if (skeleton->events[v].get_deficiency() < -std::numeric_limits<double>::epsilon()) {
      if (!expansion(v)) throw std::runtime_error("Unable to expand singleton spacetime complex!");
      hyphantic_ops += 'E';
      regularization(false);
      s << "  <Operation>E," << v << "</Operation>" << std::endl;
      s.close();
      return;
    }
  }
  s.close();
  std::sort(candidates.begin(),candidates.end(),SYNARMOSMA::pair_predicate_dbl);

  if (weaving == Hyphansis::dynamic) {
    dynamic_hyphansis(candidates);
  }
  else {
    musical_hyphansis(candidates);
  }
}

bool Spacetime::event_fusion(int n1,int n2)
{
  int i,j,l,k,m,n,im1,in1;
  std::set<int> S,vx,duplicate;
  std::set<int>::reverse_iterator it;
  std::vector<int> mindex;
  const int ulimit = skeleton->dimension();
  std::vector<int>* mutation = new std::vector<int>[ulimit+1];

  // Do the actual event swap...
  for(i=1; i<=ulimit; ++i) {
    m = (signed) skeleton->simplices[i].size();
    for(j=0; j<m; ++j) {
      if (skeleton->simplices[i][j].exchange(n1,n2)) mutation[i].push_back(j);
    }
  }

  // Place simplices in the right dimension...
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
          skeleton->simplices[i][l].activate();
        }
      }
    }
    for(it=duplicate.rbegin(); it!=duplicate.rend(); ++it) {
      in1 = *it;
      skeleton->simplices[i].erase(skeleton->simplices[i].begin() + in1);
    }
    duplicate.clear();
  }
  skeleton->events[n1].activate();
  skeleton->events[n2].deactivate();
  delete[] mutation;

  if (!skeleton->active_event(n2)) {
    if (!skeleton->events[n2].zero_energy()) {
      skeleton->events[n1].increment_energy(skeleton->events[n2].get_energy());
      skeleton->events[n2].nullify_energy();
    }
  }

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
  skeleton->compute_entourages();
  skeleton->compute_neighbours();
  return true;
}

bool Spacetime::event_twist()
{
  // This method will fuse two 0-simplices with each other, so as to twist the
  // complex's topology, creating non-orientability and torsion groups in the
  // simplicial homology.
  int i,j,n1,n2,delta,nc,tolerance,n = (signed) skeleton->events.size();
  std::vector<std::string> ekeys;
  SYNARMOSMA::hash_map::const_iterator qt;
  std::set<int> parents;
  std::set<int>::const_iterator it;
  std::vector<int> candidates,used;
  std::vector<std::pair<int,int> > slist;
  static bool first = true;

#ifdef DEBUG
  assert(skeleton->consistent());
#endif

  tolerance = (first) ? 2 : 1;

  for(i=0; i<n; ++i) {
    if (skeleton->active_event(i)) candidates.push_back(i);
  }
  nc = (signed) candidates.size();
  if (nc < 3) return false;
#ifdef VERBOSE
  std::cout << "Event twist with " << nc << " candidates." << std::endl;
#endif
  for(i=0; i<nc; ++i) {
    used.push_back(0);
  }
  // Grab a random pair (v1,v2) of elements from the vector candidates, set v1
  // to be inactive and then go through all d-simplices, d > 0, setting every
  // occurrence of v1 to v2.
  for(i=0; i<nc; ++i) {
    for(j=1+i; j<nc; ++j) {
      delta = skeleton->combinatorial_distance(candidates[i],candidates[j]);
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
  skeleton->events[n2].deactivate();
  if (!skeleton->active_event(n2)) {
    if (!skeleton->events[n2].zero_energy()) {
      skeleton->events[n1].increment_energy(skeleton->events[n2].get_energy());
      skeleton->events[n2].nullify_energy();
    }
  }
  event_fusion(n1,n2);
  if (first) first = false;
  return true;
}

bool Spacetime::event_deletion(int n)
{
  if (!skeleton->active_event(n)) return false;
  int i,vx[2],ne = (signed) skeleton->simplices[1].size();
  skeleton->events[n].deactivate();
  for(i=0; i<ne; ++i) {
    if (!skeleton->active_simplex(1,i)) continue;
    skeleton->simplices[1][i].get_vertices(vx);
    if (vx[0] == n || vx[1] == n) skeleton->simplex_deletion(1,i);
  }
  return true;  
}

int Spacetime::event_addition(const std::vector<double>& xc)
{
  int n = (signed) skeleton->events.size();
  Event vt;

  geometry->vertex_addition(xc);
  skeleton->events.push_back(vt);

  return n;
}

int Spacetime::event_addition(const std::set<int>& antecedents)
{
  int n = (signed) skeleton->events.size();
  Event vt;

  geometry->vertex_addition(antecedents);
  skeleton->events.push_back(vt);

  return n;
}

int Spacetime::event_addition(int base)
{
  int n = (signed) skeleton->events.size();
  Event vt;

  if (skeleton->events[base].zero_energy()) {
    geometry->vertex_addition(base);
  }
  else {
    geometry->vertex_addition(base,1.0/(1.0 + skeleton->events[base].get_energy()));
  }
  skeleton->events.push_back(vt);

  return n;
}
