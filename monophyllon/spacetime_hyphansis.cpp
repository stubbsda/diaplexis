#include "spacetime.h"

using namespace DIAPLEXIS;

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
  std::string output = "NULL";
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
      break;
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
  std::string output = "NULL";
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
      break;
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
  //if (skeleton->RND->drandom() < 0.1) return 'G';
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

int Spacetime::select_vertex(const std::vector<int>& candidates,double intensity) const
{
  if (candidates.empty()) return -1;
  // The closer the intensity is to unity, the more we should try to choose an element 
  // of candidates close to the beginning
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
  tdeficit = std::abs(skeleton->events[vcandidates[0]].deficiency - skeleton->events[vcandidates[n-1]].deficiency);
  for(i=0; i<n; ++i) {
    output = vcandidates[i];
    cdeficit = std::abs(skeleton->events[output].deficiency - skeleton->events[vcandidates[0]].deficiency);
    if (intensity <= cdeficit/tdeficit) break;
  }
  return output;
}

void Spacetime::musical_hyphansis(const std::vector<std::pair<int,double> >& candidates)
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
    neutral_vertices.push_back(v);
    if (skeleton->events[v].deficiency < -std::numeric_limits<double>::epsilon()) {
      m_vertices.push_back(v);
    }
    else if (skeleton->events[v].deficiency > std::numeric_limits<double>::epsilon()) {
      x_vertices.push_back(v);
    }
  }

  if (!m_keys.empty()) m_width = double(m_keys[0] - m_keys.back());
  if (!x_keys.empty()) x_width = double(x_keys.back() - x_keys[0]);
  for(i=0; i<opcount; ++i) {
    j = key_list[i];
    if (j > 40) {
      v = select_vertex(m_vertices,double(j - m_keys.back())/m_width);
    }
    else if (j < 40) {
      v = select_vertex(x_vertices,double(j - x_keys[0])/x_width);
    }
    else {
      // Any active vertex will do for a neutral operator
      v = skeleton->RND->irandom(neutral_vertices);
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
      success = fission(v,pvalues[0]);
      opstring << "," << pvalues[0]; 
    }
    else if (op == "Um") {
      success = fusion_m(v);
    }
    else if (op == "Om") {
      success = foliation_m(v);
    }
    else if (op == "E") {
      success = expansion(v,pvalues[0]);
      opstring << "," << pvalues[0]; 
    }
    else if (op == "I") {
      success = inflation(v,pvalues[0]);
      opstring << "," << pvalues[0]; 
    }
    else if (op == "P") {
      success = perforation(v,0);
      opstring << ",0"; 
    }
    else if (op == "V") {
      success = circumvolution(v);
    }
    else if (op == "D") {
      success = deflation(v);
    }
    else if (op == "Ux") {
      success = fusion_x(v,pvalues[0]);
      opstring << "," << pvalues[0]; 
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
      success = skeleton->reduction(v);
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
      opstring << "," << pvalues[0]; 
    }
    else if (op == "A") {
      success = amputation(v,10.0);
      opstring << ",10.0"; 
    }
    else if (op == "G") {
      success = germination(v);
    }
    else if (op == "T") {
      success = skeleton->edge_parity_mutation(v);
    }
    if (success) {
      s << "  <Operation>" << opstring.str() << "</Operation>" << std::endl;
      regularization(false);
      hyphantic_ops += op;
    }
    opstring.str("");
  }

  // We're done, so close the hyphantic log file and return
  s.close(); 
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
#ifdef VERBOSE
    if (!skeleton->events[i].zero_energy()) nze++;
    if (skeleton->events[i].deficiency > std::numeric_limits<double>::epsilon()) npos++;
    if (skeleton->events[i].deficiency < -std::numeric_limits<double>::epsilon()) nneg++;
#endif
    // Eliminate vertices whose structural deficiency is close to zero...
    alpha = std::abs(skeleton->events[i].deficiency);
    if (alpha < std::numeric_limits<double>::epsilon()) continue;
    candidates.push_back(std::pair<int,double>(i,alpha));
  }
#ifdef VERBOSE
  std::cout << "There are " << nze << " vertices with positive energy or " << 100.0*double(nze)/double(skeleton->cardinality(0)) << " percent of the total." << std::endl;
  std::cout << "There are " << npos << " positive vertices and " << nneg << " negative vertices in the spacetime complex." << std::endl;
#endif
  if (candidates.empty()) {
    s.close();
    return;
  }

  if (candidates.size() == 1) {
    v = candidates[0].first;
    if (skeleton->events[v].deficiency < -std::numeric_limits<double>::epsilon()) {
      assert(expansion(v));
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

void Spacetime::dynamic_hyphansis(const std::vector<std::pair<int,double> >& candidates)
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
    if (!skeleton->active_event(v)) continue;
    alpha = skeleton->events[v].deficiency;
    if (alpha < -std::numeric_limits<double>::epsilon()) {
      implication(op);
      opstring << op << "," << v;
      if (op == "F") {
        success = fission(v,0.4);
        opstring << ",0.4";
      }
      else if (op == "Um") {
        success = fusion_m(v);
      }
      else if (op == "Om") {
        success = foliation_m(v);
      }
      else if (op == "E") {
        success = expansion(v,0.15);
        opstring << ",0.15";
      }
      else if (op == "I") {
        success = inflation(v,0.25);
        opstring << ",0.25";
      }
      else if (op == "P") {
        success = perforation(v,0);
        opstring << ",0";
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
          opstring << "D," << v;
        }
        else {
          op = "R";
          success = skeleton->reduction(v);
          opstring << "R," << v;
        }
      }
      else {
        explication(op);
        opstring << op << "," << v;
        if (op == "C") {
          success = correction(v);
        }
        else if (op == "N") {
          success = contraction(v,1.2);
          opstring << ",1.2";
        }
        else if (op == "Ux") {
          success = fusion_x(v,0.5);
          opstring << "0.5";
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
          opstring << ",10.0";
        }
        if (!success && op == "R") {
          opstring.str("");
          if (skeleton->RND->irandom(2) == 0) {
            op = "N";
            success = contraction(v,1.2);
            opstring << "N," << v << ",1.2"; 
          }
          else {
            op = "Ux";
            success = fusion_x(v,0.5);
            opstring << "Ux," << v << ",0.5";
          }
        }
        if (!success) {
          op = "Sg";
          success = compensation_g(v);
          opstring.str("");
          opstring << "Sg," << v; 
        }
      }
    }
    if (success) {
      s << "  <Operation>" << opstring.str() << "</Operation>" << std::endl;
      regularization(false);
      hyphantic_ops += op;
      nsuccess++;
      opstring.str("");
    }
    if (double(nsuccess)/nactive > 0.1) break;
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
  s.close();
}

bool Spacetime::interplication(int centre,double size,int D)
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
  assert(vertex_deletion(centre));
  for(i=0; i<(signed) skeleton->events.size(); ++i) {
    if (!skeleton->active_event(i)) continue;
    l = geometry->get_squared_distance(centre,i,false);
    if (l < size) {
      assert(vertex_deletion(i));
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
      xc.push_back(cvertex[j] + width*(skeleton->RND->drandom() - 0.5));
    }
    if (!geometry->get_uniform()) {
      for(j=g_dim; j<D; ++j) {
        xc.push_back(-dm + dm/2.0*skeleton->RND->drandom()); 
      }
    }
    S.insert(vertex_addition(xc));
    xc.clear();
  }
  skeleton->simplex_addition(S,vx_delta);

  // Now add a set of lower-dimensional simplices to this knot and also tie it 
  // in to the existing spacetime complex...
  base = S;
  kvertex = S;
  S.clear();
  for(i=D-1; i>=3; --i) {
    // As the dimension shrinks, the number of simplices should grow
    l = double(i);
    m = int(skeleton->RND->nrandom(dm - l,2.0));
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
        if (skeleton->RND->drandom() < (0.5*l/dm)) {
          for(k=0; k<g_dim; ++k) {
            xc.push_back(cvertex[k] + width*(skeleton->RND->drandom() - 0.5));
          }
          if (!geometry->get_uniform()) {
            for(k=g_dim; k<i; ++k) {
              alpha = skeleton->RND->nrandom(-l,dm);
              if (alpha < -dm) alpha = -dm + 0.1 + 1.5*skeleton->RND->drandom();
              if (alpha > -0.1) alpha = -0.1 - skeleton->RND->drandom();
              xc.push_back(alpha); 
            }
          }
          q = vertex_addition(xc);
          xc.clear();
          nvertex.insert(q);
          kvertex.insert(q);
        }
        else {
          its = 0;
          do {
            its++;
            q = skeleton->RND->irandom(base);
            // Make sure it isn't already in S,
            if (S.count(q) == 0) break;
          } while(its < 2*((signed) base.size()));
        }
        if (q == -1) {
          // Create a new vertex from scratch
          for(k=0; k<g_dim; ++k) {
            xc.push_back(cvertex[k] + width*(skeleton->RND->drandom() - 0.5));
          }
          if (!geometry->get_uniform()) {
            for(k=g_dim; k<i; ++k) {
              alpha = skeleton->RND->nrandom(-l,dm);
              if (alpha < -dm) alpha = -dm + 0.1 + 1.5*skeleton->RND->drandom();
              if (alpha > -0.1) alpha = -0.1 - skeleton->RND->drandom();
              xc.push_back(alpha); 
            }
          }
          q = vertex_addition(xc);
          xc.clear();
          nvertex.insert(q);
          kvertex.insert(q);          
        }
        S.insert(q);
      }
#ifdef VERBOSE
      std::cout << "Created " << i << "-simplex with " << nvertex.size() << " new vertices" << std::endl;
#endif
      skeleton->simplex_addition(S,vx_delta);
      S.clear();
      nc++;
    } while(nc < m);
    base = nvertex;
    nvertex.clear();
  }
  // Finally the 2-simplexes to bridge the knot with the ambient spacetime complex
  d = bvertex.size()*(bvertex.size() - 1)/2;
  nbridge = int((0.1 + 0.15*skeleton->RND->drandom())*d);
#ifdef VERBOSE
  std::cout << "Adding " << nbridge << " 2-simplices" << std::endl;
#endif
  d = 0;
  do {
    q = skeleton->RND->irandom(bvertex);
    if (skeleton->RND->drandom() < 0.33) {
      do { 
        p = skeleton->RND->irandom(bvertex);
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
    if (skeleton->RND->drandom() < 0.5) {
      S.insert(ambient[0].first);
    }
    else {
      S.insert(ambient[1].first);
    }
    if (skeleton->simplex_addition(S,vx_delta)) d++;
    S.clear();
  } while(d < nbridge);

  regularization(true);
  return true;
}

bool Spacetime::stellar_deletion(int base)
{
  // If this vertex is already part of a d-simplex, d >= 2, then
  // this operation is pointless...
  if (skeleton->vertex_dimension(base) >= 2) return false;
  int m,vx[2];
  std::set<int> nset;
  std::set<int>::const_iterator it;

  for(it=skeleton->events[base].entourage.begin(); it!=skeleton->events[base].entourage.end(); ++it) {
    if (skeleton->active_simplex(1,*it)) {
      skeleton->simplices[1][*it].get_vertices(vx);
      m = (vx[0] == base) ? vx[1] : vx[0];
      nset.insert(m);
    }
  }
  if (nset.size() != 3) return false;    
  vertex_deletion(base);
  skeleton->simplex_addition(nset,vx_delta);
  regularization(true);
  return true;
}

bool Spacetime::stellar_addition(int base)
{
  if (skeleton->vertex_dimension(base) < 2) return false;
  // Take one of these 2-simplices and eliminate it in favour of a 
  // new vertex
  int i,m = (signed) skeleton->simplices[2].size();
  unsigned int l;
  std::vector<int> sx;
  std::vector<double> xc,xtemp;
  std::set<int> candidates,S;
  SYNARMOSMA::hash_map::const_iterator qt;

  for(i=0; i<m; ++i) {
    if (!skeleton->active_simplex(2,i)) continue;
    if (skeleton->simplices[2][i].contains(base)) candidates.insert(i);
  }
  // Choose a 2-simplex at random, get its three vertices and then 
  // eliminate each of the three edges in succession...
  m = skeleton->RND->irandom(candidates);
  skeleton->simplices[2][m].get_vertices(sx);

  S.clear();
  S.insert(sx[0]); S.insert(sx[1]);
  qt = skeleton->index_table[1].find(S);
  skeleton->simplex_deletion(1,qt->second);

  S.clear();
  S.insert(sx[0]); S.insert(sx[2]);
  qt = skeleton->index_table[1].find(S);
  skeleton->simplex_deletion(1,qt->second);

  S.clear();
  S.insert(sx[1]); S.insert(sx[2]);
  qt = skeleton->index_table[1].find(S);
  skeleton->simplex_deletion(1,qt->second);    
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
  m = vertex_addition(xc);
  for(i=0; i<3; ++i) {
    S.clear();
    S.insert(m);
    S.insert(sx[i]);
    skeleton->simplex_addition(S,vx_delta);
  }
  regularization(true);
  return true;    
}

bool Spacetime::germination(int base)
{
  // This method constructs new neighbour vertices w_i for the vertex base which are
  // unit distance from v and orthogonal to v's existing edges, if possible.
  if (skeleton->events[base].boundary) return false;

  bool good,modified = false;
  double a,b,d_min,delta;
  std::set<int> N,Dm2,Dm1,current,tset,free_dims;
  std::set<int>::const_iterator it,jt;
  std::vector<double> x,y,z,xc,bvector;
  SYNARMOSMA::hash_map::const_iterator qt;
  Simplex S;
  int i,j,vx[2],D1 = -1,D2 = -1,m,in1,n = -1,mi = 0;
  const int nv = (signed) skeleton->events.size();
  const int ne = (signed) skeleton->simplices[1].size();
  const int D = (signed) geometry->dimension();

  for(i=0; i<ne; ++i) {
    if (!skeleton->active_simplex(1,i)) continue;
    skeleton->simplices[1][i].get_vertices(vx);
    if (base != vx[0] && base != vx[1]) continue;
    j = (vx[0] == base) ? vx[1] : vx[0];
    N.insert(j);
  }
  for(it=N.begin(); it!=N.end(); ++it) {
    n = *it;
    if (skeleton->active_event(n)) break;
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
    if (!skeleton->active_event(in1)) {
      modified = true;
      skeleton->events[in1].active = true;
    }
    S.initialize(base,in1);
    qt = skeleton->index_table[1].find(S.vertices);
    if (qt == skeleton->index_table[1].end()) {
      modified = true;
      skeleton->simplices[1].push_back(S);
      skeleton->index_table[1][S.vertices] = (signed) skeleton->simplices[1].size() - 1;
    }
    else {
      if (!skeleton->active_simplex(1,qt->second)) {
        modified = true;
        skeleton->simplices[1][qt->second].active = true;
      }
    }
  }
  else {
    modified = true;
    m = vertex_addition(base);
    Dm1.insert(m);
    geometry->set_coordinates(m,xc);
    S.initialize(base,m);
    skeleton->simplices[1].push_back(S);
    skeleton->index_table[1][S.vertices] = (signed) skeleton->simplices[1].size() - 1;
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
    if (!skeleton->active_event(in1)) {
      modified = true;
      skeleton->events[in1].active = true;
    }
    S.initialize(base,in1);
    qt = skeleton->index_table[1].find(S.vertices);
    if (qt == skeleton->index_table[1].end()) {
      modified = true;
      skeleton->simplices[1].push_back(S);
      skeleton->index_table[1][S.vertices] = (signed) skeleton->simplices[1].size() - 1;
    }
    else {
      if (!skeleton->active_simplex(1,qt->second)) {
        modified = true;
        skeleton->simplices[1][qt->second].active = true;
      }
    }
  }
  else {
    modified = true;
    m = vertex_addition(base);
    geometry->set_coordinates(m,xc);
    S.initialize(base,m);
    skeleton->simplices[1].push_back(S);
    skeleton->index_table[1][S.vertices] = (signed) skeleton->simplices[1].size() - 1;
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
    if (!skeleton->active_event(in1)) {
      modified = true;
      skeleton->events[in1].active = true;
    }
    S.initialize(base,in1);
    qt = skeleton->index_table[1].find(S.vertices);
    if (qt == skeleton->index_table[1].end()) {
      modified = true;
      skeleton->simplices[1].push_back(S);
      skeleton->index_table[1][S.vertices] = (signed) skeleton->simplices[1].size() - 1;
    }
    else {
      skeleton->simplices[1][qt->second].active = true;
    }
  }
  else {
    modified = true;
    m = vertex_addition(base);
    geometry->set_coordinates(m,xc);
    Dm1.insert(m);
    S.initialize(base,m);
    skeleton->simplices[1].push_back(S);
    skeleton->index_table[1][S.vertices] = (signed) skeleton->simplices[1].size() - 1;
  }

  if (modified) regularization(false);

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
      if (!skeleton->active_event(in1)) {
        modified = true;
        skeleton->events[in1].active = true;
      }
      S.initialize(base,in1);
      qt = skeleton->index_table[1].find(S.vertices);
      if (qt == skeleton->index_table[1].end()) {
        modified = true;
        skeleton->simplices[1].push_back(S);
        skeleton->index_table[1][S.vertices] = (signed) skeleton->simplices[1].size() - 1;
      }
      else {
        skeleton->simplices[1][qt->second].active = true;
      }
    }
    else {
      modified = true;
      m = vertex_addition(base);
      current.insert(m);
      geometry->set_coordinates(m,xc);
      S.initialize(base,m);
      skeleton->simplices[1].push_back(S);
      skeleton->index_table[1][S.vertices] = (signed) skeleton->simplices[1].size() - 1;
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
      if (!skeleton->active_event(in1)) {
        modified = true;
        skeleton->events[in1].active = true;
      }
      S.initialize(base,in1);
      qt = skeleton->index_table[1].find(S.vertices);
      if (qt == skeleton->index_table[1].end()) {
        modified = true;
        skeleton->simplices[1].push_back(S);
        skeleton->index_table[1][S.vertices] = (signed) skeleton->simplices[1].size() - 1;
      }
      else {
        skeleton->simplices[1][qt->second].active = true;
      }
    }
    else {
      modified = true;
      m = vertex_addition(base);
      current.insert(m);
      geometry->set_coordinates(m,xc);
      S.initialize(base,m);
      skeleton->simplices[1].push_back(S);
      skeleton->index_table[1][S.vertices] = (signed) skeleton->simplices[1].size() - 1;
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

bool Spacetime::correction(int base)
{
  // This method needs to loop over all vertices in a given sheet and find those which
  // are capable of adding another vertex at a distance of (roughly) one and which is
  // orthogonal to the vertex's current set of edges.
  int i,j,n,m,in1;
  unsigned int r;
  bool active,modified = false;
  double l,d1,d2,dbest;
  SYNARMOSMA::hash_map::const_iterator qt;
  std::set<int> candidates;
  Simplex s1;
  const int nv = (signed) skeleton->events.size();

  for(i=0; i<nv; ++i) {
    if (i == base) continue;
    if (!skeleton->active_event(i)) continue;
    if (std::abs(skeleton->events[i].deficiency) < std::numeric_limits<double>::epsilon()) continue;
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
        if (skeleton->active_event(j)) {
          active = true;
          break;
        }
      }
    }
    if (active) continue;
    if (in1 >= 0) {
      if (!skeleton->active_event(in1)) {
#ifdef VERBOSE
        std::cout << "Restoring vertex " << in1 << " in orthonormal fission between " << base << "  " << i << std::endl;
#endif
        skeleton->events[in1].active = true;
        modified = true;
      }
      // Connect this vertex to v and i if necessary
      s1.initialize(base,in1);
      qt = skeleton->index_table[1].find(s1.vertices);
      if (qt == skeleton->index_table[1].end()) {
#ifdef VERBOSE
        std::cout << "Adding edge connecting " << base << " and " << in1 << " in orthonormal fission" << std::endl;
#endif
        skeleton->simplices[1].push_back(s1);
        skeleton->index_table[1][s1.vertices] = (signed) skeleton->simplices[1].size() - 1;
        modified = true;
      }
      else {
        if (!skeleton->active_simplex(1,qt->second)) {
#ifdef VERBOSE
          std::cout << "Restoring edge with key " << SYNARMOSMA::make_key(skeleton->simplices[1][qt->second].vertices) << " in orthonormal fission" << std::endl;
#endif
          skeleton->simplices[1][qt->second].active = true;
          modified = true;
        }
      }
      s1.initialize(i,in1);
      qt = skeleton->index_table[1].find(s1.vertices);
      if (qt == skeleton->index_table[1].end()) {
#ifdef VERBOSE
        std::cout << "Adding edge connecting " << i << " and " << in1 << " in orthonormal fission" << std::endl;
#endif
        skeleton->simplices[1].push_back(s1);
        skeleton->index_table[1][s1.vertices] = (signed) skeleton->simplices[1].size() - 1;
        modified = true;
      }
      else {
        if (!skeleton->active_simplex(1,qt->second)) {
#ifdef VERBOSE
          std::cout << "Restoring edge with key " << SYNARMOSMA::make_key(skeleton->simplices[1][qt->second].vertices) << " in orthonormal fission" << std::endl;
#endif
          skeleton->simplices[1][qt->second].active = true;
          modified = true;
        }
      }
      continue;
    }
    candidates.insert(i);
  }
  if (modified) regularization(true);
  if (candidates.empty()) return false;
  n = skeleton->RND->irandom(candidates);
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
    if (!skeleton->active_event(i)) continue;
    if (geometry->get_squared_distance(i,xc) <= 0.5) return false;
  }
#ifdef VERBOSE
  std::cout << "Adding vertex between " << base << " and " << n << std::endl;
#endif
  m = vertex_addition(xc);
  s1.initialize(base,m);
  skeleton->simplices[1].push_back(s1);
  skeleton->index_table[1][s1.vertices] = (signed) skeleton->simplices[1].size() - 1;
  s1.initialize(n,m);
  skeleton->simplices[1].push_back(s1);
  skeleton->index_table[1][s1.vertices] = (signed) skeleton->simplices[1].size() - 1;

  return true;
}

bool Spacetime::contraction(int base,double l)
{
  int i,u,vx[2];
  std::set<int> pool;
  const double L2 = l*l;
  const int ne = (signed) skeleton->simplices[1].size();

  for(i=0; i<ne; ++i) {
    if (!skeleton->active_simplex(1,i)) continue;
    skeleton->simplices[i][i].get_vertices(vx);
    if (vx[0] != base && vx[1] != base) continue;
    u = (vx[0] == base) ? vx[1] : vx[0];
    if (geometry->get_squared_distance(base,u,true) > L2) pool.insert(i);
  }
  if (pool.empty()) return false;
  i = skeleton->RND->irandom(pool);
#ifdef VERBOSE
  std::cout << "Deleting edge with key " << SYNARMOSMA::make_key(skeleton->simplices[1][i].vertices) << " in localization" << std::endl;
#endif
  skeleton->simplex_deletion(1,i);
  return true;
}

bool Spacetime::compensation_m(int base)
{
  int i,j,vx[2];
  double l;
  std::set<int> candidates,s1;
  SYNARMOSMA::hash_map::const_iterator qt;
  const int ne = (signed) skeleton->simplices[1].size();

  if (skeleton->vertex_dimension(base) < 2) return false;

#ifdef VERBOSE
  std::cout << "Compensation with " << skeleton->dimension() << std::endl;
#endif

  // Remove an edge from this vertex: ideally one that connects with a vertex whose degree
  // is also excessive *and* which is relatively far away (edge length >> 1)
  for(i=0; i<ne; ++i) {
    if (!skeleton->active_simplex(1,i)) continue;
    skeleton->simplices[1][i].get_vertices(vx);
    if (vx[0] != base && vx[1] != base) continue;
    j = (vx[0] == base) ? vx[1] : vx[0];
    if (skeleton->vertex_dimension(j) < 2) continue;
    if (skeleton->events[j].deficiency < -std::numeric_limits<double>::epsilon()) continue;
    l = geometry->get_squared_distance(base,j,true);
    if (l >= 0.81 && l <= 1.21) continue;
    s1.clear();
    s1.insert(base);
    s1.insert(j);
    qt = skeleton->index_table[1].find(s1);
    candidates.insert(qt->second);
  }

  if (candidates.empty()) {
#ifdef VERBOSE
    std::cout << "Unable to perform dimensional compensation, no viable candidates..." << std::endl;
#endif
    return false;
  }
  j = skeleton->RND->irandom(candidates);
#ifdef VERBOSE
  std::cout << "Deleting edge with key " << SYNARMOSMA::make_key(skeleton->simplices[1][j].vertices) << " to reduce the complex's dimensionality" << std::endl;
#endif
  skeleton->simplex_deletion(1,j);
  return true;
}

bool Spacetime::compensation_g(int base)
{
  int i,j,vx[2];
  double l;
  std::set<int> candidates,s1;
  SYNARMOSMA::hash_map::const_iterator qt;
  const int ne = (signed) skeleton->simplices[1].size();

  unsigned int sdegree = 0;
  for(i=0; i<ne; ++i) {
    if (!skeleton->active_simplex(1,i)) continue;
    skeleton->simplices[1][i].get_vertices(vx);
    if (base != vx[0] && base != vx[1]) continue;
    sdegree++;
  }
  if (sdegree == 2*geometry->dimension()) return false;
  const unsigned int idegree = 2*geometry->dimension();
  const int nv = (signed) skeleton->events.size();

  bool up = (sdegree < idegree) ? true : false;
#ifdef VERBOSE
  std::cout << "Compensation with " << base << "  " << up << "  " << sdegree - idegree << std::endl;
#endif
  if (up) {
    // Add an edge to the vertex base: ideally one that connects with a vertex whose degree is
    // also too low and which is relatively close at hand. If the only candidates are distant
    // then do nothing and return false!
    for(i=0; i<nv; ++i) {
      if (i == base) continue;
      if (!skeleton->active_event(i)) continue;
      if (skeleton->edge_exists(base,i)) continue;
      l = geometry->get_squared_distance(base,i,true);
      if (l >= 0.81 && l <= 1.21) candidates.insert(i);
    }

    if (candidates.empty()) {
#ifdef VERBOSE
      std::cout << "Failure in positive compensation..." << std::endl;
#endif
      return false;
    }
    j = skeleton->RND->irandom(candidates);
    vx_delta.insert(base);
    vx_delta.insert(j);
    s1.clear();
    s1.insert(base);
    s1.insert(j);
    qt = skeleton->index_table[1].find(s1);
    if (qt != skeleton->index_table[1].end()) {
#ifdef DEBUG
      assert(!skeleton->active_simplex(1,qt->second));
#endif
#ifdef VERBOSE
      std::cout << "Restoring edge with key " << SYNARMOSMA::make_key(skeleton->simplices[1][qt->second].vertices) << " in positive compensation" << std::endl;
#endif
      skeleton->simplices[1][qt->second].active = true;
    }
    else {
#ifdef VERBOSE
      std::cout << "Adding edge connecting " << base << " and " << j << " in positive compensation" << std::endl;
#endif
      Simplex S(base,j);
      // Now add this simplex to the spacetime complex...
      skeleton->simplices[1].push_back(S);
      skeleton->index_table[1][S.vertices] = (signed) skeleton->simplices[1].size() - 1;
    }
  }
  else {
    // Remove an edge from this vertex: ideally one that connects with a vertex whose degree
    // is also excessive *and* which is relatively far away (edge length >> 1)
    for(i=0; i<ne; ++i) {
      if (!skeleton->active_simplex(1,i)) continue;
      skeleton->simplices[1][i].get_vertices(vx);
      if (vx[0] != base && vx[1] != base) continue;
      j = (vx[0] == base) ? vx[1] : vx[0];
      l = geometry->get_squared_distance(base,j,true);
      if (l >= 0.81 && l <= 1.21) continue;
      s1.clear();
      s1.insert(base);
      s1.insert(j);
      qt = skeleton->index_table[1].find(s1);
      candidates.insert(qt->second);
    }
    if (candidates.empty()) {
#ifdef VERBOSE
      std::cout << "Failure in negative compensation..." << std::endl;
#endif
      return false;
    }
    i = skeleton->RND->irandom(candidates);
#ifdef VERBOSE
    std::cout << "Deleting edge with key " << SYNARMOSMA::make_key(skeleton->simplices[1][i].vertices) << " in negative compensation" << std::endl;
#endif
    skeleton->simplex_deletion(1,i);
  }
  return true;
}

bool Spacetime::unravel(int base)
{
  int i,j,c,vx[2],in_max = -1;
  unsigned int n1,n2;
  std::set<int>::const_iterator it;
  double vf,vf_max = 0.0;

  if (base >= 0) {
    n1 = (signed) skeleton->events[base].entourage.size();
    for(it=skeleton->events[base].entourage.begin(); it!=skeleton->events[base].entourage.end(); ++it) {
      j = *it;
      if (!skeleton->active_simplex(1,j)) continue;
      vf = double(skeleton->simplices[1][j].entourage.size());
      if (vf > vf_max) {
        vf_max = vf;
        in_max = j;
      }
    }
  }
  else {
    for(i=0; i<(signed) skeleton->simplices[1].size(); ++i) {
      if (!skeleton->active_simplex(1,i)) continue;
      vf = double(skeleton->simplices[1][i].entourage.size());
      if (vf > vf_max) {
        vf_max = vf;
        in_max = i;
      }
    }
  }
  if (in_max > -1) {
    // Eliminate the edge in_max...
    skeleton->simplex_deletion(1,in_max);
    return true;
  }
  if (base >= 0) {
    n1 = (signed) skeleton->events[base].entourage.size();
    if (n1 <= 2*geometry->dimension()) return false;
    for(it=skeleton->events[base].entourage.begin(); it!=skeleton->events[base].entourage.end(); ++it) {
      j = *it;
      if (!skeleton->active_simplex(1,j)) continue;
      skeleton->simplices[1][j].get_vertices(vx);
      c = (vx[0] == base) ? vx[1] : vx[0];
      n2 = (signed) skeleton->events[c].entourage.size();
      if (n2 <= 2*geometry->dimension()) continue;
      vf = 0.5*double(n1 + n2 - 4*geometry->dimension());
      if (vf > vf_max) {
        vf_max = vf;
        in_max = j;
      }
    }
  }
  else {
    for(i=0; i<(signed) skeleton->simplices[1].size(); ++i) {
      if (!skeleton->active_simplex(1,i)) continue;
      skeleton->simplices[1][i].get_vertices(vx);
      n1 = (signed) skeleton->events[vx[0]].entourage.size();
      if (n1 <= 2*geometry->dimension()) continue;
      n2 = (signed) skeleton->events[vx[1]].entourage.size();
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
  skeleton->simplex_deletion(1,in_max);
  return true;
}

int Spacetime::compression(double threshold,std::set<int>& vmodified)
{
  int i,n,nc,vx[2];
  std::set<int>::const_iterator it;
  double s,alpha;
  std::vector<int> candidates;
  const int nv = (signed) skeleton->events.size();
  const int ne = (signed) skeleton->simplices[1].size();

  for(i=0; i<ne; ++i) {
    if (!skeleton->active_simplex(1,i)) continue;
    alpha = std::abs(skeleton->simplices[1][i].volume);
    if (alpha > threshold) {
      s = std::exp(-(alpha - threshold));
      if (skeleton->RND->drandom() > s) candidates.push_back(i);
    }
  }
  nc = (signed) candidates.size();
  for(i=0; i<nc; ++i) {
    n = candidates[i];
    skeleton->simplex_deletion(1,n);
    // Remove the references in events[v/w].neighbours and
    // events[v/w].entourage
    skeleton->simplices[1][n].get_vertices(vx);

    it = std::find(skeleton->events[vx[0]].neighbours.begin(),skeleton->events[vx[0]].neighbours.end(),vx[1]);
    if (it != skeleton->events[vx[0]].neighbours.end()) skeleton->events[vx[0]].neighbours.erase(it);
    it = std::find(skeleton->events[vx[1]].neighbours.begin(),skeleton->events[vx[1]].neighbours.end(),vx[0]);
    if (it != skeleton->events[vx[1]].neighbours.end()) skeleton->events[vx[1]].neighbours.erase(it);

    it = std::find(skeleton->events[vx[0]].entourage.begin(),skeleton->events[vx[0]].entourage.end(),n);
    if (it != skeleton->events[vx[0]].entourage.end()) skeleton->events[vx[0]].entourage.erase(it);
    it = std::find(skeleton->events[vx[1]].entourage.begin(),skeleton->events[vx[1]].entourage.end(),n);
    if (it != skeleton->events[vx[1]].entourage.end()) skeleton->events[vx[1]].entourage.erase(it);

    vmodified.insert(vx[0]);
    vmodified.insert(vx[1]);
  }
  if (!skeleton->connected()) {
    // We need to add enough (short!) edges to keep the spacetime connected...
    int j,k,nsedge,ncomp;
    double l,tlength;
    std::set<int>::const_iterator jt,kt;
    std::set<int> linked;
    std::vector<int> component,sedge;
    std::vector<std::tuple<int,int,double> > connect;
    Simplex S;

    ncomp = skeleton->component_analysis(component);
#ifdef VERBOSE
    std::cout << "There are " << ncomp << " components in this spacetime." << std::endl;
#endif
    std::vector<int>* cvertex = new std::vector<int>[ncomp];
    for(i=0; i<nv; ++i) {
      if (!skeleton->active_event(i)) continue;
      cvertex[component[i]].push_back(i);
    }

    // Now we want to find the shortest possible edge that connects all of these
    // components together, first the brute force solution
    for(i=0; i<ncomp; ++i) {
      for(j=1+i; j<ncomp; ++j) {
        l = minimize_lengths(cvertex[i],cvertex[j],vx);
        connect.push_back(std::tuple<int,int,double>(vx[0],vx[1],l));
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
      S.initialize(j,k);
      skeleton->simplices[1].push_back(S);
      skeleton->index_table[1][S.vertices] = (signed) skeleton->simplices[1].size() - 1;
      skeleton->events[j].neighbours.insert(k);
      skeleton->events[k].neighbours.insert(j);
      vmodified.insert(k);
      vmodified.insert(j);
    }
    nc -= nsedge;
#ifdef DEBUG
    assert(skeleton->connected());
#endif
    delete[] cvertex;
  }
  return nc;
}

bool Spacetime::foliation_m(int base)
{
  int i,p,n1,n2,vx[2];
  std::set<int> candidates,s1;
  SYNARMOSMA::hash_map::iterator qt;
  const int ne = (signed) skeleton->simplices[1].size();

  for(i=0; i<ne; ++i) {
    if (!skeleton->active_simplex(1,i)) continue;
    skeleton->simplices[1][i].get_vertices(vx);
    if (vx[0] != base && vx[1] != base) continue;
    p = (vx[0] == base) ? vx[1] : vx[0];
    if (std::abs(skeleton->events[p].deficiency) < std::numeric_limits<double>::epsilon()) continue;
    candidates.insert(p);
  }
  if (candidates.size() < 2) return false;
  n1 = skeleton->RND->irandom(candidates);
  do {
    n2 = skeleton->RND->irandom(candidates);
    if (n2 != n1) break;
  } while(true);

  // Clearly adding this edge creates at least a 2-simplex (v,n1,n2) which should be
  // added to the spacetime complex!
  std::set<int> leaf;
  leaf.insert(base);
  leaf.insert(n1);
  leaf.insert(n2);
  Simplex S(leaf);
  qt = skeleton->index_table[2].find(S.vertices);
  if (qt == skeleton->index_table[2].end()) {
    skeleton->simplices[2].push_back(S);
    skeleton->index_table[2][S.vertices] = (signed) skeleton->simplices[2].size() - 1;
    return true;
  }
  if (skeleton->active_simplex(2,qt->second)) return false;
  skeleton->simplices[2][qt->second].active = true;
  s1.insert(n1); s1.insert(n2);
  qt = skeleton->index_table[1].find(s1);
  skeleton->simplices[1][qt->second].active = true;
  vx_delta.insert(n1);
  vx_delta.insert(n2);
  return true;
}

bool Spacetime::foliation_x(int base)
{
  int i,p,n1,n2,vx[2];
  std::set<int> candidates,S;
  SYNARMOSMA::hash_map::iterator qt;
  const int ne = (signed) skeleton->simplices[1].size();

  for(i=0; i<ne; ++i) {
    if (!skeleton_active_simplex(1,i)) continue;
    skeleton->simplices[1][i].get_vertices(vx);
    if (vx[0] != base && vx[1] != base) continue;
    p = (vx[0] == base) ? vx[1] : vx[0];
    if (std::abs(skeleton->events[p].deficiency) < std::numeric_limits<double>::epsilon()) continue;
    candidates.insert(p);
  }
  if (candidates.size() < 2) return false;
  n1 = skeleton->RND->irandom(candidates);
  do {
    n2 = skeleton->RND->irandom(candidates);
    if (n2 != n1) break;
  } while(true);
  S.insert(n1); S.insert(n2);
  qt = skeleton->index_table[1].find(S);
  if (qt == skeleton->index_table[1].end()) return false;
  if (!skeleton->active_simplex(1,qt->second)) return false;
  skeleton->simplex_deletion(1,qt->second);
  return true;
}

bool Spacetime::amputation(int base,double tolerance)
{
  int i,j,n,p,vx[2];
  std::set<int> candidates;
  const int ne = (signed) skeleton->simplices[1].size();
  const int ulimit = skeleton->dimension();

  if (skeleton->events[base].deficiency > tolerance) candidates.insert(base);
  for(i=0; i<ne; ++i) {
    if (!skeleton->active_simplex(1,i)) continue;
    skeleton->simplices[1][i].get_vertices(vx);
    if (base != vx[0] && base != vx[1]) continue;
    n = (vx[0] == base) ? vx[1] : vx[0];
    if (skeleton->events[n].deficiency > tolerance) candidates.insert(n);
  }
  
  if (candidates.empty()) return false;
  n = skeleton->RND->irandom(candidates);
 
#ifdef VERBOSE
  std::cout << "Amputating vertex " << n << " and all its dependent simplices" << std::endl;
#endif
  // Delete this vertex and all its edges...
  skeleton->events[n].active = false;
  vx_delta.insert(n);
  for(i=ulimit; i>=1; --i) {
    p = (signed) skeleton->simplices[i].size();
    for(j=0; j<p; ++j) {
      if (!skeleton->active_simplex(i,j)) continue;
      if (skeleton->simplices[i][j].contains(n)) skeleton->simplex_deletion(i,j);
    }
  }
  if (!skeleton->active_event(n)) {
    if (!skeleton->events[n].zero_energy()) {
      skeleton->events[base].increment_energy(events[n].get_energy());
      skeleton->events[n].nullify_energy();
    }
  }
  
  return true;
}

bool Spacetime::fusion_x(int base,double tolerance)
{
  int i,u;
  std::set<int> candidates;
  const int nv = (signed) skeleton->events.size();

  for(i=0; i<nv; ++i) {
    if (i == base) continue;
    if (!skeleton->active_event(i)) continue;
    if (std::abs(skeleton->events[i].deficiency) < std::numeric_limits<double>::epsilon()) continue;
    if (geometry->get_squared_distance(base,i,true) > tolerance) continue;
    candidates.insert(i);
  }
  if (candidates.empty()) return false;
  u = skeleton->RND->irandom(candidates);
#ifdef VERBOSE
  std::cout << "Fusing deficient vertices: " << u << " => " << base << std::endl;
#endif
  vertex_fusion(base,u);
  return true;
}

bool Spacetime::fusion_m(int base)
{
  int i,u,vx[2];
  std::set<int> candidates;
  const int ne = (signed) skeleton->simplices[1].size();

  for(i=0; i<ne; ++i) {
    if (!skeleton->active_simplex(1,i)) continue;
    skeleton->simplices[1][i].get_vertices(vx);
    if (vx[0] != base && vx[1] != base) continue;
    u = (vx[0] == base) ? vx[1] : vx[0];
    candidates.insert(u);
  }
  if (candidates.empty()) return false;
  u = skeleton->RND->irandom(candidates);
#ifdef VERBOSE
  std::cout << "Fusing vertices: " << u << " => " << base << std::endl;
#endif
  vertex_fusion(base,u);
  return true;
}

bool Spacetime::fission(int base,double density)
{
  int i,p,q,n,vx[2];
  Simplex S;
  const int ne = (signed) skeleton->simplices[1].size();

  if (base >= 0) {
    n = 0;

    p = vertex_addition(base);

    for(i=0; i<ne; ++i) {
      if (!skeleton->active_simplex(1,i)) continue;
      skeleton->simplices[1][i].get_vertices(vx);
      if (vx[0] != base && vx[1] != base) continue;
      q = (vx[0] == base) ? vx[1] : vx[0];
      if (skeleton->RND->drandom() < density) {
        S.initialize(q,p);
        skeleton->simplices[1].push_back(S);
        skeleton->index_table[1][S.vertices] = (signed) skeleton->simplices[1].size() - 1;
        n++;
      }
    }
    S.initialize(base,p);
    skeleton->simplices[1].push_back(S);
    skeleton->index_table[1][S.vertices] = (signed) skeleton->simplices[1].size() - 1;
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

    if (skeleton->dimension() < 1) return false;

    if (skeleton->RND->drandom() < 0.5) {
      // The simplest case: we just need to add a new vertex, clone the antecedent
      // vertex's one-dimensional entourage and lastly create a 1-simplex joining
      // the two vertices.
      int d;

      do {
        p = skeleton->RND->irandom(skeleton->events.size());
        if (skeleton->active_event(p)) break;
      } while(true);
      q = vertex_addition(p);
      nsimplex.insert(p);
      nsimplex.insert(q);

      for(i=0; i<ne; ++i) {
        if (!skeleton->active_simplex(1,i)) continue;
        skeleton->simplices[1][i].get_vertices(vx);
        if (vx[0] != p && vx[1] != p) continue;
        n = (vx[0] == p) ? vx[1] : vx[0];
        nsimplex.insert(n);
      }
      d = (signed) nsimplex.size() - 1;
      if (d >= Complex::ND) {
        std::set<int> s1;
        d = skeleton->RND->irandom(2,Complex::ND);
        do {
          n = skeleton->RND->irandom(nsimplex);
          s1.insert(n);
        } while((signed) s1.size() < (1+d));
        nsimplex = s1;
      }
      S.initialize(nsimplex);
#ifdef VERBOSE
      std::cout << "Vertex fission on " << p << " with new vertex " << q << " and simplex dimension " << d << std::endl;
#endif
      skeleton->simplices[d].push_back(S);
      skeleton->index_table[d][S.vertices] = skeleton->simplices[d].size() - 1;
    }
    else {
      // The plot thickens: we must create a new 1-simplex out of an existing one, perhaps by forming
      // a pair of 2-simplices: (w1,u_new,v_new) and (w2,u_new,v_new)
      std::set<int> antecedent;

      do {
        n = skeleton->RND->irandom(skeleton->simplices[1].size());
        if (skeleton->active_simplex(1,n)) break;
      } while(true);
      skeleton->simplices[1][n].get_vertices(vx);
      antecedent.insert(vx[0]);
      antecedent.insert(vx[1]);
      p = vertex_addition(antecedent);
      antecedent.insert(p);
      q = vertex_addition(antecedent);

      nsimplex.insert(vx[0]);
      nsimplex.insert(p);
      nsimplex.insert(q);
      S.initialize(nsimplex);
      skeleton->simplices[2].push_back(S);
      skeleton->index_table[2][S.vertices] = (signed) skeleton->simplices[2].size() - 1;
      nsimplex.clear();

      nsimplex.insert(vx[1]);
      nsimplex.insert(p);
      nsimplex.insert(q);
      S.initialize(nsimplex);
      skeleton->simplices[2].push_back(S);
      skeleton->index_table[2][S.vertices] = (signed) skeleton->simplices[2].size() - 1;
    }
    return true;
  }
}

void Spacetime::regularization(bool minimal)
{
  if (!minimal) skeleton->simplicial_implication();

  skeleton->compute_neighbours();
  if (skeleton->connected()) {
    skeleton->compute_entourages();
    return;
  }

  int i,j,k,v1,v2,n1,n2,nc,loc = -1;
  unsigned int max = 0;
  double l,mdelta;
  std::string sx;
  std::set<int> colours;
  std::vector<int> component;
  Simplex S;
  SYNARMOSMA::hash_map::const_iterator qt;
  const int nv = (signed) skeleton->events.size();

  nc = skeleton->component_analysis(component);
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

  for(i=0; i<nv; ++i) {
    if (!skeleton->active_event(i)) continue;
    cvertex[component[i]].push_back(i);
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
      if (skeleton->events[cvertex[i][0]].zero_energy()) {
        // Eliminate this vertex altogether...
#ifdef VERBOSE
        std::cout << "Eliminating isolated vertex " << cvertex[i][0] << std::endl;
#endif
        skeleton->events[cvertex[i][0]].active = false;
        cvertex[i].erase(cvertex[i].begin());
      }
    }
  }
  n1 = (signed) cvertex[loc].size();

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
    S = Simplex(v1,v2);
    qt = skeleton->index_table[1].find(S.vertices);
    if (qt == skeleton->index_table[1].end()) {
#ifdef VERBOSE
      std::cout << "Adding edge connecting " << v1 << " and " << v2 << " to maintain the complex's connectedness" << std::endl;
#endif
      skeleton->simplices[1].push_back(S);
      skeleton->index_table[1][S.vertices] = (signed) skeleton->simplices[1].size() - 1;
    }
    else {
#ifdef DEBUG
      assert(!skeleton->active_simplex(1,qt->second));
#endif
#ifdef VERBOSE
      std::cout << "Restoring the simplex with key " << SYNARMOSMA::make_key(skeleton->simplices[1][qt->second].vertices) << " to maintain the complex's connectedness" << std::endl;
#endif
      skeleton->simplices[1][qt->second].active = true;
    }
    // Need to ensure that v1 and v2 are also coloured by sheet...
    skeleton->events[v1].active = true;
    skeleton->events[v2].active = true;
  }
  skeleton->compute_neighbours();
#ifdef DEBUG
  assert(skeleton->connected());
#endif
  skeleton->compute_entourages();
  delete[] cvertex;
}

bool Spacetime::circumvolution()
{
  // This method attempts a circumvolution using boundary edges, it is
  // normally only called when the global energy reaches a critical threshold.
  int i,j,k,l,u[2],w[2];
  std::set<int> edge_set,S;
  std::set<int>::const_iterator it,jt;
  SYNARMOSMA::hash_map::const_iterator qt;
  const int D = skeleton->dimension();
  int vx[D];

  for(i=0; i<(signed) skeleton->simplices[D].size(); ++i) {
    if (!skeleton->active_simplex(D,i)) continue;
    for(j=0; j<=D; ++j) {
      qt = skeleton->index_table[D-1].find(skeleton->simplices[D][i].faces[j]);
      if (skeleton->simplices[D-1][qt->second].entourage.size() == 1) {
        skeleton->simplices[D-1][qt->second].get_vertices(vx);
        for(k=0; k<D; ++k) {
          for(l=k+1; l<D; ++l) {
            S.clear();
            S.insert(vx[k]);
            S.insert(vx[l]);
            qt = skeleton->index_table[1].find(S);
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

      skeleton->simplices[1][n1].get_vertices(w);
      skeleton->simplices[1][n2].get_vertices(u);

      d = skeleton->combinatorial_distance(w[0],u[0]);
      if (d <= 2) continue;

      d = skeleton->combinatorial_distance(w[1],u[0]);
      if (d <= 2) continue;

      d = skeleton->combinatorial_distance(w[0],u[1]);
      if (d <= 2) continue;

      d = skeleton->combinatorial_distance(w[1],u[1]);
      if (d <= 2) continue;

      candidates.push_back(n1);
      candidates.push_back(n2);
    }
  }
  if (candidates.empty()) return false;
  np = candidates.size()/2;
  i = skeleton->RND->irandom(np);
  n1 = candidates[2*i];
  n2 = candidates[2*i+1];
  skeleton->simplices[1][n1].get_vertices(w);
  skeleton->simplices[1][n2].get_vertices(u);
  vx_delta.insert(u[0]);
  vx_delta.insert(u[1]);
  if (skeleton->RND->drandom() < 0.5) {
    vertex_fusion(w[0],u[0]);
    vertex_fusion(w[1],u[1]);
  }
  else {
    vertex_fusion(w[0],u[1]);
    vertex_fusion(w[1],u[0]);
  }
  return true;
}

bool Spacetime::circumvolution(int base)
{
  // This method seeks to fuse together two d-simplices, one of
  // which contains the vertex base
  int i,d,nd,s1,s2;
  Simplex S;
  std::set<int> candidates;
  std::vector<int> order;

  if (base >= 0) {
    nd = skeleton->vertex_dimension(base);
    if (nd < 1) return false;
    d = skeleton->RND->irandom(1,nd);
    nd = (signed) skeleton->simplices[d].size();
    for(i=0; i<nd; ++i) {
      if (!skeleton->active_simplex(d,i)) continue;
      if (skeleton->simplices[d][i].contains(base)) candidates.insert(i);
    }
    s1 = skeleton->RND->irandom(candidates);
    candidates.clear();
    for(i=0; i<nd; ++i) {
      if (!skeleton->active_simplex(d,i)) continue;
      if (i == s1) continue;
      S = skeleton->simplices[d][s1] ^ skeleton->simplices[d][i];
      if (S.empty()) candidates.insert(i);
    }
    if (candidates.empty()) return false;
    // Now find which of these candidates is the closest to the d-simplex s1, where the
    // distance between two simplices is defined to be the maximum inter-vertex distance
    int j,vx1[1+d],vx2[1+d];
    double D1,delta;
    std::set<int> slist;
    std::set<int>::const_iterator it;

    skeleton->simplices[d][s1].get_vertices(vx1);
    for(it=candidates.begin(); it!=candidates.end(); ++it) {
      skeleton->simplices[d][*it].get_vertices(vx2);
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
    s2 = skeleton->RND->irandom(slist);
  }
  else {
    do {
      d = skeleton->RND->irandom(1,skeleton->dimension());
      for(i=0; i<(signed) skeleton->simplices[d].size(); ++i) {
        if (skeleton->active_simplex(d,i)) candidates.insert(i);
      }
      if (candidates.size() > 3) break;
    } while(true);
    do {
      s1 = skeleton->RND->irandom(candidates);
      s2 = skeleton->RND->irandom(candidates);
      if (s1 == s2) continue;
      // Now check to make sure that the intersection of
      // these two d-simplices is null...
      S = skeleton->simplices[d][s1] ^ skeleton->simplices[d][s2];
      if (S.empty()) break;
    } while(true);
  }
#ifdef VERBOSE
  std::cout << "Circumvolving the " << d << "-simplices: " << SYNARMOSMA::make_key(skeleton->simplices[d][s2].vertices) << " => " << SYNARMOSMA::make_key(skeleton->simplices[d][s1].vertices) << std::endl;
#endif
  int v1[d+1],v2[d+1];
  for(i=0; i<=d; ++i) {
    order.push_back(i);
  }
  std::random_shuffle(order.begin(),order.end());
  skeleton->simplices[d][s1].get_vertices(v1);
  skeleton->simplices[d][s2].get_vertices(v2);

  for(i=0; i<=d; ++i) {
    vx_delta.insert(v2[order[i]]);
    vertex_fusion(v1[i],v2[order[i]]);
  }
  return true;
}

bool Spacetime::expansion(int base,double creativity)
{
  // Create an entirely new d-simplex
  int i,d,k,novum = 0;
  std::set<int> vx;

  if (base >= 0) {
    int n1 = skeleton->vertex_dimension(base);
    if (n1 == Complex::ND) return false;
    int u,vtx[2],its = 0;
    std::set<int> M;
    double tau = std::abs(skeleton->events[base].deficiency) + 25.0;
    const int ne = (signed) skeleton->simplices[1].size();

    d = int(double(Complex::ND-1-n1)*1.0/(1.0 + std::exp(tau)) + double(1 + n1));

    for(i=0; i<ne; ++i) {
      if (!skeleton->active_simplex(1,i)) continue;
      skeleton->simplices[1][i].get_vertices(vtx);
      if (vtx[0] != base && vtx[1] != base) continue;
      u = (vtx[0] == base) ? vtx[1] : vtx[0];
      if (skeleton->events[u].deficiency < -std::numeric_limits<double>::epsilon()) M.insert(u);
    }
    if (M.empty()) creativity = 1.0;
    vx.insert(base);
    do {
      its++;
      if (its == Complex::ND) creativity = 1.0;
      if (skeleton->RND->drandom() < creativity) {
        // Create a new vertex...
        k = vertex_addition(base);
        novum++;
      }
      else {
        k = skeleton->RND->irandom(M);
      }
      vx.insert(k);
    } while((signed) vx.size() < (1+d));
  }
  else {
    int j,nv;
    bool success;
    std::set<int>::const_iterator it;

    i = skeleton->dimension();
    d = (i <= 2) ? 3 : skeleton->RND->irandom(2,i);
    if (double(d+1)/double(skeleton->events.size()) > 0.8) creativity = 1.0;
    do {
      if (skeleton->RND->drandom() < creativity) {
        // Create a new vertex...
        k = vertex_addition(vx);
        novum++;
      }
      else {
        // Grab an existing vertex...
        j = 0;
        nv = (signed) skeleton->events.size();
        success = false;
        do {
          k = skeleton->RND->irandom(nv);
          it = std::find(vx.begin(),vx.end(),k);
          if (it == vx.end()) {
            skeleton->events[k].active = true;
            vx_delta.insert(k);
            success = true;
            break;
          }
          ++j;
        } while(j < nv);
        if (success == false) {
          k = vertex_addition(vx);
          novum++;
        }
      }
      vx.insert(k);
    } while((signed) vx.size() < (1+d));
  }
#ifdef VERBOSE
  std::cout << "Created a " << d << "-simplex with " << novum << " new vertices." << std::endl;
#endif
  skeleton->simplex_addition(vx,vx_delta);
  return true;
}

bool Spacetime::expansion(int base)
{
  int n = skeleton->vertex_dimension(base);
  if (n == Complex::ND) return false;
  int i,u,d,m,vtx[2];
  std::set<int> vx,s,N;
  std::set<int>::const_iterator it;
  SYNARMOSMA::hash_map::const_iterator qt;
  Simplex S;
  const int ne = (signed) skeleton->simplices[1].size();
  double tau = std::abs(skeleton->events[base].deficiency) + 25.0;

  d = int(double(Complex::ND-1-n)*1.0/(1.0 + std::exp(tau)) + double(1 + n));

  vx.insert(base);
  for(i=0; i<d; ++i) {
    // Create a new vertex...
    m = vertex_addition(base);
    vx.insert(m);
  }
#ifdef VERBOSE
  std::cout << "Created a " << d << "-simplex with base " << base << std::endl;
#endif
  skeleton->simplex_addition(vx,vx_delta);
  for(i=0; i<ne; ++i) {
    if (!skeleton->active_simplex(1,i)) continue;
    skeleton->simplices[1][i].get_vertices(vtx);
    if (vtx[0] != base && vtx[1] != base) continue;
    u = (vtx[0] == base) ? vtx[1] : vtx[0];
    if (skeleton->events[u].deficiency < -std::numeric_limits<double>::epsilon()) N.insert(u);
  }

  for(it=N.begin(); it!=N.end(); ++it) {
    n = *it;
    m = skeleton->RND->irandom(vx);
    s.clear();
    s.insert(base);
    s.insert(n);
    qt = skeleton->index_table[1].find(s);
    skeleton->simplices[1][qt->second].active = false;
    S.initialize(m,n);
    qt = skeleton->index_table[1].find(S.vertices);
    if (qt == skeleton->index_table[1].end()) {
      skeleton->simplices[1].push_back(S);
      skeleton->index_table[1][S.vertices] = (signed) skeleton->simplices[1].size() - 1;
    }
    else {
      skeleton->simplices[1][qt->second].active = true;
    }
  }
  return true;
}

void Spacetime::vertex_fusion(int n1,int n2)
{
  int i,j,l,k,m,n,im1,in1;
  std::set<int> vx,duplicate;
  std::set<int>::reverse_iterator it;
  std::vector<int> mindex;
  const int ulimit = skeleton->dimension();
  std::vector<int>* mutation = new std::vector<int>[ulimit+1];

  // Do the actual vertex swap...
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
      vx = skeleton->simplices[i][in1].vertices;
      for(l=0; l<n; ++l) {
        if (l == in1) continue;
        if (skeleton->simplices[i][l].vertices == vx) {
          duplicate.insert(in1);
          skeleton->simplices[i][l].active = true;
        }
      }
    }
    for(it=duplicate.rbegin(); it!=duplicate.rend(); ++it) {
      in1 = *it;
      skeleton->simplices[i].erase(skeleton->simplices[i].begin() + in1);
    }
    duplicate.clear();
  }
  skeleton->events[n1].active = true;
  skeleton->events[n2].active = false;
  delete[] mutation;

  if (!skeleton->active_event(n2)) {
    if (!skeleton->events[n2].zero_energy()) {
      skeleton->events[n1].increment_energy(skeleton->events[n2].get_energy());
      skeleton->events[n2].nullify_energy();
    }
  }

  // Then recalculate the skeleton->index_table hash map..
  for(i=1; i<=ulimit; ++i) {
    skeleton->index_table[i].clear();
    m = (signed) skeleton->simplices[i].size();
    for(j=0; j<m; ++j) {
      skeleton->simplices[i][j].entourage.clear();
      skeleton->index_table[i][skeleton->simplices[i][j].vertices] = j;
    }
  }
  for(i=0; i<(signed) skeleton->events.size(); ++i) {
    skeleton->events[i].entourage.clear();
  }

  // Then recalculate the entourages...
  skeleton->compute_entourages();
  skeleton->compute_neighbours();
}

void Spacetime::superposition_fusion(std::set<int>& vmodified)
{
  int i,j,m,nf,v1,v2,nfail = 0,nfused = 0;
  double delta,pfusion = 0.0,alpha = 0.1;
  std::vector<int> modified;
  std::vector<std::pair<int,int> > candidates;
  const double na = double(skeleton->cardinality(0));
  const int nv = (signed) skeleton->events.size();
  const int ulimit = skeleton->dimension();

  // Vertex fusion - if two vertices are close enough together, they
  // should coalesce.
  for(i=0; i<nv; ++i) {
    modified.push_back(0);
    if (!skeleton->active_event(i)) continue;
    if (!skeleton->events[i].zero_energy()) continue;
    for(j=1+i; j<nv; ++j) {
      if (!skeleton->active_event(j)) continue;
      if (!skeleton->events[j].zero_energy()) continue;
      delta = geometry->get_squared_distance(i,j,false);
      if (delta < alpha) candidates.push_back(std::pair<int,int>(i,j));
    }
  }
  if (candidates.empty()) return;

  // Now we need to carry out the fusions...
  nf = (signed) candidates.size();
  do {
    i = skeleton->RND->irandom(nf);
    v1 = candidates[i].first;
    v2 = candidates[i].second;
    if (modified[v1] == 1 || modified[v2] == 1) {
      nfail++;
      continue;
    }
    modified[v1] = 1;
    modified[v2] = 1;
#ifdef VERBOSE
    std::cout << "Fusing vertices " << v2 << " => " << v1 << " via superposition" << std::endl;
#endif
    vertex_fusion(v1,v2);
    vmodified.insert(v1);
    nfused++;
    pfusion = double(nfused)/na;
  } while(pfusion < 0.05 && nfused < nf && nfail < 20);

  // Then recalculate the skeleton->index_table hash map..
  for(i=1; i<=ulimit; ++i) {
    skeleton->index_table[i].clear();
    m = (signed) skeleton->simplices[i].size();
    for(j=0; j<m; ++j) {
      skeleton->simplices[i][j].entourage.clear();
      skeleton->index_table[i][skeleton->simplices[i][j].vertices] = j;
    }
  }
  for(i=0; i<nv; ++i) {
    skeleton->events[i].entourage.clear();
  }

  // Then recalculate the entourages...
  skeleton->compute_entourages();
  skeleton->compute_neighbours();
}

void Spacetime::superposition_fission(std::set<int>& vmodified)
{
  int i,n,nc;
  std::set<int>::const_iterator it;
  Event vx;
  Simplex S;
  const int nv = (signed) skeleton->events.size();
  const int fmax = int(1.05*double(nv));

  // Finally, the opposite possibility - that a given vertex might undergo
  // spontaneous fission, creating a new vertex in its immediate vicinity...
  nc = nv;
  do {
    n = skeleton->RND->irandom(nc);
    if (!skeleton->active_event(n)) continue;
    if (skeleton->RND->drandom() > 0.01) continue;
    vx = skeleton->events[n];
    for(it=vx.neighbours.begin(); it!=vx.neighbours.end(); ++it) {
      i = *it;
      skeleton->events[i].neighbours.insert(nc);
      vmodified.insert(i);
      S.initialize(nc,i);
      skeleton->simplices[1].push_back(S);
      skeleton->index_table[1][S.vertices] = (signed) skeleton->simplices[1].size() - 1;
      S.clear();
    }
    vmodified.insert(nc);
    geometry->vertex_addition(n);
    skeleton->events.push_back(vx);
    nc++;
#ifdef VERBOSE
    std::cout << "Created vertex " << nc - 1 << " from " << n << std::endl;
#endif
  } while(nc < fmax);
}

bool Spacetime::vertex_twist()
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
  std::stringstream s;
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
  skeleton->events[n2].active = false;
  vx_delta.insert(n2);
  if (!skeleton->active_event(n2)) {
    if (!skeleton->events[n2].zero_energy()) {
      skeleton->events[n1].increment_energy(events[n2].get_energy());
      skeleton->events[n2].nullify_energy();
    }
  }
  vertex_fusion(n1,n2);
  if (first) first = false;
  return true;
}

bool Spacetime::perforation(int base,int d)
{
  int i,j,k,n,nd;
  bool good,found;
  std::set<int> vx,candidates;

  if (base >= 0) {
    // This call of the perforation operator is localized
    int n1 = skeleton->vertex_dimension(base);
    if (n1 < 2) return false;
    d = skeleton->RND->irandom(2,n1);
    nd = (signed) skeleton->simplices[d].size();
    for(i=0; i<nd; ++i) {
      if (!skeleton->active_simplex(d,i)) continue;
      if (!skeleton->simplices[d][i].contains(base)) continue;
      // We need to check now if each of this simplex's 1+d faces is
      // also the face of another d-simplex
      good = true;
      for(j=0; j<1+d; ++j) {
        found = false;
        vx = skeleton->simplices[d][i].faces[j];
        for(k=0; k<nd; ++k) {
          if (k == i || !skeleton->active_simplex(d,k)) continue;
          if (skeleton->simplices[d][k].face(vx)) {
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
    if (d < 2 || d >= Complex::ND) return false;
    if (skeleton->simplices[d].empty()) return false;
    nd = (signed) skeleton->simplices[d].size();
    for(i=0; i<nd; ++i) {
      if (!skeleton->active_simplex(d,i)) continue;
      // We need to check now if each of this simplex's 1+d faces is
      // also the face of another d-simplex
      good = true;
      for(j=0; j<1+d; ++j) {
        found = false;
        vx = skeleton->simplices[d][i].faces[j];
        for(k=0; k<nd; ++k) {
          if (k == i || !skeleton->active_simplex(d,k)) continue;
          if (skeleton->simplices[d][k].face(vx)) {
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
  n = skeleton->RND->irandom(candidates);
  skeleton->simplex_deletion(d,n);
  return true;
}

bool Spacetime::deflation(int base)
{
  int d = skeleton->vertex_dimension(base);
  if (d < 2) return false;
  int i,n,dw = 1;
  std::set<int> candidates;
  if (d > 2) dw = skeleton->RND->irandom(1,d);
  const int m = (signed) skeleton->simplices[dw].size();
  for(i=0; i<m; ++i) {
    if (!skeleton->active_simplex(dw,i)) continue;
    if (skeleton->simplices[dw][i].contains(base)) candidates.insert(i);
  }
  n = skeleton->RND->irandom(candidates);
#ifdef VERBOSE
  std::cout << "Deflating a " << d << "-simplex with base " << base << " to a " << dw << "-simplex." << std::endl;
#endif
  skeleton->simplex_deletion(dw,n);
  return true;
}

bool Spacetime::vertex_deletion(int n)
{
  if (!skeleton->active_event(n)) return false;
  int i,vx[2],ne = (signed) skeleton->simplices[1].size();
  skeleton->events[n].active = false;
  for(i=0; i<ne; ++i) {
    if (!simplices[1][i].active) continue;
    skeleton->simplices[1][i].get_vertices(vx);
    if (vx[0] == n || vx[1] == n) skeleton->simplex_deletion(1,i);
  }
  return true;  
}

int Spacetime::vertex_addition(const std::vector<double>& xc)
{
  int n = (signed) skeleton->events.size();
  Event vt;

  geometry->vertex_addition(xc);
  skeleton->events.push_back(vt);

  return n;
}

int Spacetime::vertex_addition(const std::set<int>& antecedents)
{
  int n = (signed) skeleton->events.size();
  Event vt;

  geometry->vertex_addition(antecedents);
  skeleton->events.push_back(vt);

  return n;
}

int Spacetime::vertex_addition(int base)
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

bool Spacetime::inflation(int base,double creativity)
{
  // Performs an inflation on the vertex base with colour sheet
  int i,k,n1,na,delta,its = 0;
  Event vt;
  bool success;
  std::set<int>::const_iterator it;
  SYNARMOSMA::hash_map::const_iterator qt;
  std::set<int> vx,candidates;

  if (base >= 0) {
    n1 = skeleton->vertex_dimension(base);
    if (n1 == Complex::ND) return false;
    int j,vtx[2];
    std::set<int>::const_iterator jt;
    std::set<int> M,N;
    double tau;
    const int ne = (signed) skeleton->simplices[1].size();

    // We want to inflate this to a d-simplex, where d is between
    // 1+n1 and ND - the greater the magnitude of v's deficiency,
    // the higher the dimension d should be...
    tau = std::abs(skeleton->events[base].deficiency) + 25.0;
    delta = int(double(Complex::ND-1-n1)*1.0/(1.0 + std::exp(tau)) + double(1 + n1));

    for(i=0; i<(signed) skeleton->simplices[n1].size(); ++i) {
      if (!skeleton->active_simplex(n1,i)) continue;
      if (skeleton->simplices[n1][i].contains(base)) candidates.insert(i);
    }
    if (candidates.empty()) return false;

    for(i=0; i<ne; ++i) {
      if (!skeleton->active_simplex(1,i)) continue;
      skeleton->simplices[1][i].get_vertices(vtx);
      if (base != vtx[0] && base != vtx[1]) continue;
      j = (vtx[0] == base) ? vtx[1] : vtx[0];
      N.insert(j);
    }
    vx = skeleton->simplices[n1][skeleton->RND->irandom(candidates)].vertices;
    for(it=N.begin(); it!=N.end(); ++it) {
      jt = std::find(vx.begin(),vx.end(),*it);
      if (jt == vx.end()) {
        if (skeleton->events[*it].deficiency < -std::numeric_limits<double>::epsilon()) M.insert(*it);
      }
    }
    if (M.empty()) creativity = 1.0;
    do {
      its++;
      if (its == Complex::ND) creativity = 1.0;
      if (skeleton->RND->drandom() < creativity) {
        // Create a new vertex...
        k = vertex_addition(base);
      }
      else {
        k = skeleton->RND->irandom(M);
      }
      vx.insert(k);
    } while((signed) vx.size() < (1 + delta));
#ifdef VERBOSE
    std::cout << "Inflated a " << n1 << "-simplex into a " << delta << "-simplex based on vertex " << base << std::endl;
#endif
  }
  else {
    if (skeleton->dimension() == 0) {
      n1 = 0;
    }
    else if (skeleton->dimension() < (signed) geometry->dimension()) {
      n1 = 1;
    }
    else {
      n1 = 1 + skeleton->RND->irandom(skeleton->dimension());
    }
    for(i=0; i<(signed) skeleton->events.size(); ++i) {
      if (skeleton->active_event(i)) candidates.insert(i);
    }
    na = (signed) candidates.size();
    if (n1 == 0) {
      delta = skeleton->RND->irandom(2,2*geometry->dimension());
      vx.insert(skeleton->RND->irandom(candidates));
    }
    else {
      candidates.clear();
      for(i=0; i<(signed) skeleton->simplices[n1].size(); ++i) {
        if (skeleton->active_simplex(n1,i)) candidates.insert(i);
      }
      if (candidates.empty()) return false;
      vx = skeleton->simplices[n1][skeleton->RND->irandom(candidates)].vertices;
      delta = n1 + skeleton->RND->irandom(1,geometry->dimension());
    }
    if (double(delta)/double(na) > 0.25) creativity = 1.0;
    if (double(vx.size())/double(na) > 0.9) creativity = 1.0;
    do {
      if (skeleton->RND->drandom() < creativity) {
        // Create a new vertex...
        k = vertex_addition(vx);
      }
      else {
        // Grab an existing vertex...
        its = 0;
        success = false;
        do {
          // We should perhaps alter this to favour vertices that are
          // few hops away
          k = skeleton->RND->irandom(skeleton->events.size());
          if (!skeleton->active_event(k)) continue;
          it = std::find(vx.begin(),vx.end(),k);
          if (it == vx.end()) {
            success = true;
            break;
          }
          its++;
        } while(its < (signed) skeleton->events.size());
        if (success == false) {
          k = (signed) skeleton->events.size();
          skeleton->events.push_back(vt);
        }
      }
      vx.insert(k);
    } while((signed) vx.size() < (1+delta));
#ifdef VERBOSE
    std::cout << "Inflated a " << n1 << "-simplex into a " << delta << "-simplex" << std::endl;
#endif
  }
  skeleton->simplex_addition(vx,vx_delta);
  return true;
}
