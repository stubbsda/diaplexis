#include "sheet.h"

using namespace DIAPLEXIS;

Sheet::Sheet()
{
  H = new SYNARMOSMA::Homology(SYNARMOSMA::Homology::Field::mod2,SYNARMOSMA::Homology::Method::native);
  pi1 = new SYNARMOSMA::Homotopy;
}

Sheet::Sheet(int n,SYNARMOSMA::Homology::Field f,SYNARMOSMA::Homology::Method m)
{
  H = new SYNARMOSMA::Homology(f,m);
  pi1 = new SYNARMOSMA::Homotopy;

  index = n;
}

Sheet::Sheet(int n,int p,SYNARMOSMA::Homology::Field f,SYNARMOSMA::Homology::Method m)
{
  H = new SYNARMOSMA::Homology(f,m);
  pi1 = new SYNARMOSMA::Homotopy;

  index = n;
  parent = p;
}

Sheet::Sheet(const Sheet& source)
{
  index = source.index;
  parent = source.parent;
  voice = source.voice;
  sleep = source.sleep;
  fertility = source.fertility;
  drowsiness = source.drowsiness;
  score_allocated = source.score_allocated;
  hyphantic_ops = source.hyphantic_ops;
  H = new SYNARMOSMA::Homology(*source.H);
  pi1 = new SYNARMOSMA::Homotopy(*source.pi1);
  pseudomanifold = source.pseudomanifold;
  boundary = source.boundary;
  orientable = source.orientable;
  if (score_allocated > 0) {
    hyphantic_notes = new std::vector<int>[score_allocated];
    for(int i=0; i<score_allocated; ++i) {
      hyphantic_notes[i] = source.hyphantic_notes[i];
    }
  }
}

Sheet& Sheet::operator =(const Sheet& source)
{
  if (this == &source) return *this;

  delete H;
  delete pi1;
  if (score_allocated > 0) delete[] hyphantic_notes;

  index = source.index;
  parent = source.parent;
  voice = source.voice;
  sleep = source.sleep;
  fertility = source.fertility;
  drowsiness = source.drowsiness;
  score_allocated = source.score_allocated;
  hyphantic_ops = source.hyphantic_ops;
  H = new SYNARMOSMA::Homology(*source.H);
  pi1 = new SYNARMOSMA::Homotopy(*source.pi1);
  pseudomanifold = source.pseudomanifold;
  boundary = source.boundary;
  orientable = source.orientable;

  if (score_allocated > 0) {
    hyphantic_notes = new std::vector<int>[score_allocated];
    for(int i=0; i<score_allocated; ++i) {
      hyphantic_notes[i] = source.hyphantic_notes[i];
    }    
  }

  return *this;
}

Sheet::~Sheet()
{
  delete H;
  delete pi1;
  if (score_allocated > 0) delete[] hyphantic_notes;
}

void Sheet::clear()
{
  hyphantic_ops = "";
  sleep = 0;
  parent = -1;
  voice = -1;
  index = -1;
  fertility = 0.0;
  drowsiness = 0.0;
  H->clear();
  pi1->clear();
  pseudomanifold = false;
  boundary = false;
  orientable = false;
  if (score_allocated > 0) delete[] hyphantic_notes;
  score_allocated = -1;  
}

void Sheet::set_topology(const SYNARMOSMA::Homology* K,const SYNARMOSMA::Homotopy* p,bool pm,bool bd,bool orient)
{
  delete H;
  delete pi1;

  H = new SYNARMOSMA::Homology(*K);
  pi1 = new SYNARMOSMA::Homotopy(*p); 
  pseudomanifold = pm;
  boundary = bd;
  orientable = orient;
}

void Sheet::parse_music_score(int max_iter,std::string& filename)
{
  int i,n,v,its,nsilent = 0;
  std::string line;
  std::vector<std::string> elements;
  std::ifstream mscore(filename);
  const char delimiter = '/';

  score_allocated = 1 + max_iter;
  hyphantic_notes = new std::vector<int>[score_allocated];

  // Now read the measure that corresponds to this iteration and sheet...
  while(mscore.good()) {
    getline(mscore,line);
    // If the line is empty or doesn't contain a forward slash, ignore it...
    if (line.empty()) continue;
    if (line.find(delimiter) == std::string::npos) continue;
    // Tokenize the line at the forward slash...
    SYNARMOSMA::split(line,delimiter,elements);
    its = std::stoi(elements[0]) - 1;
    if (its > max_iter) continue;
    // So this is a line for this relaxation step, check if it is the right sheet/voice...
    v = std::stoi(elements[1]);
    if (v != voice) continue;
    n = std::stoi(elements[2]);
    hyphantic_notes[its].push_back(n);
  }
  // Check for iterations during which there is no hyphansis...
  for(i=0; i<=max_iter; ++i) {
    if (hyphantic_notes[i].empty()) nsilent++;
  }
  if (nsilent > 0) std::cout << "Warning: For sheet " << index << ", there are " << nsilent << " relaxation steps during which there will be no hyphantic operations!" << std::endl;
  
  // Close the score file
  mscore.close();
}

int Sheet::serialize(std::ofstream& s) const
{
  int i,count = 0,n = (signed) hyphantic_ops.size();

  s.write((char*)(&index),sizeof(int)); count += sizeof(int);
  s.write((char*)(&parent),sizeof(int)); count += sizeof(int);
  s.write((char*)(&voice),sizeof(int)); count += sizeof(int);
  s.write((char*)(&score_allocated),sizeof(int)); count += sizeof(int);
  s.write((char*)(&sleep),sizeof(int)); count += sizeof(int);
  s.write((char*)(&fertility),sizeof(double)); count += sizeof(double);
  s.write((char*)(&drowsiness),sizeof(double)); count += sizeof(double);
  s.write((char*)(&n),sizeof(int)); count += sizeof(int);
  for(i=0; i<n; ++i) {
    s.write((char*)(&hyphantic_ops[i]),sizeof(char)); count += sizeof(char);
  }
  // Now the algebraic properties...
  count += H->serialize(s);
  count += pi1->serialize(s);
  s.write((char*)(&pseudomanifold),sizeof(bool)); count += sizeof(bool);
  s.write((char*)(&boundary),sizeof(bool)); count += sizeof(bool);
  s.write((char*)(&orientable),sizeof(bool)); count += sizeof(bool);
  if (score_allocated > 0) {
    int j,m;
    
    for(i=0; i<score_allocated; ++i) {
      n = (signed) hyphantic_notes[i].size();
      s.write((char*)(&n),sizeof(int));
      for(j=0; j<n; ++j) {
        m = hyphantic_notes[i][j];
        s.write((char*)(&m),sizeof(int));
      }
    }

  }
  return count;
}

int Sheet::deserialize(std::ifstream& s) 
{
  int i,n,count = 0;
  char c;

  clear();

  s.read((char*)(&index),sizeof(int)); count += sizeof(int);
  s.read((char*)(&parent),sizeof(int)); count += sizeof(int);
  s.read((char*)(&voice),sizeof(int)); count += sizeof(int);
  s.read((char*)(&score_allocated),sizeof(int)); count += sizeof(int);
  s.read((char*)(&sleep),sizeof(int)); count += sizeof(int);
  s.read((char*)(&fertility),sizeof(double)); count += sizeof(double);
  s.read((char*)(&drowsiness),sizeof(double)); count += sizeof(double);
  s.read((char*)(&n),sizeof(int)); count += sizeof(int);
  for(i=0; i<n; ++i) {
    s.read((char*)(&c),sizeof(char)); count += sizeof(char);
    hyphantic_ops += c;
  }
  // Now the algebraic properties...
  count += H->deserialize(s);
  count += pi1->deserialize(s);
  s.read((char*)(&pseudomanifold),sizeof(bool)); count += sizeof(bool);
  s.read((char*)(&boundary),sizeof(bool)); count += sizeof(bool);
  s.read((char*)(&orientable),sizeof(bool)); count += sizeof(bool);
  if (score_allocated > 0) {
    int j,q;
    std::vector<int> bar;

    hyphantic_notes = new std::vector<int>[score_allocated];

    for(i=0; i<score_allocated; ++i) {
      s.read((char*)(&n),sizeof(int));
      for(j=0; j<n; ++j) {
        s.read((char*)(&q),sizeof(int));
        bar.push_back(q);
      }
      hyphantic_notes[i] = bar;
      bar.clear();
    }
  }
  return count;
}


