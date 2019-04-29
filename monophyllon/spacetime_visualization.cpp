#include "spacetime.h"

using namespace DIAPLEXIS;

void Spacetime::get_energy_extrema(std::pair<double,double>& output) const
{
  unsigned int i;
  double alpha,u_ex = 0.0,l_ex = 0.0;
  const unsigned int nv = skeleton->events.size();

  for(i=0; i<nv; ++i) {
    if (!skeleton->active_event(i)) continue;
    u_ex = skeleton->events[i].get_energy();
    break;
  }
  l_ex = u_ex;
  for(i=0; i<nv; ++i) {
    if (!skeleton->active_event(i)) continue;
    alpha = skeleton->events[i].get_energy();
    if (u_ex < alpha) u_ex = alpha;
    if (l_ex > alpha) l_ex = alpha;
  }
  output.first = u_ex;
  output.second = l_ex;
}

void Spacetime::get_deficiency_extrema(std::pair<double,double>& output) const
{
  unsigned int i;
  double t,u_ex = 0.0,l_ex = 0.0;
  const unsigned int nv = skeleton->events.size();

  for(i=0; i<nv; ++i) {
    if (!skeleton->active_event(i)) continue;
    u_ex = skeleton->events[i].get_deficiency();
    break;
  }
  l_ex = u_ex;
  for(i=0; i<nv; ++i) {
    if (!skeleton->active_event(i)) continue;
    t = skeleton->events[i].get_deficiency();
    if (u_ex < t) u_ex = t;
    if (l_ex > t) l_ex = t;
  }
  output.first = u_ex;
  output.second = l_ex;
}

void Spacetime::compute_colours(std::vector<unsigned char>& chi,bool use_energy) const
{
  int i,vx[2];
  unsigned char out[3];
  double x_max,x_min,delta,theta;
  std::pair<double,double> pr;
  const int nv = (signed) skeleton->events.size();
  const int ne = (signed) skeleton->simplices[1].size();
  double xvalue[nv];

  // Here we will make use the energy or deficiency values of the vertices to
  // set the vertex colour according to a thermal palette.
  if (use_energy) {
    get_energy_extrema(pr);
  }
  else {
    get_deficiency_extrema(pr);
  }
  x_max = pr.first; x_min = pr.second;
  delta = x_max - x_min;
  if (delta < 0.05) {
    // In this case just paint everything black...
    for(i=0; i<nv; ++i) {
      if (!skeleton->active_event(i)) {
        chi.push_back(255);
        chi.push_back(255);
        chi.push_back(255);
        continue;
      }
      chi.push_back(0);
      chi.push_back(0);
      chi.push_back(0);
    }
    for(i=0; i<ne; ++i) {
      if (!skeleton->active_simplex(1,i)) {
        chi.push_back(255);
        chi.push_back(255);
        chi.push_back(255);
        continue;
      }
      chi.push_back(0);
      chi.push_back(0);
      chi.push_back(0);
    }
    return;
  }

  for(i=0; i<nv; ++i) {
    if (!skeleton->active_event(i)) {
      chi.push_back(255);
      chi.push_back(255);
      chi.push_back(255);
      continue;
    }
    theta = (M_PI/2.0)*(xvalue[i] - x_min)/delta;
    SYNARMOSMA::RGB_intensity(theta,out);
    xvalue[i] = theta;
    chi.push_back(out[0]);
    chi.push_back(out[1]);
    chi.push_back(out[2]);
  }

  for(i=0; i<ne; ++i) {
    if (!skeleton->active_simplex(1,i)) {
      chi.push_back(255);
      chi.push_back(255);
      chi.push_back(255);
      continue;
    }
    skeleton->simplices[1][i].get_vertices(vx);
    theta = 0.5*(xvalue[vx[0]] + xvalue[vx[1]]);
    SYNARMOSMA::RGB_intensity(theta,out);
    chi.push_back(out[0]);
    chi.push_back(out[1]);
    chi.push_back(out[2]);
  }
}

void Spacetime::export_visual_data(std::vector<float>& vcoords,std::vector<int>& evertex,std::pair<int,int>& vdata) const
{
  int i,j = 0,D,vx[2];
  std::vector<double> x;
  const int nv = (signed) skeleton->events.size();
  const int ne = (signed) skeleton->simplices[1].size();
  int offset[nv];

  D = geometry->compute_coordinates(x);
  vcoords.clear();
  evertex.clear();

  vdata.first = skeleton->cardinality(0);
  vdata.second = skeleton->cardinality(1);

  for(i=0; i<nv; ++i) {
    offset[i] = -1;
    if (!skeleton->active_event(i)) continue;
    offset[i] = j; j++;
  }

  if (D == 1) {
    for(i=0; i<nv; ++i) {
      if (!skeleton->active_event(i)) continue;
      vcoords.push_back(float(x[offset[i]]));
      vcoords.push_back(0.0);
      vcoords.push_back(0.0);
    }
  }
  else if (D == 2) {
    for(i=0; i<nv; ++i) {
      if (!skeleton->active_event(i)) continue;
      vcoords.push_back(float(x[2*offset[i]]));
      vcoords.push_back(float(x[2*offset[i]+1]));
      vcoords.push_back(0.0);
    }
  }
  else if (D >= 3) {
    for(i=0; i<nv; ++i) {
      if (!skeleton->active_event(i)) continue;
      vcoords.push_back(float(x[D*offset[i]]));
      vcoords.push_back(float(x[D*offset[i]+1]));
      vcoords.push_back(float(x[D*offset[i]+2]));
    }
  }
  // Now the neighbour table...
  for(i=0; i<ne; ++i) {
    if (!skeleton->active_simplex(1,i)) continue;
    skeleton->simplices[1][i].get_vertices(vx);
    evertex.push_back(offset[vx[0]]);
    evertex.push_back(offset[vx[1]]);
  }
}

void Spacetime::export_visual_data(std::vector<float>& colours,std::vector<float>& vcoords,std::vector<int>& evertex,std::pair<int,int>& vdata,bool use_energy) const
{
  int i,j = 0,D,vx[2];
  std::vector<double> x;
  std::vector<unsigned char> chi;
  const int nv = (signed) skeleton->events.size();
  const int ne = (signed) skeleton->simplices[1].size();
  int offset[nv];

  compute_colours(chi,use_energy);
  D = geometry->compute_coordinates(x);

  colours.clear();
  vcoords.clear();
  evertex.clear();

  vdata.first = nv;
  vdata.second = ne;

  for(i=0; i<nv; ++i) {
    offset[i] = -1;
    if (!skeleton->active_event(i)) continue;
    offset[i] = j; j++;
  }

  for(i=0; i<3*nv; ++i) {
    colours.push_back(float(chi[i])/255.0);
  }
  for(i=0; i<3*ne; ++i) {
    colours.push_back(float(chi[3*nv+i])/255.0);
  }
  if (D == 1) {
    for(i=0; i<nv; ++i) {
      if (!skeleton->active_event(i)) {
        vcoords.push_back(0.0);
        vcoords.push_back(0.0);
        vcoords.push_back(0.0);
        continue;
      }
      vcoords.push_back(float(x[offset[i]]));
      vcoords.push_back(0.0);
      vcoords.push_back(0.0);
    }
  }
  else if (D == 2) {
    for(i=0; i<nv; ++i) {
      if (!skeleton->active_event(i)) {
        vcoords.push_back(0.0);
        vcoords.push_back(0.0);
        vcoords.push_back(0.0);
        continue;
      }
      vcoords.push_back(float(x[2*offset[i]]));
      vcoords.push_back(float(x[2*offset[i]+1]));
      vcoords.push_back(0.0);
    }
  }
  else if (D >= 3) {
    for(i=0; i<nv; ++i) {
      if (!skeleton->active_event(i)) {
        vcoords.push_back(0.0);
        vcoords.push_back(0.0);
        vcoords.push_back(0.0);
        continue;
      }
      vcoords.push_back(float(x[D*offset[i]]));
      vcoords.push_back(float(x[D*offset[i]+1]));
      vcoords.push_back(float(x[D*offset[i]+2]));
    }
  }
  for(i=0; i<ne; ++i) {
    if (!skeleton->active_simplex(1,i)) {
      vx[0] = 0;
      vx[1] = 0;
    }
    else {
      skeleton->simplices[1][i].get_vertices(vx);
    }
    evertex.push_back(vx[0]);
    evertex.push_back(vx[1]);
  }
}

