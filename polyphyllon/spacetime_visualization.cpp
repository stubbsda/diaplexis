#include "spacetime.h"

using namespace DIAPLEXIS;

void Spacetime::get_energy_extrema(std::pair<double,double>& output) const
{
  int i;
  double alpha,u_ex = 0.0,l_ex = 0.0;
  const int nv = (signed) skeleton->events.size();

  for(i=0; i<nv; ++i) {
    if (!skeleton->events[i].active()) continue;
    u_ex = skeleton->events[i].get_energy();
    break;
  }
  l_ex = u_ex;
  for(i=0; i<nv; ++i) {
    if (!skeleton->events[i].active()) continue;
    alpha = skeleton->events[i].get_energy();
    if (u_ex < alpha) u_ex = alpha;
    if (l_ex > alpha) l_ex = alpha;
  }
  output.first = u_ex;
  output.second = l_ex;
}

void Spacetime::get_deficiency_extrema(std::pair<double,double>& output) const
{
  int i;
  double alpha,u_ex = 0.0,l_ex = 0.0;
  const int nv = (signed) skeleton->events.size();

  for(i=0; i<nv; ++i) {
    if (!skeleton->events[i].active()) continue;
    u_ex = skeleton->events[i].get_deficiency();
    break;
  }
  l_ex = u_ex;
  for(i=0; i<nv; ++i) {
    if (!skeleton->events[i].active()) continue;
    alpha = skeleton->events[i].get_deficiency();
    if (u_ex < alpha) u_ex = alpha;
    if (l_ex > alpha) l_ex = alpha;
  }
  output.first = u_ex;
  output.second = l_ex;
}

void Spacetime::compute_colours(std::vector<unsigned char>& chi,bool use_sheets,bool use_energy) const
{
  int i,j;
  const int nv = (signed) skeleton->events.size();
  const int ne = (signed) skeleton->simplices[1].size();

  if (use_sheets) {
    if (codex.size() == 1) {
      // A monocosmos, so everything is black...
      for(i=0; i<nv; ++i) {
        if (!skeleton->events[i].active()) {
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
        if (!skeleton->simplices[1][i].active()) {
          chi.push_back(255);
          chi.push_back(255);
          chi.push_back(255);
          continue;
        }
        chi.push_back(0);
        chi.push_back(0);
        chi.push_back(0);
      }
    }
    else if (codex.size() == 2) {
      // Very special case: 2^2 - 1 = 3, the number of elements of
      // an RGB colour vector
      for(i=0; i<nv; ++i) {
        if (!skeleton->events[i].active()) {
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
        if (!skeleton->simplices[1][i].active()) {
          chi.push_back(255);
          chi.push_back(255);
          chi.push_back(255);
          continue;
        }
        if (skeleton->simplices[1][i].active(0) && skeleton->simplices[1][i].active(1)) {
          chi.push_back(0);
          chi.push_back(0);
          chi.push_back(255);
        }
        else {
          if (skeleton->simplices[1][i].active(0)) {
            chi.push_back(255);
            chi.push_back(0);
            chi.push_back(0);
          }
          else {
            chi.push_back(0);
            chi.push_back(255);
            chi.push_back(0);
          }
        }
      }
    }
    else {
      // The generic case: construct a basis for the colour space...
      const int nt = (signed) codex.size();
      int p,q,divisors[3],level[3],mlevel[3];
      float cbasis[3*nt],cvalues[3];

      for(i=0; i<nv; ++i) {
        if (!skeleton->events[i].active()) {
          chi.push_back(255);
          chi.push_back(255);
          chi.push_back(255);
          continue;
        }
        chi.push_back(0);
        chi.push_back(0);
        chi.push_back(0);
      }

      p = nt/3;
      q = nt%3;
      for(i=0; i<q; ++i) {
        mlevel[i] = p + 2;
      }
      for(i=q; i<3; ++i) {
        mlevel[i] = p + 1;
      }
      level[0] = 1;
      level[1] = 1;
      level[2] = 1;
      for(i=0; i<3*nt; ++i) {
        cbasis[i] = 0.0;
      }
      for(i=0; i<nt; ++i) {
        for(j=0; j<3; ++j) {
          if (i%3 == j) {
            cbasis[3*i+j] = float(level[j])/float(mlevel[j]);
            level[j] += 1;
          }
        }
      }
      for(i=0; i<ne; ++i) {
        if (!skeleton->simplices[1][i].active()) {
          chi.push_back(255);
          chi.push_back(255);
          chi.push_back(255);
          continue;
        }
        cvalues[0] = 0.0;
        cvalues[1] = 0.0;
        cvalues[2] = 0.0;
        divisors[0] = 0;
        divisors[1] = 0;
        divisors[2] = 0;
        for(j=0; j<nt; ++j) {
          if (skeleton->simplices[1][i].active(j)) {
            cvalues[0] += cbasis[3*j];
            cvalues[1] += cbasis[3*j+1];
            cvalues[2] += cbasis[3*j+2];
            divisors[j%3] += 1;
          }
        }
        for(j=0; j<3; ++j) {
          if (divisors[j] > 1) cvalues[j] = cvalues[j]/float(divisors[j]);
          chi.push_back(int(255.0*cvalues[j]));
        }
      }
    }
  }
  else {
    // Here rather than using the sheet ubiquity for colouring everything
    // we will make use the energy or deficiency values of the vertices to
    // set the vertex colour according to a thermal palette.
    int vx[2];
    unsigned char out[3];
    double x_max,x_min,delta,theta,xvalue[nv];
    std::pair<double,double> pr;

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
        if (!skeleton->events[i].active()) {
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
        if (!skeleton->simplices[1][i].active()) {
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
      if (!skeleton->events[i].active()) {
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
      if (!skeleton->simplices[1][i].active()) {
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
}

void Spacetime::export_visual_data(std::vector<float>& vcoords,std::vector<int>& evertex,int* vdata,int sheet) const
{
  int i,j = 0,D,vx[2];
  std::vector<double> x;
  const int nv = (signed) skeleton->events.size();
  const int ne = (signed) skeleton->simplices[1].size();
  int offset[nv];

  D = geometry->compute_coordinates(x);
  vcoords.clear();
  evertex.clear();

  vdata[0] = skeleton->cardinality(0,sheet);
  vdata[1] = skeleton->cardinality(1,sheet);

  if (sheet == -1) {
    for(i=0; i<nv; ++i) {
      offset[i] = -1;
      if (!skeleton->events[i].active()) continue;
      offset[i] = j; j++;
    }
    if (D == 1) {
      for(i=0; i<nv; ++i) {
        if (!skeleton->events[i].active()) continue;
        vcoords.push_back(float(x[offset[i]]));
        vcoords.push_back(0.0);
        vcoords.push_back(0.0);
      }
    }
    else if (D == 2) {
      for(i=0; i<nv; ++i) {
        if (!skeleton->events[i].active()) continue;
        vcoords.push_back(float(x[2*offset[i]]));
        vcoords.push_back(float(x[2*offset[i]+1]));
        vcoords.push_back(0.0);
      }
    }
    else if (D >= 3) {
      for(i=0; i<nv; ++i) {
        if (!skeleton->events[i].active()) continue;
        vcoords.push_back(float(x[D*offset[i]]));
        vcoords.push_back(float(x[D*offset[i]+1]));
        vcoords.push_back(float(x[D*offset[i]+2]));
      }
    }
    // Now the neighbour table...
    for(i=0; i<ne; ++i) {
      if (!skeleton->simplices[1][i].active()) continue;
      skeleton->simplices[1][i].get_vertices(vx);
      evertex.push_back(offset[vx[0]]);
      evertex.push_back(offset[vx[1]]);
    }
  }
  else {
    for(i=0; i<nv; ++i) {
      offset[i] = -1;
      if (!skeleton->events[i].active(sheet)) continue;
      offset[i] = j; j++;
    }
    if (D == 1) {
      for(i=0; i<nv; ++i) {
        if (!skeleton->events[i].active(sheet)) continue;
        vcoords.push_back(float(x[offset[i]]));
        vcoords.push_back(0.0);
        vcoords.push_back(0.0);
      }
    }
    else if (D == 2) {
      for(i=0; i<nv; ++i) {
        if (!skeleton->events[i].active(sheet)) continue;
        vcoords.push_back(float(x[2*offset[i]]));
        vcoords.push_back(float(x[2*offset[i]+1]));
        vcoords.push_back(0.0);
      }
    }
    else if (D >= 3) {
      for(i=0; i<nv; ++i) {
        if (!skeleton->events[i].active(sheet)) continue;
        vcoords.push_back(float(x[D*offset[i]]));
        vcoords.push_back(float(x[D*offset[i]+1]));
        vcoords.push_back(float(x[D*offset[i]+2]));
      }
    }
    // Now the neighbour table...
    for(i=0; i<ne; ++i) {
      if (!skeleton->simplices[1][i].active(sheet)) continue;
      skeleton->simplices[1][i].get_vertices(vx);
      evertex.push_back(offset[vx[0]]);
      evertex.push_back(offset[vx[1]]);
    }
  }
}

void Spacetime::export_visual_data(std::vector<float>& colours,std::vector<float>& vcoords,std::vector<int>& evertex,int* vdata,bool use_energy) const
{
  int i,j = 0,D,vx[2];
  std::vector<double> x;
  std::vector<unsigned char> chi;
  const int nv = (signed) skeleton->events.size();
  const int ne = (signed) skeleton->simplices[1].size();
  int offset[nv];

  compute_colours(chi,false,use_energy);
  D = geometry->compute_coordinates(x);

  colours.clear();
  vcoords.clear();
  evertex.clear();

  vdata[0] = nv;
  vdata[1] = ne;

  for(i=0; i<nv; ++i) {
    offset[i] = -1;
    if (!skeleton->events[i].active()) continue;
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
      if (!skeleton->events[i].active()) {
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
      if (!skeleton->events[i].active()) {
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
      if (!skeleton->events[i].active()) {
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
    if (!skeleton->simplices[1][i].active()) {
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

