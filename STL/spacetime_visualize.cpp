#include "spacetime.h"

void Spacetime::compute_colours(std::vector<unsigned char>& chi,bool use_sheets,bool use_energy) const
{
  int i,j;
  const int nv = (signed) events.size();
  const int ne = (signed) simplices[1].size();

  if (use_sheets) {
    if (codex.size() == 1) {
      // A monocosmos, so everything is black...
      for(i=0; i<nv; ++i) {
        if (ghost(events[i].ubiquity)) {
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
        if (ghost(simplices[1][i].ubiquity)) {
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
        if (ghost(events[i].ubiquity)) {
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
        if (ghost(simplices[1][i].ubiquity)) {
          chi.push_back(255);
          chi.push_back(255);
          chi.push_back(255);
          continue;
        }
        if (simplices[1][i].ubiquity[0] == 1 && simplices[1][i].ubiquity[1] == 1) {
          chi.push_back(0);
          chi.push_back(0);
          chi.push_back(255);
        }
        else {
          if (simplices[1][i].ubiquity[0] == 1) {
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
        if (ghost(events[i].ubiquity)) {
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
        if (ghost(simplices[1][i].ubiquity)) {
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
          if (simplices[1][i].ubiquity[j] == 1) {
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
    double x,x_max,x_min,delta,theta,xvalue[nv];

    for(i=0; i<nv; ++i) {
      if (ghost(events[i].ubiquity)) continue;
      x_min = (use_energy) ? events[i].energy : events[i].deficiency;
      x_max = x_min;
      break;
    }

    if (use_energy) {
      for(i=0; i<nv; ++i) {
        if (ghost(events[i].ubiquity)) continue;
        x = events[i].energy;
        xvalue[i] = x;
        if (x > x_max) x_max = x;
        if (x < x_min) x_min = x;
      }      
    }
    else {
      for(i=0; i<nv; ++i) {
        if (ghost(events[i].ubiquity)) continue;
        x = events[i].deficiency;
        xvalue[i] = x;
        if (x > x_max) x_max = x;
        if (x < x_min) x_min = x;
      }
    }
    delta = x_max - x_min;
    if (delta < 0.05) {
      // In this case just paint everything black...
      for(i=0; i<nv; ++i) {
        if (ghost(events[i].ubiquity)) {
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
        if (ghost(simplices[1][i].ubiquity)) {
          chi.push_back(255);
          chi.push_back(255);
          chi.push_back(255);
          continue;
          continue;
        }
        chi.push_back(0);
        chi.push_back(0);
        chi.push_back(0);
      }
      return;
    }

    for(i=0; i<nv; ++i) {
      if (ghost(events[i].ubiquity)) {
        chi.push_back(255);
        chi.push_back(255);
        chi.push_back(255);
        continue;
      }
      theta = (M_PI/2.0)*(xvalue[i] - x_min)/delta;
      RGB_intensity(theta,out);
      xvalue[i] = theta;
      chi.push_back(out[0]);
      chi.push_back(out[1]);
      chi.push_back(out[2]);
    }

    for(i=0; i<ne; ++i) {
      if (ghost(simplices[1][i].ubiquity)) {
        chi.push_back(255);
        chi.push_back(255);
        chi.push_back(255);
        continue;
      }
      simplices[1][i].get_vertices(vx);      
      theta = 0.5*(xvalue[vx[0]] + xvalue[vx[1]]);
      RGB_intensity(theta,out);
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
  const int nv = (signed) events.size();
  const int ne = (signed) simplices[1].size();
  int offset[nv];

  D = geometry->compute_coordinates(x);
  vcoords.clear();
  evertex.clear();

  if (sheet == -1) {
    for(i=0; i<nv; ++i) {
      offset[i] = -1;
      if (ghost(events[i].ubiquity)) continue;
      offset[i] = j; j++;
    }
    vdata[0] = j;
    if (D == 1) {
      for(i=0; i<nv; ++i) {
        if (ghost(events[i].ubiquity)) continue;
        vcoords.push_back(float(x[i]));
        vcoords.push_back(0.0);
        vcoords.push_back(0.0);
      }
    }
    else if (D == 2) {
      for(i=0; i<nv; ++i) {
        if (ghost(events[i].ubiquity)) continue;
        vcoords.push_back(float(x[2*i]));
        vcoords.push_back(float(x[2*i+1]));
        vcoords.push_back(0.0);
      }
    }
    else if (D >= 3) {
      for(i=0; i<nv; ++i) {
        if (ghost(events[i].ubiquity)) continue;
        vcoords.push_back(float(x[D*i]));
        vcoords.push_back(float(x[D*i+1]));
        vcoords.push_back(float(x[D*i+2]));
      }
    }
    // Now the neighbour table...
    j = 0;
    for(i=0; i<ne; ++i) {
      if (ghost(simplices[1][i].ubiquity)) continue;
      simplices[1][i].get_vertices(vx);
      evertex.push_back(offset[vx[0]]);
      evertex.push_back(offset[vx[1]]);
      ++j;
   }
   vdata[1] = j;
  }
  else {
    for(i=0; i<nv; ++i) {
      offset[i] = -1;
      if (events[i].ubiquity[sheet] == 0) continue;
      offset[i] = j; j++;
    }
    vdata[0] = j;
    if (D == 1) {
      for(i=0; i<nv; ++i) {
        if (events[i].ubiquity[sheet] == 0) continue;
        vcoords.push_back(float(x[i]));
        vcoords.push_back(0.0);
        vcoords.push_back(0.0);
      }
    }
    else if (D == 2) {
      for(i=0; i<nv; ++i) {
        if (events[i].ubiquity[sheet] == 0) continue;
        vcoords.push_back(float(x[2*i]));
        vcoords.push_back(float(x[2*i+1]));
        vcoords.push_back(0.0);
      }
    }
    else if (D >= 3) {
      for(i=0; i<nv; ++i) {
        if (events[i].ubiquity[sheet] == 0) continue;
        vcoords.push_back(float(x[D*i]));
        vcoords.push_back(float(x[D*i+1]));
        vcoords.push_back(float(x[D*i+2]));
      }
    }
    // Now the neighbour table...
    j = 0;
    for(i=0; i<ne; ++i) {
      if (simplices[1][i].ubiquity[sheet] == 0) continue;
      simplices[1][i].get_vertices(vx);
      evertex.push_back(offset[vx[0]]);
      evertex.push_back(offset[vx[1]]);
      ++j;
    }
    vdata[1] = j;
  }
}

void Spacetime::export_visual_data(std::vector<float>& colours,std::vector<float>& vcoords,std::vector<int>& evertex,int* vdata,bool use_energy) const
{
  int i,D,vx[2];
  std::vector<double> x;
  std::vector<unsigned char> chi;

  const int nv = (signed) events.size();
  const int ne = (signed) simplices[1].size();

  compute_colours(chi,false,use_energy);
  D = geometry->compute_coordinates(x);
  colours.clear();
  vcoords.clear();
  evertex.clear();
  vdata[0] = nv;
  vdata[1] = ne;
  for(i=0; i<3*nv; ++i) {
    colours.push_back(float(chi[i])/255.0);
  }
  for(i=0; i<3*ne; ++i) {
    colours.push_back(float(chi[3*nv+i])/255.0);
  }
  if (D == 1) {
    for(i=0; i<nv; ++i) {
      vcoords.push_back(float(x[Geometry::background_dimension*i]));
      vcoords.push_back(0.0);
      vcoords.push_back(0.0);
    }
  }
  else if (D == 2) {
    for(i=0; i<nv; ++i) {
      vcoords.push_back(float(x[Geometry::background_dimension*i]));
      vcoords.push_back(float(x[Geometry::background_dimension*i+1]));
      vcoords.push_back(0.0);
    }
  }
  else if (D >= 3) {
    for(i=0; i<nv; ++i) {
      vcoords.push_back(float(x[Geometry::background_dimension*i]));
      vcoords.push_back(float(x[Geometry::background_dimension*i+1]));
      vcoords.push_back(float(x[Geometry::background_dimension*i+2]));
    }
  }
  for(i=0; i<ne; ++i) {
    if (ghost(simplices[1][i].ubiquity)) {
      vx[0] = 0;
      vx[1] = 0;
    }
    else {
      simplices[1][i].get_vertices(vx);
    }
    evertex.push_back(vx[0]);
    evertex.push_back(vx[1]);
  }
}

