/////////////////////////////////////////////////////////////
//                                                         //
// Copyright (c) 2003-2017 by The University of Queensland //
// Centre for Geoscience Computing                         //
// http://earth.uq.edu.au/centre-geoscience-computing      //
//                                                         //
// Primary Business: Brisbane, Queensland, Australia       //
// Licensed under the Open Software License version 3.0    //
// http://www.apache.org/licenses/LICENSE-2.0              //
//                                                         //
/////////////////////////////////////////////////////////////

#include "pfcc.h"
#include <cstdlib>

// --- STL includes ---
#include <fstream>
#include <vector>
#include <map>

using std::ifstream;
using std::map;
using std::vector;

using std::make_pair;
/*!
  convert pfc (vvf) file to geo

  \param infilename the input file
  \param outfilename the output file
  \param bbx_min the minimum corner of the bounding box
  \param bbx_max the maximum corner of the bounding box
  \param cx circular boundary cond. in x-dir
  \param cy circular boundary cond. in y-dir
  \param cz circular boundary cond. in z-dir
  \param scale scaling factor
*/
void pfc_convert(const string& infilename,const string& outfilename,const Vec3& bbx_min,const Vec3& bbx_max,int cx,int cy,int cz,double scale)
{
  string ins;
  string::size_type idx;
  map<string,int> colmap;

  ifstream infile(infilename.c_str()); 

  // read 1st line -> dummy,nr.of particles
  infile >> ins;
  idx=ins.find(",");
  string s2=ins.substr(idx+1);
  int npart=atoi(s2.c_str());
  cout << "npart= " << npart << endl;
  
  // read 2nd line (formatting) -> throw away
  infile >> ins;

  //read 3rd line (colum description)
  infile >> ins;
  // split into strings 
  int count=0;
  string::size_type idx_old=0;
  do {
    idx=ins.find(",",idx_old);
    s2=ins.substr(idx_old,idx-idx_old);
    colmap.insert(make_pair(s2,count));
    count++;
    idx_old=idx+1;
  } while (idx!=string::npos);
  // exctact locations of the interesting columns
  int id_col=colmap["EVENT_LABEL"];
  int posx_col=colmap["PFC_BALL_LOCX"];
  int posy_col=colmap["PFC_BALL_LOCY"];
  int posz_col=colmap["PFC_BALL_LOCZ"];
  int rad_col=colmap["PFC_BALL_RADIUS"];
  int tag_col=colmap["PFC_BALL_ATTR"];

  // write geo header
  ofstream outfile(outfilename.c_str());
  outfile << "LSMGeometry 1.2" << endl;
  outfile << "BoundingBox " << bbx_min << " " << bbx_max << endl;
  outfile << "PeriodicBoundaries " << cx << " " << cy << " " << cz << endl;
  outfile << "Dimensions 3D" << endl;
  outfile << "BeginParticles" << endl;
  outfile << "Simple" << endl;
  outfile << npart << endl;
  vector<string> temp(count);
  for(int ipart=0;ipart<npart;ipart++){
    // read pfc Line
    for(int j=0;j<count;j++){
      infile >> temp[j];
    }
    int pid=atoi((temp[id_col]).c_str());
    Vec3 pos=Vec3(atof((temp[posx_col]).c_str()),atof((temp[posy_col]).c_str()),atof((temp[posz_col]).c_str()));
    double rad=atof((temp[rad_col]).c_str());
    int tag=atoi((temp[tag_col]).c_str());
    // write geo line
    outfile << pos << " " << rad << " " << pid << " " << tag << endl;
  }
  // write geo footer
  outfile << "EndParticles" << endl;
  outfile << "BeginConnect" << endl;
  outfile << "0" << endl;
  outfile << "EndConnect" << endl;
  // close files
  outfile.close();
  infile.close();
}

/*!
  read bounding box file
*/
pair<Vec3,Vec3> read_bbx(const string& filename)
{
  Vec3 min,max;
  double xmin=0.0,ymin=0.0,zmin=0.0;
  double xmax=0.0,ymax=0.0,zmax=0.0;
  double tmp=0.0;
  string id,dummy;

  ifstream infile(filename.c_str());

  for(int i=0;i<6;i++){
    infile >> id >> dummy >> tmp;
    if(id=="wx1"){
      xmax=tmp;
    } else if(id=="wx2"){
      xmin=tmp;
    } else if(id=="wy1"){
      ymax=tmp;
    } else if(id=="wy2"){
      ymin=tmp;
    } else if(id=="wz1"){
      zmax=tmp;
    } else if(id=="wz2"){
      zmin=tmp;
    }
  }
  infile.close();

  if (xmin>xmax) {
    tmp=xmin;xmin=xmax;xmax=tmp;
  }
  if (ymin>ymax) {
    tmp=ymin;ymin=ymax;ymax=tmp;
  }
  if (zmin>zmax) {
    tmp=zmin;zmin=zmax;zmax=tmp;
  }

  min=Vec3(xmin,ymin,zmin);
  max=Vec3(xmax,ymax,zmax);

  return make_pair(min,max);
}
