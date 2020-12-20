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

#include "fracwriter.h"

// --- IO includes ---
#include <fstream>
#include <iostream>
#include <sstream>

using std::ostringstream;
using std::ofstream;
using std::endl;

fwdata::fwdata(const FracFrame::fdata& fd,int t)
{
  pos=fd.pos;
  normal=fd.normal;
  size=fd.size;
  id1=fd.id1;
  id2=fd.id2;
  time=t;
  tag=fd.tag;
  ptag1=fd.ptag1;
  ptag2=fd.ptag2;
  dist=fd.dist;
}

/*!
  constructor
*/
FracWriter::FracWriter()
{
  with_plane=false;
}

/*!
  add point data to the writer
*/
void FracWriter::addData(const vector<FracFrame::fdata>& newdata,int t)
{
  for(vector<FracFrame::fdata>::const_iterator iter=newdata.begin();
      iter!=newdata.end();
      iter++){
    m_data.push_back(fwdata(*iter,t));
    m_nbrk_map[iter->id1]++;
    m_nbrk_map[iter->id2]++;
  }
  std::cerr << "added " << newdata.size() << " fractures" << std::endl;
}

/*!
  add a plane to the data writer
*/
void FracWriter::addPlane(const Plane3D& P)
{
  with_plane=true;
  Vec3 center=P.getPos();
  double maxdist=0.0;
  for(vector<fwdata>::iterator iter=m_data.begin();
      iter!=m_data.end();
      iter++){
    double dist=(center-(iter->pos)).norm();
    maxdist=(dist > maxdist) ? dist : maxdist;
  }
  Vec3 v1=maxdist*P.GetV();
  Vec3 v2=maxdist*P.GetU();
  m_c1=center-v1-v2;
  m_c2=center+v1-v2;
  m_c3=center+v1+v2;
  m_c4=center-v1+v2;
}
 
/*!
  write data int VTK-XML file

  \param filename the filename
*/
void FracWriter::write(const string& filename)
{
  // open file
  ofstream outfile(filename.c_str());

  // write VTK-XML header
  outfile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n";
  outfile << "<UnstructuredGrid>\n";
  outfile << "<Piece NumberOfPoints=\"" << m_data.size()<< "\" NumberOfCells=\"0\">\n";

  // write positions
  outfile << "<Points>\n";
  outfile << "<DataArray NumberOfComponents=\"3\" type=\"Float64\" format=\"ascii\">\n";
  for(vector<fwdata>::iterator iter=m_data.begin();
      iter!=m_data.end();
      iter++){
    outfile << iter->pos << endl;
  }
  outfile << "</DataArray>\n";
  outfile << "</Points>\n";
  
  // write times
  outfile << "<PointData Scalars=\"time\">\n";
  outfile << "<DataArray type=\"Int32\" Name=\"time\" NumberOfComponents=\"1\" format=\"ascii\">\n";
  for(vector<fwdata>::iterator iter=m_data.begin();
      iter!=m_data.end();
      iter++){
    outfile << iter->time << endl;
  }
  outfile << "</DataArray>\n";
  // write radii
  outfile << "<DataArray type=\"Float32\" Name=\"radius\" NumberOfComponents=\"1\" format=\"ascii\">\n";
  for(vector<fwdata>::iterator iter=m_data.begin();
      iter!=m_data.end();
      iter++){
    outfile << iter->size << endl;
  }
  outfile << "</DataArray>\n";
 // write particle distances
  outfile << "<DataArray type=\"Float32\" Name=\"distance\" NumberOfComponents=\"1\" format=\"ascii\">\n";
  for(vector<fwdata>::iterator iter=m_data.begin();
      iter!=m_data.end();
      iter++){
    outfile << iter->dist << endl;
  }
  outfile << "</DataArray>\n";
  // write normals
  outfile << "<DataArray type=\"Float32\" Name=\"normal\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for(vector<fwdata>::iterator iter=m_data.begin();
      iter!=m_data.end();
      iter++){
    outfile << iter->normal.X() << " " << iter->normal.Y() << " " << iter->normal.Z()  << endl;
  }
  outfile << "</DataArray>\n";
  
  // write (ex-)bond tags
  outfile << "<DataArray type=\"Int32\" Name=\"tag\" NumberOfComponents=\"1\" format=\"ascii\">\n";
  for(vector<fwdata>::iterator iter=m_data.begin();
      iter!=m_data.end();
      iter++){
    outfile << iter->tag  << endl;
  }
  outfile << "</DataArray>\n";
 
   // write particle tags of (ex-) bonds
  outfile << "<DataArray type=\"Int32\" Name=\"particle tag 1\" NumberOfComponents=\"1\" format=\"ascii\">\n";
  for(vector<fwdata>::iterator iter=m_data.begin();
      iter!=m_data.end();
      iter++){
    outfile << iter->ptag1  << endl;
  }
  outfile << "</DataArray>\n";
 
  outfile << "<DataArray type=\"Int32\" Name=\"particle tag 2\" NumberOfComponents=\"1\" format=\"ascii\">\n";
  for(vector<fwdata>::iterator iter=m_data.begin();
      iter!=m_data.end();
      iter++){
    outfile << iter->ptag2  << endl;
  }
  outfile << "</DataArray>\n";
  outfile << "</PointData>\n";
  
  // write empty cell block
  outfile << "<Cells>\n";
  outfile << "<DataArray type=\"Int32\" NumberOfComponents=\"1\" Name=\"connectivity\" format=\"ascii\">\n";
  outfile << "</DataArray>";
  outfile << "<DataArray type=\"Int32\" NumberOfComponents=\"1\" Name=\"offsets\" format=\"ascii\">\n";
  outfile << "</DataArray>\n";
  outfile << "<DataArray type=\"UInt8\" NumberOfComponents=\"1\" Name=\"types\" format=\"ascii\">\n";
  outfile << "</DataArray>\n";  
  outfile << "</Cells>\n";

  // write footer
  outfile << "</Piece>\n";
  outfile << "</UnstructuredGrid>\n";
  outfile << "</VTKFile>\n";


  // close file
  outfile.close();

  // if we have a plane, setup filename and write it
  if(with_plane){
    ostringstream planefilename;

    planefilename << "plane_" << filename;
    writePlane(planefilename.str());
  }
}


/*!
  write data as text file

  \param filename the filename
*/
void FracWriter::writeText(const string& filename)
{
  // open file
  ofstream outfile(filename.c_str());

  // write positions & times
  for(vector<fwdata>::iterator iter=m_data.begin();
      iter!=m_data.end();
      iter++){
    outfile << iter->pos << " " << iter->size << " " << iter->time << " " << iter->normal << " " << iter->id1 << " " << iter->id2 << " " << iter->tag << " " << iter->ptag1 << " " << iter->ptag2 << endl;
  }

  // close file
  outfile.close();
}

/*!
  write data as text file

  \param filename the filename
*/
void FracWriter::writeParticleList(const string& filename)
{
  // open file
  ofstream outfile(filename.c_str());

  // write positions & times
  for(map<int,int>::iterator iter=m_nbrk_map.begin();
      iter!=m_nbrk_map.end();
      iter++){
    outfile << iter->first << " " << iter->second << endl;
  }

  // close file
  outfile.close();
}

/*!
  write fracture distribution profile (y-direction)

  \param ymin minimum y-coordinate
  \param ymax maximum y-coordinate
  \nbin number of bins in the distribution
  \param filename the filename
*/								
void FracWriter::writeProfile(double ymin,double ymax,int nbin,const string& filename)
{
  // open file
  ofstream outfile(filename.c_str());
  
  vector<double> prof=vector<double>(nbin,0.0);

  // calculate profile
  for(vector<fwdata>::iterator iter=m_data.begin();
      iter!=m_data.end();
      iter++){
    int bnr=int(floor(nbin*((iter->pos.Y()-ymin)/(ymax-ymin))));
    if((bnr>=0) && (bnr<nbin)){
      prof[bnr]+=(iter->size)*(iter->size);
    }
  }

  // write profile
  for(int i=0;i<nbin;i++){
    outfile << ymin+(double(i)+0.5)*((ymax-ymin)/double(nbin)) << " " << prof[i] << std::endl;
  }
  // close file
  outfile.close();
}

/*!
  write plane data int VTK-XML file
*/
void FracWriter::writePlane(const string& filename)
{
  // open file
  ofstream outfile(filename.c_str());

  // write VTK-XML header
  outfile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n";
  outfile << "<UnstructuredGrid>\n";
  outfile << "<Piece NumberOfPoints=\"4\" NumberOfCells=\"1\">\n";

  // write positions
  outfile << "<Points>\n";
  outfile << "<DataArray NumberOfComponents=\"3\" type=\"Float64\" format=\"ascii\">\n";
  outfile << m_c1 << endl;
  outfile << m_c2 << endl;
  outfile << m_c3 << endl;
  outfile << m_c4 << endl;
  outfile << "</DataArray>\n";
  outfile << "</Points>\n";
  
  // write cells
  outfile << "<Cells>\n";
  outfile << "<DataArray type=\"Int32\" Name=\"connectivity\">\n";
  outfile << "0 1 2 3\n";
  outfile << "</DataArray>\n";
  outfile << "<DataArray type=\"Int32\" Name=\"offsets\">\n";
  outfile << "4\n";
  outfile << "</DataArray>\n";
  outfile << "<DataArray type=\"UInt8\" Name=\"types\">\n";
  outfile << "9\n";
  outfile << "</DataArray>\n";
  outfile << "</Cells>\n";

  // write dummy cell data
  outfile << "<CellData>\n";
  outfile << "<DataArray NumberOfComponents=\"1\" type=\"Float32\" Name=\"dummy\">\n";
  outfile << "1.0\n";
  outfile << "</DataArray>\n";
  outfile << "</CellData>\n";
  outfile << "</Piece>\n";
  outfile << "</UnstructuredGrid>\n";
  outfile << "</VTKFile>\n";
  // close file
  outfile.close();
 
}
