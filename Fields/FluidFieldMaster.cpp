/////////////////////////////////////////////////////////////
//                                                         //
// Copyright (c) 2003-2014 by The University of Queensland //
// Centre for Geoscience Computing                         //
// http://earth.uq.edu.au/centre-geoscience-computing      //
//                                                         //
// Primary Business: Brisbane, Queensland, Australia       //
// Licensed under the Open Software License version 3.0    //
// http://www.apache.org/licenses/LICENSE-2.0          //
//                                                         //
/////////////////////////////////////////////////////////////

// --- TML includes ---
#include "tml/comm/comm.h"

//--- project includes ---
#include "FieldMaster.h"
#include "Foundation/vec3.h"
#include "field_const.h"
#include "FluidFieldMaster.h"

// --- TML includes ---
#include "tml/comm/comm.h"

//--- IO includes ---
#include <iostream>
#include <fstream>
#include <algorithm>

//--- STL inculdes ---
#include <string>
#include <map>


using std::map;
using std::multimap;
using std::cout;
using std::endl;
using std::ofstream;
using std::string;

//==== SCALAR PFM ===

/*!
  Constructor. Set up the Master and broadcast parameters to the slaves.

  \param comm the communicator
  \param fieldname the name of the field to be saved
  \param filename the name of the file to be saved into or the base for the generation of the filenames if the saving format requires multiple files
  \param savetype the way to save data, currently supported are DX,SUM
  \param t0 the first timestep to be saved
  \param tend the last timestep to be saved
  \param dt the interval between timesteps to be saving
*/
ScalarFluidFieldMaster::ScalarFluidFieldMaster(TML_Comm* comm,const string& fieldname,const string& filename,const string& savetype,int t0,int tend,int dt)
  :AFieldMaster(comm,fieldname,filename,savetype,t0,tend,dt)
{
  m_comm->broadcast_cont(fieldname);
  m_comm->broadcast(m_id);
}


void ScalarFluidFieldMaster::collect()
{
  // send field id to slave
  m_comm->broadcast(m_id);

  if((m_write_type==WRITE_TYPE_SUM)||(m_write_type==WRITE_TYPE_MAX)) {
    collectSum();
  } else {
    collectFull();
  }
}

/*!
  collect the full set of data, i.e. value, and posn of fluid cells
*/
void ScalarFluidFieldMaster::collectFull()
{
  multimap<int,pair<Vec3,double> > temp_mm;

  // send type of collect to slaves
  int coll_type=1;
  m_comm->broadcast(coll_type);

  // get data from slaves
  m_comm->gather(temp_mm);

  // collate receive data from multimap into single map
  for(multimap<int,pair<Vec3,double> >::iterator iter=temp_mm.begin();
      iter!=temp_mm.end();
      iter++){
    m_save_vector.push_back(iter->second);
  }
}


/*!
  collect sum of values only
*/
void ScalarFluidFieldMaster::collectSum()
{

  multimap<int,double> temp_mm;

  // send type of collect to slaves
  m_comm->broadcast(m_write_type);

  // get data from slaves
  m_comm->gather(temp_mm);

  // collate receive data from multimap into single map
  for(multimap<int,double>::iterator iter=temp_mm.begin();
      iter!=temp_mm.end();
      iter++){
    m_sum_vector.push_back(iter->second);
  }
}


/*!
  sum data and write them out into a single continuous file

  \warning n
*/
void ScalarFluidFieldMaster::writeAsSUM()
{

  // sum data
  double sum_data=0.0;
  for(vector<double>::iterator iter=m_sum_vector.begin();
      iter!=m_sum_vector.end();
      iter++){
    sum_data+=*iter;
  }
  // open file for appending
  ofstream out_file(m_file_name.c_str(),ios::app);
  // write data
  out_file << sum_data << endl;
  // close file
  out_file.close();
  m_save_vector.erase(m_save_vector.begin(),m_save_vector.end());
  m_sum_vector.erase(m_sum_vector.begin(),m_sum_vector.end());
}

/*!
  find max datum and write it out into a single continuous file
*/
void ScalarFluidFieldMaster::writeAsMAX()
{

  // get max data
  double max_data=*m_sum_vector.begin();
  for(vector<double>::iterator iter=m_sum_vector.begin();
      iter!=m_sum_vector.end();
      iter++){
    max_data=(*iter > max_data) ? *iter : max_data;
  }
  // open file for appending
  ofstream out_file(m_file_name.c_str(),ios::app);
  // write data
  out_file << max_data << endl;
  // close file
  out_file.close();
  m_save_vector.erase(m_save_vector.begin(),m_save_vector.end());
  m_sum_vector.erase(m_sum_vector.begin(),m_sum_vector.end());
}

/*!
  write data as a raw series of values, one row of values per time step,
  all timesteps into the same file
*/
void ScalarFluidFieldMaster::writeAsRAW_SERIES()
{

  // open file
  ofstream out_file(m_file_name.c_str(),ios::app);

  // write data
  for(vector<pair<Vec3,double> >::iterator iter=m_save_vector.begin();
      iter!=m_save_vector.end();
      iter++){
    out_file << iter->first  << " " << iter->second << endl;
  }

  out_file << endl;
  out_file.close();
  m_save_vector.erase(m_save_vector.begin(),m_save_vector.end());
  m_sum_vector.erase(m_sum_vector.begin(),m_sum_vector.end());
}

/*!
  write data as raw position, value, one time step per file
*/

void ScalarFluidFieldMaster::writeAsRAW()
{

  //generate filename
  string fn=makeFilename();

  // open file
  ofstream out_file(fn.c_str());

  // write data

  for(vector<pair<Vec3,double> >::iterator iter=m_save_vector.begin();
      iter!=m_save_vector.end();
      iter++){
    out_file << iter->first  << " " << iter->second << endl;
  }

  // close file
  out_file.close();

  // clean up
  m_save_vector.erase(m_save_vector.begin(),m_save_vector.end());
  m_sum_vector.erase(m_sum_vector.begin(),m_sum_vector.end());
}


bool sortOnZ(const pair<Vec3,double> a, const pair<Vec3,double> b) {
  Vec3 a_index=a.first;
  Vec3 b_index=b.first;
  return a_index.Z() < b_index.Z();
}

bool sortOnY(const pair<Vec3,double> a, const pair<Vec3,double> b) {
  Vec3 a_index=a.first;
  Vec3 b_index=b.first;
  return a_index.Y() < b_index.Y();
}

bool sortOnX(const pair<Vec3,double> a, const pair<Vec3,double> b) {
  Vec3 a_index=a.first;
  Vec3 b_index=b.first;
  return a_index.X() < b_index.X();
}

vector<pair<Vec3,double> > sortVector(vector<pair<Vec3,double> > unsorted) 
{
  vector<pair<Vec3,double> > sorted = unsorted;
  std::stable_sort(sorted.begin(), sorted.end(), sortOnX);
  std::stable_sort(sorted.begin(), sorted.end(), sortOnY);
  std::stable_sort(sorted.begin(), sorted.end(), sortOnZ);
  return sorted;
}


/*!
  write data as VTI file for paraview, one time step per file
*/

void ScalarFluidFieldMaster::writeAsVTI()
{

  //generate filename
  string fn=makeFilename();

  // open file
  ofstream out_file(fn.c_str());

  vector<pair<Vec3,double> > sorted;
  sorted=sortVector(m_save_vector);
  double minX,minY,minZ,maxX,maxY,maxZ;
  vector<pair<Vec3,double> >::iterator iter_begin=sorted.begin();
  Vec3 pos_begin=iter_begin->first;
  minX=pos_begin.X();minY=pos_begin.Y();minZ=pos_begin.Z();
  vector<pair<Vec3,double> >::iterator iter_end=sorted.end()-1;
  Vec3 pos_end=iter_end->first;
  maxX=pos_end.X();maxY=pos_end.Y();maxZ=pos_end.Z();

  int z_count=0;
  for(vector<pair<Vec3,double> >::iterator iter0=sorted.begin();iter0!=sorted.end();iter0++){
    Vec3 pos=iter0->first;
    if(pos.Z()==minZ) {z_count++;} 
    else {break;};
  }
  int nz=sorted.size()/z_count; 

  int y_count=0;
  for(vector<pair<Vec3,double> >::iterator iter0=sorted.begin();iter0!=sorted.end();iter0++){
    Vec3 pos=iter0->first;
    if(pos.Y()==minY) {y_count++;} 
    else {break;};
  }
  int ny=z_count/y_count;

  int x_count=0;
  for(vector<pair<Vec3,double> >::iterator iter0=sorted.begin();iter0!=sorted.end();iter0++){
    Vec3 pos=iter0->first;
    if(pos.X()==minX) {x_count++;} 
    else {break;};
  }
  int nx=y_count/x_count;

  double xside,yside,zside;
  if(nx==1){
    xside=maxX-minX;
  }else{
    xside=(maxX-minX)/double(nx-1);
  }

  if(ny==1){
    yside=maxY-minY;
  }else{
    yside=(maxY-minY)/double(ny-1);
  }
  if(nz==1){
    zside=maxZ-minZ;
  }else{
    zside=(maxZ-minZ)/double(nz-1);
  }

  out_file<<"<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">"<< endl;
  out_file<<"<ImageData WholeExtent=\"0 "<<nx<<" 0 "<<ny<<" 0 " <<nz<<"\" Origin=\""<<minX-0.5*xside<<" "<<minY-0.5*yside<<" "<<minZ-0.5*zside<<"\" Spacing=\""\
          <<xside<<" "<<yside<<" "<<zside<<"\">" << endl;
  out_file << "<Piece Extent=\"0 "<<nx<<" 0 "<<ny<<" 0 "<<nz<<"\">" << endl;
  out_file << "<CellData>" << endl;
  out_file << "<DataArray type=\"Float32\" Name=\"fluid scalar data\" format=\"ascii\">" << endl;

  vector<pair<Vec3,double> >::iterator iter=sorted.begin();
  for (int z=0;z<nz;z++){
    for (int y=0;y<ny;y++){
      for (int x=0;x<nx;x++){
        out_file << iter->second << " ";
        iter++;
      };
      out_file<<endl;
    };
    out_file<<endl;  
  };  

  out_file << "</DataArray>" << endl <<endl;
  out_file << "</CellData>" << endl;
  out_file << "</Piece>" << endl;
  out_file << "</ImageData>" << endl;
  out_file << "</VTKFile>" << endl;

  // close file
  out_file.close();

  // clean up
  m_save_vector.erase(m_save_vector.begin(),m_save_vector.end());
  m_sum_vector.erase(m_sum_vector.begin(),m_sum_vector.end());
}


// === VECTOR PFM ===
/*!
  Constructor. Set up the Master and broadcast parameters to the slaves.

  \param comm the communicator
  \param fieldname the name of the field to be saved
  \param filename the name of the file to be saved into or the base for the generation of the filenames if the saving format requires multiple files
  \param savetype the way to save data, currently supported are DX,SUM
  \param t0 the first timestep to be saved
  \param tend the last timestep to be saved
  \param dt the interval between timesteps to be saving
*/
VectorFluidFieldMaster::VectorFluidFieldMaster(TML_Comm* comm,const string& fieldname,const string& filename,const string& savetype,int t0,int tend,int dt)
  :AFieldMaster(comm,fieldname,filename,savetype,t0,tend,dt)
{
  m_comm->broadcast_cont(fieldname);
  m_comm->broadcast(m_id);
}

void VectorFluidFieldMaster::collect()
{
  multimap<int,pair<Vec3,Vec3> > temp_mm;

  // send field id to slaves
  m_comm->broadcast(m_id);

  //receive data from fields
  m_comm->gather(temp_mm);

  // collate receive data from multimap into single map
  for(multimap<int,pair<Vec3,Vec3> >::iterator iter=temp_mm.begin();
      iter!=temp_mm.end();
      iter++){
    m_save_vector.push_back(iter->second);
  }
}


/*!
  sum data and write them out into a single continuous file

  \warning untested
*/
void VectorFluidFieldMaster::writeAsSUM()
{
  Vec3 sum_data=Vec3(0.0,0.0,0.0);
  for(vector<pair<Vec3,Vec3> >::iterator iter=m_save_vector.begin();
      iter!=m_save_vector.end();
      iter++){
    sum_data+=iter->second;
  }
  // open file for appending
  ofstream out_file(m_file_name.c_str(),ios::app);
  // write data
  out_file << sum_data << endl;
  // close file
  out_file.close();
  m_save_vector.erase(m_save_vector.begin(),m_save_vector.end());
}

/*!
  get max datum and write it out into a single continuous file

  \warning untested
*/
void VectorFluidFieldMaster::writeAsMAX()
{
  Vec3 max_data=(m_save_vector.begin())->second;
  for(vector<pair<Vec3,Vec3> >::iterator iter=m_save_vector.begin();
      iter!=m_save_vector.end();
      iter++){
    max_data=(iter->second.norm2() > max_data.norm2()) ? iter->second : max_data;
  }
  // open file for appending
  ofstream out_file(m_file_name.c_str(),ios::app);
  // write data
  out_file << max_data << endl;
  // close file
  out_file.close();
  m_save_vector.erase(m_save_vector.begin(),m_save_vector.end());
}


/*!
  write data as a raw series of values, one row of values per time step,
  all timesteps into the same file
*/
void VectorFluidFieldMaster::writeAsRAW_SERIES()
{
  // open file
  ofstream out_file(m_file_name.c_str(),ios::app);

  // write data
  for(vector<pair<Vec3,Vec3> >::iterator iter=m_save_vector.begin();
      iter!=m_save_vector.end();
      iter++){
    out_file << iter->first  << " " << iter->second << endl;
  }
  out_file << endl;
  out_file.close();
  m_save_vector.erase(m_save_vector.begin(),m_save_vector.end());
}

/*!
  write data as raw position, value pairs, one time step per file
*/
void VectorFluidFieldMaster::writeAsRAW()
{
  //generate filename
  string fn=makeFilename();

  // open file
  ofstream out_file(fn.c_str());

  for(vector<pair<Vec3,Vec3> >::iterator iter=m_save_vector.begin();
      iter!=m_save_vector.end();
      iter++){
    out_file << iter->first  << " " << iter->second << endl;
  }

  // close file
  out_file.close();

  // clean up
  m_save_vector.erase(m_save_vector.begin(),m_save_vector.end());
}

/*!
  write data as VTU file for paraview, position, value, one time step per file
*/
void VectorFluidFieldMaster::writeAsVTU()
{
  //generate filename
  string fn=makeFilename();

  // open file
  ofstream out_file(fn.c_str());

  int n=0;
  for(vector<pair<Vec3,Vec3> >::iterator iter=m_save_vector.begin();
      iter!=m_save_vector.end();
      iter++){
    n++;
  }
  double xside,yside,zside;
  Vec3 number0;
  vector<pair<Vec3,Vec3> >::iterator iter0=m_save_vector.begin();
  Vec3 pos0=iter0->first;

  for(iter0=m_save_vector.begin();iter0!=m_save_vector.end();iter0++){
    number0=iter0->first;
    if(number0.X()!=pos0.X()){xside=fabs(number0.X()-pos0.X());break;};
  }
  for(iter0=m_save_vector.begin();iter0!=m_save_vector.end();iter0++){
    number0=iter0->first;
    if(number0.Y()!=pos0.Y()){yside=fabs(number0.Y()-pos0.Y());break;};
  }
  for(iter0=m_save_vector.begin();iter0!=m_save_vector.end();iter0++){
    number0=iter0->first;
    if(number0.Z()!=pos0.Z()){zside=fabs(number0.Z()-pos0.Z());break;};
  }

  double cellsize=(xside+yside+zside)/3.0;
  out_file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">" << endl;
  out_file << "<UnstructuredGrid>" << endl;
  out_file << "<Piece NumberOfPoints=\"" << n << "\" NumberOfCells=\"" << 0 << "\">" << endl;

  //write positions
  out_file << "<Points>" << endl;
  out_file << "<DataArray NumberOfComponents=\"3\" type=\"Float64\" format=\"ascii\">" << endl;
  for(vector<pair<Vec3,Vec3> >::iterator iter=m_save_vector.begin();
      iter!=m_save_vector.end();
      iter++){
    out_file << iter->first<< endl;
  }
  out_file << "</DataArray>" << endl;
  out_file << "</Points>" << endl;

  //write cell size
  out_file << "<PointData Scalars=\"cellsize\">" << endl;
  out_file << "<DataArray type=\"Float64\" Name=\"cellsize\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
  for(vector<pair<Vec3,Vec3> >::iterator iter=m_save_vector.begin();
      iter!=m_save_vector.end();
      iter++){
    out_file << cellsize << endl;
  }
  out_file << "</DataArray>" << endl;

  //write vector value
  out_file << "<DataArray type=\"Float64\" Name=\"vector value\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;
  for(vector<pair<Vec3,Vec3> >::iterator iter=m_save_vector.begin();
      iter!=m_save_vector.end();
      iter++){
    out_file << iter->second << endl;
  }

  out_file << "</DataArray>" << endl;
  out_file << "</PointData>" << endl;
  out_file << "<Cells>" << endl;
  out_file << "<DataArray type=\"Int32\" NumberOfComponents=\"1\" Name=\"connectivity\" format=\"ascii\">" << endl;
  out_file << "</DataArray>" << endl;
  out_file << "<DataArray type=\"Int32\" NumberOfComponents=\"1\" Name=\"offsets\" format=\"ascii\">" << endl;
  out_file << "</DataArray>" << endl;
  out_file << "<DataArray type=\"UInt8\" NumberOfComponents=\"1\" Name=\"types\" format=\"ascii\">" << endl;
  out_file << "</DataArray>" << endl;
  out_file << "</Cells>" << endl;
  out_file << "</Piece>" << endl;
  out_file << "</UnstructuredGrid>" << endl;
  out_file << "</VTKFile>" << endl;

  // close file
  out_file.close();

  // clean up
  m_save_vector.erase(m_save_vector.begin(),m_save_vector.end());
}


