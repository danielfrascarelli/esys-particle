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

#include "DataExtractor.h"
#include "SnapFileHelp.h"

// --- Project includes ---
#include "ntable/src/handle.h"

// --- std. includes ---
#include <cmath>

using std::atan;
using std::sqrt;

// --- IO includes ---
#include <fstream>

using std::ifstream;
using std::ofstream;

/*!
  constructor

  \param x nr. of grid points in x-direction
  \param y nr. of grid points in y-direction
  \param z nr. of grid points in z-direction
  \param range grid spacing
  \param p0_global minimal corner (origin) of the global search space
*/
DataExtractor::DataExtractor(int x,int y,int z,double range,const Vec3& p0_global,const Vec3& pmax_global)://modified (fluid contents)
  m_data(x,y,z,range,0.0,p0_global,pmax_global,0,0,0) //modified (fluid contents)
{
  
}

/*!
  read snapshot

  \param filename the name of the "metadata" file, usually *_0.txt
*/
void DataExtractor::read(const string& infilename)
{
  ifstream headerfile(infilename.c_str());
  int version=get_version(infilename);
  vector<string> filenames=get_filenames(infilename,version);

  // get main files
  for(vector<string>::iterator iter=filenames.begin();
      iter!=filenames.end();
      iter++){
    cout << *iter << endl;
    ifstream datafile(iter->c_str());

    // get particles
    int npart;
    Vec3 pos;
    Vec3 oldpos;
    Vec3 initpos;
    Vec3 force;
    Vec3 vel;
    Vec3 angvel;
    Vec3 circ_shift;
    double rad;
    double mass;
    double q1,q2,q3,q4;
    int id;
    int tag;
    datafile >> npart;
    if(version<2){
      for(int i=0;i<npart;i++){
	// read data
	datafile >> pos  >> rad >> id >> tag >> mass >> initpos >> oldpos >> vel >> force >> q1 >> q2 >> q3 >> q4 >> angvel;
	DataParticle dp=DataParticle(pos,initpos,rad,id);
	m_data.insert(dp);
      }
    } else {
      for(int i=0;i<npart;i++){
	// read data
	datafile >> pos  >> rad >> id >> tag >> mass >> initpos >> oldpos >> vel >> force >> circ_shift >> q1 >> q2 >> q3 >> q4 >> angvel;
	DataParticle dp=DataParticle(pos,initpos,rad,id);
	m_data.insert(dp);	
      }
    }
    datafile.close();
  }
  std::cout << "inserted " << m_data.size() << "  particles " << std::endl;
  // m_data.build();
}

/*!
  write tensor data as unstructured grid VTK-XML file

  \param filename the name of the output file
  \param dataname the name of the data in the VTK file
*/
void DataExtractor::writeTensorDataVtk(const string& filename,const string& dataname)
{

}

/*!
  write scalar data as unstructured grid VTK-XML file

  \param filename the name of the output file
  \param dataname the name of the data in the VTK file
*/
void DataExtractor::writeScalarDataVtk(const string& filename,const string& dataname)
{
  ofstream vtkfile(filename.c_str());
  // write the file
  // write header 
  vtkfile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n";
  vtkfile << "<UnstructuredGrid>\n";
  vtkfile << "<Piece NumberOfPoints=\"" << m_data.size() << "\" NumberOfCells=\"0\">\n";
  // write particle pos
  vtkfile << "<Points>\n";
  vtkfile << "<DataArray NumberOfComponents=\"3\" type=\"Float64\" format=\"ascii\">\n";
  for(NeighborTable<DataParticle>::iterator iter=m_data.begin();
      iter!=m_data.end();
      iter++){
    vtkfile << iter->getPos() << "  ";
  }
  vtkfile << "</DataArray>\n";
  vtkfile << "</Points>\n";
  
  // --- write particle data ---
  // radius
  vtkfile << "<PointData Scalars=\"radius\">\n";
  vtkfile << "<DataArray type=\"Float64\" Name=\"radius\" NumberOfComponents=\"1\" format=\"ascii\">\n";
  for(NeighborTable<DataParticle>::iterator iter=m_data.begin();
      iter!=m_data.end();
      iter++){
    vtkfile << iter->getRad() << "  ";
  }
  vtkfile << "</DataArray>\n";
  vtkfile << "</PointData>\n";
  
  // data
  vtkfile << "<PointData Scalars=\"" << dataname << "\">\n";
  vtkfile << "<DataArray type=\"Float64\" Name=\"radius\" NumberOfComponents=\"1\" format=\"ascii\">\n";
  for(NeighborTable<DataParticle>::iterator iter=m_data.begin();
      iter!=m_data.end();
      iter++){
    vtkfile << iter->getScalarData() << "  ";
  }
  vtkfile << "</DataArray>\n";
  vtkfile << "</PointData>\n";

  // write empty cell block
  vtkfile << "<Cells>\n";
  vtkfile << "<DataArray type=\"Int32\" NumberOfComponents=\"1\" Name=\"connectivity\" format=\"ascii\">\n";
  vtkfile << "</DataArray>\n";
  vtkfile << "<DataArray type=\"UInt8\" NumberOfComponents=\"1\" Name=\"types\" format=\"ascii\">\n";
  vtkfile << "</DataArray>\n";
  vtkfile << "</Cells>\n";
  
  // write footer
  vtkfile << "</Piece>\n";
  vtkfile << "</UnstructuredGrid>\n";
  vtkfile << "</VTKFile>\n";
  
  // close file
  vtkfile.close();
}



/*!
  extract best fit strain tensors from the data and write the result
  into the tensor data member of the particles

  \param rad the search radius for the neighbour particles
*/
void DataExtractor::StrainToTensorData(double rad)
{
  // for each particle
  for(NeighborTable<DataParticle>::iterator iter=m_data.begin();
      iter!=m_data.end();
      iter++){
    // get list of neighbours
    T_Handle<NeighborTable<DataParticle>::particlelist> plist=m_data.getParticlesNearPoint(iter->getPos());
    std::cout << "pos: " << iter->getPos() << "  list size : " << plist->size() << std::endl; 
    if(plist->size()>3){
      // init sums
      Matrix3 M;
      Vec3 Su,Sv,Sw;
      // for each neighbour
      for(NeighborTable<DataParticle>::particlelist::iterator n_iter=plist->begin();
	  n_iter!=plist->end();
	  n_iter++){
	// get relative position & displacement
	Vec3 disp=(*n_iter)->getDisplacement()-iter->getDisplacement();
	Vec3 rpos=(*n_iter)->getPos()-iter->getPos();
	// --- udate sums
	// M
	M(0,0)+=rpos.X()*rpos.X();
	M(0,1)+=rpos.X()*rpos.Y();
	M(0,2)+=rpos.X()*rpos.Z();
	M(1,1)+=rpos.Y()*rpos.Y();
	M(1,2)+=rpos.Y()*rpos.Z();
	M(2,2)+=rpos.Z()*rpos.Z();
	// Su,v,w
	Su+=disp.X()*rpos;
	Sv+=disp.Y()*rpos;
	Sw+=disp.Z()*rpos;
      }
      // make M symmetric
      M(1,0)=M(0,1);
      M(2,0)=M(0,2);
      M(2,1)=M(1,2);
      // solve equations to get displacement gradient tensor components
      Vec3 Du=M.solve(Su);
      Vec3 Dv=M.solve(Sv);
      Vec3 Dw=M.solve(Sw);
      // make deformation gradient tensor
      Matrix3 A;
      A(0,0)=Du.X()+1; A(0,1)=Du.Y();     A(0,2)=Du.Z();
      A(1,0)=Dv.X();      A(1,1)=Dv.Y()+1; A(1,2)=Dv.Z();
      A(2,0)=Dw.X();     A(2,1)=Dw.Y();     A(2,2)=Dw.Z()+1;
      // right Cauchy-Green deformation tensor
      Matrix3 C=A.trans()*A;
      // write to tensor data member  of the particle
      iter->setTensorData(C);
    }
  }
}

/*!
  write maximum shear strain to scalar data member. Needs strain tensor calculation done before
*/
void DataExtractor::MaxShearToScalarData()
{
    for(NeighborTable<DataParticle>::iterator iter=m_data.begin();
      iter!=m_data.end();
      iter++){
      // get eigenvalues of tensor data
      Vec3 D1,D2,D3;
      double e1,e2,e3;
      (iter->getTensorData()).eigen(D1,D2,D3,e1,e2,e3);
      // set scalar data to max shear strain
      iter->setScalarData(atan(sqrt(e3/e1)));
    }
}
