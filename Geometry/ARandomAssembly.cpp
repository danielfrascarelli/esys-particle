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

#include "ARandomAssembly.h"

//-- system includes --
#include <cstdlib>
using std::rand;

//-- STL includes --
#include <map>
#include <utility>

using std::map;
using std::pair;

// --- IO includes ---
#include <fstream>

using std::ofstream;

double ARandomAssembly::m_small_value=1e-7;

/*!
  helper function, return random value between min & max
*/
double ARandomAssembly::m_random(double imin,double imax)
{
  return imin+((imax-imin)*((double)(rand())/(double)(RAND_MAX)));
} 

/*!
  get the list of neighbors of a particle

  \param Po the Particle
*/
vector<SimpleParticle> ARandomAssembly::getNeighborList(const SimpleParticle& Po)
{
  const vector<SimpleParticle> *NLP;

  NLP=m_snt->getNeighbors(Po.getPos());

  return *NLP;
}

/*!
  Get the n closest neighbors of a particle (sorted)

  \param Po the particle
  \param n max nr. or neighbours returned
*/
vector<SimpleParticle> ARandomAssembly::getClosestNeighbors(const SimpleParticle& Po, int n)
{
  vector<SimpleParticle> CL;
  map<double,const SimpleParticle*> pmap;
  const vector<SimpleParticle> *NLP;
  double max_dist = 0.0;

  // get neighbour list pointer
  NLP=m_snt->getNeighbors(Po.getPos());

  for(vector<SimpleParticle>::const_iterator iter=NLP->begin();
      iter!=NLP->end();
      iter++){
    double dist=(Po.getPos()-iter->getPos()).norm()-iter->getRad();
    if(pmap.size()<4){ // less than 4 in pvec -> insert
      pmap.insert(make_pair(dist,&(*iter)));
      max_dist=(pmap.rbegin())->first;
    } else if(dist<max_dist){ // closer than 4th -> insert
      pmap.erase(max_dist);
      pmap.insert(make_pair(dist,&(*iter)));
      max_dist=(pmap.rbegin())->first;
    }
  }
  for(map<double,const SimpleParticle*>::iterator iter=pmap.begin();
      iter!=pmap.end();
      iter++){
    CL.push_back(*(iter->second));
  }
  return CL;
}

/*!
  Get the 3 clostest neighbors of a particle (sorted)

  \param Po the particle
  \param NL the list of neighbors

  \todo Current implementation is lazy (NlogN), implement cN
*/
vector<SimpleParticle> ARandomAssembly::get3ClosestNeighbors(const SimpleParticle& Po,const vector<SimpleParticle>& NL)
{
  vector<SimpleParticle> CL;
  vector<SimpleParticle>::const_iterator iter;

  if(NL.size()<2){ // 0 or 1 neighbor -> no need to sort
    CL=NL;
  } else if (NL.size()==2) { // 2 neighbors 
    double dist1=(Po.getPos()-NL[0].getPos()).norm()-NL[0].getRad();
    double dist2=(Po.getPos()-NL[1].getPos()).norm()-NL[1].getRad();
    if(dist1<dist2){
      CL.push_back(NL[0]);
      CL.push_back(NL[1]);
    } else {
      CL.push_back(NL[1]);
      CL.push_back(NL[0]);
    }
  } else { // 3 or more neighbors
    map<double,SimpleParticle> nmap;

    for(iter=NL.begin();iter!=NL.end();iter++){
      double dist=(Po.getPos()-iter->getPos()).norm()-iter->getRad();
      nmap.insert(pair<double,SimpleParticle>(dist,*iter));
    }
    map<double,SimpleParticle>::iterator m_iter=nmap.begin();
    CL.push_back(m_iter->second);
    m_iter++;
    CL.push_back(m_iter->second);
    m_iter++;
    CL.push_back(m_iter->second);
  } 
  return CL;
}



/*!
  get particle closest to a particle (on surface separation)
  
  \param Po the particle
  \param NL the list of neighbors
*/
SimpleParticle ARandomAssembly::getClosestParticle(const SimpleParticle& Po, const vector<SimpleParticle>& NL)
{
  SimpleParticle CP=*(NL.begin());

  double dist=(Po.getPos()-CP.getPos()).norm()-CP.getRad();

  for(vector<SimpleParticle>::const_iterator citer=NL.begin();citer!=NL.end();citer++){
    double ndist=(Po.getPos()-citer->getPos()).norm()-citer->getRad();
    if(ndist<dist){
      CP=*citer;
      dist=ndist;
    }
  }
  return CP;
}

void ARandomAssembly::writeToVtkFile(const string& filename)
{
  std::cout << "FaultedBlock2D::writeToVtkFile( " << filename << " )" << std::endl;

  // open output file
  ofstream vtkfile(filename.c_str());
  // write header 
  vtkfile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n";
  vtkfile << "<UnstructuredGrid>\n";
  vtkfile << "<Piece NumberOfPoints=\"" << m_bpart.size() << "\" NumberOfCells=\"" << m_iset.size() << "\">\n";

  // write particle pos
  vtkfile << "<Points>\n";
  vtkfile << "<DataArray NumberOfComponents=\"3\" type=\"Float64\" format=\"ascii\">\n";
  for(vector<SimpleParticle>::const_iterator iter=m_bpart.begin();
      iter!=m_bpart.end();
      iter++){
    vtkfile << iter->getPos() << std::endl;
  }
  vtkfile << "</DataArray>\n";
  vtkfile << "</Points>\n";

   // --- write particle data ---
  // radius
  vtkfile << "<PointData Scalars=\"radius\">\n";
  vtkfile << "<DataArray type=\"Float64\" Name=\"radius\" NumberOfComponents=\"1\" format=\"ascii\">\n";
   for(vector<SimpleParticle>::const_iterator iter=m_bpart.begin();
      iter!=m_bpart.end();
      iter++){
     vtkfile << iter->getRad() << std::endl;
   }
  vtkfile << "</DataArray>\n";
  // tag
  vtkfile << "<DataArray type=\"Int32\" Name=\"particleTag\" NumberOfComponents=\"1\" format=\"ascii\">\n";
   for(vector<SimpleParticle>::const_iterator iter=m_bpart.begin();
      iter!=m_bpart.end();
      iter++){
     vtkfile << iter->getTag() << std::endl;
   }
  vtkfile << "</DataArray>\n";
  // id 
  vtkfile << "<DataArray type=\"Int32\" Name=\"Id\" NumberOfComponents=\"1\" format=\"ascii\">\n";
  for(vector<SimpleParticle>::const_iterator iter=m_bpart.begin();
      iter!=m_bpart.end();
      iter++){
     vtkfile << iter->getID() << std::endl;
   }
  vtkfile << "</DataArray>\n";
  vtkfile << "</PointData>\n";

  // wite bonds
  vtkfile << "<Cells>\n";
  vtkfile << "<DataArray type=\"Int32\" NumberOfComponents=\"1\" Name=\"connectivity\" format=\"ascii\">\n";
  for(set<BasicInteraction,BILess>::const_iterator iter=m_iset.begin();
      iter!=m_iset.end();
      iter++){
    vtkfile << iter->first() << " " << iter->second() << std::endl;
  }
  vtkfile << "</DataArray>";
  // offsets
  vtkfile << "<DataArray type=\"Int32\" NumberOfComponents=\"1\" Name=\"offsets\" format=\"ascii\">\n";
  for (size_t i = 1; i < m_iset.size()*2; i+=2) vtkfile << i+1 << "\n";
  vtkfile << "</DataArray>\n";
  // element type
  vtkfile << "<DataArray type=\"UInt8\" NumberOfComponents=\"1\" Name=\"types\" format=\"ascii\">\n";
  const int CELL_LINE_TYPE = 3;
  for (size_t i = 0; i < m_iset.size(); i++) vtkfile << CELL_LINE_TYPE << "\n";
  vtkfile << "</DataArray>\n";  
  vtkfile << "</Cells>\n";
  // bond data
  vtkfile << "<CellData>\n";
  vtkfile << "<DataArray type=\"Int32\" Name=\"tag\" NumberOfComponents=\"1\" format=\"ascii\">\n";
    for(set<BasicInteraction,BILess>::const_iterator iter=m_iset.begin();
      iter!=m_iset.end();
      iter++){
      vtkfile << iter->getTag() << std::endl;
  }
  vtkfile << "</DataArray>\n";  
  vtkfile << "</CellData>\n";
  // write footer
  vtkfile << "</Piece>\n";
  vtkfile << "</UnstructuredGrid>\n";
  vtkfile << "</VTKFile>\n";

  // close file
  vtkfile.close();
}
