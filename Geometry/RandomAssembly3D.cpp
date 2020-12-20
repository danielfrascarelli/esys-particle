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

#include "Foundation/console.h"
#include "Geometry/RandomAssembly3D.h"

//-- project includes --
#include "Geometry/Sphere3d.h"

//-- STL includes --
#include <map>
#include <utility>

using std::map;
using std::pair;

/*!
  get closest plane to a particle

  \param Po the particle
*/
Plane3D ARandomAssembly3D::getClosestPlane(const SimpleParticle& Po)
{
  //cout << "getClosestPlane : " << Po.getPos() << endl;
  Plane3D PL=*(Borders.begin());

  Vec3 PoPos=Po.getPos();
  double dist=(PL.sep(PoPos));
  //cout << "plane: " << PL.GetO() << PL.GetN() << dist << endl;
  for(vector<Plane3D>::iterator iter=Borders.begin();iter!=Borders.end();iter++){
    double ndist=iter->sep(PoPos);
    //cout << "plane: " << iter->GetO() << iter->GetN() << ndist << endl;
    if(ndist<dist){
      PL=*iter;
      dist=ndist;
    }
  }
  //cout << "closest plane: " << PL.GetO() << PL.GetN() << dist << endl;

  return PL;
}

/*!
  Find a fit for a sphere using the list of neigbors

  \param Po the particle to fit
  \param NL the list of neighbors
*/
bool ARandomAssembly3D::findAFit(SimpleParticle& Po, const vector<SimpleParticle>& NL)
{
  //cout << "findAFit - 4 particles\n";
  bool find_a_fit ;
  Vec3 M;
  double r;
  int id=Po.getID();

  if(NL.size()<4){
    find_a_fit=false;
    cout << "less than 4 neighbors" << endl; // can't happen
  } else {
     Vec3 Pos1=NL[0].getPos();
     Vec3 Pos2=NL[1].getPos();
     Vec3 Pos3=NL[2].getPos();
     Vec3 Pos4=NL[3].getPos();
     double r1=NL[0].getRad();
     double r2=NL[1].getRad();
     double r3=NL[2].getRad();
     double r4=NL[3].getRad();

     find_a_fit=Sphere3D::FillIn(Pos1,Pos2,Pos3,Pos4,r1,r2,r3,r4,M,r);
     Po=SimpleParticle(M,r,id);
     //cout << "found " << M << " , " << r << endl; 
  }
   
  return find_a_fit ;
}

/*!
  Find a fit for a sphere using the list of neigbor list and a plane

  \param Po the particle to fit
  \param NL the list of neighbors
  \param L the Plane3D
*/
bool ARandomAssembly3D::findAFit(SimpleParticle& Po, const vector<SimpleParticle>& NL, const Plane3D& L)
{
  //cout << "findAFit - 3 particles, 1 plane \n";
  bool find_a_fit ;
  Vec3 M;
  double r;
  int id=Po.getID();

  if(NL.size()<3){
    find_a_fit=false;
    cout << "less than 3 neighbors" << endl; // can't happen
  } else {
      Vec3 Pos1=NL[0].getPos();
      Vec3 Pos2=NL[1].getPos();
      Vec3 Pos3=NL[2].getPos();
      double r1=NL[0].getRad();
      double r2=NL[1].getRad();
      double r3=NL[2].getRad();
      Vec3 WallO=L.GetO();
      Vec3 WallD=L.GetW();
      
      find_a_fit=Sphere3D::FillInWP(Pos1,Pos2,Pos3,WallO,WallD,r1,r2,r3,M,r);
      
      Po=SimpleParticle(M,r,id);
      //if(find_a_fit) cout << "found " << M << " , " << r << endl; 
  }
  return find_a_fit;
}

/*!
  check if Po is within the Space and is not crossing any boundary or 
  overlapping with other particles.

  \param Po the particle
*/
bool ARandomAssembly3D::checkAFit(const SimpleParticle& Po)
{
  bool fail=false;
  
  // check vs. radius
  if((Po.getRad()<m_rmin) || (Po.getRad()>m_rmax)) {
    fail=true;
  }
  // check vs. borders
  double px=Po.getPos().X();
  double py=Po.getPos().Y();
  double pz=Po.getPos().Z();
  //double r=Po.getRad();
  if((px<m_xmin-m_small_value) || (px>m_xmax+m_small_value) || 
     (py<m_ymin-m_small_value) || (py>m_ymax+m_small_value) ||
     (pz<m_zmin-m_small_value) || (pz>m_zmax+m_small_value)){
    fail=true;
    //cout << "Fail : outside" << endl;
  }
  // check vs. all neighbors
  if(!fail){
    vector<SimpleParticle> NL=getNeighborList(Po); // get the NL here because Po may have moved during fit
    vector<SimpleParticle>::const_iterator iter=NL.begin();
    while(!fail && iter!=NL.end()){
      double dist=(Po.getPos()-iter->getPos()).norm()+m_small_value;;
      if(dist < (Po.getRad()+iter->getRad())){
        fail=true;
        //cout << "Fail : particle collision" << endl;
      }
      iter++;
    }
  }
  // check vs. closest plane
  if(!fail) {
    Plane3D L=getClosestPlane(Po);
    fail=(Po.getRad()-L.sep(Po.getPos()))>m_small_value;
    //cout << "Fail: collision with plane" << endl;
  }

  return !fail ;
}

/*!
  Fill the space in the skeleton after it has been seeded

  \param tries the number of tries
*/
void ARandomAssembly3D::fillSpace(int tries)
{
  int countfail=0;
  int countfound=0;
  int countwithplane=0;
  int trywithplane=0;
  //bool fail,findfit;

  while(countfail<tries){
    bool findfit=false;
    bool foundwithplane=false;
    Vec3 P=getAPoint();
    double r=m_random(m_rmin,m_rmax);
    SimpleParticle Po=SimpleParticle(P,r,getNParts());
    vector<SimpleParticle> T4=getClosestNeighbors(Po,4);    // get list of closest neighbors (max 4)
    Plane3D L=getClosestPlane(Po);                            // get closest plane/line
    if(T4.size()>3){ // at least 4 neighbors 
      SimpleParticle Pi=T4[0];
      double ndist=(Po.getPos()-Pi.getPos()).norm();
      if( ndist==0.0){
        findfit=false;
      } else {
        if( ndist < Pi.getRad()){ // if Po inside Pi -> move Po to the surface of Pi
          Vec3 npos=Pi.getPos()+((Po.getPos()-Pi.getPos())*(Pi.getRad()/ndist));
          Po.moveTo(npos);
        }
        Vec3 PoPos=Po.getPos();
        double dist_p=L.sep(PoPos);
        double dist_3=(PoPos-T4[3].getPos()).norm()-T4[3].getRad();
        if(dist_p>dist_3){  // 4th particle closer than plane -> fit 4 particles
          findfit=findAFit(Po,T4);
        } else { // plane closer than 4th particle -> fit 3 particles + plane
          findfit=findAFit(Po,T4,L);
          foundwithplane=findfit;
          trywithplane++;
        }
      }
    } else if(T4.size()==3) { // 3 neighbors  -> try  3 particles + plane
      SimpleParticle Pi=T4[0];
      double ndist=(Po.getPos()-Pi.getPos()).norm();
      if( ndist==0.0){
        findfit=false;
            } else {
        if( ndist < Pi.getRad()){ // if Po inside Pi -> move Po to the surface of Pi
          Vec3 npos=Pi.getPos()+((Po.getPos()-Pi.getPos())*(Pi.getRad()/ndist));
          Po.moveTo(npos);
        }
        findfit=findAFit(Po,T4,L);
        foundwithplane=findfit;
        trywithplane++;
      }
    } 
    if(findfit){ // found something, check
      findfit=checkAFit(Po);
      foundwithplane=foundwithplane & findfit;
    } 
    if(findfit){  // found & checked -> insert
      insertParticle(Po);
      if(countfail*10>tries){
        cout << "found particle " << Po.getID() << " after " << countfail << " tries" << endl;
      }
      countfail=0;
      countfound++;
      if(foundwithplane) countwithplane++;
    } else {
      countfail++;
    }
  }
  
  console.Info() << "inserted " << countfound << " Particles" << "\n";
  console.Info() << "found " << trywithplane << " with 3 Particles and 1 Plane, accepted " << countwithplane << "\n";
}
