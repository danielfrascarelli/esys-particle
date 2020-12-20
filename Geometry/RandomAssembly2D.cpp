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

#include "Geometry/RandomAssembly2D.h"

//-- project includes --
#include "Foundation/vec3.h"
#include "Geometry/SimpleParticle.h"

//-- STL includes --
#include <map>
#include <utility>

using std::map;
using std::pair;

//-- IO includes --
#include <iostream>

using std::cout;
using std::endl;


/*!
  get closest plane to a particle

  \param Po the particle
*/
Line *ARandomAssembly2D::getClosestPlane(const SimpleParticle& Po)
{
  //cout << "getClosestPlane : " << Po.getPos() << endl;
  vector<Line>::iterator PL=Borders.begin();

  Vec3 PoPos=Po.getPos();
  double dist=PL->sep(PoPos);
  //cout << "plane: " << PL.GetO() << PL.GetN() << dist << endl;
  for(vector<Line>::iterator iter=Borders.begin();iter!=Borders.end();iter++){
    double ndist=iter->sep(PoPos);
    //cout << "plane: " << iter->GetO() << iter->GetN() << ndist << endl;
    if(ndist<dist){
      PL=iter;
      dist=ndist;
    }
  }
  //cout << "closest plane: " << PL.GetO() << PL.GetN() << dist << endl;

  return &(*PL);
}


/*!
  Find a fit for a sphere using the list of neigbor list and a plane

  \param Po the particle to fit
  \param NL the list of neighbors
  \param PL the Plane

  \todo check for at least 2 particles
*/
bool ARandomAssembly2D::findAFit(SimpleParticle& Po, const vector<SimpleParticle>& NL, const Line& L)
{
  bool find_a_fit ;
  Vec3 M;
  double r;
  int id=Po.getID();

  Vec3 Pos1=NL[0].getPos();
  Vec3 Pos2=NL[1].getPos();
  Vec3 WallO=L.GetO();
  Vec3 WallD=L.GetU();

  find_a_fit=Sphere2D::FillInWP(Pos1,Pos2,WallO,WallD,NL[0].getRad(),NL[1].getRad(),M,r);

  Po=SimpleParticle(M,r,id);
  
  //cout << "fit with plane : " << L.GetO() << endl;

  return find_a_fit ;
}

/*!
  Find a fit for a sphere using the list of neigbors

  \param Po the particle to fit
  \param NL the list of neighbors

  \todo check for at least 3 particles
*/
bool ARandomAssembly2D::findAFit(SimpleParticle& Po, const vector<SimpleParticle>& NL)
{
  bool find_a_fit ;
  Vec3 M;
  double r;
  int id=Po.getID();

  Vec3 Pos1=NL[0].getPos();
  Vec3 Pos2=NL[1].getPos();
  Vec3 Pos3=NL[2].getPos();
  find_a_fit=Sphere2D::FillIn(Pos1,Pos2,Pos3,NL[0].getRad(),NL[1].getRad(),NL[2].getRad(),M,r);

  Po=SimpleParticle(M,r,id);

  return find_a_fit ;
}


/*!
  check if Po is within the Space and is not crossing any boundary or 
  overlapping with other particles.

  \param Po the particle
*/
bool ARandomAssembly2D::checkAFit(const SimpleParticle& Po)
{
  bool fail=false;
 
  // check vs. radius
  if((Po.getRad()<m_rmin) || (Po.getRad()>m_rmax)) {
    fail=true;
  }
  // check vs. borders
  double px=Po.getPos().X();
  double py=Po.getPos().Y();
  // double r=Po.getRad(); unused.
  if((px<m_xmin-m_small_value) || (px-m_small_value>m_xmax) || (py<m_ymin-m_small_value) || (py-m_small_value>m_ymax)){
    fail=true;
    //cout << "Fail : outside" << endl;
  }
  // check vs. all neighbors
  if(!fail){
    vector<SimpleParticle> NL=getNeighborList(Po); // get the NL here because Po may have moved during fit
    vector<SimpleParticle>::const_iterator iter=NL.begin();
    while(!fail && iter!=NL.end()){
      double dist=(Po.getPos()-iter->getPos()).norm()+m_small_value;;
      if(dist<(Po.getRad()+iter->getRad())){
        fail=true;
        //cout << "Fail : particle collision" << endl;
      }
      iter++;
    }
  }
  // check vs. closest plane
  if(!fail){
    Line *L=getClosestPlane(Po);
    fail=(Po.getRad()-L->sep(Po.getPos()))>m_small_value;
  }
  
  return !fail ;
}


    
/*!
  Fill the space in the skeleton after it has been seeded

  \param tries the number of tries
*/
void ARandomAssembly2D::fillSpace(int tries)
{
  int countfail=0;
  int countfound=0;
  //bool fail;
  
  while(countfail<tries){
    bool findfit=false;
    Vec3 P=getAPoint();
    double r=m_random(m_rmin,m_rmax);
    SimpleParticle Po=SimpleParticle(P,r,getNParts());
    vector<SimpleParticle> T3=getClosestNeighbors(Po,3); // get closest neighbors (max 3) 
    Line* L=getClosestPlane(Po);                            // get closest plane/line
    if(T3.size()>2){ // 3 neighbors, closest 
      SimpleParticle Pi=T3[0];
      double ndist=(Po.getPos()-Pi.getPos()).norm();
      if( ndist==0.0){
        findfit=false;
      } else {
        if( ndist < Pi.getRad()){ // if Po inside Pi -> move Po to the surface of Pi
          Vec3 npos=Pi.getPos()+((Po.getPos()-Pi.getPos())*(Pi.getRad()/ndist));
          Po.moveTo(npos);
        }
        Vec3 PoPos=Po.getPos();
        double dist_p=L->sep(PoPos);
        double dist_3=(PoPos-T3[2].getPos()).norm()-T3[2].getRad();
        if(dist_p>dist_3){  // 3rd particle closer than plane -> fit 3 particles
          findfit=findAFit(Po,T3);
        } else { // plane closer than 3rd particle -> fit 2 particles + plane
          findfit=findAFit(Po,T3,*L);
        }
      }
    } else if(T3.size()==2) { // 2 neighbors  -> try  2 particles + plane
      SimpleParticle Pi=T3[0];
      double ndist=(Po.getPos()-Pi.getPos()).norm();
      if( ndist==0.0){
        findfit=false;
      } else {
        if( ndist < Pi.getRad()){ // if Po inside Pi -> move Po to the surface of Pi
          Vec3 npos=Pi.getPos()+((Po.getPos()-Pi.getPos())*(Pi.getRad()/ndist));
          Po.moveTo(npos);
        }
        findfit=findAFit(Po,T3,*L);
      }
    } 
    if(findfit){ // found something, check
      findfit=checkAFit(Po);
    } 
    if(findfit){  // found & checked -> insert
      insertParticle(Po);
      //cout << "found particle " << Po.getID() << " after " << countfail << " tries" << endl;
      countfail=0;
      countfound++;
    } else {
      countfail++;
    }
  }
  
  console.Info() << "inserted " << countfound << " Particles" << "\n";
}


