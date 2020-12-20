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

#include "Geometry/SimpleNTable3D.h"

/*!
  Return the grid index of a position.

  \param pos the position
  \warning does not check if pos is within space
*/
int CSimple3DNTable::index(const Vec3& pos) const
{
  int ix=int((pos.X()-m_p0.X())/m_dim);
  int iy=int((pos.Y()-m_p0.Y())/m_dim);
  int iz=int((pos.Z()-m_p0.Z())/m_dim);
  int idx=ix+m_xsize*iy+iz*m_xsize*m_ysize;

  return idx; 
}

/*!
  Get all indices to which a particle at a given position will be added. 

  \param pos the position
  \warning does not check if pos is within space
*/
vector<int> CSimple3DNTable::allidx(const Vec3& pos) const
{
  vector<int> res;

  int ix=int((pos.X()-m_p0.X())/m_dim);
  int iy=int((pos.Y()-m_p0.Y())/m_dim);
  int iz=int((pos.Z()-m_p0.Z())/m_dim);
  int idx=ix+m_xsize*iy+iz*m_xsize*m_ysize;

  res.push_back(idx);
  // iz-layer
  if(ix>0){
    res.push_back(idx-1);
    if(iy>0){
      res.push_back(idx-m_xsize-1);
    }
    if(iy<m_ysize-1){
      res.push_back(idx+m_xsize-1);
    }
  }
  if(ix<m_xsize-1){
    res.push_back(idx+1);
    if(iy>0){
      res.push_back(idx-m_xsize+1);
    }
    if(iy<m_ysize-1){
      res.push_back(idx+m_xsize+1);
    }
  }
  if(iy>0){
    res.push_back(idx-m_xsize);
  }
  if(iy<m_ysize-1){
    res.push_back(idx+m_xsize);
  }
  // iz-1 layer
  if(iz>0){  
    int idb=idx-m_xsize*m_ysize;
    res.push_back(idb);
    if(ix>0){
      res.push_back(idb-1);
      if(iy>0){
        res.push_back(idb-m_xsize-1);
      }
      if(iy<m_ysize-1){
        res.push_back(idb+m_xsize-1);
      }
    }
    if(ix<m_xsize-1){
      res.push_back(idb+1);
      if(iy>0){
        res.push_back(idb-m_xsize+1);
      }
      if(iy<m_ysize-1){
        res.push_back(idb+m_xsize+1);
      }
    }
    if(iy>0){
      res.push_back(idb-m_xsize);
    }
    if(iy<m_ysize-1){
      res.push_back(idb+m_xsize);
    }
  }
  // iz+1 layer
  if(iz<m_zsize-1){  
    int idb=idx+m_xsize*m_ysize;
    res.push_back(idb);
    if(ix>0){
      res.push_back(idb-1);
      if(iy>0){
        res.push_back(idb-m_xsize-1);
      }
      if(iy<m_ysize-1){
        res.push_back(idb+m_xsize-1);
      }
    }
    if(ix<m_xsize-1){
      res.push_back(idb+1);
      if(iy>0){
        res.push_back(idb-m_xsize+1);
      }
      if(iy<m_ysize-1){
        res.push_back(idb+m_xsize+1);
      }
    }
    if(iy>0){
      res.push_back(idb-m_xsize);
    }
    if(iy<m_ysize-1){
      res.push_back(idb+m_xsize);
    }
  }
  return res;
}

/*!
  insert circular images of the particle

  \param cbp the particle
*/
void CSimple3DNTable::insertParticleCircular(SimpleParticle cbp)
{
  if(m_xcirc){
    if (int((cbp.getPos().X()-m_p0.X())/m_dim)==1) {
      cbp.moveTo(cbp.getPos()+m_xshift);
      const vector<int> idx=allidx(cbp.getPos());
      for(vector<int>::const_iterator citer=idx.begin();citer!=idx.end();citer++){
        m_data[*citer].push_back(cbp);
      }
    } else if (int((cbp.getPos().X()-m_p0.X())/m_dim)==m_xsize-2) {
      cbp.moveTo(cbp.getPos()-m_xshift);
      const vector<int> idx=allidx(cbp.getPos());
      for(vector<int>::const_iterator citer=idx.begin();citer!=idx.end();citer++){
        m_data[*citer].push_back(cbp);
      }
    }
  }
}


/*!
  Constructor 

  \param pos position of the (xmin,ymin,zmin) point
  \param dim size of the space
  \param r grid spacing
*/
CSimple3DNTable::CSimple3DNTable(const Vec3& pos,const Vec3& dim ,double r ,bool xcirc,bool ycirc,bool zcirc)
{
  //cout << "CSimple3DNTable(" << pos << " , " << dim << " , " << r << ")" << endl;
  m_xsize=int(ceil(dim.X()/r));
  m_ysize=int(ceil(dim.Y()/r));
  m_zsize=int(ceil(dim.Z()/r));
  m_xcirc=xcirc;
  m_ycirc=ycirc;
  m_zcirc=zcirc;
  m_p0=pos;
  m_dim=r;
  // pad in circular directions
  if (m_xcirc) {
    m_xsize+=2;
    m_p0-=Vec3(r,0.0,0.0);
    m_xshift=Vec3(dim.X(),0.0,0.0);
  }
  if (m_ycirc) {
    m_ysize+=2;
    m_p0-=Vec3(0.0,r,0.0);
    m_yshift=Vec3(0.0,dim.Y(),0.0);
  }
  if (m_zcirc) {
    m_zsize+=2;
    m_p0-=Vec3(0.0,0.0,r);
    m_yshift=Vec3(0.0,0.0,dim.Z());
  }
  //cout << "NT size : " << m_xsize << " , " << m_ysize << " , " << m_zsize << endl;
  m_data=new vector<SimpleParticle>[m_xsize*m_ysize*m_zsize];
}

/*!
  Put all interactions into a set

  \param iset the set into which to put them
  \param dmax max distance for the creation of an interaction
*/
void CSimple3DNTable::getInteractions(set<BasicInteraction,BILess>& iset, double dmax)
{
  for(int i=0;i<m_xsize;i++){
    for(int j=0;j<m_ysize;j++){
      for(int k=0;k<m_zsize;k++){
        int idx=i+m_xsize*j+k*m_xsize*m_ysize;
        //cout << "-- " << i << " , " << j << " , " << idx << endl;
        if(m_data[idx].size()>=2){
          for(vector<SimpleParticle>::iterator iter=m_data[idx].begin();
              iter!=m_data[idx].end()-1;
              iter++)
          {
            for (
              vector<SimpleParticle>::iterator iter2=iter+1;
              iter2!=m_data[idx].end();
              iter2++)
            {
              // cout << "Pair: " << iter->getID() << " - " << iter2 ->getID() << endl;
              if((iter->getPos()-iter2->getPos()).norm() < dmax*(iter->getRad()+iter2->getRad()))
              {
                iset.insert(BasicInteraction(iter->getID(),iter2->getID()));
              }
            }
          }
        }
      }
    }
  }
}

// debugging only 
void CSimple3DNTable::print()
{ 
  for(int i=0;i<m_xsize;i++){
    for(int j=0;j<m_ysize;j++){
      for(int k=0;k<m_zsize;k++){
        int idx=i+m_xsize*j+k*m_xsize*m_ysize;
        //cout << "-- " << i << " , " << j << " , " << " , " << j << " , " << idx << endl;
        for(vector<SimpleParticle>::iterator iter=m_data[idx].begin();
            iter!=m_data[idx].end();
            iter++){
          cout << iter->getPos() << " , " << iter->getRad() << endl;
        }
      }
    }
  }
}
