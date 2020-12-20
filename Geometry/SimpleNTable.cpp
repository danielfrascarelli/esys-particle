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


#include "Geometry/SimpleNTable.h"

//------ ASimpleNTable member functions ----

ASimpleNTable::ASimpleNTable()
  : m_data(NULL),
    m_p0(Vec3::ZERO),
    m_dim(0.0),
    m_numInsertedParticles(0)  
{
}

ASimpleNTable::~ASimpleNTable()
{
  if (m_data != NULL) delete [] m_data;
}

/*!
  get all particles near a given position

  \param pos the position
*/
const vector<SimpleParticle>* ASimpleNTable::getNeighbors(const Vec3& pos) const
{
  return &(m_data[index(pos)]);
}

/*!
  Add particle to all neighorlists it belongs to

  \param cbp the particle 
*/
void ASimpleNTable::insertParticle(SimpleParticle cbp)
{
  //cout << "ASimpleNTable::insertParticle(" << cbp << ")" << endl;
  const vector<int> idx=allidx(cbp.getPos());
  if (idx.size() > 0) {
    m_numInsertedParticles++;
  }
  for(vector<int>::const_iterator citer=idx.begin();citer!=idx.end();citer++){
    m_data[*citer].push_back(cbp);
  }
  insertParticleCircular(cbp);
}

int ASimpleNTable::getNumInsertedParticles() const
{
  return m_numInsertedParticles;
}

/*!
  get particle closest to given position

  \param pos the position

  \warning doesn't check if position is in space
*/
int ASimpleNTable::getClosestParticleID(const Vec3& pos) const
{
  vector<SimpleParticle>* parts=&(m_data[index(pos)]);
  int id=-1;

  double dist=m_dim;
  for(vector<SimpleParticle>::iterator iter=parts->begin();
      iter!=parts->end();
      iter++){
    double ndist=(pos-iter->getPos()).norm();
    if(ndist<dist){
      dist=ndist;
      id=iter->getID();
    }
  }
  
  return id;
}


//------ CSimple2DNTable member functions ----
/*!
  Constructor 

  \param pos position of the (xmin,ymin) point
  \param dim size of the space
  \param r grid spacing
*/
CSimple2DNTable::CSimple2DNTable(const Vec3& pos,const Vec3& dim,double r,bool xcirc,bool ycirc)
{
  //cout << "CSimple2DNTable(" << pos << " , " << dim << " , " << r << ")" << endl;
  m_xsize=int(ceil(dim.X()/r));
  m_ysize=int(ceil(dim.Y()/r));
  m_xcirc=xcirc;
  m_ycirc=ycirc;
  m_p0=pos;
  m_dim=r;
  // pad in circular directions
  if(m_xcirc) {
    m_xsize+=2;
    m_p0-=Vec3(r,0.0,0.0);
    m_xshift=Vec3(dim.X(),0.0,0.0);
  }
  if(m_ycirc) {
    m_ysize+=2;
    m_p0-=Vec3(0.0,r,0.0);
    m_yshift=Vec3(0.0,dim.Y(),0.0);
  }
  //cout << "NT size : " << m_xsize << " , " << m_ysize << endl;
  m_data=new vector<SimpleParticle>[m_xsize*m_ysize];
}

/*!
  insert circular images of the particle

  \param cbp the particle
*/
void CSimple2DNTable::insertParticleCircular(SimpleParticle cbp)
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
  Return the grid index of a position.

  \param pos the position
  \warning does not check if pos is within space
*/
int CSimple2DNTable::index(const Vec3& pos) const
{
  int ix=int((pos.X()-m_p0.X())/m_dim);
  int iy=int((pos.Y()-m_p0.Y())/m_dim);
  int idx=ix+m_xsize*iy;

  return idx;
}

/*!
  Get all indices to which a particle at a given position will be added. 

  \param pos the position
  \warning does not check if pos is within space
*/
vector<int> CSimple2DNTable::allidx(const Vec3& pos) const
{
  vector<int> res;

  int ix=int((pos.X()-m_p0.X())/m_dim);
  int iy=int((pos.Y()-m_p0.Y())/m_dim);
  int idx=ix+m_xsize*iy;
  res.push_back(idx);
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
  return res;
}



/*!
  Put all interactions into a set

  \param iset the set into which to put them
*/  
void CSimple2DNTable::getInteractions(set<BasicInteraction,BILess>& iset,double dmax)
{
  for(int i=0;i<m_xsize;i++){
    for(int j=0;j<m_ysize;j++){
      int idx=i+m_xsize*j;
      //cout << "-- " << i << " , " << j << " , " << idx << endl;
      if(m_data[idx].size()>=2){
        for(vector<SimpleParticle>::iterator iter=m_data[idx].begin();
            iter!=m_data[idx].end()-1;
            iter++){
          for(vector<SimpleParticle>::iterator iter2=iter+1;
              iter2!=m_data[idx].end();
              iter2++){
            //	    cout << "Pair: " << iter->getID() << " - " << iter2 ->getID() << endl;
            if((iter->getPos()-iter2->getPos()).norm()<dmax*(iter->getRad()+iter2->getRad())){
              iset.insert(BasicInteraction(iter->getID(),iter2 ->getID()));
            }
          }
        }
      }
    }
  }
}

// debugging only 
void CSimple2DNTable::print()
{
  for(int i=0;i<m_xsize;i++){
    for(int j=0;j<m_ysize;j++){
      int idx=i+m_xsize*j;
      cout << "-- " << i << " , " << j << " , " << idx << endl;
      for(
        vector<SimpleParticle>::iterator iter=m_data[idx].begin();
        iter!=m_data[idx].end();
        iter++)
      {
        cout << iter->getPos() << " , " << iter->getRad() << endl;
      }
    }
  }
}
