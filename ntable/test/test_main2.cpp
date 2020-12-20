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

// --- project includes ---
#include "ntable.h"
#include "nt_slab.h"
#include "nt_block.h"
#include "BasicParticle.h"
#include "vec3.h"

//--- STL includes ---
#include <vector>

using std::vector;

//--- IO includes ---
#include <iostream>

using std::cout;
using std::endl;

int main(int argc,char **argv)
{
  NeighborTable<CBasicParticle> NT(3,3,4,1.0,Vec3(0.0,0.0,0.0),0,0,0);

  NT.insert(CBasicParticle(1,Vec3(0.50,1.50,2.50),1.0));
  NT.insert(CBasicParticle(2,Vec3(0.54,0.82,2.34),1.0));
  NT.insert(CBasicParticle(3,Vec3(0.46,2.00,0.14),1.0));
  NT.insert(CBasicParticle(4,Vec3(4.00,2.00,1.00),1.0)); // outside
  NT.insert(CBasicParticle(5,Vec3(1.08,1.38,3.19),1.0));
  NT.insert(CBasicParticle(6,Vec3(0.49,0.31,2.49),1.0));
  NT.insert(CBasicParticle(7,Vec3(1.26,1.49,1.40),1.0));

  cout << NT << endl;

  // test slabs and slab iterators
  NTSlab<CBasicParticle> NTS1=NT.xy_slab(2);
  cout << "xy_slab(2), size:" << NTS1.size() << endl;
  for(NTSlab<CBasicParticle>::iterator iter=NTS1.begin();
      iter!=NTS1.end();
      iter++){
    cout << iter->getID() << " ";
  }
  cout << endl;
  cout << "xy_slab(2) backwards" << endl;
  for(NTSlab<CBasicParticle>::iterator iter=NTS1.rbegin();
      iter!=NTS1.rend();
      iter--){
    cout << iter->getID() << " ";
  }
  cout << endl;
  NTSlab<CBasicParticle> NTS2=NT.xz_slab(1);
  cout << "xz_slab(1)" << endl;
  for(NTSlab<CBasicParticle>::iterator iter=NTS2.begin();
      iter!=NTS2.end();
      iter++){
    cout << iter->getID() << " ";
  }
  cout << endl;
  NTSlab<CBasicParticle> NTS3=NT.yz_slab(1);
  cout << "yz_slab(1)" << endl;
  for(NTSlab<CBasicParticle>::iterator iter=NTS3.begin();
      iter!=NTS3.end();
      iter++){
    cout << iter->getID() << " ";
  }
  cout << endl;
  
  // test insert
  NTS1.insert(NTS1.end(),CBasicParticle(8,Vec3(1.29,1.47,1.42),1.0));
  cout << NT << endl;
  
  // test nblist
  T_Handle<NeighborTable<CBasicParticle>::pairlist> plist=NT.getFullList();

  cout << "nr. of pairs found: " << plist->size() << endl;

  for(NeighborTable<CBasicParticle>::pairlist::iterator iter=plist->begin();
      iter!=plist->end();
      iter++){
    cout << iter->first->getID() << "," << iter->second->getID() << endl;
  }
  plist->erase(plist->begin(),plist->end());

  cout << "erased" << endl;

  // test access by id
  CBasicParticle *P1=NT.ptr_by_id(2);
  cout << "P1 : " << P1->getPos() << endl;

  // test block
  NTBlock<CBasicParticle> NTB=NT.block(1,1,1,2,1,3);
  cout << NTB << endl;
  
  // test block iter
  cout << "block iter" << endl;
  for(NTBlock<CBasicParticle>::iterator iter=NTB.begin();
      iter!=NTB.end();
      iter++){
    cout << iter->getID() << " ";
  }
  cout << endl;

  // inner block
  NTBlock<CBasicParticle> IB=NT.inner();
  cout << IB << endl;
 
  // test particles along plane
  cout << "particles at plane (1,0,0)(1,0,0)" << endl;
  Vec3 Orig=Vec3(1.0,0.0,0.0);
  Vec3 Normal=Vec3(1.0,0.0,0.0);
  T_Handle<NeighborTable<CBasicParticle>::particlelist> parlist=NT.getParticlesAtPlane(Orig,Normal);

  for(NeighborTable<CBasicParticle>::particlelist::iterator iter=parlist->begin();
      iter!=parlist->end();
      iter++){
    cout << (*iter)->getID() << " ";
  }
  cout << endl;

  // test getAllParticles
  parlist=NT.getAllParticles();
  cout << "all particles " << endl;
  for(NeighborTable<CBasicParticle>::particlelist::iterator iter=parlist->begin();
      iter!=parlist->end();
      iter++){
    cout << (*iter)->getID() << " ";
  }
  cout << endl;

  // test getNearest
  Vec3 pos=Vec3(1.0,1.0,1.0);
  cout << "particle closest to " << pos << " " << (NT.getNearestPtr(pos))->getID() << endl;  

  return 0;
}
