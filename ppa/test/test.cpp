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

/* 
   test for parallel particle array
*/

//--- project includes ---
#include "pp_array.h"
#include "BasicParticle.h"
#include "comm_world.h"
#include "vec3.h"

//--- IO includes ---
#include <iostream>
using std::cout;

int main(int argc,char **argv)
{
  ParallelParticleArray<CBasicParticle> *ppa;
 
  MPI_Init(&argc,&argv);

  TML_CommWorld worldcomm;
  Vec3 min=Vec3(0.0,0.0,0.0);
  Vec3 max=Vec3(8.0,8.0,8.0);
  double range=2.0;

  ppa=new ParallelParticleArray<CBasicParticle>(&worldcomm,3,min,max,range);
  // cout << *ppa << endl;
  ppa->rebuild();
  
  ppa->forAllParticles(&CBasicParticle::setRad,1.0);
  ppa->forParticle(2,&CBasicParticle::setRad,1.0);
  vector<int> idv;
  ppa->forAllParticlesGet(idv,&CBasicParticle::getID);
 
  cout << "IDs " << idv.size() << endl;
  for(vector<int>::iterator iter=idv.begin();
      iter!=idv.end();
      iter++){
    cout << *iter << " ";
  }
  cout << endl;

  delete ppa;

  MPI_Finalize();
 
  return 0;
}
