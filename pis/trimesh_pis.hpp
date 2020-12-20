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

/*!
  constructor

  \param mesh_p
  \param ppa_p
*/
template <class ParticleType> 
TriMesh_PIS<ParticleType>::TriMesh_PIS(TriMesh* mesh_p,ParallelParticleArray<ParticleType>* ppa_p):
  AParallelInteractionStorage(ppa_p)
{
  m_mesh=mesh_p;
}

template <class ParticleType> 
TriMesh_PIS<ParticleType>::~TriMesh_PIS()
{}

/*!
  add excluding IG

  \param exig_p
*/
template<class ParticleType>
void TriMesh_PIS<ParticleType>::addExIG(AParallelInteractionStorage* exig_p)
{
  std::cerr << "Setting an exclusing in a Mesh2D interaction group is not supported" << std::endl;
}


template <class ParticleType> 
AFieldSlave* TriMesh_PIS<ParticleType>::generateNewScalarFieldSlave(TML_Comm*,const string&,int,int,int,int)
{
  AFieldSlave* new_fs=NULL;

  return new_fs;
}

template <class ParticleType> 
AFieldSlave* TriMesh_PIS<ParticleType>::generateNewVectorFieldSlave(TML_Comm*,const string&,int,int,int,int)
{
  AFieldSlave* new_fs=NULL;

  return new_fs;
}
