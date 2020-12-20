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


#include "Model/BTriMeshIP.h"

BTriMeshIP::BTriMeshIP() : k(0.0), brk(0.0), m_name(), m_meshName()
{
}

BTriMeshIP::BTriMeshIP(
  const std::string &interactionName,
  const std::string &meshName,
  double normalK,
  double breakDistance
) : k(normalK), brk(breakDistance), m_name(interactionName), m_meshName(meshName)
{
}
  
BTriMeshIP::~BTriMeshIP()
{
}

void BTriMeshIP::setMeshName(const std::string &meshName)
{
  m_meshName = meshName;
}

const std::string &BTriMeshIP::getMeshName() const
{
  return m_meshName;
}

void BTriMeshIP::setName(const std::string &name)
{
  m_name = name;
}

const std::string &BTriMeshIP::getName() const
{
  return m_name;
}

std::string BTriMeshIP::getTypeString() const
{
  return "Bonded";
}

MeshBuildPrms::MeshBuildPrms()
{
}

MeshBuildPrms::~MeshBuildPrms()
{
}

MeshTagBuildPrms::MeshTagBuildPrms() : MeshBuildPrms(), m_tag(-1), m_mask(0)
{
}

MeshTagBuildPrms::MeshTagBuildPrms(int tag, int mask) : MeshBuildPrms(), m_tag(tag), m_mask(mask)
{
}

std::string MeshTagBuildPrms::getTypeString() const
{
  return "BuildByTag";
}

MeshGapBuildPrms::MeshGapBuildPrms() : MeshBuildPrms(), m_maxGap(0.0)
{
}

MeshGapBuildPrms::MeshGapBuildPrms(double maxGap) : MeshBuildPrms(), m_maxGap(maxGap)
{
}

std::string MeshGapBuildPrms::getTypeString() const
{
  return "BuildByGap";
}

