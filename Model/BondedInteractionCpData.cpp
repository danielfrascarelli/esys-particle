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

#include "Model/BondedInteraction.h"
#include "Model/BondedInteractionCpData.h"
#include "Model/RotBondedInteraction.h"
#include "Model/RotThermBondedInteraction.h"
#include "Model/BrittleBeamSC.h"
#include "Model/BrittleBeamDZC.h"

BondedInteractionCpData::BondedInteractionCpData() : m_p1Id(-1), m_p2Id(-1), m_tag(-1)
{
}

BondedInteractionCpData::BondedInteractionCpData(const CBondedInteraction &bondedInteraction)
  : m_p1Id(bondedInteraction.getPairID().first),
    m_p2Id(bondedInteraction.getPairID().second),
    m_tag(bondedInteraction.getTag())
{
}

BondedInteractionCpData::BondedInteractionCpData(const CRotBondedInteraction &bondedInteraction)
  : m_p1Id(bondedInteraction.getPairID().first),
    m_p2Id(bondedInteraction.getPairID().second),
    m_tag(bondedInteraction.getTag())
{
}

BondedInteractionCpData::BondedInteractionCpData(
  const CRotThermBondedInteraction &bondedInteraction
)
  : m_p1Id(bondedInteraction.getPairID().first),
    m_p2Id(bondedInteraction.getPairID().second),
    m_tag(bondedInteraction.getTag())
{
}

BondedInteractionCpData::BondedInteractionCpData(const BrittleBeamSCInteraction &bondedInteraction)
  : m_p1Id(bondedInteraction.getPairID().first),
    m_p2Id(bondedInteraction.getPairID().second),
    m_tag(bondedInteraction.getTag())
{}

BondedInteractionCpData::BondedInteractionCpData(const BrittleBeamDZCInteraction &bondedInteraction)
  : m_p1Id(bondedInteraction.getPairID().first),
    m_p2Id(bondedInteraction.getPairID().second),
    m_tag(bondedInteraction.getTag())
{}

BondedInteractionCpData::BondedInteractionCpData(
  ParticleId particle1Id,
  ParticleId particle2Id,
  InteractionTag interactionTag
)
  : m_p1Id(particle1Id),
    m_p2Id(particle2Id),
    m_tag(interactionTag)
{
}

void BondedInteractionCpData::set(const CBondedInteraction &bondedInteraction)
{
  m_p1Id = bondedInteraction.getPairID().first;
  m_p2Id = bondedInteraction.getPairID().second;
  m_tag  = bondedInteraction.getTag();
}

void BondedInteractionCpData::set(ParticleId particle1Id, ParticleId particle2Id, InteractionTag interactionTag)
{
  m_p1Id = particle1Id;
  m_p2Id = particle2Id;
  m_tag  = interactionTag;
}

BondedInteractionCpData::ParticleId BondedInteractionCpData::getP1Id() const
{
  return m_p1Id;
}

BondedInteractionCpData::ParticleId BondedInteractionCpData::getP2Id() const
{
  return m_p2Id;
}

BondedInteractionCpData::InteractionTag BondedInteractionCpData::getTag() const
{
  return m_tag;
}

void BondedInteractionCpData::saveCheckPointData(std::ostream &oStream)
{
  const char delim = ' ';
  oStream
    << m_p1Id << delim
    << m_p2Id << delim
    << m_tag;
}

void BondedInteractionCpData::loadCheckPointData(std::istream &iStream)
{
  iStream
    >> m_p1Id
    >> m_p2Id
    >> m_tag;
}
