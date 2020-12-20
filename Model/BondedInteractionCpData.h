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

#ifndef __BONDEDINTERACTIONCPDATA_H
#define __BONDEDINTERACTIONCPDATA_H

#include "Parallel/CheckPointable.h"

class CBondedInteraction;
class CRotBondedInteraction;
class CRotThermBondedInteraction;
class BrittleBeamSCInteraction;
class BrittleBeamDZCInteraction;

/**
 * Helper class for checkpointing BondedInteraction data.
 */
class BondedInteractionCpData : public esys::lsm::CheckPointable
{
public:
  typedef int ParticleId;
  typedef int InteractionTag;

  BondedInteractionCpData();

  virtual ~BondedInteractionCpData()
  {
  }

  BondedInteractionCpData(const CBondedInteraction &bondedInteraction);
  BondedInteractionCpData(const CRotBondedInteraction &bondedInteraction);
  BondedInteractionCpData(const CRotThermBondedInteraction &bondedInteraction);
  BondedInteractionCpData(const BrittleBeamSCInteraction &bondedInteraction);
  BondedInteractionCpData(const BrittleBeamDZCInteraction &bondedInteraction);

  BondedInteractionCpData(
    ParticleId particle1Id,
    ParticleId particle2Id,
    InteractionTag interactionTag
  );

  void set(const CBondedInteraction &bondedInteraction);

  void set(ParticleId particle1Id, ParticleId particle2Id, InteractionTag interactionTag);

  ParticleId getP1Id() const;

  ParticleId getP2Id() const;

  InteractionTag getTag() const;

  virtual void saveCheckPointData(std::ostream &oStream);

  virtual void loadCheckPointData(std::istream &iStream);

private:
  ParticleId     m_p1Id;
  ParticleId     m_p2Id;
  InteractionTag m_tag;
};

#endif //__BONDEDINTERACTIONCPDATA_H
