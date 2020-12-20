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


#ifndef ESYS_LSMPARTICLEFIELDSAVERPRMSPY_H
#define ESYS_LSMPARTICLEFIELDSAVERPRMSPY_H

#include "Python/esys/lsm/FieldSaverPrmsPy.h"

//--- STL includes ---
#include <string>

//--- Boost includes ---
#include <boost/python.hpp>

namespace esys
{
  namespace lsm
  {
    class ParticleFieldSaverPrmsPy : public FieldSaverPrmsPy
    {
    public:
      ParticleFieldSaverPrmsPy(
        const std::string &fieldName,
        const std::string &fileName,
        const std::string &fileFormat,
        int beginTimeStep,
        int endTimeStep,
        int timeStepIncr
      );

    private:
    };

    class ParticleScalarFieldSaverPrmsPy : public ParticleFieldSaverPrmsPy
    {
    public:
      ParticleScalarFieldSaverPrmsPy(
        const std::string &fieldName,
        const std::string &fileName,
        const std::string &fileFormat,
        int beginTimeStep,
        int endTimeStep,
        int timeStepIncr
      );
    };

    class ParticleVectorFieldSaverPrmsPy : public ParticleFieldSaverPrmsPy
    {
    public:
      ParticleVectorFieldSaverPrmsPy(
        const std::string &fieldName,
        const std::string &fileName,
        const std::string &fileFormat,
        int beginTimeStep,
        int endTimeStep,
        int timeStepIncr
      );
    };

    class TaggedParticleScalarFieldSaverPrmsPy : public ParticleScalarFieldSaverPrmsPy
    {
    private:
      int m_tag;
      int m_mask;
      
    public:
      TaggedParticleScalarFieldSaverPrmsPy(
        const std::string &fieldName,
        const std::string &fileName,
        const std::string &fileFormat,
        int beginTimeStep,
        int endTimeStep,
        int timeStepIncr,
	      int tag,
        int mask
      );

      int getTag() const {return m_tag;}
      int getMask() const {return m_mask;}
    };

    class TaggedParticleVectorFieldSaverPrmsPy : public ParticleVectorFieldSaverPrmsPy
    {
    private:
      int m_tag;
      int m_mask;
      
    public:
      TaggedParticleVectorFieldSaverPrmsPy(
        const std::string &fieldName,
        const std::string &fileName,
        const std::string &fileFormat,
        int beginTimeStep,
        int endTimeStep,
        int timeStepIncr,
	      int tag,
	      int mask
      );

      int getTag() const {return m_tag;}
      int getMask() const {return m_mask;}
    };

    void exportParticleFieldSaverPrms();
  } // namespace lsm
} // namespace esys

#endif
