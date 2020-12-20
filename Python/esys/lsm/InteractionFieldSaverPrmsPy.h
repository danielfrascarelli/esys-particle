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


#ifndef ESYS_LSMINTERACTIONFIELDSAVERPRMSPY_H
#define ESYS_LSMINTERACTIONFIELDSAVERPRMSPY_H

#include "Python/esys/lsm/FieldSaverPrmsPy.h"

//--- STL includes ---
#include <string>

//--- Boost includes ---
#include <boost/python.hpp>

namespace esys
{
  namespace lsm
  {
    class InteractionFieldSaverPrmsPy : public FieldSaverPrmsPy
    {
    public:
      InteractionFieldSaverPrmsPy(
        const std::string &interactionName,
        const std::string &fieldName,
        const std::string &fileName,
        const std::string &fileFormat,
        int beginTimeStep,
        int endTimeStep,
        int timeStepIncr
      );

      const std::string &getInteractionName() const
      {
        return m_interactionName;
      }
    private:
      std::string m_interactionName;
    };

    class InteractionScalarFieldSaverPrmsPy : public InteractionFieldSaverPrmsPy
    {
    public:
      InteractionScalarFieldSaverPrmsPy(
        const std::string &interactionName,
        const std::string &fieldName,
        const std::string &fileName,
        const std::string &fileFormat,
        int beginTimeStep,
        int endTimeStep,
        int timeStepIncr
      );
    };

    class CheckedInteractionScalarFieldSaverPrmsPy : public InteractionFieldSaverPrmsPy
    {
    public:
      CheckedInteractionScalarFieldSaverPrmsPy(
        const std::string &interactionName,
        const std::string &fieldName,
        const std::string &fileName,
        const std::string &fileFormat,
        int beginTimeStep,
        int endTimeStep,
        int timeStepIncr
      );
    };

    class TaggedInteractionScalarFieldSaverPrmsPy : public InteractionScalarFieldSaverPrmsPy
    {
    public:
      TaggedInteractionScalarFieldSaverPrmsPy(
        const std::string &interactionName,
        const std::string &fieldName,
        const std::string &fileName,
        const std::string &fileFormat,
        int beginTimeStep,
        int endTimeStep,
        int timeStepIncr,
	int tag,
	int mask
      );

      int getTag() const {return m_tag;};
      int getMask() const {return m_mask;};

    protected:
      int m_tag;
      int m_mask;
    };

    class InteractionVectorFieldSaverPrmsPy : public InteractionFieldSaverPrmsPy
    {
    public:
      InteractionVectorFieldSaverPrmsPy(
        const std::string &interactionName,
        const std::string &fieldName,
        const std::string &fileName,
        const std::string &fileFormat,
        int beginTimeStep,
        int endTimeStep,
        int timeStepIncr
      );
    };

    class CheckedInteractionVectorFieldSaverPrmsPy : public InteractionFieldSaverPrmsPy
    {
    public:
      CheckedInteractionVectorFieldSaverPrmsPy(
        const std::string &interactionName,
        const std::string &fieldName,
        const std::string &fileName,
        const std::string &fileFormat,
        int beginTimeStep,
        int endTimeStep,
        int timeStepIncr
      );
    };
    
    void exportInteractionFieldSaverPrms();
  } // namespace lsm
} // namespace esys

#endif
