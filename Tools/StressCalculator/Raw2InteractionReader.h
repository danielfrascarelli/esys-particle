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

#ifndef ESYS_LSMRAW2INTERACTIONREADER_H
#define ESYS_LSMRAW2INTERACTIONREADER_H

#include "Foundation/vec3.h"
#include "Foundation/StringUtil.h"
#include "Tools/StressCalculator/Contact.h"

#include <iostream>

namespace esys
{
  namespace lsm
  {
    class Raw2InteractionReader
    {
    public:
      typedef Contact InteractionData;

      Raw2InteractionReader(std::istream &iStream)\
        : m_pIStream(&iStream)
      {
      }

      bool hasNext() const
      {
        return ((m_pIStream->peek()) != std::istream::traits_type::eof());
      }

      InteractionData next()
      {
        InteractionData data;
        data.read(*m_pIStream);
        return data;
      }

      private:
        Raw2InteractionReader();
        Raw2InteractionReader(const Raw2InteractionReader &reader);
        Raw2InteractionReader &operator=(const Raw2InteractionReader &reader);

        std::istream *m_pIStream;
    };
  }
}

#endif
