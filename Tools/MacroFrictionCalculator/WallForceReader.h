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


#ifndef ESYS_LSMWALLFORCEREADER_H
#define ESYS_LSMWALLFORCEREADER_H

#include "Foundation/vec3.h"
#include "Foundation/StringUtil.h"

#include <vector>
#include <iostream>
#include <stdexcept>

namespace esys
{
  namespace lsm
  {
    class WallForcesRecord
    {
    public:
      typedef int WallId;
      class WallForce
      {
      public:
        WallForce() : m_id(-1), m_force()
        {
        }
        
        WallForce(WallId id, const Vec3 &force) : m_id(id), m_force(force)
        {
        }

        WallForce(
          StringUtil::StringVector::const_iterator begin,
          StringUtil::StringVector::const_iterator end
        )
          : m_id(-1), m_force()
        {
          if (end-begin >= 4) {
            m_id = StringUtil::to<int>(*begin);
            for (int i = 0; i < 3; i++) {
              m_force[i] = StringUtil::to<double>(*(begin + i + 1));
            }
          }
        }
                
        const WallId &getId() const
        {
          return m_id;
        }
        
        const Vec3 getForce() const
        {
          return m_force;
        }
      private:
        WallId m_id;
        Vec3   m_force;
      };

      WallForcesRecord(const std::string &line)
      {
        parseLine(line);
      }

      void parseLine(const std::string &line)
      {
        StringUtil::StringVector elemVector = StringUtil::splitStrings(StringUtil::trim(line), " ");
        const int elemsPerWallForce = 4;
        m_wallForceVector.clear();
        m_wallForceVector.reserve(elemVector.size()/elemsPerWallForce);
        if ((elemVector.size() % elemsPerWallForce) == 0)
        {
          for (
            StringUtil::StringVector::const_iterator it = elemVector.begin();
            it != elemVector.end();
            it += elemsPerWallForce
          )
          {
            const WallForce wallForce(it, it + elemsPerWallForce);
            if (wallForce.getId() >= static_cast<int>(m_wallForceVector.size()))
            {
              m_wallForceVector.insert(m_wallForceVector.end(), wallForce.getId()-m_wallForceVector.size()+1, WallForce());
            }
            m_wallForceVector[wallForce.getId()] = wallForce;
          }
        }
        else
        {
          std::stringstream msg;
          msg 
            << "Record '"
            << line
            << "' does not contain a number of elements divisible by "
            << elemsPerWallForce;
          throw std::runtime_error(msg.str());
        }
      }
      
      const WallForce &get(WallId id) const
      {
        return m_wallForceVector[id];
      }

    private:
      typedef std::vector<WallForce> WallForceVector;
      WallForceVector m_wallForceVector;
    };

    class WallForceReader
    {
    public:
      typedef std::pair<Vec3,Vec3>     WallForcePair;
      typedef WallForcesRecord::WallId WallId;
      

      WallForceReader(WallId wallId1, WallId wallId2, std::istream &iStream)
        : m_pIStream(&iStream),
          m_lineBuffer(),
          m_wallId1(wallId1),
          m_wallId2(wallId2)
      {
      }

      WallForceReader(WallId wallId1, WallId wallId2)
        : m_pIStream(NULL),
          m_lineBuffer(),
          m_wallId1(wallId1),
          m_wallId2(wallId2)
      {
      }
      
      void setStream(std::istream &iStream)
      {
        m_pIStream = &iStream;
      }

      bool hasNext() const
      {
        return ((m_pIStream != NULL) && ((m_pIStream->peek()) != std::istream::traits_type::eof()));
      }

      WallForcePair next()
      {
        m_pIStream->getline(m_lineBuffer, MAX_LINE_SIZE);
        WallForcesRecord record(m_lineBuffer);
        return WallForcePair(record.get(m_wallId1).getForce(), record.get(m_wallId2).getForce());
      }

    private:
      static const int MAX_LINE_SIZE = 4096;
    
      std::istream *m_pIStream;
      char         m_lineBuffer[MAX_LINE_SIZE];
      WallId       m_wallId1;
      WallId       m_wallId2;
    };
  }
}

#endif
