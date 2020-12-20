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


#ifndef ESYS_LSMCONTACTCOLLECTION_H
#define ESYS_LSMCONTACTCOLLECTION_H

#include "Foundation/vec3.h"
#include "Foundation/StlIterator.h"
#include "Tools/StressCalculator/Contact.h"
#include "Tools/StressCalculator/Vec3Comparer.h"

#include <vector>
#include <map>

namespace esys
{
  namespace lsm
  {
    class ContactCollection
    {
    public:
      typedef esys::lsm::Contact                          Contact;
      typedef esys::lsm::Vec3XyzComparer                  Vec3Comparer;
      typedef std::vector<Contact>                        ContactVector;
      typedef ForwardIterator<ContactVector>              ContactIterator;
      typedef std::map<Vec3, ContactVector, Vec3Comparer> ContactMap;

      class ContactIteratorIterator
      {
      public:
        typedef ForwardIterator<ContactMap> ContactMapIterator;
        ContactIteratorIterator(ContactMap &map)
          : m_contactMapIterator(map)
        {
        }
        
        bool hasNext() const
        {
          return m_contactMapIterator.hasNext();
        }

        ContactIterator next()
        {
          return ContactIterator(m_contactMapIterator.next().second);
        }

        ContactIterator current() const
        {
          return ContactIterator(m_contactMapIterator.current().second);
        }

      private:
        ContactMapIterator m_contactMapIterator;
      };
      

      ContactCollection()
        : m_contactMap()
      {
      }

      void addContact(const Contact &contact)
      {
        ContactMap::iterator it = m_contactMap.find(contact.getCentrePos1());
        if (it == m_contactMap.end()) {
          it = (m_contactMap.insert(ContactMap::value_type(contact.getCentrePos1(), ContactVector()))).first;
        }
        it->second.push_back(contact);
      }

      void addContactWithReverse(const Contact &contact)
      {
        addContact(contact);
        addContact(
          Contact(
            contact.getParticle2(),
            contact.getParticle1(),
            contact.getForcePos(),
            contact.getForce()*-1.0
          )
        );
      }

      template <typename TmplInteractionData>
      void addInteraction(const TmplInteractionData &interactionData)
      {
        if (interactionData.getForce() != Vec3::ZERO) {
          addContactWithReverse(
            Contact(
              interactionData.getParticle1(),
              interactionData.getParticle2(),
              interactionData.getForcePos(),
              interactionData.getForce()
            )
          );
        }
      }

      template <typename TmplIterator>
      void addInteractions(TmplIterator &iterator)
      {
        while (iterator.hasNext())
        {
          addInteraction(iterator.next());
        }
      }

      ContactVector getContactVector(const Vec3 &centrePt) const
      {
        ContactMap::const_iterator it = m_contactMap.find(centrePt);
        if (it != m_contactMap.end()) {
          return it->second;
        }
        return ContactVector();
      }
      
      ContactIteratorIterator getContactIteratorIterator()
      {
        return ContactIteratorIterator(m_contactMap);
      }

    private:
      ContactMap m_contactMap;
    };
  }
}

#endif
