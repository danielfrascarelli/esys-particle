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


#ifndef __ETRIMESHIP_H
#define __ETRIMESHIP_H

class ETriMeshIP
{
public:

  ETriMeshIP() : k(0.0), m_name(), m_meshName()
  {
  }

  virtual ~ETriMeshIP()
  {
  }
  
  ETriMeshIP(
    const std::string &interactionName,
    const std::string &meshName,
    double normalK
  )
    : k(normalK),
      m_name(interactionName),
      m_meshName(meshName)
  {
  }

  void setName(const std::string &name)
  {
    m_name = name;
  }

  const std::string &getName() const
  {
    return m_name;
  }

  void setMeshName(const std::string &name)
  {
    m_meshName = name;
  }

  const std::string &getMeshName() const
  {
    return m_meshName;
  }

  virtual std::string getTypeString() const
  {
    return "Elastic";
  }

public:
  double k;
private:
  std::string m_name;
  std::string m_meshName;
};

#endif // __ETRIMESHIP_H
