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

#ifndef __BMESH2DIP_H
#define __BMESH2DIP_H

class BMesh2DIP
{
public:
  double k;
  double brk;

  BMesh2DIP() : k(0.0), brk(0.0), m_name(), m_meshName()
  {
  }

  virtual ~BMesh2DIP()
  {
  }

  BMesh2DIP(
    const std::string& interactionName,
    const std::string& meshName,
    double normalK,
    double breakDistance
  ) : 
      k(normalK),
      brk(breakDistance),
      m_name(interactionName),
      m_meshName(meshName)
  {
  }

  void setMeshName(const std::string &meshName)
  {
    m_meshName = meshName;
  }

  const std::string &getMeshName() const
  {
    return m_meshName;
  }

  void setName(const std::string &name)
  {
    m_name = name;
  }

  const std::string &getName() const
  {
    return m_name;
  }

  virtual std::string getTypeString() const
  {
    return "Bonded";
  }

private:
  std::string m_name;
  std::string m_meshName;
};


#endif //__BMESH2DIP_H
