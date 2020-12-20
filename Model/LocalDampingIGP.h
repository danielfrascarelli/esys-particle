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

#ifndef __LOCALDAMPING_IGP_H
#define __LOCALDAMPING_IGP_H

// -- project includes --
#include "Model/IGParam.h"
#include "Foundation/vec3.h"
// -- STL includes --
#include <string>

using std::string;

/*!
  \brief Interaction group parameters for CLocalDampingGroup
*/
class CLocalDampingIGP : public AIGParam
{
protected:
  string m_type; // type of damping (rot/lin)
  double m_visc; //!< damping coefficient
  double m_dt;   //!< time step

public:
  CLocalDampingIGP();
  CLocalDampingIGP(const string& type,
    const string &name,
    double viscosity,
    double dt
  );

  virtual void  packInto(CVarMPIBuffer*) const;
  void setType(const string& type){m_type=type;}
  void setVisc(double v){m_visc=v;}
  double getVisc()const{return m_visc;}
  void setTimeStep(double t){m_dt=t;}
  void setTimeStepSize(double t){setTimeStep(t);}
  double getTimeStep()const{return m_dt;}

  virtual std::string getTypeString() const {return m_type;}
};

CLocalDampingIGP* extractLocalDampingIGP(AMPIBuffer*);

#endif //__LOCALDAMPING_IGP_H
