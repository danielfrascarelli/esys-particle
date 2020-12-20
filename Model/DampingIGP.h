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

#ifndef __DAMPING_IGP_H
#define __DAMPING_IGP_H

// -- project includes --
#include "Model/IGParam.h"
#include "Foundation/vec3.h"
// -- STL includes --
#include <string>

using std::string;
 
/*!
  \brief Interaction group parameters for CDampingGroup
*/
class CDampingIGP : public AIGParam
{
protected:
  string m_type; // type of damping (rot/lin)
  Vec3 m_vref; //!< reference velocity
  double m_visc; //!< artificial viscosity
  double m_dt;   //!< time step
  int m_max_iter; //!< max nr. of iterations

public:
  CDampingIGP();
  CDampingIGP(const string& type,
    const string &name,
    double viscosity,
    double dt,
    int maxIteractions,
    const Vec3 &refVelocity = Vec3::ZERO
  );
  
  virtual void  packInto(CVarMPIBuffer*) const;
  void setType(const string& type){m_type=type;}
  void setVRef(const Vec3 V){m_vref=V;}
  Vec3 getVRef()const{return m_vref;}
  void setVisc(double v){m_visc=v;}
  double getVisc()const{return m_visc;}
  void setTimeStep(double t){m_dt=t;}
  void setTimeStepSize(double t){setTimeStep(t);}
  double getTimeStep()const{return m_dt;}
  void setMaxIter(int mi){m_max_iter=mi;}
  int getMaxIter()const {return m_max_iter;}
  
  virtual std::string getTypeString() const {return m_type;}
};

CDampingIGP* extractDampingIGP(AMPIBuffer*);

#endif //__DAMPING_IGP_H
