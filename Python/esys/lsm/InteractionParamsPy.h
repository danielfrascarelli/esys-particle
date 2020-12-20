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

#ifndef ESYS_LSMINTERACTIONPARAMSPY_H
#define ESYS_LSMINTERACTIONPARAMSPY_H

/*---------------------------------------------------
 *
 * wrapper classes for interaction parameters
 *
 *--------------------------------------------------*/

//--- project includes ---
#include "Model/BondedInteraction.h"
#include "Model/CappedBondedInteraction.h"
#include "Model/ElasticInteraction.h"
#include "Model/RotElasticInteraction.h"
#include "Model/FrictionInteraction.h"
#include "Model/SpringDashpotFrictionInteraction.h"
#include "Model/RotBondedInteraction.h"
#include "Model/VWFrictionInteraction.h"
#include "Model/RotFricInteraction.h"
#include "Model/RotThermBondedInteraction.h"
#include "Model/RotThermFricInteraction.h"
#include "Model/RotThermElasticInteraction.h"
#include "Model/DampingIGP.h"
#include "Model/LocalDampingIGP.h"
#include "Model/ABCDampingIGP.h"
#include "Model/BodyForceGroup.h"
#include "Model/HertzianElasticInteraction.h"
#include "Model/HertzianViscoElasticFrictionInteraction.h"
#include "Model/HertzianViscoElasticInteraction.h"
#include "Model/HertzMindlinInteraction.h"
#include "Model/HertzMindlinViscoInteraction.h"
#include "Model/LinearDashpotInteraction.h"
#include "Model/BrittleBeamSC.h"
#include "Model/BrittleBeamDZC.h"

//--- STL includes ---
#include <string>

namespace esys
{
  namespace lsm
  {
    class Vec3Py;

    /**
     * Base class for python-exposed interaction parameter classes.
     * Purely artificial base class to aid documentation.
     */
    class InteractionPrmsPy
    {
    public:
      InteractionPrmsPy();
    };

    /*!
      \class DampingPrmsPy
      \brief wrapper for CDampingIGP 
    */
    class DampingPrmsPy : public CDampingIGP
    {
    public:
      DampingPrmsPy(
         const std::string &type,
         const std::string &name,
         double viscosity,
         int maxIterations
      );
    };

    class LinDampingPrmsPy : public DampingPrmsPy
    {
    public:
      LinDampingPrmsPy(
         const std::string &name,
         double viscosity,
         int maxIterations
      );
    };

    class RotDampingPrmsPy : public DampingPrmsPy
    {
    public:
      RotDampingPrmsPy(
         const std::string &name,
         double viscosity,
         int maxIterations
      );
    };

    /*!
      \class LocalDampingPrmsPy
      \brief wrapper for CLocalDampingIGP 
    */
    class LocalDampingPrmsPy : public CLocalDampingIGP
    {
    public:
      LocalDampingPrmsPy(
         const std::string &name,
         double viscosity
      );
    };

    /*!
      \class RotLocalDampingPrmsPy
      \brief wrapper for CLocalDampingIGP 
    */
    class RotLocalDampingPrmsPy : public CLocalDampingIGP
    {
    public:
      RotLocalDampingPrmsPy(
         const std::string &name,
         double viscosity
      );
    };

    /*!
      \class ABCDampingPrmsPy
      \brief wrapper for ABCDampingIGP 
    */
    class ABCDampingPrmsPy : public ABCDampingIGP
    {
    public:
      ABCDampingPrmsPy(
         const std::string &type,
         const std::string &name,
         double viscosity,
         int maxIterations,
         const Vec3& vref,
         const Vec3& pos,
         const Vec3& normal,
         double c1
      );
    };

    /*!
      \class NRotBondPrmsPy
      \brief wrapper for CBondedIGP
    */
    class NRotBondPrmsPy : public CBondedIGP
    {
    public:
      NRotBondPrmsPy(const std::string&,double,double,int);
      NRotBondPrmsPy(const std::string&,double,double,int,bool);
    };

    /*!
      \class CappedNRotBondPrmsPy
      \brief wrapper for CCappedBondedIGP
    */
    class CappedNRotBondPrmsPy : public CCappedBondedIGP
    {
    public:
      CappedNRotBondPrmsPy(const std::string&,double,double,double,int);
    };

  /*!
      \class NRotShortBondPrmsPy
      \brief wrapper for CBondedIGP (used in construction of short bonded IG)
    */
    class NRotShortBondPrmsPy : public CBondedIGP
    {
    public:
      NRotShortBondPrmsPy(const std::string&,double,double,int);
    };

    /*!
      \class NRotElasticPrmsPy
      \brief wrapper for CElasticIGP
    */
    class NRotElasticPrmsPy : public CElasticIGP
    {
    public:
      NRotElasticPrmsPy(
        const  std::string &name,
        double normalK
      );

      NRotElasticPrmsPy(
        const  std::string &name,
        double normalK,
        bool   scaling
      );
    };

    /*!
      \class HertzianElasticPrmsPy
      \brief wrapper for CHertzianElasticIGP
    */
    class HertzianElasticPrmsPy : public CHertzianElasticIGP
    {
    public:
      HertzianElasticPrmsPy(const std::string&,double,double);
    };
 
    /*!
      \class HertzianViscoElasticFrictionPrmsPy
      \brief wrapper for CHertzianViscoElasticFrictionIGP
    */
    class HertzianViscoElasticFrictionPrmsPy
    : public CHertzianViscoElasticFrictionIGP
    {
    public:
      HertzianViscoElasticFrictionPrmsPy(
        const std::string&,
        double,
        double,
        double,
        double,
        double
      );
    };

    /*!
      \class HertzianViscoElasticPrmsPy
      \brief wrapper for CHertzianViscoElasticIGP
    */
    class HertzianViscoElasticPrmsPy : public CHertzianViscoElasticIGP
    {
    public:
      HertzianViscoElasticPrmsPy(const std::string&,double,double,double);
    };
 
    /*!
      \class HertzMindlinPrmsPy
      \brief wrapper for CHertzMindlinIGP
    */
    class HertzMindlinPrmsPy
    : public CHertzMindlinIGP
    {
    public:
      HertzMindlinPrmsPy(
        const std::string&,
        double,
        double,
        double
      );
    };

    /*!
      \class HertzMindlinViscoPrmsPy
      \brief wrapper for CHertzMindlinViscoIGP
    */
    class HertzMindlinViscoPrmsPy
    : public CHertzMindlinViscoIGP
    {
    public:
      HertzMindlinViscoPrmsPy(
        const std::string&,
        double,
        double,
        double,
        double
      );
    };

    /*!
      \class LinearDashpotPrmsPy
      \brief wrapper for CLinearDashpotIGP
    */
    class LinearDashpotPrmsPy : public CLinearDashpotIGP
    {
    public:
      LinearDashpotPrmsPy(const std::string&,double,double);
    };

    /*!
      \class NRotFrictionPrmsPy
      \brief wrapper for CFrictionIGP
    */
    class NRotFrictionPrmsPy : public CFrictionIGP
    {
    public: 
      NRotFrictionPrmsPy(
        const std::string &name,
        double normalK,
        double dynamicMu,
        double shearK,
        bool scaling
      );

      NRotFrictionPrmsPy(
        const std::string &name,
        double normalK,
        double dynamicMu,
        double shearK
      );
    };

    /*!
      \class SpringDashpotFrictionPrmsPy
      \brief wrapper for SpringDashpotFrictionIGP
    */
    class SpringDashpotFrictionPrmsPy : public CSpringDashpotFrictionIGP
    {
    public:
      SpringDashpotFrictionPrmsPy(
        const std::string &name,
        double youngModulus,
        double poissonRatio,
        double cor,
        double dynamicMu
      );
    };

    /*!
      \class BrittleBeamPrmsPy
      \brief wrapper for CRotBondedIGP
    */
    class BrittleBeamPrmsPy : public CRotBondedIGP
    {
    private:
      std::string m_name; // unused ??
    public:
      BrittleBeamPrmsPy(
        const  std::string &name,
        double youngsModulus,
        double poissonsRatio,
        double cohesion,
        double tanAngle,
        int    aTag
      );

      BrittleBeamPrmsPy(
        const  std::string &name,
        double youngsModulus,
        double poissonsRatio,
        double cohesion,
        double tanAngle,
        int    aTag,
        bool   meanR_scaling
      );

      BrittleBeamPrmsPy(
        const  std::string &name,
        double youngsModulus,
        double poissonsRatio,
        double cohesion,
        double tanAngle,
        int    aTag,
        bool   meanR_scaling,
        double truncated,
        double beta1,
        double beta2
      );
    };

    
    /*!
      \class BrittleBeamSCPrmsPy
      \brief wrapper for BrittleBeamSCIGP
    */
    class BrittleBeamSCPrmsPy : public BrittleBeamSCIGP
    {
    private:
      std::string m_name; // unused ??
    public:
      BrittleBeamSCPrmsPy(
        const  std::string &name,
        double youngsModulus,
        double poissonsRatio,
        double tensileStrength,
        double frictionAngle,
	double tensileCutoff,
        int    aTag
      );
   };

    /*!
      \class BrittleBeamDZCPrmsPy
      \brief wrapper for BrittleBeamDZCIGP
    */
    class BrittleBeamDZCPrmsPy : public BrittleBeamDZCIGP
    {
    private:
      std::string m_name; // unused ??
    public:
      BrittleBeamDZCPrmsPy(
        const  std::string &name,
        double youngsModulus,
        double poissonsRatio,
        double tensileStrength,
        double frictionAngle,
    double tensileCutoff,
    double compressCutoff,
    double beta1,
    double beta2,
        int    aTag
      );
   };
    
    /*!
      \class RotBondPrmsPy
      \brief wrapper for CRotBondedIGP
    */
    class RotBondPrmsPy : public CRotBondedIGP
    {
    private:
      std::string m_name; // unused ??
    public:
      RotBondPrmsPy(
        const  std::string &name,
        double normalK,
        double shearK,
        double torsionK,
        double bendingK,
        double normalBrkForce,
        double shearBrkForce,
        double torsionBrkForce,
        double bendingBrkForce,
        int    aTag
      );

      RotBondPrmsPy(
        const  std::string &name,
        double normalK,
        double shearK,
        double torsionK,
        double bendingK,
        double normalBrkForce,
        double shearBrkForce,
        double torsionBrkForce,
        double bendingBrkForce,
        int    aTag,
        bool   scaling
      );

      RotBondPrmsPy(
        const  std::string &name,
        double normalK,
        double shearK,
        double torsionK,
        double bendingK,
        double normalBrkForce,
        double shearBrkForce,
        double torsionBrkForce,
        double bendingBrkForce,
        int    aTag,
        bool   scaling,
        bool   meanR_scaling
      );

      RotBondPrmsPy(
        const  std::string &name,
        double normalK,
        double shearK,
        double torsionK,
        double bendingK,
        double normalBrkForce,
        double shearBrkForce,
        double torsionBrkForce,
        double bendingBrkForce,
        int    aTag,
        bool   scaling,
        bool   meanR_scaling,
        double truncated
      );
    };

    /*!
      \class FrictionPrmsPy
      \brief wrapper for CRotFrictionIGP
    */
    class FrictionPrmsPy : public CRotFrictionIGP
    {
    public:
      FrictionPrmsPy(
        const  std::string &name,
        double youngsModulus,
        double poissonsRatio,
        double dynamicMu,    // sliding frictional coefficient
        double staticMu,     // max static frictional coefficient
        bool   rigid,
        bool   meanR_scaling
      );

      FrictionPrmsPy(
        const  std::string &name,
        double youngsModulus,
        double poissonsRatio,
        double dynamicMu,    // sliding frictional coefficient
        double staticMu     // max static frictional coefficient
      );

      FrictionPrmsPy(
        const  std::string &name,
        double youngsModulus,
        double poissonsRatio,
        double dynamicMu,    // sliding frictional coefficient
        double staticMu,     // max static frictional coefficient
        bool   rigid
      );
    };
    
    /*!
      \class RotFrictionPrmsPy
      \brief wrapper for CRotFrictionIGP
    */
    class RotFrictionPrmsPy : public CRotFrictionIGP
    {
    public:
      RotFrictionPrmsPy(
        const  std::string &name,
        double normalK,
        double dynamicMu,    // sliding frictional coefficient
        double staticMu,     // max static frictional coefficient
        double shearK
      );

      RotFrictionPrmsPy(
        const  std::string &name,
        double normalK,
        double dynamicMu,    // sliding frictional coefficient
        double staticMu,     // max static frictional coefficient
        double shearK,
        bool   scaling
      );

      RotFrictionPrmsPy(
        const  std::string &name,
        double normalK,
        double dynamicMu,    // sliding frictional coefficient
        double staticMu,     // max static frictional coefficient
        double shearK,
        bool   scaling,
        bool   rigid
      );

      RotFrictionPrmsPy(
        const  std::string &name,
        double normalK,
        double dynamicMu,    // sliding frictional coefficient
        double staticMu,     // max static frictional coefficient
        double shearK,
        bool   scaling,
        bool   rigid,
        bool   meanR_scaling
      );
    };

    /*!
      \class RotElasticPrmsPy
      \brief wrapper for CRotElasticIGP
    */
    class RotElasticPrmsPy : public CRotElasticIGP
    {
    public:
      RotElasticPrmsPy(
        const  std::string& name,
        double normalK
      );

      RotElasticPrmsPy(
        const  std::string& name,
        double normalK,
        bool   scaling
      );
    };

    /*!
      \class RotThermElasticPrmsPy
      \brief wrapper for CRotThermElasticIGP
    */
    class RotThermElasticPrmsPy : public CRotThermElasticIGP
    {
    public:
      RotThermElasticPrmsPy(
        const  std::string &name,
        double normalK,
        double diffusivity
      );
    };

    /*!
      \class RotThermFrictionPrmsPy
      \brief wrapper for CRotThermFrictionIGP
    */
    class RotThermFrictionPrmsPy : public CRotThermFrictionIGP
    {
    public:
      RotThermFrictionPrmsPy(
        const std::string &name,
        double normalK,
        double dynamicMu,    // sliding frictional coefficient
        double staticMu,     // max static frictional coefficient
        double shearK,
        double diffusivity
      );
    };

    /*!
      \class RotThermBondPrmsPy
      \brief wrapper for CRotThermBondedIGP
    */
    class RotThermBondPrmsPy : public CRotThermBondedIGP
    {
    public:
      static const int INVALID_BOND_TAG;
      RotThermBondPrmsPy(
        const std::string &name,
        double normalK,
        double shearK,
        double torsionK,
        double bendingK,
        double normalBrkForce,
        double shearBrkForce,
        double torsionBrkForce,
        double bendingBrkForce,
        double diffusivity,
        int    aTag
      );
    };
    
    /*!
      \class GravityParamsPy 
      \brief wrapper for GravityIGP
    */
    class GravityPrmsPy : public GravityIGP
    {
    public:
      GravityPrmsPy(
        const std::string &name,
        const Vec3Py &acceleration
      );
    };

    /*!
      \class BuoyancyParamsPy 
      \brief wrapper for BuoyancyIGP
    */
    class BuoyancyPrmsPy : public BuoyancyIGP
    {
    public:
      BuoyancyPrmsPy(
        const std::string &name,
        const Vec3Py &acceleration,
        const double &fluidDensity,
        const double &fluidHeight
      );
    };

    /*!
      \class VWFrictionPrmsPy
      \brief wrapper for VWFrictionIGP
    */
    class VWFrictionPrmsPy : public VWFrictionIGP
    {
    public:
      VWFrictionPrmsPy(
        const  std::string &name,
        double normalK,
        double dynamicMu,
        double shearK,
        double alpha
      );
    };

    void exportInteractionPrms();
    
  } // namespace lsm
} // namespace esys

#endif // ESYS_LSMINTERACTIONPARAMSPY_H
