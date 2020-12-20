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

#include <mpi.h>
#include <boost/version.hpp>
#include <boost/python.hpp>
#include "Python/esys/lsm/InteractionParamsPy.h"
#include "Python/esys/lsm/util/Vec3Py.h"

namespace esys
{
  namespace lsm
  {
    InteractionPrmsPy::InteractionPrmsPy()
    {
    }

    /*!
      constructor for DampingPrmsPy

      \param type the type of damping, "Damping" or "RotDamping"
      \param name 
      \param viscosity the damping constant
      \param dt the time step
      \param maxIteration the maximum number of iterations
    */
    DampingPrmsPy::DampingPrmsPy(
      const std::string &type,
      const std::string &name,
      double viscosity,
      int maxIterations
    ) :
      CDampingIGP(type,name, viscosity, 0.0, maxIterations)
    {
    }

    LinDampingPrmsPy::LinDampingPrmsPy(
      const std::string &name,
      double viscosity,
      int maxIterations
    )
      : DampingPrmsPy("Damping", name, viscosity, maxIterations)
    {
    }

    RotDampingPrmsPy::RotDampingPrmsPy(
      const std::string &name,
      double viscosity,
      int maxIterations
    )
      : DampingPrmsPy("RotDamping", name, viscosity,  maxIterations)
    {
    }

    LocalDampingPrmsPy::LocalDampingPrmsPy(
      const std::string &name,
      double viscosity
    )
      : CLocalDampingIGP ("LocalDamping",name,viscosity,0.0)
    {
    }

    RotLocalDampingPrmsPy::RotLocalDampingPrmsPy(
      const std::string &name,
      double viscosity
    )
      : CLocalDampingIGP ("RotLocalDamping",name,viscosity,0.0)
    {
    }

    /*!
      constructor for ABCDampingPrmsPy

      \param name 
      \param viscosity the damping constant
      \param dt the time step
      \param maxIteration the max. nr of iterations
    */
    ABCDampingPrmsPy::ABCDampingPrmsPy(
      const std::string &type,
      const std::string &name,
      double viscosity,
      int maxIterations,
      const Vec3& vref,
      const Vec3& pos,
      const Vec3& normal,
      double c1
    )
      : ABCDampingIGP(type,name, viscosity, 0.0, maxIterations,vref,pos,normal,c1)
    {
    }

    /*!
      constructor with scaling parameter for NRotBondPrmsPy
    
      \param name
      \param normalK
      \param breakDistance
      \param aTag
      \param scaling scaling of normal stiffness with particle size
    */
    NRotBondPrmsPy::NRotBondPrmsPy(
      const  std::string &name,
      double normalK,
      double breakDistance,
      int    aTag,
      bool   scaling
    )
      : CBondedIGP(name, aTag, normalK, breakDistance,scaling)
    {
    }

    /*!
      Constructor without scaling parameter for NRotBondPrmsPy. Scaling is set 
      to "true" by default and a warning message is given.
    
      \param name
      \param normalK
      \param breakDistance
      \param aTag
    */
    NRotBondPrmsPy::NRotBondPrmsPy(
      const  std::string &name,
      double normalK,
      double breakDistance,
      int    aTag
    )
      : CBondedIGP(name, aTag, normalK, breakDistance,true)
    {
      std::cerr << "\n--- WARNING ---\n";
      std::cerr << "ESyS-Particle 2.0 by default scales elastic stiffness";
      std::cerr << " according to particle dimensions.\n";
      std::cerr << "To disable scaling set \"scaling=False\" in the";
      std::cerr << " NRotBondedPrms argument list.\n";
      std::cerr << "To remove this warning set \"scaling=True\".\n";
      std::cerr << "For more information about scaling";
      std::cerr << " refer to the ESyS-Particle Tutorial available at:\n";
      std::cerr << "https://wiki.geocomp.uq.edu.au/index.php/ESyS-Particle\n" << endl;
    }
    
    /*!
      constructor for CappedNRotBondPrmsPy
    
      \param name
      \param normalK
      \param breakDistance
      \param maxForce
      \param aTag
    */
    CappedNRotBondPrmsPy::CappedNRotBondPrmsPy(
      const std::string &name,
      double normalK,
      double breakDistance,
      double maxForce,
      int aTag
    )
      : CCappedBondedIGP(name, aTag, normalK, breakDistance,maxForce)
    {
    }
    
    /*!
      constructor for NRotShortBondPrmsPy
    
      \param name
      \param normalK
      \param breakDistance
      \param aTag
    */
    NRotShortBondPrmsPy::NRotShortBondPrmsPy(
      const std::string &name,
      double normalK,
      double breakDistance,
      int aTag
    )
      : CBondedIGP(name, aTag, normalK, breakDistance)
    {
    }
    
   /*!
      constructor for NRotElasticPrmsPy
    
      \param name
      \param normalK
      \param scaling scaling of normal stiffness with particle size
   */
    NRotElasticPrmsPy::NRotElasticPrmsPy(
      const  std::string &name,
      double normalK,
      bool   scaling
    ) :
      CElasticIGP(name, normalK, scaling)
    {
    }

    /*!
      Constructor without scaling parameter for NRotElasticPrmsPy. Scaling is set 
      to "true" by default and a warning message is given.
    
      \param name
      \param normalK
   */
    NRotElasticPrmsPy::NRotElasticPrmsPy(
      const std::string &name,
      double normalK
    ) :
      CElasticIGP(name, normalK)
    {
      std::cerr << "\n--- WARNING ---\n";
      std::cerr << "ESyS-Particle 2.0 by default scales elastic stiffness";
      std::cerr << " according to particle dimensions.\n";
      std::cerr << "To disable scaling set \"scaling=False\" in the";
      std::cerr << " NRotElasticPrms argument list.\n";
      std::cerr << "To remove this warning set \"scaling=True\".\n";
      std::cerr << "For more information about scaling";
      std::cerr << " refer to the ESyS-Particle Tutorial available at:\n";
      std::cerr << "https://wiki.geocomp.uq.edu.au/index.php/ESyS-Particle\n" << endl;
    } 

    /*!
      constructor for HertzianElasticPrmsPy
    
      \param name
      \param E Youngs Modulus
      \param nu poisson ratio
   */
    HertzianElasticPrmsPy::HertzianElasticPrmsPy(
      const std::string &name,
      double E,
      double nu
    )
      : CHertzianElasticIGP(name, E, nu)
    {}

    /*!
      constructor for HertzianViscoElasticFrictionPrmsPy
    
      \param name
      \param A Damping Constant
      \param E Youngs Modulus
      \param nu poisson ratio
      \param dynamicMu friction coefficient
      \param shearK spring constant
   */
    HertzianViscoElasticFrictionPrmsPy::HertzianViscoElasticFrictionPrmsPy(
      const std::string &name,
      double A,
      double E,
      double nu,
      double dynamicMu,
      double shearK
    )
      : CHertzianViscoElasticFrictionIGP(name, A, E, nu, dynamicMu, shearK, 0.0)
    {}

    /*!
      constructor for HertzianViscoElasticPrmsPy
    
      \param name
      \param A Damping Constant
      \param E Youngs Modulus
      \param nu poisson ratio
   */
    HertzianViscoElasticPrmsPy::HertzianViscoElasticPrmsPy(
      const std::string &name,
      double A,
      double E,
      double nu
    )
      : CHertzianViscoElasticIGP(name, A, E, nu)
    {}

    /*!
      constructor for HertzMindlinPrmsPy

      \param name
      \param E Youngs Modulus
      \param nu poisson ratio
      \param dynamicMu friction coefficient
   */
    HertzMindlinPrmsPy::HertzMindlinPrmsPy(
      const std::string &name,
      double E,
      double nu,
      double dynamicMu
    )
      : CHertzMindlinIGP(name, E, nu, dynamicMu, 0.0)
    {}

    /*!
      constructor for HertzMindlinViscoPrmsPy

      \param name
      \param E Youngs Modulus
      \param nu poisson ratio
      \param dynamicMu friction coefficient
      \param cor restitution coefficient
   */
    HertzMindlinViscoPrmsPy::HertzMindlinViscoPrmsPy(
      const std::string &name,
      double E,
      double nu,
      double dynamicMu,
      double cor
    )
      : CHertzMindlinViscoIGP(name, E, nu, dynamicMu, cor, 0.0)
    {}

    /*!
      constructor for LinearDashpotPrmsPy
    
      \param name
      \param damp damping coefficient ("viscosity")
      \param cutoff interaction range, relative to particle radii
   */
    LinearDashpotPrmsPy::LinearDashpotPrmsPy(
      const std::string &name,
      double damp,
      double cutoff
    ) :
      CLinearDashpotIGP(name, damp, cutoff)
    {}
    
    /*!
      constructor for NRotFrictionPrmsPy
    
      \param name
      \param normalK
      \param dynamicMu
      \param shearK
      \param scaling scaling of normal stiffness with particle size
    */
    NRotFrictionPrmsPy::NRotFrictionPrmsPy(
      const std::string &name,
      double normalK,
      double dynamicMu,
      double shearK,
      bool scaling
    )
      :CFrictionIGP(name, normalK, dynamicMu, shearK, 0.0,scaling)
    {
    }

    /*!
      Constructor without scaling parameter for NRotFrictionPrmsPy. Scaling is set 
      to "true" by default and a warning message is given.
    
      \param name
      \param normalK
      \param dynamicMu
      \param shearK
    */
    NRotFrictionPrmsPy::NRotFrictionPrmsPy(
      const std::string &name,
      double normalK,
      double dynamicMu,
      double shearK
      ) : CFrictionIGP(name, normalK, dynamicMu, shearK, 0.0,false)
    {
      std::cerr << "\n--- WARNING ---\n";
      std::cerr << "ESyS-Particle 2.0 by default scales elastic stiffness";
      std::cerr << " according to particle dimensions.\n";
      std::cerr << "To disable scaling set \"scaling=False\" in the";
      std::cerr << " NRotFrictionPrms argument list.\n";
      std::cerr << "To remove this warning set \"scaling=True\".\n";
      std::cerr << "For more information about scaling";
      std::cerr << " refer to the ESyS-Particle Tutorial available at:\n";
      std::cerr << "https://wiki.geocomp.uq.edu.au/index.php/ESyS-Particle\n" << endl;
    }     

    /*!
      constructor for SpringDashpotFrictionPrmsPy
    
      \param name
      \param youngModulus
      \param poissonRatio
      \param cor
      \param dynamicMu
    */
    SpringDashpotFrictionPrmsPy::SpringDashpotFrictionPrmsPy(
      const std::string &name,
      double youngModulus,
      double poissonRatio,
      double cor,
      double dynamicMu
    )
      :CSpringDashpotFrictionIGP(name, youngModulus, poissonRatio, cor, dynamicMu, 0.0)
    {
    }

    /*!
      Constructor without scaling parameter for RotBondPrmsPy. Scaling is set 
      to "true" by default and a warning message is given.
    
      \param name
      \param normalK
      \param shearK
      \param torsionK
      \param bendingK
      \param normalBrkForce  
      \param shearBrkForce
      \param torsionBrkForce
      \param bendingBrkForce
      \param aTag
    */
    RotBondPrmsPy::RotBondPrmsPy(
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
    )
      : CRotBondedIGP(
          name,normalK,shearK,torsionK,bendingK,normalBrkForce,shearBrkForce,
          torsionBrkForce,bendingBrkForce,aTag,true,true,1.0,1.0,1.0)
    {
      std::cerr << "\n--- WARNING ---\n";
      std::cerr << "ESyS-Particle 2.x by default scales elastic stiffness";
      std::cerr << " according to particle dimensions.\n";
      std::cerr << "To disable scaling set \"scaling=False\" in the";
      std::cerr << " RotBondPrms argument list.\n";
      std::cerr << "To remove this warning set \"scaling=True\".\n";
      std::cerr << "For more information about scaling";
      std::cerr << " refer to the ESyS-Particle Tutorial available at:\n";
      std::cerr << "https://wiki.geocomp.uq.edu.au/index.php/ESyS-Particle\n" << endl;
    }
    
    /*!
      constructor for RotBondPrmsPy
    
      \param name
      \param normalK
      \param shearK
      \param torsionK
      \param bendingK
      \param normalBrkForce  
      \param shearBrkForce
      \param torsionBrkForce
      \param bendingBrkForce
      \param aTag
      \param scaling scaling of normal stiffness with particle size
    */
    RotBondPrmsPy::RotBondPrmsPy(
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
    )
      : CRotBondedIGP(
          name,normalK,shearK,torsionK,bendingK,normalBrkForce,shearBrkForce,
          torsionBrkForce,bendingBrkForce,aTag,scaling,true,1.0,1.0,1.0)
    {
    }

    /*!
      constructor for RotBondPrmsPy
    
      \param name
      \param normalK
      \param shearK
      \param torsionK
      \param bendingK
      \param normalBrkForce  
      \param shearBrkForce
      \param torsionBrkForce
      \param bendingBrkForce
      \param aTag
      \param scaling scaling of normal stiffness with particle size
      \param meanR_scaling scales stiffness using mean of two particles sizes instead of the minimum particle size
    */
    RotBondPrmsPy::RotBondPrmsPy(
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
    )
      : CRotBondedIGP(
          name,normalK,shearK,torsionK,bendingK,normalBrkForce,shearBrkForce,
          torsionBrkForce,bendingBrkForce,aTag,scaling,meanR_scaling,1.0,1.0,1.0)
    {
    }

    /*!
      constructor for RotBondPrmsPy
    
      \param name
      \param normalK
      \param shearK
      \param torsionK
      \param bendingK
      \param normalBrkForce  
      \param shearBrkForce
      \param torsionBrkForce
      \param bendingBrkForce
      \param aTag
      \param scaling scaling of normal stiffness with particle size
      \param meanR_scaling scales stiffness using mean of two particles sizes instead of the minimum particle size
      \param truncated a factor determining amount by which tensile strength is truncated in the failure criterion.
      \param beta1 a factor suppressing bending in failure criterion
      \param beta2 a factor suppressing torsion in failure criterion
    */
    RotBondPrmsPy::RotBondPrmsPy(
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
    )
      : CRotBondedIGP(
          name,normalK,shearK,torsionK,bendingK,normalBrkForce,shearBrkForce,
          torsionBrkForce,bendingBrkForce,aTag,scaling,meanR_scaling,truncated,1.0,1.0)
    {
    }

    /*!
      constructor for RotThermBondPrmsPy
    
      \param name
      \param normalK
      \param shearK
      \param torsionK
      \param bendingK
      \param normalBrkForce  
      \param shearBrkForce
      \param torsionBrkForce
      \param bendingBrkForce
      \param diffusivity
      \param aTag
    */
    RotThermBondPrmsPy::RotThermBondPrmsPy(
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
      int aTag
    )
      : CRotThermBondedIGP(
          name,normalK,shearK,torsionK,bendingK,
          normalBrkForce,shearBrkForce,torsionBrkForce,bendingBrkForce,
          diffusivity,aTag
        )
    {
      std::cerr << " -- WARNING : EXPERIMENTAL FEATURE --- " << std::endl;
      std::cerr << " RotThermBond interactions are not fully tested yet " << std::endl;
      std::cerr << " Use at your own risk  " << std::endl;
    }

    /*!
      constructor for BrittleBeamPrmsPy
    
      \param name
      \param youngsModulus
      \param poissonsRatio
      \param cohesion
      \param tanAngle
      \param aTag
    */
    BrittleBeamPrmsPy::BrittleBeamPrmsPy(
      const  std::string &name,
      double youngsModulus,
      double poissonsRatio,
      double cohesion,
      double tanAngle,
      int    aTag
    )
      : CRotBondedIGP(
          name,youngsModulus,poissonsRatio,cohesion,tanAngle,aTag,true,1.0,1.0,1.0)
    {
    }

    /*!
      constructor for BrittleBeamPrmsPy
    
      \param name
      \param youngsModulus
      \param poissonsRatio
      \param cohesion
      \param tanAngle
      \param aTag
      \param meanR_scaling
    */
    BrittleBeamPrmsPy::BrittleBeamPrmsPy(
      const  std::string &name,
      double youngsModulus,
      double poissonsRatio,
      double cohesion,
      double tanAngle,
      int    aTag,
      bool   meanR_scaling
    )
      : CRotBondedIGP(
          name,youngsModulus,poissonsRatio,cohesion,tanAngle,aTag,meanR_scaling,1.0,1.0,1.0)
    {
    }

    /*!
      constructor for BrittleBeamPrmsPy
    
      \param name
      \param youngsModulus
      \param poissonsRatio
      \param cohesion
      \param tanAngle
      \param aTag
      \param meanR_scaling
      \param truncated
      \param beta1
      \param beta2
    */
    BrittleBeamPrmsPy::BrittleBeamPrmsPy(
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
    )
      : CRotBondedIGP(
          name,youngsModulus,poissonsRatio,cohesion,tanAngle,aTag,meanR_scaling,truncated,beta1,beta2)
    {
    }

    /*!
      constructor for BrittleBeamSCPrmsPy
    
      \param name
      \param youngsModulus
      \param poissonsRatio
      \param cohesion
      \param frictionAngle
      \param tensileCutoff
      \param aTag
    */
    BrittleBeamSCPrmsPy::BrittleBeamSCPrmsPy(
        const  std::string &name,
        double youngsModulus,
        double poissonsRatio,
        double cohesion,
        double frictionAngle,
        double tensileCutoff,
        int aTag ) : 
	BrittleBeamSCIGP(name, youngsModulus, poissonsRatio, cohesion, frictionAngle, tensileCutoff, aTag)
    {
        std::cerr << " -- WARNING : EXPERIMENTAL FEATURE --- " << std::endl;
        std::cerr << "  BrittleBeamSC interactions are not fully tested yet " << std::endl;
        std::cerr << " Use at your own risk  " << std::endl;    
    }
    
    // --- END ROTBONDED ---
    
   /*!
      constructor for BrittleBeamDZCPrmsPy
      
      \param name
      \param youngsModulus
      \param poissonsRatio
      \param cohesion
      \param frictionAngle
      \param tensileCutoff
      \param compressCutoff
      \param beta1
      \param beta2
      \param aTag
    */
    BrittleBeamDZCPrmsPy::BrittleBeamDZCPrmsPy(
        const  std::string &name,
        double youngsModulus,
        double poissonsRatio,
        double cohesion,
        double frictionAngle,
        double tensileCutoff,
        double compressCutoff,
        double beta1,
        double beta2,
        int aTag ) :
    BrittleBeamDZCIGP(name, youngsModulus, poissonsRatio, cohesion, frictionAngle, tensileCutoff, compressCutoff, beta1, beta2, aTag)
    {
        std::cerr << " -- WARNING : EXPERIMENTAL FEATURE --- " << std::endl;
        std::cerr << "  BrittleBeamDZC interactions are not fully tested yet " << std::endl;
        std::cerr << " Use at your own risk  " << std::endl;
    }
    
    // --- END ROTBONDED ---
    
    FrictionPrmsPy::FrictionPrmsPy(
      const  std::string &name,
      double youngsModulus,
      double poissonsRatio,
      double dynamicMu,    // sliding frictional coefficient
      double staticMu,     // max static frictional coefficient
      bool   rigid,
      bool   meanR_scaling
    ) : CRotFrictionIGP(name, youngsModulus, poissonsRatio, dynamicMu, staticMu, 0.0, rigid, meanR_scaling)
    {
    }

    FrictionPrmsPy::FrictionPrmsPy(
      const  std::string &name,
      double youngsModulus,
      double poissonsRatio,
      double dynamicMu,    // sliding frictional coefficient
      double staticMu,     // max static frictional coefficient
      bool   rigid
    ) : CRotFrictionIGP(name, youngsModulus, poissonsRatio, dynamicMu, staticMu, 0.0, rigid, true)
    {
    }

    FrictionPrmsPy::FrictionPrmsPy(
      const  std::string &name,
      double youngsModulus,
      double poissonsRatio,
      double dynamicMu,    // sliding frictional coefficient
      double staticMu     // max static frictional coefficient
    ) : CRotFrictionIGP(name, youngsModulus, poissonsRatio, dynamicMu, staticMu, 0.0, false, true)
    {
    }
    
    /*!
      Constructor for RotFrictionPrmsPy
      \param name
      \param normalK
      \param dynamicMu
      \param staticMu
      \param shearK
    */
    RotFrictionPrmsPy::RotFrictionPrmsPy(
      const std::string &name,
      double normalK,
      double dynamicMu,    // sliding frictional coefficient
      double staticMu,     // max static frictional coefficient
      double shearK
    ) : CRotFrictionIGP(name, normalK, dynamicMu, staticMu, shearK, 0.0, true, false, true)
    {
      std::cerr << "\n--- WARNING ---\n";
      std::cerr << "ESyS-Particle 2.0 by default scales elastic stiffness";
      std::cerr << " according to particle dimensions.\n";
      std::cerr << "To disable scaling set \"scaling=False\" in the";
      std::cerr << " RotFrictionPrms argument list.\n";
      std::cerr << "To remove this warning set \"scaling=True\".\n";
      std::cerr << "For more information about scaling";
      std::cerr << " refer to the ESyS-Particle Tutorial available at:\n";
      std::cerr << "https://wiki.geocomp.uq.edu.au/index.php/ESyS-Particle\n" << endl;
    }
    
   /*!
      Constructor for RotFrictionPrmsPy
      \param name
      \param normalK
      \param dynamicMu
      \param staticMu
      \param shearK
      \param scaling scaling of normal stiffness with particle size
    */
    RotFrictionPrmsPy::RotFrictionPrmsPy(
      const std::string &name,
      double normalK,
      double dynamicMu,    // sliding frictional coefficient
      double staticMu,     // max static frictional coefficient
      double shearK,
      bool   scaling
    ) : CRotFrictionIGP(name, normalK, dynamicMu, staticMu, shearK, 0.0, scaling, false, true)
    {
    }

    RotFrictionPrmsPy::RotFrictionPrmsPy(
      const std::string &name,
      double normalK,
      double dynamicMu,    // sliding frictional coefficient
      double staticMu,     // max static frictional coefficient
      double shearK,
      bool   scaling,
      bool   rigid
    ) : CRotFrictionIGP(name, normalK, dynamicMu, staticMu, shearK, 0.0, scaling, rigid, true)
    {
    }

    RotFrictionPrmsPy::RotFrictionPrmsPy(
      const std::string &name,
      double normalK,
      double dynamicMu,    // sliding frictional coefficient
      double staticMu,     // max static frictional coefficient
      double shearK,
      bool   scaling,
      bool   rigid,
      bool   meanR_scaling
    ) : CRotFrictionIGP(name, normalK, dynamicMu, staticMu, shearK, 0.0, scaling, rigid, meanR_scaling)
    {
    }
    
   /*!
      Constructor for RotThermFrictionPrmsPy
      \param name
      \param normalK
      \param dynamicMu
      \param staticMu
      \param shearK
      \param diffusivity
    */
    RotThermFrictionPrmsPy::RotThermFrictionPrmsPy(
      const std::string &name,
      double normalK,
      double dynamicMu,    // sliding frictional coefficient
      double staticMu,     // max static frictional coefficient
      double shearK,
      double diffusivity
    ) : CRotThermFrictionIGP(
          name, normalK, dynamicMu, staticMu, shearK, diffusivity, 0.0
        )
    {
      std::cerr << " -- WARNING : EXPERIMENTAL FEATURE --- " << std::endl;
      std::cerr << " RotThermFriction interactions are not fully tested yet " << std::endl;
      std::cerr << " Use at your own risk  " << std::endl;
    }
    
    /*!
      Constructor for RotElasticPrmsPy
      \param name
      \param normalK
    */
    RotElasticPrmsPy::RotElasticPrmsPy(
      const std::string &name,
      double normalK
    ) :
      CRotElasticIGP(name, normalK, true)
    {
      std::cerr << "\n--- WARNING ---\n";
      std::cerr << "ESyS-Particle 2.0 by default scales elastic stiffness";
      std::cerr << " according to particle dimensions.\n";
      std::cerr << "To disable scaling set \"scaling=False\" in the";
      std::cerr << " RotFrictionPrms argument list.\n";
      std::cerr << "To remove this warning set \"scaling=True\".\n";
      std::cerr << "For more information about scaling";
      std::cerr << " refer to the ESyS-Particle Tutorial available at:\n";
      std::cerr << "https://wiki.geocomp.uq.edu.au/index.php/ESyS-Particle\n" << endl;
    }
    
   /*!
      Constructor for RotElasticPrmsPy
      \param name
      \param normalK
      \param scaling scaling of normal stiffness with particle size
    */
    RotElasticPrmsPy::RotElasticPrmsPy(
      const std::string &name,
      double normalK,
      bool   scaling
    ) : CRotElasticIGP(name, normalK, scaling)
    {
    }

    /*!
      constructor for RotThermElasticPrmsPy
    
      \param name
      \param normalK
      \param diffusivity
   */
    RotThermElasticPrmsPy::RotThermElasticPrmsPy(
      const std::string &name,
      double normalK,
      double diffusivity
    )
      : CRotThermElasticIGP(name, normalK, diffusivity)
    {  
      std::cerr << " -- WARNING : EXPERIMENTAL FEATURE --- " << std::endl;
      std::cerr << " RotThermElastic interactions are not fully tested yet " << std::endl;
      std::cerr << " Use at your own risk  " << std::endl;
    }

    /*!
      constructor for GravityPrmsPy

      \param name
      \param accel the graviational acceleration
    */
    GravityPrmsPy::GravityPrmsPy(const std::string &name,const Vec3Py& accel) 
      : GravityIGP(name, accel)
    {
    }

    /*!
      constructor for BuoyancyPrmsPy

      \param name
      \param accel the graviational acceleration
    */
    BuoyancyPrmsPy::BuoyancyPrmsPy(const std::string &name,const Vec3Py& accel, const double& fluidDensity, const double& fluidHeight)
      : BuoyancyIGP(name, accel, fluidDensity, fluidHeight)
    {
    }

    /*!
      constructor for VWFrictionPrmsPy
    
      \param name
      \param normalK
      \param dynamicMu
      \param shearK
      \param alpha
    */
    VWFrictionPrmsPy::VWFrictionPrmsPy(
      const std::string &name,
      double normalK,
      double dynamicMu,
      double shearK,
      double alpha
    )
      : VWFrictionIGP(name, normalK, dynamicMu, shearK, 0.0, alpha)
    {
    }


    using boost::python::arg;
    
    /*!
      export the interfaces to Python via boost
    */
    void exportInteractionPrms()
    {
      // Disable autogeneration of C++ signatures (Boost 1.34.0 and higher)
      // for Epydoc which stumbles over indentation in the automatically generated strings.
      boost::python::docstring_options no_autogen(true,false);

      boost::python::class_<InteractionPrmsPy>(
        "InteractionPrms",
        "Base class for interaction parameters.",
        boost::python::init<>()
      )
      ;

      boost::python::class_<DampingPrmsPy,boost::python::bases<InteractionPrmsPy> >(
        "DampingPrms",
        "Viscous damping parameters.",
        boost::python::init<const std::string&,const std::string&,double,int>(
          (
            arg("type"),
            arg("name"),
            arg("viscosity"),
            arg("maxIterations")=100
          ),
          "Parameters defining damping of individual particle motions or rotations\n"
          "@type type: string\n"
          "@kwarg type: Type of damping: 'Damping', 'RotDamping'\n"
          "@type name: string\n"
          "@kwarg name: Name of the interaction\n"
          "@type viscosity: float\n"
          "@kwarg viscosity: Viscosity coefficient for damping\n"
          "@type maxIterations: int\n"
          "@kwarg maxIterations: maximum number of viscosity iterations\n"
        )
      );
      boost::python::class_<LinDampingPrmsPy,boost::python::bases<DampingPrmsPy> >(
        "LinDampingPrms",
        "Linear velocity damping parameters.",
        boost::python::init<const std::string&,double,int>(
          (
            arg("name"),
            arg("viscosity"),
            arg("maxIterations")=100
          ),
          "Parameters defining damping of individual particle motions\n"
          "@type name: string\n"
          "@kwarg name: Name of the interaction\n"
          "@type viscosity: float\n"
          "@kwarg viscosity: Viscosity coefficient for damping\n"
          "@type maxIterations: int\n"
          "@kwarg maxIterations: maximum number of viscosity iterations\n"
        )
      );

      boost::python::class_<RotDampingPrmsPy,boost::python::bases<DampingPrmsPy> >(
        "RotDampingPrms",
        "Rotational/angular velocity damping parameters.",
        boost::python::init<const std::string&,double,int>(
          (
            arg("name"),
            arg("viscosity"),
            arg("maxIterations")=100
          ),
          "Parameters defining damping of individual particle rotations\n"
          "@type name: string\n"
          "@kwarg name: Name of the interaction\n"
          "@type viscosity: float\n"
          "@kwarg viscosity: Viscosity coefficient for damping\n"
          "@type maxIterations: int\n"
          "@kwarg maxIterations: maximum number of viscosity iterations\n"
        )
      );

      boost::python::class_<LocalDampingPrmsPy,boost::python::bases<InteractionPrmsPy> >(
        "LocalDampingPrms",
        "Local damping parameters for translational hysteretic damping.",
        boost::python::init<const std::string&,double>(
          (
            arg("name"),
            arg("viscosity")
          ),
          "Parameters defining hysteretic damping of individual particle motions\n"
          "@type name: string\n"
          "@kwarg name: Name of the interaction\n"
          "@type viscosity: float\n"
          "@kwarg viscosity: Damping coefficient\n"
        )
       )
       ;

      boost::python::class_<RotLocalDampingPrmsPy,boost::python::bases<InteractionPrmsPy> >(
        "RotLocalDampingPrms",
        "Local damping parameters for rotational hysteretic damping.",
        boost::python::init<const std::string&,double>(
          (
            arg("name"),
            arg("viscosity")
          ),
          "Parameters defining hysteretic damping of individual particle rotations\n"
          "@type name: string\n"
          "@kwarg name: Name of the interaction\n"
          "@type viscosity: float\n"
          "@kwarg viscosity: Damping coefficient\n"
        )
       )
       ;
      
      boost::python::class_<ABCDampingPrmsPy,boost::python::bases<InteractionPrmsPy> >(
        "ABCDampingPrms",
        "Class defining Viscous damping parameters for absorbing boundary conditions.",
        boost::python::init<const std::string&,const std::string&,double,int,const Vec3Py&,const Vec3Py,const Vec3Py&,double>(
          (
            arg("type"),
            arg("name"),
            arg("viscosity"),
            arg("maxIterations"),
            arg("Vref"),
            arg("pos"),
            arg("normal"),
            arg("c1")
          ),
          "Parameters defining absorbing boundary conditions.\n"
          "@type type: string\n"
          "@kwarg type: Type of damping: 'Damping', 'RotDamping'\n"
          "@type name: string\n"
          "@kwarg name: Name of the interaction\n"
          "@type viscosity: float\n"
          "@kwarg viscosity: Viscosity coefficient for damping\n"
          "@type maxIterations: int\n"
          "@kwarg maxIterations: maximum number of viscosity iterations\n"
          "@type Vref: vec3\n"
          "@kwarg Vref: reference velocity\n"
          "@type pos: vec3\n"
          "@kwarg pos: position of origin of absorbing boundary plane\n"
          "@type normal: vec3\n"
          "@kwarg normal: normal vector of absorbing boundary plane\n"
          "@type c1: float\n"
          "@kwarg c1: exponent defining rate of decay of damping with distance from absorbing boundary plane\n"
        )
       )
       ;
      boost::python::class_<NRotBondPrmsPy,boost::python::bases<InteractionPrmsPy> >(
        "NRotBondPrms",
        "Parameters for linear elastic bonded interactions with a specified breaking distance (or strain).",
        boost::python::init<const std::string &, double, double, int, bool>(
          (
            arg("name"),
            arg("normalK"),
            arg("breakDistance"),
            arg("tag"),
            arg("scaling")
          ),
          "Parameters for breakable linear elastic bonded interactions.\n"
          "@type name: string\n"
          "@kwarg name: Name assigned to this group of interactions.\n"
          "@type normalK: float\n"
          "@kwarg normalK: spring constant used when calculating linear"
          " elastic normal force.\n"
          "@type breakDistance: float\n"
          "@kwarg breakDistance: When particles are separated by this distance"
          " the bond breaks.\n"
          "@type tag: int\n"
          "@kwarg tag: Connections which are tagged with C{tag}"
          " will be created with these parameters.\n"
          "@type scaling: bool\n"
          "@kwarg scaling: When True (default), normal stiffness is scaled"
          " with with particle size.\n"
        )
      )
      .def(boost::python::init<const std::string &,double,double,int>(
        (
          arg("name"),
          arg("normalK"),
          arg("breakDistance"),
          arg("tag")
        )
      ))
      .def(
        "getName",
        &NRotBondPrmsPy::getName,
        boost::python::return_value_policy<boost::python::copy_const_reference>(),
        "@rtype: string\n"
        "@return: Name assigned to this group of interactions.\n"
      )
      ;

      boost::python::class_<CappedNRotBondPrmsPy,boost::python::bases<InteractionPrmsPy> >(
        "CappedNRotBondPrms",
        "Parameters defining linear elastic bonded interactions with a capped separation distance.",
        boost::python::init<const std::string &, double, double,double, int>(
          (
            arg("name"),
            arg("normalK"),
            arg("breakDistance"),
            arg("maxForce"),
            arg("tag")
          ),
          "Parameters for capped linear elastic bonded interactions.\n"
          "@type name: string\n"
          "@kwarg name: Name assigned to this group of interactions.\n"
          "@type normalK: float\n"
          "@kwarg normalK: spring constant used when calculating linear"
          " elastic normal force.\n"
          "@type breakDistance: float\n"
          "@kwarg breakDistance: When particles are separated by this distance"
          " the bond breaks.\n"
          "@type tag: int\n"
          "@kwarg tag: Connections which are tagged with C{tag}"
          " will be created with these parameters.\n"
        )
      )
      .def(
        "getName",
        &CappedNRotBondPrmsPy::getName,
        boost::python::return_value_policy<boost::python::copy_const_reference>(),
        "@rtype: string\n"
        "@return: Name assigned to this group of interactions.\n"
      )
      ;

      boost::python::class_<NRotShortBondPrmsPy,boost::python::bases<InteractionPrmsPy> >(
        "NRotShortBondPrms",
        "Parameters for short linear elastic bond interactions. A bonded interaction where the equilibrium distance is determined by the initial separation of the particles, thus permitting overlapping particles.",
        boost::python::init<const std::string &, double, double, int>(
          (
            arg("name"),
            arg("normalK"),
            arg("breakDistance"),
            arg("tag")
          ),
          "Parameters for short linear elastic bonded interactions.\n"
          "@type name: string\n"
          "@kwarg name: Name assigned to this group of interactions.\n"
          "@type normalK: float\n"
          "@kwarg normalK: spring constant used when calculating linear"
          " elastic normal force.\n"
          "@type breakDistance: float\n"
          "@kwarg breakDistance: When particles are separated by this distance"
          " the bond breaks.\n"
          "@type tag: int\n"
          "@kwarg tag: Connections which are tagged with C{tag}"
          " will be created with these parameters.\n"
        )
      )
      .def(
        "getName",
        &NRotShortBondPrmsPy::getName,
        boost::python::return_value_policy<boost::python::copy_const_reference>(),
        "@rtype: string\n"
        "@return: Name assigned to this group of interactions.\n"
      )
      ;

      boost::python::class_<NRotElasticPrmsPy,boost::python::bases<InteractionPrmsPy> >(
        "NRotElasticPrms",
        "Parameters for linear elastic contact interactions.",
        boost::python::init<const std::string &, double, bool>(
          (
            arg("name"),
            arg("normalK"),
            arg("scaling")
          ),
          "Parameters for linear elastic contact interactions.\n"
          "@type name: string\n"
          "@kwarg name: Name assigned to this group of interactions.\n"
          "@type normalK: float\n"
          "@kwarg normalK: spring constant used when calculating linear"
          " elastic normal force.\n"
          "@type scaling: bool\n"
          "@kwarg scaling: When True (default), normal stiffness is scaled"
          " with particle size.\n"
        )
      )
      .def(boost::python::init<const std::string &, double>(
        (
          arg("name"),
          arg("normalK")
        )
      ))
      .def(
        "getName",
        &NRotElasticPrmsPy::getName,
        boost::python::return_value_policy<boost::python::copy_const_reference>()
      )
      ;

      boost::python::class_<HertzianElasticPrmsPy,boost::python::bases<InteractionPrmsPy> >(
        "HertzianElasticPrms",
        "Parameters for Hertzian elastic contact interactions.",
        boost::python::init<const std::string &, double, double>(
          (
            arg("name"),
            arg("E"),
            arg("nu")
          ),
          "Parameters for Hertzian elastic contact interactions.\n"
          "@type name: string\n"
          "@kwarg name: Name assigned to this group of interactions.\n"
          "@type E: float\n"
          "@kwarg E: Young's modulus used for force calculation.\n"
          "@type nu: float\n"
          "@kwarg nu: poisson ratio used for force calculation.\n"
        )
      )
      .def(
        "getName",
        &HertzianElasticPrmsPy::getName,
        boost::python::return_value_policy<boost::python::copy_const_reference>()
      )
      ;

      boost::python::class_<
        HertzianViscoElasticFrictionPrmsPy,
        boost::python::bases<InteractionPrmsPy>                  
      >(
        "HertzianViscoElasticFrictionPrms",
        "Parameters for Hertzian viscoelastic contact interactions with friction.",
        boost::python::init<
          const std::string &,
          double,             
          double,
          double,
          double,
          double
        >(
          (
            arg("name"),
            arg("A"),
            arg("E"),
            arg("nu"),
            arg("dynamicMu"),
            arg("shearK")
          ),
          "Parameters for Hertzian viscoelastic and frictional interactions.\n"
          "@type name: string\n"
          "@kwarg name: Name assigned to this group of interactions.\n"
          "@type A: float\n"
          "@kwarg A: Damping constant used for force calculation.\n"
          "@type E: float\n"
          "@kwarg E: Young's modulus used for force calculation.\n"
          "@type nu: float\n"
          "@kwarg nu: poisson ratio used for force calculation.\n"
          "@type dynamicMu: float\n"
          "@kwarg dynamicMu: friction coefficient"
          " when contact is sliding.\n"
          "@type shearK: float\n"
          "@kwarg shearK: spring constant used when calculating linear"
          " viscoelastic shear force.\n"
        )
      )
      .def(
        "getName",
        &HertzianViscoElasticFrictionPrmsPy::getName,
        boost::python::return_value_policy<boost::python::copy_const_reference>()
      )
      ;

      boost::python::class_<
        HertzianViscoElasticPrmsPy,
        boost::python::bases<InteractionPrmsPy> 
      >(
        "HertzianViscoElasticPrms",
        "Parameters for Hertzian viscoelastic contact interactions.",
        boost::python::init<const std::string &, double, double, double>(
          (
            arg("name"),
            arg("A"),
            arg("E"),
            arg("nu")
          ),
          "Parameters for Hertzian viscoelastic contact interactions.\n"
          "@type name: string\n"
          "@kwarg name: Name assigned to this group of interactions.\n"
          "@type A: float\n"
          "@kwarg A: Damping constant used for force calculation.\n"
          "@type E: float\n"
          "@kwarg E: Young's modulus used for force calculation.\n"
          "@type nu: float\n"
          "@kwarg nu: poisson ratio used for force calculation.\n"
        )
      )
      .def(
        "getName",
        &HertzianViscoElasticPrmsPy::getName,
        boost::python::return_value_policy<boost::python::copy_const_reference>()
      )
      ;

      boost::python::class_<
        HertzMindlinPrmsPy,
        boost::python::bases<InteractionPrmsPy>
      >(
        "HertzMindlinPrms",
        "Parameters for HertzMindlin interactions.",
        boost::python::init<
          const std::string &,
          double,
          double,
          double
        >(
          (
            arg("name"),
            arg("youngsModulus"),
            arg("poissonsRatio"),
            arg("dynamicMu")
          ),
          "Parameters for HertzMindlin interactions.\n"
          "@type name: string\n"
          "@kwarg name: Name assigned to this group of interactions.\n"
          "@type youngsModulus: float\n"
          "@kwarg youngsModulus: Young's modulus used for force calculation.\n"
          "@type poissonsRatio: float\n"
          "@kwarg poissonsRatio: Poisson's ratio used for force calculation.\n"
          "@type dynamicMu: float\n"
          "@kwarg dynamicMu: friction coefficient"
          " when contact is sliding.\n"
        )
      )
      .def(
        "getName",
        &HertzMindlinPrmsPy::getName,
        boost::python::return_value_policy<boost::python::copy_const_reference>()
      )
      ;

      boost::python::class_<
        HertzMindlinViscoPrmsPy,
        boost::python::bases<InteractionPrmsPy>
      >(
        "HertzMindlinViscoPrms",
        "Parameters for HertzMindlin interactions with damping.",
        boost::python::init<
          const std::string &,
          double,
          double,
          double,
          double
        >(
          (
            arg("name"),
            arg("youngsModulus"),
            arg("poissonsRatio"),
            arg("dynamicMu"),
            arg("restitution")
          ),
          "Parameters for HertzMindlin interactions with damping.\n"
          "@type name: string\n"
          "@kwarg name: Name assigned to this group of interactions.\n"
          "@type youngsModulus: float\n"
          "@kwarg youngsModulus: Young's modulus used for force calculation.\n"
          "@type poissonsRatio: float\n"
          "@kwarg poissonsRatio: Poisson's ratio used for force calculation.\n"
          "@type dynamicMu: float\n"
          "@kwarg dynamicMu: friction coefficient"
          " when contact is sliding.\n"
          "@type restitution: float\n"
          "@kwarg restitution: restitution coefficient.\n"
        )
      )
      .def(
        "getName",
        &HertzMindlinViscoPrmsPy::getName,
        boost::python::return_value_policy<boost::python::copy_const_reference>()
      )
      ;

      boost::python::class_<
        LinearDashpotPrmsPy,
        boost::python::bases<InteractionPrmsPy>
      >(
        "LinearDashpotPrms",
        "Parameters for linear dashpot interactions. This interaction group can be used in parallel with C{NRotElasticPrms} to define spring-dashpot interactions.",
        boost::python::init<const std::string &, double, double>(
          (
            arg("name"),
            arg("damp"),
            arg("cutoff")
          ),
          "Define parameters for linear dashpot interactions.\n"
          "@type name: string\n"
          "@kwarg name: Name assigned to this group of interactions.\n"
          "@type damp: float\n"
          "@kwarg damp: damping coefficient (viscosity) used for force"
          " calculation.\n"
          "@type cutoff: float\n"
          "@kwarg cutoff: interaction range relative to particle radii.\n"
        )
      )
      .def(
        "getName",
        &LinearDashpotPrmsPy::getName,
        boost::python::return_value_policy<boost::python::copy_const_reference>()
      )
      ;

      boost::python::class_<NRotFrictionPrmsPy,boost::python::bases<InteractionPrmsPy> >(
        "NRotFrictionPrms",
        "Parameters for non-rotational friction, a Coulomb frictional"
        " law with shear-stiffness. Forces are applied at particle centres,"
        " not at the contact point.",
        boost::python::init<const std::string &,double,double,double,bool>(
          (
            arg("name"),
            arg("normalK"),
            arg("dynamicMu"),
            arg("shearK"),
            arg("scaling")
          ),
          "Define parameters for linear elastic frictional contact interactions.\n"
          "@type name: string\n"
          "@kwarg name: Name assigned to this group of interactions.\n"
          "@type normalK: float\n"
          "@kwarg normalK: spring constant used when calculating linear"
          " elastic normal force.\n"
          "@type dynamicMu: float\n"
          "@kwarg dynamicMu: friction coefficient"
          " when contact is sliding.\n"
          "@type shearK: float\n"
          "@kwarg shearK: spring constant used when calculating linear"
          " elastic shear force.\n"
          "@type scaling: bool\n"
          "@kwarg scaling: When True (default), normal stiffness is scaled"
          " with particle size.\n"
        )
      )
      .def(boost::python::init<const std::string &,double,double,double>(
        (
          arg("name"),
          arg("normalK"),
          arg("dynamicMu"),
          arg("shearK")
        )
      ))
      .def(
        "getName",
        &NRotFrictionPrmsPy::getName,
        boost::python::return_value_policy<boost::python::copy_const_reference>(),
        "@rtype: string\n"
        "@return: Name assigned to this group of interactions."
      )
      ;

      boost::python::class_<SpringDashpotFrictionPrmsPy,boost::python::bases<InteractionPrmsPy> >(
        "SpringDashpotFrictionPrms",
        "Parameters for non-rotational spring-dashpot friction, a Coulomb frictional"
        " law with shear-stiffness. Forces are applied at particle centres,"
        " not at the contact point.",
        boost::python::init<const std::string &,double,double,double,double>(
          (
            arg("name"),
            arg("youngsModulus"),
            arg("poissonsRatio"),
            arg("restitution"),
            arg("dynamicMu")
          ),
          "Define parameters for linear elastic frictional contact interactions.\n"
          "@type name: string\n"
          "@kwarg name: Name assigned to this group of interactions.\n"
          "@type youngsModulus: float\n"
          "@kwarg youngsModulus: elastic modulus used when calculating linear"
          " elastic normal force.\n"
          "@type poissonsRatio: float\n"
          "@kwarg poissonsRatio: Poisson's ratio used when calculating linear"
          " elastic shear force.\n"
          "@type restitution: float\n"
          "@kwarg restitution: Coefficient of Restitution in normal direction"
          " (0 < restitution < 1). \n"
          "@type dynamicMu: float\n"
          "@kwarg dynamicMu: friction coefficient"
          " when contact is sliding.\n"
        )
      )
      .def(
        "getName",
        &SpringDashpotFrictionPrmsPy::getName,
        boost::python::return_value_policy<boost::python::copy_const_reference>(),
        "@rtype: string\n"
        "@return: Name assigned to this group of interactions."
      )
      ;
      
      boost::python::class_<VWFrictionPrmsPy,boost::python::bases<InteractionPrmsPy> >(
        "VWFrictionPrms",
        "Parameters for velocity weakening friction. Forces are applied at particle centres,"
        " not at the contact point.",
        boost::python::init<
          const std::string &,
          double,
          double,
          double,
          double
        >(
          (
            arg("name"),
            arg("normalK"),
            arg("staticMu"),
            arg("shearK"),
            arg("alpha")
          ),
          "Define parameters for velocity weakening frictional contact interactions.\n"
          "@type name: string\n"
          "@kwarg name: Name assigned to this group of interactions.\n"
          "@type normalK: float\n"
          "@kwarg normalK: spring constant used when calculating linear"
          " elastic normal force.\n"
          "@type dynamicMu: float\n"
          "@kwarg dynamicMu: friction coefficient"
          " when contact is sliding.\n"
          "@type shearK: float\n"
          "@kwarg shearK: spring constant used when calculating linear"
          " elastic shear force.\n"
          "@type alpha: float\n"
          "@kwarg alpha: coefficient for the amount of weakening\n"
        )
      )
      .def(
        "getName",
        &VWFrictionPrmsPy::getName,
        boost::python::return_value_policy<boost::python::copy_const_reference>(),
        "@rtype: string\n"
        "@return: Name assigned to this group of interactions.\n"
      )
      ;

      boost::python::class_<BrittleBeamPrmsPy,boost::python::bases<InteractionPrmsPy> >(
        "BrittleBeamPrms",
        "Parameters for rotational bonded interactions based on elastic beam"
        " theory with a Mohr-Coulomb failure criterion.\n",
        boost::python::init<
          const std::string &,
          double,
          double,
          double,
          double,
          int
        >(
          (
            arg("name"),
            arg("youngsModulus"),
            arg("poissonsRatio"),
            arg("cohesion"),
            arg("tanAngle"),
            arg("tag")
          )
        )
      )
      .def(
        boost::python::init<
          const std::string &,
          double,
          double,
          double,
          double,
          int,
          bool
        >(
          (
            arg("name"),
            arg("youngsModulus"),
            arg("poissonsRatio"),
            arg("cohesion"),
            arg("tanAngle"),
            arg("tag"),
            arg("meanR_scaling")
          )
        )
      )
      .def(
        boost::python::init<
          const std::string &,
          double,
          double,
          double,
          double,
          int,
          bool,
          double,
          double,
          double
        >(
          (
            arg("name"),
            arg("youngsModulus"),
            arg("poissonsRatio"),
            arg("cohesion"),
            arg("tanAngle"),
            arg("tag"),
            arg("meanR_scaling"),
            arg("truncated"),
            arg("beta1"),
            arg("beta2")
          ),
          "Parameters defining elastic-brittle beam interactions between bonded particles.\n"
          "@type name: string\n"
          "@kwarg name: Name assigned to this group of interactions.\n"
          "@type youngsModulus: float\n"
          "@kwarg youngsModulus: Youngs Modulus for bonds (stress units).\n"
          "@type poissonsRatio: float\n"
          "@kwarg poissonsRatio: Poisson's ratio for bonds (dimensionless).\n"
          "@type cohesion: float\n"
          "@kwarg cohesion: Mohr-Coulomb cohesion factor (stress units).\n"
          "@type tanAngle: float\n"
          "@kwarg tanAngle: tan(angle of internal friction) for Mohr-Coulomb failure criterion (dimensionless).\n"
          "@type tag: int\n"
          "@kwarg tag: Connections which are tagged with C{tag}"
          " will be created with these parameters.\n"
          "@type meanR_scaling: bool\n"
          "@kwarg meanR_scaling: use mean particle radius instead of "
          "min particle radius as the radius of bonds\n"
          "@type truncated: double\n"
          "@kwarg truncated: factor by which to truncate tensile strength\n"
          "@type beta1: double\n"
          "@kwarg beta1: factor by which to suppress bending in failure criterion\n"
          "@type beta2: double\n"
          "@kwarg beta2: factor by which to suppress torsion in failure criterion\n"
        )
      )
      .def(
        "getName",
        &BrittleBeamPrmsPy::getName,
        boost::python::return_value_policy<boost::python::copy_const_reference>(),
        "@rtype: string\n"
        "@return: Name assigned to this group of interactions."
      )
      ;

      boost::python::class_<BrittleBeamSCPrmsPy,boost::python::bases<InteractionPrmsPy> >(
        "BrittleBeamSCPrms",
        "Parameters for rotational bonded interactions based on elastic beam"
        " theory with a particle stress based failure criterion.\n",
        boost::python::init<const std::string &,double,double,double,double,double,int
        >(
          (
            arg("name"),
            arg("youngsModulus"),
            arg("poissonsRatio"),
            arg("cohesion"),
            arg("frictionAngle"),
            arg("tensileCutoff"),
            arg("tag")
          )
	  ,
          "Parameters defining elastic-brittle beam interactions with a particle-stress based MC failure criterion between bonded particles.\n"
          "@type name: string\n"
          "@kwarg name: Name assigned to this group of interactions.\n"
          "@type youngsModulus: float\n"
          "@kwarg youngsModulus: Youngs Modulus for bonds (stress units).\n"
          "@type poissonsRatio: float\n"
          "@kwarg poissonsRatio: Poisson's ratio for bonds (dimensionless).\n"
          "@type cohesion: float\n"
          "@kwarg cohesion: particle scale cohesion (stress units).\n"
          "@type frictionAngle: float\n"
          "@kwarg frictionAngle: angle of internal friction for Mohr-Coulomb failure criterion in degrees (dimensionless).\n"
          "@type tensileCutoff: float\n"
          "@kwarg tensileCutoff: tensile cut-off (stress units).\n"
          "@type tag: int\n"
          "@kwarg tag: Connections which are tagged with C{tag}"
          " will be created with these parameters.\n"
        )
      )
             
      .def(
        "getName",
        &BrittleBeamSCPrmsPy::getName,
        boost::python::return_value_policy<boost::python::copy_const_reference>(),
        "@rtype: string\n"
        "@return: Name assigned to this group of interactions."
      )
      ;

      boost::python::class_<BrittleBeamDZCPrmsPy,boost::python::bases<InteractionPrmsPy> >(
        "BrittleBeamDZCPrms",
        "Parameters for rotational bonded interactions based on elastic beam"
        " theory with Ding and Zhang (2014) failure criterion.\n",
        boost::python::init<const std::string &,double,double,double,double,double,double,double,double,int
        >(
          (
            arg("name"),
            arg("youngsModulus"),
            arg("poissonsRatio"),
            arg("cohesion"),
            arg("frictionAngle"),
            arg("tensileCutoff"),
            arg("compressCutoff"),
            arg("beta1"),
            arg("beta2"),
            arg("tag")
          )
      ,
          "Parameters defining elastic-brittle beam interactions with a particle-stress based MC failure criterion between bonded particles.\n"
          "@type name: string\n"
          "@kwarg name: Name assigned to this group of interactions.\n"
          "@type youngsModulus: float\n"
          "@kwarg youngsModulus: Youngs Modulus for bonds (stress units).\n"
          "@type poissonsRatio: float\n"
          "@kwarg poissonsRatio: Poisson's ratio for bonds (dimensionless).\n"
          "@type cohesion: float\n"
          "@kwarg cohesion: particle scale cohesion (stress units).\n"
          "@type frictionAngle: float\n"
          "@kwarg frictionAngle: angle of internal friction for Mohr-Coulomb failure criterion in degrees (dimensionless).\n"
          "@type tensileCutoff: float\n"
          "@kwarg tensileCutoff: tensile cut-off (stress units).\n"
          "@type compressCutoff: float\n"
          "@kwarg compressCutoff: shear stress cut-off under compression (stress units).\n"
          "@type beta1: float\n"
          "@kwarg beta1: supression factor for bending moments (between 0 and 1).\n"
          "@type beta2: float\n"
          "@kwarg beta2: supression factor for torsion moments (between 0 and 1).\n"
          "@type tag: int\n"
          "@kwarg tag: Connections which are tagged with C{tag}"
          " will be created with these parameters.\n"
        )
      )

      .def(
        "getName",
        &BrittleBeamDZCPrmsPy::getName,
        boost::python::return_value_policy<boost::python::copy_const_reference>(),
        "@rtype: string\n"
        "@return: Name assigned to this group of interactions."
      )
      ;
      
      boost::python::class_<FrictionPrmsPy,boost::python::bases<InteractionPrmsPy> >(
        "FrictionPrms",
        "Parameters for rotational friction, a Coulomb frictional"
        " law with shear stiffness. Forces are applied at the contact"
        " point, which generates moments. Normal and shear stiffnesses are"
        " specified via a contact Young's Modulus and Poisson's Ratio.",
        boost::python::init<
          const std::string &,
          double,
          double,
          double,
          double
        >(
          (
            arg("name"),
            arg("youngsModulus"),
            arg("poissonsRatio"),
            arg("dynamicMu"),    // sliding frictional coefficient
            arg("staticMu")     // max static frictional coefficient
          )
        )
      )
      .def(
        boost::python::init<
          const std::string &,
          double,
          double,
          double,
          double,
          bool
        >(
          (
            arg("name"),
            arg("youngsModulus"),
            arg("poissonsRatio"),
            arg("dynamicMu"),    // sliding frictional coefficient
            arg("staticMu"),     // max static frictional coefficient
            arg("rigid")
          )
        )
      )
      .def(
        boost::python::init<
          const std::string &,
          double,
          double,
          double,
          double,
          bool,
          bool
        >(
          (
            arg("name"),
            arg("youngsModulus"),
            arg("poissonsRatio"),
            arg("dynamicMu"),    // sliding frictional coefficient
            arg("staticMu"),    // max static frictional coefficient
            arg("rigid"),
            arg("meanR_scaling")
          ),
          "Parameters defining beam-like frictional interactions between rotational particles.\n"
          "@type name: string\n"
          "@kwarg name: Name assigned to this group of interactions.\n"
          "@type youngsModulus: float\n"
          "@kwarg youngsModulus: elastic contact modulus"
          " used when calculating linear elastic normal force.\n"
          "@type poissonsRatio: float\n"
          "@kwarg poissonsRatio: contact Poisson's Ratio"
          " used when calculating linear elastic shear force.\n"
          "@type dynamicMu: float\n"
          "@kwarg dynamicMu: friction coefficient"
          " when contact is sliding.\n"
          "@type staticMu: float\n"
          "@kwarg staticMu: friction coefficient which governs the"
          " transition from static contact to dynamic frictional sliding.\n"
          "@type rigid: bool\n"
          "@kwarg rigid: When True (default is False), rigid body rotations of"
          "  touching particle-pairs are taken into account.\n"
          "@type meanR_scaling: bool\n"
          "@kwarg meanR_scaling: use mean particle"
          " radius instead of min particle radius as the radius of bonds.\n"
        )
      )
      .def(
        "getName",
        &FrictionPrmsPy::getName,
        boost::python::return_value_policy<
          boost::python::copy_const_reference
        >(),
        "@rtype: string\n"
        "@return: Name assigned to this group of interactions.\n"
      )
      ;

      boost::python::class_<RotElasticPrmsPy,boost::python::bases<InteractionPrmsPy> >(
        "RotElasticPrms",
        "Parameters defining elastic contacts between rotational particles.",
        boost::python::init<
          const std::string &,
          double,
          bool
        >(
          (
            arg("name"),
            arg("normalK"),
            arg("scaling")
          ),
          "Parameters for rotational elastic contact interactions.\n"
          "@type name: string\n"
          "@kwarg name: Name assigned to this group of interactions.\n"
          "@type normalK: float\n"
          "@kwarg normalK: elastic contact modulus if C{scaling} is True"
          " (default) or spring constant if C{scaling} is False, used when"
          " calculating linear elastic normal force.\n"
          "@type scaling: bool\n"
          "@kwarg scaling: When True (default), elastic stiffnesses are scaled"
          " according to particle dimensions.\n"
        )
      )
      .def(
        boost::python::init<
          const std::string &,
          double
        >(
          (
            arg("name"),
            arg("normalK")
          )
        )
      )
      .def(
        "getName",
        &RotElasticPrmsPy::getName,
        boost::python::return_value_policy<
          boost::python::copy_const_reference
        >(),
        "@rtype: string\n"
        "@return: Name assigned to this group of interactions.\n"
      )
      ;

      boost::python::class_<RotFrictionPrmsPy,boost::python::bases<InteractionPrmsPy> >(
        "RotFrictionPrms",
        "Parameters for rotational friction, a Coulomb frictional"
        " law with shear stiffness. Forces are applied at the contact"
        " point, which generates moments.",
        boost::python::init<
          const std::string &,
          double,
          double,
          double,
          double
        >(
          (
            arg("name"),
            arg("normalK"),
            arg("dynamicMu"),    // sliding frictional coefficient
            arg("staticMu"),     // max static frictional coefficient
            arg("shearK")
          )
        )
      )
      .def(
        boost::python::init<
          const std::string &,
          double,
          double,
          double,
          double,
          bool
        >(
          (
            arg("name"),
            arg("normalK"),
            arg("dynamicMu"),    // sliding frictional coefficient
            arg("staticMu"),     // max static frictional coefficient
            arg("shearK"),
            arg("scaling")
          )
        )
      )
      .def(
        boost::python::init<
          const std::string &,
          double,
          double,
          double,
          double,
          bool,
          bool
        >(
          (
            arg("name"),
            arg("normalK"),
            arg("dynamicMu"),    // sliding frictional coefficient
            arg("staticMu"),     // max static frictional coefficient
            arg("shearK"),
            arg("scaling"),
            arg("rigid")
          )
        )
      )
      .def(
        boost::python::init<
          const std::string &,
          double,
          double,
          double,
          double,
          bool,
          bool,
          bool
        >(
          (
            arg("name"),
            arg("normalK"),
            arg("dynamicMu"),    // sliding frictional coefficient
            arg("staticMu"),     // max static frictional coefficient
            arg("shearK"),
            arg("scaling"),
            arg("rigid"),
            arg("meanR_scaling")
          ),
          "Parameters defining frictional interactions between rotational particles.\n"
          "@type name: string\n"
          "@kwarg name: Name assigned to this group of interactions.\n"
          "@type normalK: float\n"
          "@kwarg normalK: elastic contact modulus if C{scaling} is True"
          " (default) or spring constant if C{scaling} is False, used when"
          " calculating linear elastic normal force.\n"
          "@type dynamicMu: float\n"
          "@kwarg dynamicMu: friction coefficient"
          " when contact is sliding.\n"
          "@type staticMu: float\n"
          "@kwarg staticMu: friction coefficient which governs the"
          " transition from static contact to dynamic frictional sliding.\n"
          "@type shearK: float\n"
          "@kwarg shearK: elastic contact modulus if C{scaling} is True"
          " (default) or spring constant if C{scaling} is False, used when"
          " calculating linear elastic shear force.\n"
          "@type scaling: bool\n"
          "@kwarg scaling: When True (default), elastic stiffnesses are scaled"
          " according to particle dimensions.\n"
          "@type rigid: bool\n"
          "@kwarg rigid: When True (default is False), rigid body rotations of"
          " touching particle-pairs are taken into account.\n"
          "@type meanR_scaling: bool\n"
          "@kwarg meanR_scaling: use mean particle radius instead of "
          "min particle radius as the radius of bonds\n"
        )
      )
      .def(
        "getName",
        &RotFrictionPrmsPy::getName,
        boost::python::return_value_policy<
          boost::python::copy_const_reference
        >(),
        "@rtype: string\n"
        "@return: Name assigned to this group of interactions.\n"
      )
      ;

      boost::python::class_<RotBondPrmsPy,boost::python::bases<InteractionPrmsPy> >(
        "RotBondPrms",
        "Parameters for rotational bonded interactions. Parameters include"
        " spring constants for normal, shear, torsion and bending forces.\n",
        boost::python::init<
          const std::string &,
          double,
          double,
          double,
          double,
          double,
          double,
          double,
          double,
          int
        >(
          (
            arg("name"),
            arg("normalK"),
            arg("shearK"),
            arg("torsionK"),
            arg("bendingK"),
            arg("normalBrkForce"),
            arg("shearBrkForce"),
            arg("torsionBrkForce"),
            arg("bendingBrkForce"),
            arg("tag")
          )
        )
      )
      .def(
        boost::python::init<
          const std::string &,
          double,
          double,
          double,
          double,
          double,
          double,
          double,
          double,
          int,
          bool
        >(
          (
            arg("name"),
            arg("normalK"),
            arg("shearK"),
            arg("torsionK"),
            arg("bendingK"),
            arg("normalBrkForce"),
            arg("shearBrkForce"),
            arg("torsionBrkForce"),
            arg("bendingBrkForce"),
            arg("tag"),
            arg("scaling")
          )
        )
      )
      .def(
        boost::python::init<
          const std::string &,
          double,
          double,
          double,
          double,
          double,
          double,
          double,
          double,
          int,
          bool,
          bool
        >(
          (
            arg("name"),
            arg("normalK"),
            arg("shearK"),
            arg("torsionK"),
            arg("bendingK"),
            arg("normalBrkForce"),
            arg("shearBrkForce"),
            arg("torsionBrkForce"),
            arg("bendingBrkForce"),
            arg("tag"),
            arg("scaling"),
            arg("meanR_scaling")
          )
        )
      )
      .def(
        boost::python::init<
          const std::string &,
          double,
          double,
          double,
          double,
          double,
          double,
          double,
          double,
          int,
          bool,
          bool,
          double
        >(
          (
            arg("name"),
            arg("normalK"),
            arg("shearK"),
            arg("torsionK"),
            arg("bendingK"),
            arg("normalBrkForce"),
            arg("shearBrkForce"),
            arg("torsionBrkForce"),
            arg("bendingBrkForce"),
            arg("tag"),
            arg("scaling"),
            arg("meanR_scaling"),
            arg("truncated")
          ),
          "Parameters defining brittle-elastic interactions between rotational particles.\n"
          "@type name: string\n"
          "@kwarg name: Name assigned to this group of interactions.\n"
          "@type normalK: float\n"
          "@kwarg normalK: elastic contact modulus if C{scaling} is True"
          " (default) or spring constant if C{scaling} is False, used when"
          " calculating linear elastic normal force.\n"
          "@type shearK: float\n"
          "@kwarg shearK: elastic contact modulus if C{scaling} is True"
          " (default) or spring constant if C{scaling} is False, used when"
          " calculating linear elastic shear force.\n"
          "@type torsionK: float\n"
          "@kwarg torsionK: elastic contact modulus if C{scaling} is True"
          " (default) or spring constant if C{scaling} is False, used when"
          " calculating elastic torsion force.\n"
          "@type bendingK: float\n"
          "@kwarg bendingK: elastic contact modulus if C{scaling} is True"
          " (default) or spring constant if C{scaling} is False, used when"
          " calculating elastic bending force.\n"
          "@type normalBrkForce: float\n"
          "@kwarg normalBrkForce: A breaking stress if C{scaling} is True"
          " (default) or a breaking force if C{scaling} is False.\n"
          "@type shearBrkForce: float\n"
          "@kwarg shearBrkForce: A breaking stress if C{scaling} is True"
          " (default) or a breaking force if C{scaling} is False.\n"
          "@type torsionBrkForce: float\n"
          "@kwarg torsionBrkForce: A breaking stress if C{scaling} is True"
          " (default) or a breaking force if C{scaling} is False.\n"
          "@type bendingBrkForce: float\n"
          "@kwarg bendingBrkForce: A breaking stress if C{scaling} is True"
          " (default) or a breaking force if C{scaling} is False.\n"
          "@type tag: int\n"
          "@kwarg tag: Connections which are tagged with C{tag}"
          " will be created with these parameters.\n"
          "@type scaling: bool\n"
          "@kwarg scaling: When True (default), elastic stiffnesses and"
          " breaking forces are scaled according to particle dimensions.\n"
          "@type meanR_scaling: bool\n"
          "@kwarg meanR_scaling: use mean particle radius instead of "
          "min particle radius as the radius of bonds\n"
          "@type truncated: double\n"
          "@kwarg truncated: factor by which to truncate tensile strength\n"
        )
      )
      .def(
        "getName",
        &RotBondPrmsPy::getName,
        boost::python::return_value_policy<boost::python::copy_const_reference>(),
        "@rtype: string\n"
        "@return: Name assigned to this group of interactions.\n"
      )
      ;

      boost::python::class_<RotThermElasticPrmsPy,boost::python::bases<InteractionPrmsPy> >(
        "RotThermalElasticPrms",
        "EXPERIMENTAL "
        "Parameters for linear elastic contact interactions"
        " with heat transfer.\n",
        boost::python::init<const std::string &, double, double>(
          (
            arg("name"),
            arg("normalK"),
            arg("diffusivity")
          ),
          "Parameters for thermal-elastic contact interactions.\n"
          "@type name: string\n"
          "@kwarg name: Name assigned to this group of interactions.\n"
          "@type normalK: float\n"
          "@kwarg normalK: spring constant used when calculating linear"
          " elastic normal force.\n"
          "@type diffusivity: float\n"
          "@kwarg diffusivity: Thermal diffusivity.\n"
        )
      )
      .def(
        "getName",
        &RotThermElasticPrmsPy::getName,
        boost::python::return_value_policy<boost::python::copy_const_reference>()
      )
      ;

      boost::python::class_<RotThermFrictionPrmsPy,boost::python::bases<InteractionPrmsPy> >(
        "RotThermalFrictionPrms",
        "EXPERIMENTAL "
        "Parameters for rotational friction, a Coulomb frictional"
        " force law with shear stiffness and heat generation due to friction."
        " Forces are applied at the contact point, which generates moments.\n",
        boost::python::init<
          const std::string &,
          double,
          double,
          double,
          double,
          double
        >(
          (
            arg("name"),
            arg("normalK"),
            arg("dynamicMu"),    // sliding frictional coefficient
            arg("staticMu"),     // max static frictional coefficient
            arg("shearK"),
            arg("diffusivity")
          ),
          "Parameters for thermal-frictional contact interactions.\n"
          "@type name: string\n"
          "@kwarg name: Name assigned to this group of interactions.\n"
          "@type normalK: float\n"
          "@kwarg normalK: spring constant used when calculating linear"
          " elastic normal force.\n"
          "@type dynamicMu: float\n"
          "@kwarg dynamicMu: friction coefficient"
          " when contact is sliding.\n"
          "@type staticMu: float\n"
          "@kwarg staticMu: friction coefficient which governs the"
          " transition from static contact to dynamic frictional sliding.\n"
          "@type shearK: float\n"
          "@kwarg shearK: spring constant used when calculating linear"
          " elastic shear force.\n"
          "@type diffusivity: float\n"
          "@kwarg diffusivity: Thermal diffusivity.\n"
        )
      )
      .def(
        "getName",
        &RotThermFrictionPrmsPy::getName,
        boost::python::return_value_policy<
          boost::python::copy_const_reference
        >(),
        "@rtype: string\n"
        "@return: Name assigned to this group of interactions.\n"
      )
      ;

      boost::python::class_<RotThermBondPrmsPy,boost::python::bases<InteractionPrmsPy> >(
        "RotThermalBondPrms",
        "EXPERIMENTAL "
        "Parameters for rotational, thermal bonded interactions.\n",
        boost::python::init<
          const std::string &,
          double,
          double,
          double,
          double,
          double,
          double,
          double,
          double,
          double,
          int
        >(
          (
            arg("name"),
            arg("normalK"),
            arg("shearK"),
            arg("torsionK"),
            arg("bendingK"),
            arg("normalBrkForce"),
            arg("shearBrkForce"),
            arg("torsionBrkForce"),
            arg("bendingBrkForce"),
            arg("diffusivity"),
            arg("tag")
          ),
          "Parameters for thermal-brittle-elastic bonded interactions.\n"
          "@type name: string\n"
          "@kwarg name: Name assigned to this group of interactions.\n"
          "@type normalK: float\n"
          "@kwarg normalK: spring constant used when calculating linear"
          " elastic normal force.\n"
          "@type shearK: float\n"
          "@kwarg shearK: spring constant used when calculating linear"
          " elastic shear force.\n"
          "@type torsionK: float\n"
          "@kwarg torsionK: spring constant used when calculating"
          " elastic torsion force.\n"
          "@type bendingK: float\n"
          "@kwarg bendingK: spring constant used when calculating"
          " elastic bending force.\n"
          "@type normalBrkForce: float\n"
          "@kwarg normalBrkForce: When the normal force between particles"
          " exceeds this amount, the bond breaks."
          "@type shearBrkForce: float\n"
          "@kwarg shearBrkForce: When the shear force between particles"
          " exceeds this amount, the bond breaks."
          "@type torsionBrkForce: float\n"
          "@kwarg torsionBrkForce: When the torsion force between particles"
          " exceeds this amount, the bond breaks."
          "@type bendingBrkForce: float\n"
          "@kwarg bendingBrkForce: When the bending force between particles"
          " exceeds this amount, the bond breaks."
          "@type diffusivity: float\n"
          "@kwarg diffusivity: Thermal diffusivity.\n"
          "@type tag: int\n"
          "@kwarg tag: Connections which are tagged with C{tag}"
          " will be created with these parameters.\n"
        )
      )
      .def(
        "getName",
        &RotThermBondPrmsPy::getName,
        boost::python::return_value_policy<boost::python::copy_const_reference>(),
        "@rtype: string\n"
        "@return: Name assigned to this group of interactions."
      )
      ;

      boost::python::class_<GravityPrmsPy,boost::python::bases<InteractionPrmsPy> >(
        "GravityPrms",
        "Parameters for describing gravitational body force.",
        boost::python::init<const std::string &, const Vec3Py&>(
          (
            arg("name"),
            arg("acceleration")
          ),
          "Gravitational like body force applied to all particles.\n"
          "@type name: string\n"
          "@kwarg name: name of this interaction.\n"
          "@type acceleration: L{Vec3<esys.lsm.util.FoundationPy.Vec3>}\n"
          "@kwarg acceleration: Acceleration vector (magnitude and direction).\n"
        )
      )
      .def(
        "getName",
        &GravityPrmsPy::getName,
        boost::python::return_value_policy<boost::python::copy_const_reference>(),
        "@rtype: string\n"
        "@return: Name assigned to this group of interactions.\n"
      )
      ;

      boost::python::class_<BuoyancyPrmsPy,boost::python::bases<InteractionPrmsPy> >(
        "BuoyancyPrms",
        "Parameters for describing simple buoyancy body forces.",
        boost::python::init<const std::string &, const Vec3Py&, const double&, const double&>(
          (
            arg("name"),
            arg("acceleration"),
            arg("fluidDensity"),
            arg("fluidHeight")
          ),
          "Simple Buoyancy-like body force applied to all (submerged) particles.\n"
          "@type name: string\n"
          "@kwarg name: name of this interaction.\n"
          "@type acceleration: L{Vec3<esys.lsm.util.FoundationPy.Vec3>}\n"
          "@kwarg acceleration: Acceleration vector (magnitude and direction).\n"
          "@type fluidDensity: float\n"
          "@kwarg fluidDensity: density of the fluid.\n"
          "@type fluidHeight: float\n"
          "@kwarg fluidHeight: height of the fluid in direction of buoyancy force.\n"
        )
      )
      .def(
        "getName",
        &GravityPrmsPy::getName,
        boost::python::return_value_policy<boost::python::copy_const_reference>(),
        "@rtype: string\n"
        "@return: Name assigned to this group of interactions.\n"
      )
      ;
    }
  }
}
