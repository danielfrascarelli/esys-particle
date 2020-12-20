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

#ifndef __BRITTLEBEAMSC_H
#define __BRITTLEBEAMSC_H

// -- project includes --
#include "Model/ARotBondedInteraction.h"
#include "Model/RotParticle.h"
#include "Model/BondedInteractionCpData.h"
#include "Foundation/vec3.h"

// -- I/O includes --
#include <iostream>
using std::ostream;

/*!
  \struct BrittleBeamSCIGP
  \brief Interaction parameters for bonded interaction between rotational particles using the
  average stress failure criterion
  
  \author Steffen Abe, Dion Weatherley
  $Revision$
  $Date$
*/

class BrittleBeamSCIGP : public ARotBondedIGP
{
public:
    BrittleBeamSCIGP();
    BrittleBeamSCIGP(
        const  std::string &name,
        double youngsModulus,
        double poissonsRatio,
        double cohesion,
        double frictionAngle,
        double tensionCutoff,
        int    tag
    );

    BrittleBeamSCIGP(const  std::string &name,
        double kr,
        double ks,
        double kt,
        double kb,
        double cohesion,
        double frictionAngle,
	double tensionCutoff,
        int tag
    );

    virtual std::string getTypeString() const
    {
        return "BrittleBeamSC";
    }

    double cohesion;
    double tCutoff;
    double fAngle;
};

/*!
   \class BrittleBeamSCInteraction
   \brief Elastic interaction between bonded particles between rotational particles
   \author Steffen Abe, Dion Weatherley

    Force calculation is inherited from ARotBondedInteraction. Failure is calculated according
    to a truncated Mohr-Coulomb-criterion based the average stress of the two particles (Wang & Guo 2016) 

   $Revision$
   $Date$
*/
class BrittleBeamSCInteraction : public ARotBondedInteraction
{
    public: // types
        typedef BrittleBeamSCIGP ParameterType;

        /**
        * Used by PIS to save/load check-point data for objects of this type.
        */
        typedef BondedInteractionCpData CheckPointable; // Checkpointing TBD

        typedef double (BrittleBeamSCInteraction::* ScalarFieldFunction)() const; 
        typedef pair<bool,double> (BrittleBeamSCInteraction::* CheckedScalarFieldFunction)() const;
        typedef Vec3 (BrittleBeamSCInteraction::* VectorFieldFunction)() const; 

        // type & dummy implementation for parameter setting function 
        typedef void (BrittleBeamSCInteraction::* ScalarSetFunction)(double);
        static ScalarSetFunction getScalarSetFunction(const string&){return NULL;};

    protected:
        
	// pre-calculated intermediate values for failure criterion
	// cf. Eqs. 2.12ff in Wang & Guo 2016
	double m_sigma3;
	double m_cohesion;
	double m_tanAngle;
       
    public:

        BrittleBeamSCInteraction();
        BrittleBeamSCInteraction(CRotParticle*,CRotParticle*,const BrittleBeamSCIGP&);
        virtual ~BrittleBeamSCInteraction();
                                                                                    
        static ScalarFieldFunction getScalarFieldFunction(const string&);
        static CheckedScalarFieldFunction getCheckedScalarFieldFunction(const string&);
        static VectorFieldFunction getVectorFieldFunction(const string&);

        static string getType(){return "RotBonded";};

        virtual bool broken();
        double getCriterion() const;

        friend class TML_PackedMessageInterface;

        virtual void saveCheckPointData(std::ostream &oStream);
        virtual void loadCheckPointData(std::istream &iStream);

        // save/load of restart parameters
        virtual void saveRestartData(std::ostream &oStream);
        virtual void loadRestartData(std::istream &iStream);
};


#endif // __BRITTLEBEAMSC_H