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

#ifndef __BRITTLEBEAMDZC_H
#define __BRITTLEBEAMDZC_H

// -- project includes --
#include "Model/ARotBondedInteraction.h"
#include "Model/RotParticle.h"
#include "Model/BondedInteractionCpData.h"
#include "Foundation/vec3.h"

// -- I/O includes --
#include <iostream>
using std::ostream;

/*!
  \struct BrittleBeamDZCIGP
  \brief Interaction parameters for bonded interaction between rotational particles using the Ding and Zhang (2014) failure criterion
  
  \author Dion Weatherley
  $Revision$
  $Date$
*/

class BrittleBeamDZCIGP : public ARotBondedIGP
{
public:
    BrittleBeamDZCIGP();
    BrittleBeamDZCIGP(
        const  std::string &name,
        double youngsModulus,
        double poissonsRatio,
        double cohesion,
        double frictionAngle,
        double tensionCutoff,
        double compressCutoff,
        double beta1,
        double beta2,
        int    tag
    );

    BrittleBeamDZCIGP(const  std::string &name,
        double kr,
        double ks,
        double kt,
        double kb,
        double cohesion,
        double frictionAngle,
	    double tensionCutoff,
        double compressCutoff,
        double beta1,
        double beta2,
        int tag
    );

    virtual std::string getTypeString() const
    {
        return "BrittleBeamDZC";
    }

    double cohesion;
    double tCutoff;
    double fAngle;
    double cCutoff;
    double beta1;
    double beta2;
};

/*!
   \class BrittleBeamDZCInteraction
   \brief Elastic interaction between bonded particles between rotational particles
   \author Dion Weatherley

    Force calculation is inherited from ARotBondedInteraction. Failure is calculated according 
    to the modified per-beam failure criterion of Ding and Zhang (2014).

   $Revision$
   $Date$
*/
class BrittleBeamDZCInteraction : public ARotBondedInteraction
{
    public: // types
        typedef BrittleBeamDZCIGP ParameterType;

        /**
        * Used by PIS to save/load check-point data for objects of this type.
        */
        typedef BondedInteractionCpData CheckPointable; // Checkpointing TBD

        typedef double (BrittleBeamDZCInteraction::* ScalarFieldFunction)() const; 
        typedef pair<bool,double> (BrittleBeamDZCInteraction::* CheckedScalarFieldFunction)() const;
        typedef Vec3 (BrittleBeamDZCInteraction::* VectorFieldFunction)() const; 

        // type & dummy implementation for parameter setting function 
        typedef void (BrittleBeamDZCInteraction::* ScalarSetFunction)(double);
        static ScalarSetFunction getScalarSetFunction(const string&){return NULL;};

    protected:
        
	// pre-calculated values for failure criterion
	double m_tCutoff;;
	double m_cCutoff;;
	double m_cohesion;
	double m_tanAngle;
	double m_beta1;
	double m_beta2;
    double m_effR;
    double m_effA;
    double m_effI;
    double m_effJ;
       
    public:

        BrittleBeamDZCInteraction();
        BrittleBeamDZCInteraction(CRotParticle*,CRotParticle*,const BrittleBeamDZCIGP&);
        virtual ~BrittleBeamDZCInteraction();
                                                                                    
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


#endif // __BRITTLEBEAMDZC_H
