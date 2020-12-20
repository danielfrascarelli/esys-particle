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

#ifndef ESYS_LSMBONDEDTRIMESHPRMSPY_H
#define ESYS_LSMBONDEDTRIMESHPRMSPY_H

//--- Project includes ---
#include "Python/esys/lsm/MeshBuildParamsPy.h"
#include "Model/BTriMeshIP.h"

//--- STL includes ---
#include <string>

//--- Boost includes ---
#include <boost/shared_ptr.hpp>

using std::string; 

namespace esys
{
  namespace lsm
  {
    /*!
      \class NRotBondedTriMeshPrmsPy
      \brief class for bonded TriMesh interactions in python interface
    */
    class NRotBondedTriMeshPrmsPy : public BTriMeshIP
    {
    private:
      typedef boost::shared_ptr<MeshTagBuildPrmsPy> TagBuildPrmsPtr;
      typedef boost::shared_ptr<MeshGapBuildPrmsPy> GapBuildPrmsPtr;
      
      const TagBuildPrmsPtr m_tagPrmsPtr;
      const GapBuildPrmsPtr m_gapPrmsPtr;
      
    public:
      NRotBondedTriMeshPrmsPy(
        const string &name,
        const string &meshName,
        double normalK,
        double breakDistance,
        const MeshTagBuildPrmsPy &buildPrms
      );
      NRotBondedTriMeshPrmsPy(
        const string &name,
        const string &meshName,
        double normalK,
        double breakDistance,
        const MeshGapBuildPrmsPy &buildPrms
      );
      bool haveTagBuildPrms() const {return (m_tagPrmsPtr.get() != NULL);};
      bool haveGapBuildPrms() const {return (m_gapPrmsPtr.get() != NULL);};
      const MeshTagBuildPrmsPy &getTagBuildPrms() const {return *(m_tagPrmsPtr.get());};
      const MeshGapBuildPrmsPy &getGapBuildPrms() const {return *(m_gapPrmsPtr.get());};
    };
    
    void exportBondedTriMeshPrms();
  } // namespace lsm
} // namespace esys

#endif //ESYS_LSMBONDEDTRIMESHPRMSPY_H
