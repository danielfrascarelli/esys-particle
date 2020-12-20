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

#ifndef ESYS_LSMBONDEDMESH2DPARAMSPY_H
#define ESYS_LSMBONDEDMESH2DPARAMSPY_H

#include "Python/esys/lsm/MeshBuildParamsPy.h"
#include "Model/BMesh2DIP.h"
#include <string>
#include <boost/shared_ptr.hpp>

namespace esys
{
  namespace lsm
  {
    /*!
      \class NRotBondedLinMeshPrmsPy
      \brief class for bonded Mesh2D  interactions in python interface
    */
    class NRotBondedLinMeshPrmsPy : public BMesh2DIP
    {
    private:
      typedef boost::shared_ptr<MeshTagBuildPrmsPy> TagBuildPrmsPtr;
      typedef boost::shared_ptr<MeshGapBuildPrmsPy> GapBuildPrmsPtr;
      
      const TagBuildPrmsPtr m_tagPrmsPtr;
      const GapBuildPrmsPtr m_gapPrmsPtr;
      
    public: 
      NRotBondedLinMeshPrmsPy(
        const std::string &name,
        const std::string &meshName,
        double normalK,
        double breakDistance,
        const MeshTagBuildPrmsPy &buildPrms
      );
      NRotBondedLinMeshPrmsPy(
        const std::string &name,
        const std::string &meshName,
        double normalK,
        double breakDistance,
        const MeshGapBuildPrmsPy &buildPrms
      );
      bool haveTagBuildPrms() const {return (m_tagPrmsPtr.get() != NULL);}
      bool haveGapBuildPrms() const {return (m_gapPrmsPtr.get() != NULL);}
      const MeshTagBuildPrmsPy &getTagBuildPrms() const {return *(m_tagPrmsPtr.get());}
      const MeshGapBuildPrmsPy &getGapBuildPrms() const {return *(m_gapPrmsPtr.get());}
    }; // class
    
    void exportBondedMesh2dPrms();

  } // namespace lsm
} // namespace esys

#endif // ESYS_LSMBONDEDMESH2DPARAMSPY_H
