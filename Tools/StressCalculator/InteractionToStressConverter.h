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


#ifndef ESYS_LSMINTERACTIONTOSTRESSCONVERTOR_H
#define ESYS_LSMINTERACTIONTOSTRESSCONVERTOR_H

#include "Tools/StressCalculator/StressTensorCalculator.h"
#include "Tools/StressCalculator/StressTensorCollection.h"
#include "Tools/StressCalculator/VtkUnstructuredGrid.h"
#include "Tools/StressCalculator/CartesianGrid.h"
#include "Tools/StressCalculator/Vec3Comparer.h"

#include "Foundation/BoundingBox.h"

#include <iostream>
#include <iomanip>

namespace std
{
  template <>
  struct less<Vec3>
  {
    bool operator()(const Vec3 &v1, const Vec3 &v2) const
    {
      return esys::lsm::Vec3XyzComparer()(v1, v2);
    }
  };
}

namespace esys
{
  namespace lsm
  {
    typedef CartesianGrid<double>         DoubleGrid;
    typedef CartesianGrid<StressTensor *> StressTensorPtrGrid;
    typedef CartesianGrid<Tensor>         TensorGrid;

    class InteractionToStressConverter
    {
    public:
      typedef StressTensorCollection<ContactPtTensorCalculator> StressTensCollection;
      typedef StressTensCollection::StressCalculator          StressTensorCalculator;

      InteractionToStressConverter(const BoundingBox &box, double gridSpacing);

      /**
       * Reads RAW2 interaction data from the specified stream and converts
       * it to stress tensor values (unstructured grid of tensors).
       */
      void addRaw2Interactions(std::istream &iStream);

      /**
       * Writes VTK UnstructuredGrid XML format file of tensor point-data
       * and eigenvalue point-data.
       */
      void writeVtkUnstructuredXml(std::ostream &oStream);

      void writeVtkUnstructuredXmlGridInformation(std::ostream &oStream);

      /**
       * Writes VTK StructuredGrid XML format file of
       * eigenvalue point-data to the specified output-stream.
       * The irregular point data
       * is converted to a regular grid using the
       * bounding-box and grid-spacing values of this object.
       *
       * @param oStream vtx xml written to this stream.
       *
       * @see getDevRegularGrid
       */
      void writeVtkStructuredXml(std::ostream &oStream);

      StressTensorPtrGrid &getTensorIrregularGrid();
      void calcTensorIrregularGrid();

      /**
       * Converts irregular grid of average particle stress tensors (\sigma_{ij})
       * to a regular grid.
       *
       * @param bBox - regular grid is restricted to this bounding box
       *               (all irregular point data are considered even if
       *               they are outside this bounding box.
       * @param gridSpacing - distance between regular grid points
       *                      (x, y and z dimensions).
       */
      TensorGrid &getTensorRegularGrid();

      void calcTensorRegularGrid();

      /**
       * Converts irregular grid of average particle stress tensor data
       * a regular grid using the volume of intersection between a sphere
       * and the box formed by a grid-cell.
       *
       * @param bBox - regular grid is restricted to this bounding box
       *               (all irregular point data are considered even if
       *               they are outside this bounding box. But the regular
       *               grid cells are contained in the box.
       * @param gridSpacing - distance between regular grid points
       *                      (x, y and z dimensions).
       */
      DoubleGrid &getDevRegularGrid();

      void calcDevRegularGrid();

      /**
       * Writes "x y z \sigma_{max}-\sigma_{min})" records
       * to the specifed file.
       * The irregular point data is converted to a regular
       * grid specified by bBox and gridSpacing parameters.
       *
       * @see getRegularGrid
       */      
      void writeFlatStructured(std::ostream &oStream);

      /**
       * Writes "x y z r \sigma_{max}-\sigma_{min})" records
       * to the specifed file.
       * An (x,y,z) value is the centre point of the spherical
       * particle of radius r.
       *
       * @see getRegularGrid
       */
      void writeFlatUnstructured
      (
        std::ostream &oStream
      );

      /**
       * Writes open-dx format file of unstructured point data
       * (\sigma_{max}-\sigma_{min}) to the specified stream.
       *
       * @param oStream Data is written to this stream.
       */
      void writeUnstructuredDx(std::ostream &oStream);

      double getMaxRadius();

    protected:
      double getRealDevStress(const Tensor &stressTensor) const;
      double getNormDevStress(const Tensor &stressTensor) const;

      bool is3d() const;
      
    private:
      double                 m_gridSpacing;
      BoundingBox            m_bBox;
      StressTensorCalculator m_stressCalculator;
      StressTensCollection m_stressTensorCollection;
      TensorGrid             m_regTensorGrid;
      DoubleGrid             m_regDevStressGrid;
      StressTensorPtrGrid    m_irrStressTensorGrid;
    };
  }
}
#endif
