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


#include "Tools/StressCalculator/InteractionToStressConverter.h"
#include "Tools/StressCalculator/VtkStructuredGrid.h"
#include "Tools/StressCalculator/GaussianGridder.h"
#include "Tools/StressCalculator/Raw2InteractionReader.h"
#include "Tools/StressCalculator/ContactCollection.h"
#include "Tools/StressCalculator/EigenvalueCalculator.h"

#include "Geometry/SphereBoxVolCalculator.h"
#include "Geometry/CircleBoxVolCalculator.h"

namespace esys
{
  namespace lsm
  {
    template <typename TmplSphere, typename TmplBox>
    std::string getDetailsString(const TmplSphere &sphere, const TmplBox &box)
    {
      std::stringstream msg;
      msg
        << "box = ("
        << box.getMinPt() << ", " << box.getMaxPt() << ")"
        << ", sphere = (" << sphere.getCentre() << ", " << sphere.getRadius() << ")";
      return msg.str();
    }

    template <typename TmplSphere, typename TmplBox>
    void checkIntersectionVolume(double vol, const TmplSphere &sphere, const TmplBox &box)
    {
      if (false)
      {
        std::stringstream msg;
        msg 
          << "nan encountered during volume calculation: "
          << getDetailsString(sphere, box);

        throw std::runtime_error(msg.str());
      }
      else if ((vol < 0.0) && (fabs(vol) > 1.0e-6))
      {
        std::stringstream msg;
        msg 
          << "Negative intersection volume " << vol 
          << ". " << getDetailsString(sphere, box);
        throw std::runtime_error(msg.str());
      }
      else if (vol > (box.getVolume() + (box.getVolume()*1.0e-6)))
      {
        std::stringstream msg;
        msg 
          << "Volume " << vol 
          << " larger than box volume " << box.getVolume()
          << ". " << getDetailsString(sphere, box);
        throw std::runtime_error(msg.str());
      }
      else if (vol > (sphere.getVolume() + (sphere.getVolume()*1.0e-6)))
      {
        std::stringstream msg;
        msg 
          << "Volume " << vol
          << " larger than sphere volume " << sphere.getVolume()
          << ". " << getDetailsString(sphere, box);
        throw std::runtime_error(msg.str());
      }
    }

    class ThreeDIntersectionCalker : public SphereBoxVolCalculator
    {
    public:
      ThreeDIntersectionCalker(const BoundingBox &box)
        : SphereBoxVolCalculator(Box(box.getMinPt(), box.getMaxPt()))
      {
      }

      double getVolume(const Sphere &sphere)
      {
        const double vol = SphereBoxVolCalculator::getVolume(sphere);
        checkIntersectionVolume(vol, sphere, getBox());
        return vol;
      }
    };

    class TwoDIntersectionCalker  : public CircleBoxVolCalculator
    {
    public:
      TwoDIntersectionCalker(const BoundingBox &box)
        :  CircleBoxVolCalculator(Box(box.getMinPt(), box.getMaxPt()))
      {
      }

      double getVolume(const Sphere &sphere)
      {
        const double vol = CircleBoxVolCalculator::getVolume(sphere);
        checkIntersectionVolume(vol, sphere, getBox());
        return vol;
      }
    };

    class PointDataType
      : public
          vtk::DataTypeTuple<
            vtk::Float64Type,
            vtk::Float64Type,
            vtk::Matrix3Type,
            vtk::Float64Type
          >
    {
    public:
      typedef
        vtk::DataTypeTuple<
          vtk::Float64Type,
          vtk::Float64Type,
          vtk::Matrix3Type,
          vtk::Float64Type
        > Inherited;

      PointDataType()
        : Inherited(
            vtk::Float64Type("|sMax-sMin|"),
            vtk::Float64Type("Real(sMax-sMin)"),
            vtk::Matrix3Type("stressTensor"),
            vtk::Float64Type("radius")
          )
      {
      }
    };

    class PointDataTypeForGrid
      : public
          vtk::DataTypeTuple<
            vtk::Float64Type,
            vtk::Float64Type,
            vtk::Matrix3Type
          >
    {
    public:
      typedef
        vtk::DataTypeTuple<
          vtk::Float64Type,
          vtk::Float64Type,
          vtk::Matrix3Type
        > Inherited;

      PointDataTypeForGrid()
        : Inherited(
            vtk::Float64Type("|sMax-sMin|"),
            vtk::Float64Type("Real(sMax-sMin)"),
            vtk::Matrix3Type("stressTensor")
          )
      {
      }


    };

    typedef vtk::Vec3Type                                    PointType;
    typedef vtk::UnstructuredPiece<PointType, PointDataType> Piece;
    typedef vtk::UnstructuredPiece<PointType, PointDataTypeForGrid> PieceForGrid;
    typedef EigenvalueCalculator::ComplexRealImagComparer    RealImagComparer;
    typedef EigenvalueCalculator::ComplexAbsRealImagComparer AbsRealImagComparer;
    typedef EigenvalueCalculator::ComplexNormComparer        NormComparer;

    InteractionToStressConverter::InteractionToStressConverter(const BoundingBox &bBox, double gridSpacing)
      : m_gridSpacing(gridSpacing),
        m_bBox(bBox),
        m_stressCalculator(),
        m_stressTensorCollection(m_stressCalculator),
        m_regTensorGrid(bBox, gridSpacing),
        m_regDevStressGrid(bBox, gridSpacing),
        m_irrStressTensorGrid(bBox, gridSpacing)
    {
    }

    double InteractionToStressConverter::getRealDevStress(const Tensor &stressTensor) const
    {
      StressTensor::ComplexVector eigenvals = stressTensor.getEigenvalues();
      double devStress = 0.0;

      if (is3d())
      {
        std::sort(eigenvals.begin(), eigenvals.end(), RealImagComparer());
        devStress = (eigenvals[2].real() - eigenvals[0].real());
      }
      else
      {
        //
        // In 2d, make sure we only ever subtract \sigma_1 and \sigma_2,
        // so we sort by absolute value. This ensures that the zero-eigenvalue
        // (corresponding to the z-dimension) will be sorted into the
        // eigenvals[0] element.
        //
        std::sort(eigenvals.begin(), eigenvals.end(), AbsRealImagComparer());
        devStress = fabs(eigenvals[2].real() - eigenvals[1].real());
      }

      return devStress;
    }

    double InteractionToStressConverter::getNormDevStress(const Tensor &stressTensor) const
    {
      StressTensor::ComplexVector eigenvals = stressTensor.getEigenvalues();
      std::sort(eigenvals.begin(), eigenvals.end(), NormComparer());
      
      return
        (
          is3d()
          ?
          std::norm(eigenvals[2] - eigenvals[0])
          :
          std::norm(eigenvals[2] - eigenvals[1])
        );
    }

    bool InteractionToStressConverter::is3d() const
    {
      return ParticleData::is3d();
    }

    void InteractionToStressConverter::addRaw2Interactions(std::istream &iStream)
    {
      Raw2InteractionReader reader(iStream);
      ContactCollection contacts;
      contacts.addInteractions(reader);

      m_stressTensorCollection.addContactIterators(contacts.getContactIteratorIterator());
    }

    void InteractionToStressConverter::writeVtkUnstructuredXml(std::ostream &oStream)
    {
      StressTensCollection::StressTensorIterator it = m_stressTensorCollection.getIterator();
      Piece piece(PointType("points"), PointDataType());
      while (it.hasNext()) {
        StressTensCollection::StressTensorIterator::reference stressTensor = it.next();

        const double realDevStress = getRealDevStress(stressTensor);
        const double normDevStress = getNormDevStress(stressTensor);

        piece.setPoint(
          stressTensor.getPos(),
          PointDataType::DataValueTuple(
            normDevStress,
            realDevStress,
            stressTensor.getTensor(),
            stressTensor.getRad()
          )
        );
      }
      // Write the xml document header
      oStream << "<?xml version=\"1.0\"?>" << std::endl;
      vtk::UnstructuredGrid grid;
      grid.addPiece(piece);
      grid.writeXml(oStream);      
    }

    void InteractionToStressConverter::writeVtkUnstructuredXmlGridInformation(std::ostream &oStream)
    {
      TensorGrid &regularGrid = getTensorRegularGrid();
      TensorGrid::CellIterator cellIt = regularGrid.getCellIterator();
      PieceForGrid piece(PointType("points"), PointDataTypeForGrid());

      while (cellIt.hasNext())
      {
        TensorGrid::Cell::Iterator posValIt = cellIt.next().getIterator();
        while (posValIt.hasNext())
        {
          const TensorGrid::Cell::PosValuePair &pair = posValIt.next();
          const double realDevStress = getRealDevStress(pair.getValue());
          const double normDevStress = getNormDevStress(pair.getValue());

          piece.setPoint(
            pair.getPos(),
            PointDataTypeForGrid::DataValueTuple(
              normDevStress,
              realDevStress,
              pair.getValue().getTensor()
            )
          );
        }
      }
      // Write the xml document header
      oStream << "<?xml version=\"1.0\"?>" << std::endl;
      vtk::UnstructuredGrid grid;
      grid.addPiece(piece);
      grid.writeXml(oStream);
    }

    void InteractionToStressConverter::writeUnstructuredDx(std::ostream &oStream)
    {
      oStream << "points = " << m_stressTensorCollection.size() << std::endl;
      oStream << "format = ascii" << std::endl;
      oStream << "dependency = positions, positions" << std::endl;
      oStream << "interleaving = field" << std::endl;
      oStream << "field = locations, principleStressDiff" << std::endl;
      oStream << "structure = 3-vector, scalar" << std::endl;
      oStream << "type = float, float" << std::endl;
      oStream << "header = marker \"Start\\n\"" << std::endl << std::endl;
      oStream << "end" << std::endl;
      oStream << "Start" << std::endl;
      StressTensCollection::StressTensorIterator it = m_stressTensorCollection.getIterator();
      while (it.hasNext()) {
        StressTensCollection::StressTensorIterator::reference stressTensor = it.next();
        const double val = getRealDevStress(stressTensor);
        oStream << stressTensor.getPos() << " " << val << "\n";
      }
    }

    void InteractionToStressConverter::writeFlatUnstructured(std::ostream &oStream)
    {
      StressTensCollection::StressTensorIterator it = m_stressTensorCollection.getIterator();
      while (it.hasNext()) {
        StressTensCollection::StressTensorIterator::reference stressTensor = it.next();
        const double val = getRealDevStress(stressTensor);
        oStream
          << stressTensor.getPos() << " "
          << stressTensor.getRad() << " "
          << val << "\n";
      }
    }

    class StrctPointDataType
      : public
          vtk::DataTypeTuple<
            vtk::Float64Type
          >
    {
    public:
      typedef
        vtk::DataTypeTuple<
          vtk::Float64Type
        > Inherited;

      StrctPointDataType()
        : Inherited(
            vtk::Float64Type("sMax-sMin")
          )
      {
      }
    };

    typedef vtk::Vec3Type                                            StrctPointType;
    typedef vtk::StructuredPiece<StrctPointType, StrctPointDataType> StrctPiece;

    StressTensorPtrGrid &InteractionToStressConverter::getTensorIrregularGrid()
    {
      if (m_irrStressTensorGrid.size() <= 0)
      {
        calcTensorIrregularGrid();
      }
      return m_irrStressTensorGrid;
    }

    double InteractionToStressConverter::getMaxRadius()
    {
      StressTensorPtrGrid::ValueIterator it = getTensorIrregularGrid().getValueIterator();
      double maxRadius = -1.0;

      while (it.hasNext())
      {
        StressTensorPtrGrid::const_reference stressTensor = it.next();
        if (stressTensor->getRad() > maxRadius)
        {
          maxRadius = stressTensor->getRad();
        }
      }
      return maxRadius;
    }

    void InteractionToStressConverter::calcTensorIrregularGrid()
    {
      m_irrStressTensorGrid = StressTensorPtrGrid(m_bBox, m_gridSpacing);
      StressTensCollection::StressTensorIterator it = m_stressTensorCollection.getIterator();

      while (it.hasNext()) {
        StressTensCollection::StressTensorIterator::reference stressTensor = it.next();
        m_irrStressTensorGrid.insert(stressTensor.getPos(), &stressTensor);
      }
    }

    TensorGrid &InteractionToStressConverter::getTensorRegularGrid()
    {
      if (m_regTensorGrid.size() <= 0)
      {
        calcTensorRegularGrid();
      }
      return m_regTensorGrid;
    }

    template <typename TmplCellIterator, typename TmplIntsectVolCalker>
    Matrix3 getBoxTensor(
      TmplCellIterator cellIt,
      TmplIntsectVolCalker intersectCalker
    )
    {
      Matrix3 tensor;
      while (cellIt.hasNext())
      {
        typename TmplCellIterator::value_type::Iterator pairIt = cellIt.next().getIterator();
        while (pairIt.hasNext())
        {
          const StressTensor *stressTensor = pairIt.next().getValue();
          const typename TmplIntsectVolCalker::Sphere sphere(stressTensor->getPos(), stressTensor->getRad());
          const double intersectVol = intersectCalker.getVolume(sphere);
//          std::cout << "sphere: " << sphere.getCentre() << ", " << sphere.getRadius() << ", vol: " << intersectVol << std::endl;
          tensor += ((stressTensor->getTensor())*intersectVol);
        }
      }
      //
      // Return the average stress over the cell.
      //
      return tensor*(1.0/intersectCalker.getBox().getVolume());
    }

    void InteractionToStressConverter::calcTensorRegularGrid()
    {
      StressTensorPtrGrid &irregularGrid = getTensorIrregularGrid();

      m_regTensorGrid = TensorGrid(m_bBox, m_gridSpacing);

      const double maxRadius = (getMaxRadius() + m_gridSpacing);

      TensorGrid::CellIterator regCellIt = m_regTensorGrid.getCellIterator();
      while (regCellIt.hasNext())
      {
        TensorGrid::Cell &regCell = regCellIt.next();
        StressTensorPtrGrid::CellIterator irrCellIt = irregularGrid.getCellIterator(regCell.getPos(), maxRadius);

        if (is3d())
        {
          m_regTensorGrid.insert(
            regCell.getPos(),
            Tensor(
              regCell.getPos(),
              getBoxTensor(
                irrCellIt,
                ThreeDIntersectionCalker(regCell.getBox())
              )
            )
          );
        }
        else
        {
          m_regTensorGrid.insert(
            regCell.getPos(),
            Tensor(
              regCell.getPos(),
              getBoxTensor(
                irrCellIt,
                TwoDIntersectionCalker(regCell.getBox())
              )
            )
          );
        }
      }
    }

    DoubleGrid &InteractionToStressConverter::getDevRegularGrid()
    {
      if (m_regDevStressGrid.size() <= 0)
      {
        calcDevRegularGrid();
      }
      return m_regDevStressGrid;
    }

    void InteractionToStressConverter::calcDevRegularGrid()
    {
      const TensorGrid &regularTensorGrid = getTensorRegularGrid();
      m_regDevStressGrid = DoubleGrid(m_bBox, m_gridSpacing);

      TensorGrid::CellConstIterator cellIt = regularTensorGrid.getCellIterator();
      while (cellIt.hasNext())
      {
        TensorGrid::Cell::ConstIterator pairIt = cellIt.next().getIterator();
        while (pairIt.hasNext())
        {
          const Tensor &tensor = pairIt.next().getValue();
          m_regDevStressGrid.insert(
            tensor.getPos(),
            getRealDevStress(tensor)
          );
        }
      }
    }

    void InteractionToStressConverter::writeVtkStructuredXml(std::ostream &oStream)
    {
      TensorGrid &regularGrid = getTensorRegularGrid();
      StrctPiece piece(StrctPointType("points"), StrctPointDataType());
      piece.setExtent(regularGrid.getMinVecIndex(), regularGrid.getMaxVecIndex());
      TensorGrid::CellIterator cellIt = regularGrid.getCellIterator();
      while (cellIt.hasNext())
      {
        TensorGrid::Cell::Iterator posValIt = cellIt.next().getIterator();
        while (posValIt.hasNext())
        {
          const TensorGrid::Cell::PosValuePair &pair = posValIt.next();
          piece.setPoint(
            pair.getPos(),
            StrctPiece::PointData(
              getRealDevStress(pair.getValue())
              //pair.getValue().getTensor()
            )
          );
        }
      }
      // Write the xml document header
      oStream << "<?xml version=\"1.0\"?>" << std::endl;
      vtk::StructuredGrid grid;
      grid.setExtent(piece.getMinExtent(), piece.getMaxExtent());      
      grid.addPiece(piece);
      grid.writeXml(oStream);
    }
    
    void InteractionToStressConverter::writeFlatStructured(std::ostream &oStream)
    {
      DoubleGrid regular = getDevRegularGrid();
      DoubleGrid::CellIterator cellIt = regular.getCellIterator();
      while (cellIt.hasNext())
      {
        DoubleGrid::Cell::Iterator posValIt = cellIt.next().getIterator();
        while (posValIt.hasNext())
        {
          const DoubleGrid::Cell::PosValuePair &pair = posValIt.next();
          oStream << pair.getPos() << " " << pair.getValue() << "\n";
        }
      }
    }
  }
}
