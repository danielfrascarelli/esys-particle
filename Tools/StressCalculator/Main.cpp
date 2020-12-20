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


#include <iostream>
#include <fstream>
#include <stdexcept>

#include "Tools/StressCalculator/InteractionToStressConverter.h"
#include "Foundation/StringUtil.h"
#include "Foundation/BoundingBox.h"

using namespace esys::lsm;

BoundingBox getBBox(const std::string &arg)
{
  typedef std::vector<double> DoubleVector;
  DoubleVector vec = StringUtil::split<double,StringUtil::StdIStreamOp<double> >(arg, ",");
  BoundingBox bBox;
  if (vec.size() == 6)
  {
    bBox = BoundingBox(Vec3(vec[0], vec[1], vec[2]), Vec3(vec[3], vec[4], vec[5]));
  }
  else
  {
    std::stringstream msg;
    msg << "Bounding box argument " << arg << " does not contain 6 comma separated elements.";
    throw std::runtime_error(msg.str());
  }
  return bBox;
}

int main(int argc, char *argv[])
{
  try
  {
    if (argc > 4) {
      BoundingBox bBox               = getBBox(argv[1]);
      double gridSpacing             = StringUtil::to<double>(argv[2]);
      std::string inputRaw2File      = argv[3];
      std::string outputUnstrXmlGridFile = argv[4];

      std::string outputStrctFlatFile = "";
      if (argc > 5)
      {
        outputStrctFlatFile = argv[5];
      }

      std::string outputUnstrXmlFile = "";
      if (argc > 6)
      {
        outputUnstrXmlFile = argv[6];
      }

      std::string outputUnstrFlatFile = "";
      if (argc > 7)
      {
        outputUnstrFlatFile = argv[7];
      }

      ParticleData::is3d(bBox.getSizes().Z() != 0.0);
      std::cout
        << "Dim = " << (ParticleData::is3d() ? "3D" : "2D")
        << ", BBox = " << bBox
        << ", grid spacing = " << gridSpacing
        << std::endl;

      std::ifstream iStream(inputRaw2File.c_str());
      std::cerr << "Reading from file " << inputRaw2File << std::endl;
      InteractionToStressConverter converter(bBox, gridSpacing);
      converter.addRaw2Interactions(iStream);      
      std::cerr << "Writing to file " << outputUnstrXmlGridFile << std::endl;
      std::ofstream oStreamStrct(outputUnstrXmlGridFile.c_str());
      converter.writeVtkUnstructuredXmlGridInformation(oStreamStrct);
      if (!outputStrctFlatFile.empty())
      {
        std::cerr << "Writing to file " << outputStrctFlatFile << std::endl;
        std::ofstream oStreamFlat(outputStrctFlatFile.c_str());
        converter.writeFlatStructured(oStreamFlat);
      }
      if (!outputUnstrXmlFile.empty())
      {
        std::cerr << "Writing to file " << outputUnstrXmlFile << std::endl;
        std::ofstream oStreamUnstr(outputUnstrXmlFile.c_str());
        converter.writeVtkUnstructuredXml(oStreamUnstr);
      }
      if (!outputUnstrFlatFile.empty())
      {
        std::cerr << "Writing to file " << outputUnstrFlatFile << std::endl;
        std::ofstream oStreamUnstr(outputUnstrFlatFile.c_str());
        converter.writeFlatUnstructured(oStreamUnstr);
      }
    }
    else {
      std::cerr
        << "Usage: "
        << argv[0]
        << " bbox gridSpacing input_raw2_file"
        << " outputUnstrXmlGridFile [outputStrctFlatFile [outputUnstrXmlFile [outputUnstrFlatFile]]" << std::endl
        << "Converts RAW2 interaction record data (centrePt1 rad1 centrePt2 rad2 contactPt force)" << std::endl
        << "to stress values (sigma_{max}-sigma_{min})."
        << std::endl << std::endl
        << "bbox - bounding box for regular grid generation. Specified as a comma"   << std::endl
        << "       separated list of six values \"minX,minY,minZ,maxX,maxY,maxZ\""   << std::endl
        << "gridSpacing - distance between points of the regular grid."              << std::endl
        << "input_raw2_file - RAW2 interaction data output from LSM"                 << std::endl
        << "                  VectorInteractionFieldSaver."                          << std::endl
        << "outputUnstrXmlGridFile - A VTK unstructured grid is written to this"     << std::endl
        << "                         file containing dev-stress and the full tensor."<< std::endl
        << "outputStrctFlatFile - Simple ascii records (one record per line:"        << std::endl
        << "                      \"x y z \\sigma_{dev}\") are written to this file."<< std::endl
        << "outputUnstrXmlFile - A per particle VTK unstructured grid of tensor,"    << std::endl
        << "                     particle radius and dev-stress data are written"    << std::endl
        << "                     to this file."                                      << std::endl
        << "outputUnstrFlatFile - Simple ascii records (one record per line:"        << std::endl
        << "                      \"x y z \\sigma_{dev}\") are written to this file."<< std::endl        
        << std::endl << std::endl;
    }
  }
  catch (std::runtime_error &e)
  {
    std::cerr << "EXCEPTION: " << e.what() << std::endl;
    throw;
  }
  catch (const char *e)
  {
    std::cerr << "EXCEPTION: "  << std::string(e) << std::endl;
    throw;
  }
  catch (...)
  {
    std::cerr << "EXCEPTION: "  << "Unknown exception." << std::endl;
    throw;
  }
  return 0;
}
