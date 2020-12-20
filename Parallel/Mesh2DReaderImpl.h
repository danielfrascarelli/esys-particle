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

#ifndef __MESH2DREADERIMPL_H
#define __MESH2DREADERIMPL_H

// --- Project includes ---
#include "Mesh2DReader.h"

// --- STL includes ---
#include <string>

using std::string;

// --- IO includes ---
#include <iostream>

using std::istream;

namespace esys
{
  namespace lsm
  {
    /*!
      \class Mesh2DReader::Impl
      \brief implementation details for the 2d mesh reader. 
      Decouples Interface from implementation.

      \author Steffen Abe
      $Date$
      $Revision$
    */
    class Mesh2DReader::Impl
    {
    public:
      // types
      typedef std::auto_ptr<Node2DReader> NodeReaderPtr;
      typedef std::auto_ptr<Edge2DReader> EdgeReaderPtr;
      typedef std::auto_ptr<istream>    IStreamPtr;

      // variables
      NodeReaderPtr m_node_reader_ptr;
      EdgeReaderPtr m_edge_reader_ptr;
      IStreamPtr    m_istream_ptr;
      string m_file_name;

      // functions

      Impl(const string&);
      void initialise();
    };
  }
}
#endif // __MESH2DREADERIMPL_H
