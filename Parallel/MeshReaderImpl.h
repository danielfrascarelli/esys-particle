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

#ifndef __MESHREADERIMPL_H
#define __MESHREADERIMPL_H

// --- Project includes ---
#include "MeshReader.h"

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
      \class MeshReader::Impl
      \brief implementation details for the mesh reader. 
      Decouples Interface from implementation.

      \author Steffen Abe
      $Date$
      $Revision$
    */
    class MeshReader::Impl
    {
    public:
      // types
      typedef std::auto_ptr<NodeReader> NodeReaderPtr;
      typedef std::auto_ptr<TriReader>  TriReaderPtr;
      typedef std::auto_ptr<istream>    IStreamPtr;

      // variables
      NodeReaderPtr m_node_reader_ptr;
      TriReaderPtr  m_tri_reader_ptr;
      IStreamPtr    m_istream_ptr;
      string m_file_name;

      // functions

      Impl(const string&);
      void initialise();
    };
  }
}
#endif // __MESHREADERIMPL_H
