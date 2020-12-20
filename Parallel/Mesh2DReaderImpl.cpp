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

#include "Mesh2DReaderImpl.h"

// --- Project includes ---
#include "Foundation/console.h"

// --- STL includes ---
#include <stdexcept>

using std::runtime_error;

// --- IO includes ---
#include <fstream>

using std::ifstream;

namespace esys
{
  namespace lsm
  {
    /*!
      construct Mesh2DReader implementation from file
    
      \param filename the name of the file
    */
    Mesh2DReader::Impl::Impl(const string& filename)
      :m_file_name(filename)
    {
    }

    void Mesh2DReader::Impl::initialise()
    {
      // open file
      m_istream_ptr=IStreamPtr(new ifstream(m_file_name.c_str()));
      istream* the_istream=m_istream_ptr.get();

      // check if file is open
      if(!(*m_istream_ptr)){
        throw runtime_error("Can not open 2D mesh file " + m_file_name);
      }
      console.XDebug() << "Reading 2D mesh file " << m_file_name << "  open \n";

      // initialize reader with stream
      m_node_reader_ptr=NodeReaderPtr(new Node2DReader(*the_istream));
      m_edge_reader_ptr=EdgeReaderPtr(new Edge2DReader(*the_istream));
      console.XDebug() << "end Mesh2DReader::Impl::initialise\n";
    }
  } // end namespace lsm
} // end namespace esys
