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

#include "MeshReader.h"
#include "MeshReaderImpl.h"

// -- STL includes --
#include <stdexcept>

using std::runtime_error;

namespace esys
{
  namespace lsm //!Lattice Solid Model namespace.
  { 
    // ==== NodeReader functions ==========
    /*!
      construct NodeReader
      
      \param instream the stream to read from
    */
    NodeReader::NodeReader(std::istream& instream) 
      : IterativeReader<IStreamIterator<MeshNodeData> >(instream)
    {}

    /*!
      initialize node reader
    */
    void NodeReader::initialise()
    {
      string token;

      while ((!getIStream().eof()) && (token != "3D-Nodes")) {
        getIStream() >> token;
      }
      if (token != "3D-Nodes") {
        throw runtime_error("Could not find '3D-Nodes' token in stream.");
      }
      
      int nr_nodes = 0;
      getIStream() >> nr_nodes;

      setNumElements(nr_nodes);
      IterativeReader<IStreamIterator<MeshNodeData> >::initialise();
    }
    

    // ==== TriReader Functions ===========

    /*!
      construct TriReader
      
      \param instream the stream to read from
    */
    TriReader::TriReader(std::istream& instream) 
      : IterativeReader<IStreamIterator<MeshTriData> >(instream)
    {}

    /*!
      initialize triangle reader
    */
    void TriReader::initialise()
    {
      string token;

      while ((!getIStream().eof()) && (token != "Tri3")) {
        getIStream() >> token;
      }
      if (token != "Tri3") {
        throw runtime_error("Could not find 'Tri3' token in stream.");
      }
      
      int nr_nodes = 0;
      getIStream() >> nr_nodes;

      setNumElements(nr_nodes);
      IterativeReader<IStreamIterator<MeshTriData> >::initialise();
    }

    // ==== MeshReader functions ==========

    /*!
      construct a mesh reader to read from file

      \param filename the name of the file to read from
    */
    MeshReader::MeshReader(const string& filename)
      : m_impl_ptr(new Impl(filename))
    {
      m_impl_ptr->initialise();
    }

    /*!
      destroy a mesh reader. Closes all accociated files.
    */
    MeshReader::~MeshReader()
    {
      delete m_impl_ptr;
    }

    /*!
      return iterator for nodes in file/stream 
    */
    MeshReader::NodeIterator &MeshReader::getNodeIterator()
    {
      return m_impl_ptr->m_node_reader_ptr->getIterator();
    }

    /*!
      return iterator for triangles in file/stream 
    */
    MeshReader::TriIterator &MeshReader::getTriIterator()
    {   
      return m_impl_ptr->m_tri_reader_ptr->getIterator();
    }

  } // end namespace lsm
} // end namespace esys
