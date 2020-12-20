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

#include "Mesh2DReader.h"
#include "Mesh2DReaderImpl.h"

// --- Project includes ---
#include "Foundation/console.h"

// -- STL includes --
#include <stdexcept>

using std::runtime_error;

namespace esys
{
  namespace lsm //!Lattice Solid Model namespace.
  { 
    // ==== Node2DReader functions ==========
    /*!
      construct Node2DReader
      
      \param instream the stream to read from
    */
    Node2DReader::Node2DReader(std::istream& instream) 
      : IterativeReader<IStreamIterator<MeshNodeData2D> >(instream)
    {
      console.XDebug() << "Node2DReader constructor\n";
    }

    /*!
      initialize node reader
    */
    void Node2DReader::initialise()
    {
      console.XDebug() << "Node2DReader::initialise()\n";
      string token;

      while ((!getIStream().eof()) && (token != "2D-Nodes")) {
        getIStream() >> token;
      }
      if (token != "2D-Nodes") {
        throw runtime_error("Could not find '2D-Nodes' token in stream.");
      } 
      
      int nr_nodes = 0;
      getIStream() >> nr_nodes;

      setNumElements(nr_nodes);
      IterativeReader<IStreamIterator<MeshNodeData2D> >::initialise();
      console.XDebug() << "end Node2DReader::initialise()\n";
    }
    

    // ==== Edge2DReader Functions ===========

    /*!
      construct Edge2DReader
      
      \param instream the stream to read from
    */
    Edge2DReader::Edge2DReader(std::istream& instream) 
      : IterativeReader<IStreamIterator<MeshEdgeData2D> >(instream)
    {}

    /*!
      initialize edge reader
    */
    void Edge2DReader::initialise()
    {
      string token;

      while ((!getIStream().eof()) && (token != "Line2")) {
        getIStream() >> token;
      }
      if (token != "Line2") {
        throw runtime_error("Could not find 'Line2' token in stream.");
      }
      
      int nr_nodes = 0;
      getIStream() >> nr_nodes;

      setNumElements(nr_nodes);
      IterativeReader<IStreamIterator<MeshEdgeData2D> >::initialise();
    }

    // ==== Mesh2DReader functions ==========

    /*!
      construct a 2d mesh reader 

      \param filename the name of the file to read from
    */
    Mesh2DReader::Mesh2DReader(const string& filename)
      : m_impl_ptr(new Impl(filename))
    {
      m_impl_ptr->initialise();
    }

    /*!
      destroy a 2d mesh reader. Closes all accociated files.
    */
    Mesh2DReader::~Mesh2DReader()
    {
      delete m_impl_ptr;
    }

    /*!
      return iterator for nodes in file/stream 
    */
    Mesh2DReader::NodeIterator &Mesh2DReader::getNodeIterator()
    {
      console.XDebug() << "Mesh2DReader::getNodeIterator\n";
      return m_impl_ptr->m_node_reader_ptr->getIterator();
      console.XDebug() << "end Mesh2DReader::getNodeIterator\n";
    }

    /*!
      return iterator for triangles in file/stream 
    */
    Mesh2DReader::EdgeIterator &Mesh2DReader::getEdgeIterator()
    {   
      console.XDebug() << "Mesh2DReader::getEdgeIterator\n";
      return m_impl_ptr->m_edge_reader_ptr->getIterator();
      console.XDebug() << "end Mesh2DReader::getEdgeIterator\n";
    }

  } // end namespace lsm
} // end namespace esys
