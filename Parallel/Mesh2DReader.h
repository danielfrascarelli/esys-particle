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


#ifndef __MESH_2D_READER_H
#define __MESH_2D_READER_H

//-- Project includes
#include "Model/MeshData2D.h"
#include "Parallel/IterativeReader.h"

//-- STL includes --
#include <string>
using std::string;

namespace esys
{
  namespace lsm //!Lattice Solid Model namespace.
  {
    /*!
      \class Node2DReader
      \brief read a block of 2D nodes from a Finley mesh file. 

      \author Steffen Abe
      $Date$
      $Revision$
    */
    class Node2DReader : public IterativeReader<IStreamIterator<MeshNodeData2D> >
    {
    public:
      Node2DReader(std::istream&);
      virtual void initialise();
    };


    /*!
      \class Edge2DReader
      \brief read a block of edges from a Finley 2D mesh file. 
      
      \author Steffen Abe
      $Date$
      $Revision$
    */
    class Edge2DReader : public IterativeReader<IStreamIterator<MeshEdgeData2D> >
    {
    public:
      Edge2DReader(std::istream&);
      virtual void initialise();
    };

    /*!
      \class Mesh2DReader
      \brief class to read 2D meshes, or more precisely, the edges thereof, from Finley mesh format files
      \author Steffen Abe
      $Date$
      $Revision$
    */
    class Mesh2DReader
    {
    private:
      class Impl;
      Impl *m_impl_ptr; // pointer to the implementation

    public:
      // types
      typedef Node2DReader::Iterator NodeIterator;
      typedef Edge2DReader::Iterator EdgeIterator;

      // functions
      Mesh2DReader(const string&);
      ~Mesh2DReader();

      NodeIterator &getNodeIterator();
      EdgeIterator &getEdgeIterator();
    };
  } // end namespace lsm
} // end namespace esys


#endif // __MESH_READER_H
