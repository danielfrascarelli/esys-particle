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


#ifndef __MESH_READER_H
#define __MESH_READER_H

//-- Project includes
#include "Model/MeshData.h"
#include "Parallel/IterativeReader.h"

//-- STL includes --
#include <string>
using std::string;

namespace esys
{
  namespace lsm //!Lattice Solid Model namespace.
  {
    /*!
      \class NodeReader
      \brief read a block of nodes from a Finley mesh file. 

      \author Steffen Abe
      $Date$
      $Revision$
    */
    class NodeReader : public IterativeReader<IStreamIterator<MeshNodeData> >
    {
    public:
      NodeReader(std::istream&);
      virtual void initialise();
    };


    /*!
      \class TriReader
      \brief read a block of triangles from a Finley mesh file. 
      
      \author Steffen Abe
      $Date$
      $Revision$
    */
    class TriReader : public IterativeReader<IStreamIterator<MeshTriData> >
    {
    public:
      TriReader(std::istream&);
      virtual void initialise();
    };

    /*!
      \class MeshReader
      \brief class to read triangle meshes from Finley mesh format files
      \author Steffen Abe
      $Date$
      $Revision$
    */
    class MeshReader
    {
    private:
      class Impl;
      Impl *m_impl_ptr; // pointer to the implementation

    public:
      // types
      typedef NodeReader::Iterator NodeIterator;
      typedef TriReader::Iterator  TriIterator;

      // functions
      MeshReader(const string&);
      ~MeshReader();

      NodeIterator &getNodeIterator();
      TriIterator &getTriIterator();
    };
  } // end namespace lsm
} // end namespace esys


#endif // __MESH_READER_H
