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

#include "MeshData2D.h"

// --- TML includes ---
#include "tml/message/packed_message_interface.h"

// === NODE DATA ===========

/*!
  read node data from istream

  \param instream the stream to read from
*/
void MeshNodeData2D::read(istream& instream)
{
  int dummy;
  instream >> id >> dummy >> tag >> x >> y ;
}


/*!
  Pack MeshNodeData2D into a TML packed message

  \param d the data
*/
template<>
void TML_PackedMessageInterface::pack<MeshNodeData2D>(const MeshNodeData2D& d)
{
  append(d.id);
  append(d.tag);
  append(d.x);
  append(d.y);
}

/*!
  Unpack MeshNodeData2D from a TML packed message
*/
template<>
void TML_PackedMessageInterface::unpack<MeshNodeData2D>(MeshNodeData2D& d)
{
  d.id=pop_int();
  d.tag=pop_int();
  d.x=pop_double();
  d.y=pop_double();
}


// === EDGE2D DATA ===========

/*!
  read 2D edge data from istream

  \param instream the stream to read from
*/
void MeshEdgeData2D::read(istream& instream)
{
  instream >> id >> tag >> p1 >> p2 ;
}

/*!
  Pack MeshEdgeData2D into a TML packed message

  \param d the data
*/
template<>
void TML_PackedMessageInterface::pack<MeshEdgeData2D>(const MeshEdgeData2D& d)
{
  append(d.id);
  append(d.tag);
  append(d.p1);
  append(d.p2);
}

/*!
  Unpack MeshEdgeData2D from a TML packed message
*/
template<>
void TML_PackedMessageInterface::unpack<MeshEdgeData2D>(MeshEdgeData2D& d)
{
  d.id=pop_int();
  d.tag=pop_int();
  d.p1=pop_int();
  d.p2=pop_int();
}
