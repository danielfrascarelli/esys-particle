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

#include "MeshData.h"

// --- TML includes ---
#include "tml/message/packed_message_interface.h"

// === NODE DATA ===========

MeshNodeData::MeshNodeData()
  : id(-1), tag(-1), x(0.0), y(0.0), z(0.0)
{
}

MeshNodeData::MeshNodeData(int nodeId, const Vec3 &pt, int nodeTag)
  : id(nodeId), tag(nodeTag), x(pt[0]), y(pt[1]), z(pt[2])
{
}

/*!
  read node data from istream

  \param instream the stream to read from
*/
void MeshNodeData::read(std::istream& instream)
{
  int dummy;
  instream >> id >> dummy >> tag >> x >> y >> z;
}


/*!
  Pack MeshNodeData into a TML packed message

  \param d the data
*/
template<>
void TML_PackedMessageInterface::pack<MeshNodeData>(const MeshNodeData& d)
{
  append(d.id);
  append(d.tag);
  append(d.x);
  append(d.y);
  append(d.z);
}

/*!
  Unpack MeshNodeData from a TML packed message
*/
template<>
void TML_PackedMessageInterface::unpack<MeshNodeData>(MeshNodeData& d)
{
  d.id=pop_int();
  d.tag=pop_int();
  d.x=pop_double();
  d.y=pop_double();
  d.z=pop_double();
}


// === TRI DATA ===========

MeshTriData::MeshTriData()
  : id(-1), tag(-1), p1(-1), p2(-1), p3(-1)
{
}

MeshTriData::MeshTriData(int triId, int nodeId0, int nodeId1, int nodeId2, int triTag)
  : id(triId), tag(triTag), p1(nodeId0), p2(nodeId1), p3(nodeId2)
{
}

/*!
  read triangle data from istream

  \param instream the stream to read from
*/
void MeshTriData::read(std::istream& instream)
{
  instream >> id >> tag >> p1 >> p2 >> p3;
}

/*!
  Pack MeshTriData into a TML packed message

  \param d the data
*/
template<>
void TML_PackedMessageInterface::pack<MeshTriData>(const MeshTriData& d)
{
  append(d.id);
  append(d.tag);
  append(d.p1);
  append(d.p2);
  append(d.p3);
}

/*!
  Unpack MeshTriData from a TML packed message
*/
template<>
void TML_PackedMessageInterface::unpack<MeshTriData>(MeshTriData& d)
{
  d.id=pop_int();
  d.tag=pop_int();
  d.p1=pop_int();
  d.p2=pop_int();
  d.p3=pop_int();
}
