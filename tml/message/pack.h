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

#ifndef _PACK_H
#define _PACK_H

//--- project includes ---
#include "tml/message/packed_message_interface.h"


/*! \file pack.h */

/*!
\fn TML_pack(TML_Packed_Message*,const T&);
\brief templated function to pack objects into a TML_Packed_Message

Implemented as separate template function and not as member function of the object
to be packed for the following reasons: 
a) it forces redefinition in each derived class, even if there is a packing function for the base class
b) it makes packing for builtin types consistent

*/
template<typename T>
void TML_pack(TML_PackedMessageInterface*,const T&);

/*!
\fn TML_unpack(TML_Packed_Message*,T&);
\brief templated function to unpack objects from a TML_Packed_Message

Implemented as separate template function and not as member function of the object
to be packed for the following reasons: 
a) it forces redefinition in each derived class, even if there is a packing function for the base class
b) it makes packing for builtin types consistent

*/
template<typename T>
void TML_unpack(TML_PackedMessageInterface*,T&);

#endif //_PACK_H
