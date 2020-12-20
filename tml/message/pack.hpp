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

template<>
template<typename T1,typename T2> 
void TML_PackedMessageInterface::pack<pair<T1,T2> >(const pair<T1,T2>& p)
{
  pack(p.first);
  pack(p.second);
}

template<typename T1,typename T2> 
void TML_PackedMessageInterface::unpack<pair<T1,T2> >(pair<T1,T2>& p)
{
  unpack(p.first);
  unpack(p.second);
}
