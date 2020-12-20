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

template <typename TmplPackable>
void BroadcastCommand::packInto(const TmplPackable &packable)
{
  packable.packInto(&(m_varBuffer));
} 

template <typename TmplData>
void BroadcastCommand::appendTypeAndName(const TmplData &namedWithType)
{
  m_varBuffer.append(namedWithType.getTypeString().c_str());
  m_varBuffer.append(namedWithType.getName().c_str());
}
  
template <typename TmplData>
void BroadcastCommand::append(const TmplData &data)
{
  m_varBuffer.append(data);
}
