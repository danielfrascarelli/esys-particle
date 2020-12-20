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

#include "BroadCast_cmd.h"

BroadcastCommand::BroadcastCommand(const MpiRankAndComm &globalRankAndComm, int commandId) :
  m_commandId(commandId),
  m_varBuffer(globalRankAndComm.getComm(), 2000),
  m_barrier(globalRankAndComm.getComm()),
  m_cmdBuffer(globalRankAndComm.getComm(), globalRankAndComm.getRank())
{
}  

const int& BroadcastCommand::getCommandId() const
{
  return m_commandId;
}

void BroadcastCommand::broadcastCommand()
{
  m_cmdBuffer.broadcast(getCommandId());
}
  
void BroadcastCommand::broadcastBuffer()
{
  m_varBuffer.broadcast(0);  
}
  
void BroadcastCommand::wait(const std::string &msg)
{
  m_barrier.wait(msg.c_str());
}

void BroadcastCommand::broadcast()
{
  broadcastCommand();
  broadcastBuffer();
  wait("BroadcastCommand::broadcast");
}
