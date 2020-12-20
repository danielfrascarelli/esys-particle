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

#ifndef __BROADCAST_CMD_H
#define __BROADCAST_CMD_H

#include "Parallel/mpicmdbuf.h"
#include "Parallel/mpivbuf.h"
#include "Parallel/mpibarrier.h"
#include "Parallel/RankAndComm.h"

/*!
  \brief base class for broadcast commands
*/
class BroadcastCommand
{
 private:
  int            m_commandId;
  CVarMPIBuffer  m_varBuffer;
  CMPIBarrier    m_barrier;
  CMPILCmdBuffer m_cmdBuffer;

 public:
  BroadcastCommand(const MpiRankAndComm &rankAndComm, int cmdId);
  virtual ~BroadcastCommand(){}

  /**
   * Appends namedWithType.getTypeString() and namedWithType.getName()
   * strings to the data buffer.
   */
  template <typename TmplData>
  void appendTypeAndName(const TmplData &namedWithType);

  /**
   * Appends specified argument data (basic type) to the buffer.
   */
  template <typename TmplData>
  void append(const TmplData &basicTypeData);

  /**
   * Packs the specified data into the data-buffer.
   */
  template <typename TmplPackable>
  void packInto(const TmplPackable&);

  /**
   * Returns the command id of this broadcast-command.
   */
  const int& getCommandId() const;

  /**
   * Broadcasts the command (ie the command id).
   */
  void broadcastCommand();

  /**
   * Broadcasts the data buffer.
   */
  void broadcastBuffer();

  /**
   * Barrier wait.
   */
  void wait(const std::string &barrierName);

  /**
   * Broadcasts command and data buffer, then does
   * a barrier wait.
   */
  void broadcast();

};
#include "Parallel/BroadCast_cmd.hpp"

#endif // BROADCAST_CMD_H
