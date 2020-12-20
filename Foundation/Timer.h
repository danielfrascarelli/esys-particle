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

#ifndef __MPI_TIMER_H_
#define __MPI_TIMER_H_

//--- MPI includes ---
#include <mpi.h>

//--- STL includes ---
#include <string>
#include <map>

#include <boost/shared_ptr.hpp>

/*!
 * Simple timer class for storing a start and stop time.
 */
class MpiWTimer
{
public:
  MpiWTimer();

  MpiWTimer(const std::string &name);
  void setStart(const double &wTime);
  void pause(const double &wTime);
  void resume(const double &wTime);
  void setStop(const double &wTime, bool elapseIsStopMinusStart=false);
  double getTiming() const;
  const std::string &getName() const;

  bool isPaused() const;

  void zeroise();
protected:
  void isPaused(bool paused);
  
private:
  std::string m_name;
  double      m_startTime;
  double      m_stopTime;
  bool        m_isPaused;
  double      m_pauseTime;
  double      m_resumeTime;
  double      m_elapsedTime;
};

class MpiWTimers;

/**
 * Helper for appending timing data to file, caches an ofstream for a file.
 */
class TimingDataWriter
{
public:
  TimingDataWriter(const std::string &fileName, MpiWTimers &timers);

  std::ostream &getOStream();
  
  const std::string &getFileName() const;

  void writeHeader();

  void appendData();

private:
  std::string      m_fileName;
  MpiWTimers       *m_pTimers;
  bool             m_haveWrittenHeader;
  typedef boost::shared_ptr<std::ofstream> OFStreamPtr;
  OFStreamPtr      m_oFStreamPtr;
};

/*!
 * Helper class for recording various pieces of timing info.
 */
class MpiWTimers
{
public:
  MpiWTimers();
  
  void start(const std::string &name);
  void stop(const std::string &name, bool elapseIsStopMinusStart=false);
  void pause(const std::string &name);
  void resume(const std::string &name);
  void zeroise(const std::string &name);
  void zeroise();
  bool timerExists(const std::string &name) const;
  double getTiming(const std::string &name) const;

  void writeHeader(std::ostream &oStream);
  void appendData(std::ostream &oStream);
  void appendData(const std::string &fileName);
  void clear();

protected:
  //void createTimer(const std::string &timerName);

  MpiWTimer *findTimer(const std::string &timerName);
  const MpiWTimer *findTimer(const std::string &timerName) const;

  MpiWTimer &findOrCreateTimer(const std::string &timerName);

  TimingDataWriter &getWriter(const std::string &fileName);

private:
  typedef std::map<std::string, MpiWTimer> NameMpiWTimerMap;
  NameMpiWTimerMap m_timerMap;

  typedef std::map<std::string, TimingDataWriter> FileNameWriterMap;
  FileNameWriterMap m_fileNameWriterMap;
};

#endif //__MPI_TIMER_H_
