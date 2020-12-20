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

#include "Timer.h"

//--- IO includes ---
#include <fstream>
#include <iomanip>


void MpiWTimer::zeroise()
{
  m_startTime = 0.0;
  m_stopTime = 0.0;
  m_isPaused = false;
  m_pauseTime = 0.0;
  m_resumeTime = 0.0;
  m_elapsedTime = 0.0;
}

MpiWTimer::MpiWTimer()
  : m_name(),
    m_startTime(0),
    m_stopTime(0),
    m_isPaused(false),
    m_pauseTime(0),
    m_resumeTime(0),
    m_elapsedTime(0)
{
  zeroise();
}

MpiWTimer::MpiWTimer(const std::string &name)
  : m_name(name),
    m_startTime(0),
    m_stopTime(0),
    m_isPaused(false),
    m_pauseTime(0),
    m_resumeTime(0),
    m_elapsedTime(0)
{
  zeroise();
}

bool MpiWTimer::isPaused() const
{
  return m_isPaused;
}

void MpiWTimer::isPaused(bool paused)
{
  m_isPaused = paused;
}

void MpiWTimer::pause(const double &wTime)
{
  if (!isPaused()) {
    isPaused(true);
    m_pauseTime = wTime;
    m_elapsedTime += ((m_pauseTime > m_resumeTime) ? m_pauseTime - m_resumeTime : 0.0);
  }
}

void MpiWTimer::resume(const double &wTime)
{
  isPaused(false);
  m_resumeTime = wTime;
}

const std::string &MpiWTimer::getName() const
{
  return m_name;
}

void MpiWTimer::setStart(const double &wTime)
{
  m_startTime = wTime;
  m_elapsedTime = 0;
  resume(wTime);
}

void MpiWTimer::setStop(const double &wTime, bool elapseIsStopMinusStart)
{
  m_stopTime = wTime;
  pause(wTime);
  if (elapseIsStopMinusStart) {
    m_elapsedTime = m_stopTime - m_startTime;
  }
}

double MpiWTimer::getTiming() const
{
  return m_elapsedTime;
}

//======== TimingDataWriter =================
TimingDataWriter::TimingDataWriter(const std::string &fileName, MpiWTimers &timers)
  : m_fileName(fileName),
    m_pTimers(&timers),
    m_haveWrittenHeader(false),
    m_oFStreamPtr()
{
}

const std::string &TimingDataWriter::getFileName() const
{
  return m_fileName;
}

std::ostream &TimingDataWriter::getOStream()
{
  if (m_oFStreamPtr.get() == NULL) {
    m_oFStreamPtr = OFStreamPtr(new std::ofstream(getFileName().c_str()));
  }
  return *m_oFStreamPtr;
}

void TimingDataWriter::writeHeader()
{
  if (!m_haveWrittenHeader) {
    m_haveWrittenHeader = true;
    m_pTimers->writeHeader(getOStream());
  } 
}

void TimingDataWriter::appendData()
{
  writeHeader();
  m_pTimers->appendData(getOStream());
}

//======== MpiWTimers =================

MpiWTimers::MpiWTimers() : m_timerMap(), m_fileNameWriterMap()
{
}

/*
void MpiWTimers::createTimer(const std::string &timerName)
{
}
*/

MpiWTimer *MpiWTimers::findTimer(const std::string &timerName)
{
  NameMpiWTimerMap::iterator it = m_timerMap.find(timerName);
  if (it != m_timerMap.end()) {
    return &(it->second);
  }
  return NULL;
}

const MpiWTimer *MpiWTimers::findTimer(const std::string &timerName) const
{
  NameMpiWTimerMap::const_iterator it = m_timerMap.find(timerName);
  if (it != m_timerMap.end()) {
    return &(it->second);
  }
  return NULL;
}

MpiWTimer &MpiWTimers::findOrCreateTimer(const std::string &timerName)
{
  NameMpiWTimerMap::iterator it = m_timerMap.find(timerName);
  if (it == m_timerMap.end()) {
    it = (m_timerMap.insert(NameMpiWTimerMap::value_type(timerName, MpiWTimer(timerName)))).first;
    it->second.setStart(MPI_Wtime());
  }
  return it->second;
}

bool MpiWTimers::timerExists(const std::string &timerName) const
{
  return (m_timerMap.find(timerName) != m_timerMap.end());
}

void MpiWTimers::start(const std::string &timerName)
{
  findOrCreateTimer(timerName).setStart(MPI_Wtime());
}

void MpiWTimers::stop(const std::string &timerName, bool elapseIsStopMinusStart)
{
  const double time = MPI_Wtime();
  findOrCreateTimer(timerName).setStop(time, elapseIsStopMinusStart);
}

double MpiWTimers::getTiming(const std::string &timerName) const
{
  if (timerExists(timerName)) {
     findTimer(timerName)->getTiming();
  }
  return -1.0;
}

 void MpiWTimers::clear()
{
  m_timerMap.clear();
  m_fileNameWriterMap.clear();
}

void MpiWTimers::pause(const std::string &name)
{
  const double time = MPI_Wtime();
  findOrCreateTimer(name).pause(time);
}

void MpiWTimers::resume(const std::string &name)
{
  const double time = MPI_Wtime();
  findOrCreateTimer(name).resume(time);
}

void MpiWTimers::zeroise(const std::string &name)
{
  findOrCreateTimer(name).zeroise();
}

void MpiWTimers::zeroise()
{
  NameMpiWTimerMap::const_iterator it = m_timerMap.begin();
  for (; it != m_timerMap.end(); it++) 
  {
    zeroise(it->first);
  }
}

void MpiWTimers::writeHeader(std::ostream &oStream)
{
  NameMpiWTimerMap::const_iterator it = m_timerMap.begin();
  if (it != m_timerMap.end()) {
    oStream << it->first;
    it++;
    for (; it != m_timerMap.end(); it++)
    {
      oStream << "," << it->first;
    }
    oStream << std::endl;
  }
}

TimingDataWriter &MpiWTimers::getWriter(const std::string &fileName)
{
  FileNameWriterMap::iterator it = m_fileNameWriterMap.find(fileName);
  if (it == m_fileNameWriterMap.end()) {
    it = m_fileNameWriterMap.insert(
        FileNameWriterMap::value_type(
          fileName,
          TimingDataWriter(fileName, *this)
        )
      ).first;
  }
  return it->second;
}

void MpiWTimers::appendData(const std::string &fileName)
{
  getWriter(fileName).appendData();
}

void MpiWTimers::appendData(std::ostream &oStream)
{
  NameMpiWTimerMap::const_iterator it = m_timerMap.begin();
  if (it != m_timerMap.end()) {
    oStream
      << std::fixed << std::showpoint << std::setprecision(6) << std::setw(8)
      << it->second.getTiming();
    it++;
    for (; it != m_timerMap.end(); it++)
    {
      oStream
        << ", "
        << std::fixed << std::showpoint << std::setprecision(6) << std::setw(8)
        << it->second.getTiming();
    }
    oStream << "\n";
  }
}
